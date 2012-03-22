/* Copyright (c) 2012
   by authors as mentioned on:
   https://github.com/lkastner/polymake_algebra/wiki/Authors

   Project home:
   https://github.com/lkastner/polymake_algebra

   For licensing we cite the original Polymake code:

   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2, or (at your option) any
   later version: http://www.gnu.org/licenses/gpl.txt.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
*/


#include "polymake/client.h"
#include "polymake/ideal/singular.h"
#include "polymake/ListMatrix.h"
#include "polymake/Map.h"
#include "polymake/Ring.h"
#include "polymake/pair.h"

#include <libsingular.h>
#include "singular/stairc.h"

namespace polymake { namespace ideal {

//namespace {

int singular_initialized = 0;

// Mapping Polymake rings to their Singular handles.
Map<std::pair<Ring<>::id_type, Matrix<int> >, idhdl> singular_ring_map;
// Storing the handles for the Singular functions globally.
Map<std::string, idhdl> singular_function_map;

void singular_error_handler(const char* error)
{
   throw std::runtime_error(error);
}

// Initialize Singular:
void init_singular(const std::string& path) 
{
   if(singular_initialized)
   {
      cout << "Singular has already been initialized, restart polymake to retry" << endl;
      return;
   }
   
   std::string p = path+"/Singular/libsingular.so";
   char* cpath = omStrDup(p.c_str());
   siInit(cpath);
   WerrorS_callback = &singular_error_handler;
   
   sleftv arg,r1,r2;
   
   // load the singular library primdec.lib:
   // Is it really necessary to do this here?
   // Should we not only load libraries on demand?
   memset(&arg,0,sizeof(arg));
   memset(&r1,0,sizeof(r1));
   memset(&r2,0,sizeof(r2));
   arg.rtyp=STRING_CMD;
   arg.data=omStrDup("primdec.lib");
   r2.rtyp=LIB_CMD;
   int err=iiExprArith2(&r1,&r2,'(',&arg);
   if (err) printf("interpreter returns %d\n",err);
   cout << "init_singular done" << endl;
   singular_initialized = 1;
}

// This function returns the idhdl of the function to be used.
// If the handle does not exist the function is looked up and the handle
// is created.
// All handles are stored globally.
idhdl get_singular_function(std::string s) {
   if(!singular_function_map.exists(s)) {
      // now, get the procedure to call
      idhdl function=ggetid(s.c_str());
      if (function==NULL)
         throw std::runtime_error("singular function not found: "+s);
      singular_function_map[s] = function;
   }
   return singular_function_map[s];
}

// Returns the Singular equivalent for a Polymake ring.
// If the Singular ring does not exist, it is created and stored globally,
// indirectly, since it is contained in the handle.
// Also the idhdl of the ring is created here.
ring check_ring(const Ring<> r, const Matrix<int> order){
   Ring<>::id_type id = r.id();
   std::pair<Ring<>::id_type, Matrix<int> > p(id, order);
   if(!singular_ring_map.exists(p)){
      int nvars = r.n_vars();
      if(nvars == 0) 
         throw std::runtime_error("Given ring is not a polynomial ring.");
      // Create variables:
      char **n=(char**)omalloc(nvars*sizeof(char*));
      for(int i=0; i<nvars; i++)
      {
         n[i] = omStrDup(r.names()[i].c_str());
      }
      // Create Singular ring:
      ring r = rDefault(0,nvars,n);
      char* ringid = (char*) malloc(2+sizeof(unsigned int));
      sprintf(ringid,"R-%0u",id);
      // Create handle for ring:
      idhdl newRingHdl=enterid(ringid,0,RING_CMD,&IDROOT,FALSE);
      IDRING(newRingHdl)=r;
      // Store handle:
      singular_ring_map[p] = newRingHdl;

   }
   // Make it the default ring, also for the interpeter
   rSetHdl(singular_ring_map[p]);
   return currRing;
}
// If no monomial ordering is given:
ring check_ring(const Ring<> r){
   int nvars = r.n_vars();
   Matrix<int> ord = unit_matrix<int>(nvars);
   return check_ring(r, ord);
}

void check_ring(idhdl singRing) {
   rSetHdl(singRing);
}

// Convert a Singular number to a GMP rational.
Rational convert_number_to_Rational(number n, ring ring)
{
   if(rField_is_Q(ring)){
      if(SR_HDL(n) & SR_INT){
         long l = SR_TO_INT(n);
         return Rational(l, 1);
      } else {
         switch(n->s){
            case 0:
               return Rational(n->z, n->n);
            case 1:
               return Rational(n->z, n->n);
            case 3:
               return Rational(n->z);
         }
      }
   }
   throw std::runtime_error("I can has number? :P");
}

// Convert a GMP rational to a Singular number.
number convert_Rational_to_number(const Rational& r)
{
   mpz_t num, denom;
   mpz_init(num);
   mpz_set(num,numerator(r).get_rep());
   mpz_init(denom);
   mpz_set(denom,denominator(r).get_rep());

   return nlInit2gmp(num,denom);
}

poly convert_Polynomial_to_poly(const Polynomial<> mypoly)
{
   poly p = pISet(0);
   for(Entire<Polynomial<>::term_hash>::const_iterator term = entire(mypoly.get_terms()); !term.at_end(); ++term)
   {
      poly monomial = pNSet(convert_Rational_to_number(term->second)); 
      for(int k = 0; k<term->first.dim(); k++)
      {
         pSetExp(monomial,k+1,term->first[k]); 
      }
      pSetm(monomial); 
      p = pSub(p, monomial);
   }
   return p;
}

class SingularIdeal_impl : public SingularIdeal_wrap {
private:
   ::ideal singIdeal;
   idhdl singRing;

public:
  
   // Constructing singIdeal from the generators:
   SingularIdeal_impl(const Array<Polynomial<> > gens) 
   {
      if (!singular_initialized)
         throw std::runtime_error("singular not yet initialized, call init_singular(Path)");
      //cout << "creating SingularIdeal_impl from Ideal" << endl;
      // Get Singular ring:
      ring singRing = check_ring(gens[0].get_ring());
      int npoly = gens.size();
      if(!npoly)
         throw std::runtime_error("Ideal has no generators.");
      
      singIdeal = idInit(npoly,1);
      int j = 0;
      // Converting monomials as described in libsing-test2.cc.
      for(Entire<Array<Polynomial<> > >::const_iterator mypoly = entire(gens); !mypoly.at_end(); ++mypoly, ++j) {
         poly p = convert_Polynomial_to_poly(*mypoly);
         //cout << "poly: " << p_String(p,singRing,singRing) << endl;
         singIdeal->m[j]=p_Copy(p,singRing);
      }
      //cout << "DONE CREATING singular object" << endl;
   }

   SingularIdeal_impl(::ideal i)
   {
      if (!singular_initialized)
         throw std::runtime_error("singular not yet initialized, call init_singular(Path)");
      //cout << "creating SingularIdeal_impl from singular stuff" << endl;
      singIdeal=i;
   }

   ~SingularIdeal_impl() 
   {
      /* FIXME
      cout << "SingularIdeal_impl cleaning up" <<endl;
      if(singRing!=NULL) {
         if(singIdeal!=NULL)
            id_Delete(&singIdeal,singRing);
         //rKill(singRing);
      }
      cout << "SingularIdeal_impl destroyed" << endl;
      */
   }

   // Compute a groebner basis of a Polymake ideal using Singular
   void groebner(const Ring<> r, const Matrix<int> order) 
   {
      ring singRing = check_ring(r, order); 

      ::ideal res;
      res = kStd(singIdeal,NULL,testHomog,NULL);
      //cout << "DONE COMPUTING std basis" << endl;
      id_Delete(&singIdeal,singRing);
      singIdeal = res;
   }

   // Compute the dimension of an ideal.
   int dim(const Ring<> r) {
      if (!singular_initialized)
         throw std::runtime_error("singular not yet initialized, call init_singular(Path)");
      check_ring(r); 
      return scDimInt(singIdeal, NULL);
   }
   
   // Compute the radical of an ideal using primdec.lib from Singular
   SingularIdeal_wrap* radical(const Ring<> r) const {
      if (!singular_initialized)
         throw std::runtime_error("singular not yet initialized, call init_singular(Path)");

      ring singRing = check_ring(r); 
      sleftv arg;
      memset(&arg,0,sizeof(arg));

      idhdl radical=get_singular_function("radical");
      
      arg.rtyp=IDEAL_CMD;
      arg.data=(void *)idCopy(singIdeal);
      // call radical
      leftv res=iiMake_proc(radical,NULL,&arg);
      if (res==NULL) {
         errorreported = 0;
         throw std::runtime_error("radical returned an error");
      }
      return new SingularIdeal_impl((::ideal) (res->Data()));
   }

   // Converting singIdeal generators to an array of Polymake polynomials.
   Array<Polynomial<> > polynomials(const Ring<> r) const
   {
      ring singRing = check_ring(r); 
      rChangeCurrRing(singRing);

      int numgen = IDELEMS(singIdeal);
      std::vector<Polynomial<> > polys;

      int n = rVar(singRing);
      for(int j = 0; j<numgen; j++) {
         if(singIdeal->m[j] != NULL){
            ListMatrix<Vector<int> > exponents(0,n);
            poly p = singIdeal->m[j];
            std::vector<Rational> coefficients;
            while(p != NULL){
               number c = pGetCoeff(p);
               coefficients.push_back(convert_number_to_Rational(c, singRing));
               Vector<int> monomial(n);
               for(int i = 1; i<=n; i++){
                  monomial[i-1] = pGetExp(p, i);
               }
               exponents /= monomial;
               pIter(p);
            }
            polys.push_back(Polynomial<>(exponents, coefficients, r));
//            cout << "converted: " << p_String(singIdeal->m[j],singRing,singRing)<<endl;
         }
      }
      return Array<Polynomial<> >(polys);
   }
   
   friend SingularIdeal_wrap* SingularIdeal_wrap::quotient(const SingularIdeal_wrap* I, const SingularIdeal_wrap* J);
};

SingularIdeal_wrap* SingularIdeal_wrap::quotient(const SingularIdeal_wrap* I, const SingularIdeal_wrap* J){
   const ::ideal sI = static_cast<const SingularIdeal_impl*>(I)->singIdeal;
   const ::ideal sJ = static_cast<const SingularIdeal_impl*>(J)->singIdeal;
   // The first true indicates, that we receive a standard basis of I,
   // the second one that we want the output to be an ideal.
   ::ideal quot = idQuot(sI, sJ, true, true);
   return new SingularIdeal_impl(quot);
}

SingularIdeal_wrap* SingularIdeal_wrap::create(const Array<Polynomial<> > gens, Matrix<int> order) 
{
   return new SingularIdeal_impl(gens);
}

perl::Object quotient(perl::Object I, perl::Object J)
{
   Ring<> ri;
   Ring<> rj;
   I.give("RING") >> ri;
   J.give("RING") >> rj;
   if (ri.id() != rj.id())
      throw std::runtime_error("Ideals of different rings");

   check_ring(ri);
   
   // FIXME which one should be STANDARD ?
   const Array<Polynomial<> > gensI = I.give("STANDARD");
   const Array<Polynomial<> > gensJ = J.give("STANDARD | GENERATORS");
   
   SingularIdeal_impl wrapI(gensI);
   SingularIdeal_impl wrapJ(gensJ);

   SingularIdeal_wrap* quotwrap = SingularIdeal_wrap::quotient(&wrapI,&wrapJ);

   perl::Object res("Ideal");
   res.take("RING") << ri;
   res.take("GENERATORS") << quotwrap->polynomials(ri);
   delete quotwrap;
   return res;
}


UserFunction4perl("# @category Algebra"
                  "# Computes an ideal quotient via SINGULAR"
                  "# @param Ideal I"
                  "# @param Ideal J"
                  "# @return Ideal",
                  &quotient, "quotient(Ideal, Ideal)");

UserFunction4perl("# @category Other"
                  "# @param String path Path to the singular directory",
                  &init_singular, "init_singular($)");

} }


