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

#include <dlfcn.h>

#include "polymake/client.h"
#include "polymake/ideal/singular.h"
#include "polymake/ListMatrix.h"
#include "polymake/Map.h"
#include "polymake/Ring.h"
#include "polymake/pair.h"

#include <libsingular.h>
#include "singular/stairc.h"

namespace polymake { namespace ideal {

namespace {

bool singular_initialized = 0;
unsigned int ringidcounter = 0;

// Mapping Polymake rings to their Singular handles.
Map<std::pair<Ring<>::id_type, Matrix<int> >, idhdl> singular_ring_map;
// Storing the handles for the Singular functions globally.
Map<std::string, idhdl> singular_function_map;

void singular_error_handler(const char* error)
{
   throw std::runtime_error(error);
}

// Initialize Singular:
void init_singular() 
{
   if(singular_initialized)
      return;

   cout << "*** trying to determine path of libsingular ***" << endl;

   Dl_info dli;
   if (!dladdr((void*)&siInit,&dli)) {
      throw std::runtime_error("*** could not find symbol from libsingular ***");
   }
   cout << "*** found libsingular at: " << dli.dli_fname << endl;
   
   cout << "*** loading singular ***" << endl;
   char* cpath = omStrDup(dli.dli_fname);
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
   if (err) {
      printf("interpreter returns %d\n",err);
      throw std::runtime_error("*** singular loading failed ***");
   }
   cout << "*** singular done ***" << endl;
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
idhdl check_ring(const Ring<> r, const Matrix<int> order){
   init_singular();
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
      int ord_size = 2;
      int *ord = (int*)omalloc0((ord_size+1)*sizeof(int));
      int *block0 = (int*)omalloc0((ord_size+1)*sizeof(int));
      int *block1 = (int*)omalloc0((ord_size+1)*sizeof(int));
      ord[0] = ringorder_M;
      ord[1] = ringorder_c;
      block0[0] = 1;
      block1[0] = nvars;
      int **wvhdl = (int**)omalloc0((ord_size+1)*sizeof(int*));
      wvhdl[0] = (int*)omalloc0(nvars*nvars*sizeof(int));
      for(int i =0; i<nvars; i++){
         for(int j = 0; j<nvars; j++){
            wvhdl[0][i*nvars+j] = order(i,j);
         }
      }
      cout << "Monomial ordering written to Singular." << endl;
      // Create Singular ring:
      ring r = rDefault(0,nvars,n,ord_size,ord,block0,block1,wvhdl);
      char* ringid = (char*) malloc(2+sizeof(unsigned int));
      sprintf(ringid,"R_%0u",ringidcounter++); 
      // Create handle for ring:
      idhdl newRingHdl=enterid(ringid,0,RING_CMD,&IDROOT,FALSE);
      IDRING(newRingHdl)=r;
      // Store handle:
      singular_ring_map[p] = newRingHdl;

   }
   // Make it the default ring, also for the interpeter
   rSetHdl(singular_ring_map[p]);
   return singular_ring_map[p];
}

// If no monomial ordering is given:
idhdl check_ring(const Ring<> r){
   int nvars = r.n_vars();
   Matrix<int> ord = unit_matrix<int>(nvars);
   return check_ring(r, ord);
}

idhdl check_ring(idhdl singRing) {
   init_singular();
   rSetHdl(singRing);
   return singRing;
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

Polynomial<> convert_poly_to_Polynomial(const poly q, const Ring<>& r){

   int n = r.n_vars();
   poly p = pCopy(q);
   ListMatrix<Vector<int> > exponents(0,n);
   std::vector<Rational> coefficients;
   while(p != NULL){
      number c = pGetCoeff(p);
      coefficients.push_back(convert_number_to_Rational(c, currRing));
      Vector<int> monomial(n);
      for(int i = 1; i<=n; i++){
         monomial[i-1] = pGetExp(p, i);
      }
      exponents /= monomial;
      pIter(p);
   }
   return Polynomial<>(exponents, coefficients, r);
}

poly convert_Polynomial_to_poly(const Polynomial<>& mypoly)
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

   void create_singIdeal(const Array<Polynomial<> >& gens) {
      int npoly = gens.size();
      singIdeal = idInit(npoly,1);
      int j = 0;
      // Converting monomials as described in libsing-test2.cc.
      for(Entire<Array<Polynomial<> > >::const_iterator mypoly = entire(gens); !mypoly.at_end(); ++mypoly, ++j) {
         poly p = convert_Polynomial_to_poly(*mypoly);
         //cout << "poly: " << p_String(p,singRing,singRing) << endl;
         singIdeal->m[j]=p_Copy(p,currRing);
      }
   }

public:
  
   // Constructing singIdeal from the generators:
   SingularIdeal_impl(const Array<Polynomial<> >& gens, const Matrix<int>& ord)
   {
      singRing = check_ring(gens[0].get_ring(),ord);
      if(!gens.size())
         throw std::runtime_error("Ideal has no generators.");
      create_singIdeal(gens);
   }
   
   SingularIdeal_impl(const Array<Polynomial<> >& gens, idhdl r){
      singRing = check_ring(r);
      create_singIdeal(gens);
   }

   SingularIdeal_impl(::ideal i, idhdl r)
   {
      singIdeal=idCopy(i);
      singRing=r;
   }

   ~SingularIdeal_impl() 
   {
      /* This does not work yet:
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
   void groebner() 
   {
      check_ring(singRing); 

      ::ideal res;
      res = kStd(singIdeal,NULL,testHomog,NULL);
      //cout << "DONE COMPUTING std basis" << endl;
      id_Delete(&singIdeal,currRing);
      singIdeal = res;
   }

   // Compute the dimension of an ideal.
   int dim() {
      check_ring(singRing); 
      return scDimInt(singIdeal, NULL);
   }
   

   SingularIdeal_wrap* initial_ideal() const {
      check_ring(singRing); 
      ::ideal res;
      res = idHead(singIdeal);
      return new SingularIdeal_impl(res, singRing);
   }

   // Compute the radical of an ideal using primdec.lib from Singular
   SingularIdeal_wrap* radical() const {
      check_ring(singRing); 
      sleftv arg;
      memset(&arg,0,sizeof(arg));

      idhdl radical=get_singular_function("radical");
      
      arg.rtyp=IDEAL_CMD;
      arg.data=(void *)singIdeal;
      // call radical
      leftv res=iiMake_proc(radical,NULL,&arg);
      if (res==NULL) {
         errorreported = 0;
         throw std::runtime_error("radical returned an error");
      }
      return new SingularIdeal_impl((::ideal) (res->Data()), singRing);
   }

   Array<SingularIdeal_wrap*> primary_decomposition() const {
      check_ring(singRing);
      idhdl primdecSY = get_singular_function("primdecSY");
      sleftv arg;
      memset(&arg,0,sizeof(arg));
      arg.rtyp=IDEAL_CMD;
      arg.data=(void *)idCopy(singIdeal);
      // call primdecSY
      leftv res=iiMake_proc(primdecSY,NULL,&arg);
      if(res->Typ() == LIST_CMD){
         lists L = (lists)res->Data();
         Array<SingularIdeal_wrap*> result(L->nr+1);
         for(int j=0; j<=L->nr; j++){
            lists LL = (lists)L->m[j].Data();
            if(LL->m[0].Typ() == IDEAL_CMD){
               result[j] = new SingularIdeal_impl((::ideal) (LL->m[0].Data()),singRing);
            } else {
               throw std::runtime_error("Something went wrong for the primary decomposition");
            }
         }
         return result;
      } else {
         throw std::runtime_error("Something went wrong for the primary decomposition");
      }
   }

   Polynomial<> reduce(const Polynomial<>& p, const Ring<>& r) const {
      check_ring(singRing);
      poly q = kNF(singIdeal, NULL, convert_Polynomial_to_poly(p));
      return convert_poly_to_Polynomial(q,r);
   }

   // Converting singIdeal generators to an array of Polymake polynomials.
   Array<Polynomial<> > polynomials(const Ring<>& r) const
   {
      check_ring(singRing);

      int numgen = IDELEMS(singIdeal);
      std::vector<Polynomial<> > polys;

      for(int j = 0; j<numgen; j++) {
         if(singIdeal->m[j] != NULL){
            polys.push_back(convert_poly_to_Polynomial(singIdeal->m[j],r));
         }
      }
      return Array<Polynomial<> >(polys);
   }
   
   static SingularIdeal_impl* quotient(const SingularIdeal_impl* I, const SingularIdeal_impl* J);
};


SingularIdeal_impl* SingularIdeal_impl::quotient(const SingularIdeal_impl* I, const SingularIdeal_impl* J){
   // The first true indicates, that we receive a standard basis of I,
   // the second one that we want the output to be an ideal.
   ::ideal quot = idQuot(I->singIdeal, J->singIdeal, true, true);
   return new SingularIdeal_impl(quot,I->singRing);
}

} // end anonymous namespace

SingularIdeal_wrap* SingularIdeal_wrap::create(const Array<Polynomial<> >& gens, const Matrix<int>& ord) 
{
   return new SingularIdeal_impl(gens,ord);
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
   
   const Array<Polynomial<> > gensI = I.give("GROEBNER.BASIS");
   const Matrix<int> ord = I.give("GROEBNER.ORDER_MATRIX");
   const idhdl sri = check_ring(ri,ord);
   const Array<Polynomial<> > gensJ = J.give("GENERATORS");
   
   SingularIdeal_impl implI(gensI,sri);
   SingularIdeal_impl implJ(gensJ,sri);

   SingularIdeal_impl* quotimpl = SingularIdeal_impl::quotient(&implI,&implJ);

   perl::Object res("Ideal");
   res.take("RING") << ri;
   res.take("GENERATORS") << quotimpl->polynomials(ri);
   delete quotimpl;
   return res;
}


UserFunction4perl("# @category Algebra"
                  "# Computes an ideal quotient via SINGULAR"
                  "# @param Ideal I"
                  "# @param Ideal J"
                  "# @return Ideal",
                  &quotient, "quotient(Ideal, Ideal)");

} }


