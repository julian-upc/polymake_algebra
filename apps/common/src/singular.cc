#include "polymake/client.h"
#include "polymake/common/singular.h"
#include "polymake/ListMatrix.h"
#include "polymake/Map.h"
#include "polymake/Ring.h"

#include <libsingular.h>


namespace polymake { namespace common {

//namespace {

int singular_initialized = 0;

Map<Ring<>::id_type, idhdl> singular_ring_map;
Map<std::string, idhdl> singular_function_map;

void singular_error_handler(const char* error)
{
   throw std::runtime_error(error);
}

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

ring check_ring(Ring<> r){
   Ring<>::id_type id = r.id();
   if(!singular_ring_map.exists(id)){
      int nvars = r.n_vars();
      if(nvars == 0) 
         throw std::runtime_error("Given ring is not a polynomial ring.");
      char **n=(char**)omalloc(nvars*sizeof(char*));
      for(int i=0; i<nvars; i++)
      {
         n[i] = omStrDup(r.names()[i].c_str());
      }
      ring r = rDefault(0,nvars,n);
      // make it the default ring, also for the interpeter
      char* ringid = (char*) malloc(2+sizeof(unsigned int));
      sprintf(ringid,"R-%0u",id);
      idhdl newRingHdl=enterid(ringid,0,RING_CMD,&IDROOT,FALSE);
      IDRING(newRingHdl)=r;
      singular_ring_map[id] = newRingHdl;

   }
   rSetHdl(singular_ring_map[id]);
   return currRing;
}

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

number convert_Rational_to_number(const Rational& r)
{
   mpz_t num, denom;
   mpz_init(num);
   mpz_set(num,numerator(r).get_rep());
   mpz_init(denom);
   mpz_set(denom,denominator(r).get_rep());

   return nlInit2gmp(num,denom);
}

class SingularIdeal_impl : public SingularIdeal_wrap {
private:
   ideal singIdeal;

public:
   /*SingularIdeal_impl() 
   {
      if (!singular_initialized)
         throw std::runtime_error("singular not yet initialized, call init_singular(Path)");
      cout << "creating empty SingularIdeal_impl" << endl;
      singIdeal=NULL;
   }*/
   
   SingularIdeal_impl(const Array<Polynomial<> > gens) 
   {
      if (!singular_initialized)
         throw std::runtime_error("singular not yet initialized, call init_singular(Path)");
      //cout << "creating SingularIdeal_impl from Ideal" << endl;
      
      ring singRing = check_ring(gens[0].get_ring());
      int npoly = gens.size();
      if(!npoly)
         throw std::runtime_error("Ideal has no generators.");
      
      singIdeal = idInit(npoly,1); // Richtig?
      int j = 0;
      for(Entire<Array<Polynomial<> > >::const_iterator mypoly = entire(gens); !mypoly.at_end(); ++mypoly, ++j) {
         poly p = pISet(0);
         
         for(Entire<Polynomial<>::term_hash>::const_iterator term = entire(mypoly->get_terms()); !term.at_end(); ++term)
         {
            poly monomial = pNSet(convert_Rational_to_number(term->second)); // Geht das SetCoeff so?
            
            for(int k = 0; k<term->first.dim(); k++)
            {
               pSetExp(monomial,k+1,term->first[k]); // Hier war der Ring drinne, muss er wieder rein?
            }
            pSetm(monomial); // Ring noetig?
            p = pSub(p, monomial);
         }
         //cout << "poly: " << p_String(p,singRing,singRing) << endl;
         singIdeal->m[j]=p_Copy(p,singRing);
      }
      //cout << "DONE CREATING singular object" << endl;
   }

   SingularIdeal_impl(ideal i)
   {
      if (!singular_initialized)
         throw std::runtime_error("singular not yet initialized, call init_singular(Path)");
      //cout << "creating SingularIdeal_impl from singular stuff" << endl;
      singIdeal=i;
   }

   ~SingularIdeal_impl() 
   {
      cout << "SingularIdeal_impl cleaning up" <<endl;
      /*if(singRing!=NULL) {
         if(singIdeal!=NULL)
            id_Delete(&singIdeal,singRing);
         //rKill(singRing);
      }*/
      cout << "SingularIdeal_impl destroyed" << endl;
   }

   // Compute a groebner basis of a Polymake ideal using Singular
   void std(const Ring<> r) 
   {
      ring singRing = check_ring(r); 

      ideal res;
      res = kStd(singIdeal,NULL,testHomog,NULL);
      //cout << "DONE COMPUTING std basis" << endl;
      id_Delete(&singIdeal,singRing);
      singIdeal = res;
      // check if singIdeal exists
      // set singulardefaultring
      // call groebner
      // create polymake ideal maybe with singRing and singIdeal
   }

   SingularIdeal_wrap* radical(const Ring<> r) const {
      if (!singular_initialized)
         throw std::runtime_error("singular not yet initialized, call init_singular(Path)");

      ring singRing = check_ring(r); 
      // a more advanced procedure call
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
      return new SingularIdeal_impl((ideal) (res->Data()));
   }

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

};

//}

SingularIdeal_wrap* SingularIdeal_wrap::create(const Array<Polynomial<> > gens) 
{
   return new SingularIdeal_impl(gens);
}


UserFunction4perl("# @category Other"
                  "# @param String path Path to the singular directory",
                  &init_singular, "init_singular($)");

} }


