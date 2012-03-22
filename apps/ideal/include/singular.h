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
#include "polymake/Array.h"
#include "polymake/Ring.h"
#include "polymake/Polynomial.h"
#include "polymake/internal/shared_object.h"

namespace polymake { namespace ideal {

class SingularIdeal;


class SingularIdeal_wrap {
public:
   virtual ~SingularIdeal_wrap() { /* FIXME cout << "SingularIdeal destroyed" << endl; */ }

   virtual void std(const Ring<> r) = 0;
   
   virtual int dim(const Ring<> r) = 0;

   virtual SingularIdeal_wrap* radical(const Ring<> r) const = 0;

   virtual Array<Polynomial<> > polynomials(const Ring<> r) const = 0;
   
   static SingularIdeal_wrap* create(const Array<Polynomial<> > gens);

   static SingularIdeal_wrap* quotient(const SingularIdeal_wrap* I, const SingularIdeal_wrap* J);
};

class SingularIdeal {
private:
   SingularIdeal_wrap* singIdeal;

public:
   SingularIdeal(const Array<Polynomial<> > gens) {
      singIdeal = SingularIdeal_wrap::create(gens);
   }

   SingularIdeal(SingularIdeal_wrap* sI) {
      singIdeal = sI;
   }

   int dim(const Ring<> r) const  {
      return singIdeal->dim(r);
   }

   void std(const Ring<> r) const  {
      singIdeal->std(r);
   }

   SingularIdeal radical(const Ring<> r) const  {
      return SingularIdeal(singIdeal->radical(r));
   }
   
   Array<Polynomial<> > polynomials(const Ring<> r) const {
      return singIdeal->polynomials(r);
   }
   
   friend perl::Object quotient(perl::Object I, perl::Object J);
};


/*
   Ideal& operator+=(const Ideal& I) {
      if(this->get_ring() != I.get_ring()) throw std::runtime_error("Ideals of different rings.");
      append(I.size(),I.begin());
      return static_cast<Ideal&>(*this);
   }

   friend Ideal operator+(const Ideal& i1, const Ideal& i2) {
      if(i1.get_ring() != i2.get_ring()) throw std::runtime_error("Ideals of different rings.x");
      Ideal result = i1;
      result+=i2;
      return result;
   }

   Ideal groebner()
   {
      if(singObj == NULL) {
         //Ideal* writable = const_cast<Ideal*>(this);
         singObj = SingularIdeal::create(this);
      }
      SingularIdeal* basis = singObj->groebner();
      return Ideal(basis->polynomials(this->get_ring()),basis);

   }

};
*/


} }


