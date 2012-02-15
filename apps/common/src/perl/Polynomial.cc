/* Copyright (c) 1997-2010
   Ewgenij Gawrilow, Michael Joswig (Technische Universitaet Darmstadt, Germany)
   http://www.polymake.de

   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2, or (at your option) any
   later version: http://www.gnu.org/licenses/gpl.txt.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
*/

///==== this line controls the automatic file splitting: max.instances=40

#include "polymake/client.h"
#include "polymake/Polynomial.h"
#include "polymake/Rational.h"
#include "polymake/Matrix.h"
#include "polymake/IncidenceMatrix.h"
#include "polymake/Vector.h"
#include "polymake/Ring.h"
#include "polymake/SparseMatrix.h"
namespace polymake { namespace common {
///==== Automatically generated contents follow.    Please do not delete this line. ====
   template <typename T0, typename T1, typename T2, typename T3>
   FunctionInterface4perl( new_X_X_X, T0,T1,T2,T3 ) {
      perl::Value arg0(stack[1]), arg1(stack[2]), arg2(stack[3]);
      WrapperReturnNew(T0, (arg0.get<T1>(), arg1.get<T2>(), arg2.get<T3>()) );
   };

   FunctionInstance4perl(new_X_X_X, Polynomial< Rational, int >, perl::Canned< const pm::ColChain<pm::Transposed<pm::Matrix<int> > const&, pm::Matrix<int> const&> >, perl::Canned< const Vector< Rational > >, perl::Canned< const Ring< Rational, int > >);
   FunctionInstance4perl(new_X_X_X, Polynomial< Rational, int >, perl::Canned< const SparseMatrix< int, NonSymmetric > >, perl::Canned< const Vector< Rational > >, perl::Canned< const Ring< Rational, int > >);
///==== Automatically generated contents end here.  Please do not delete this line. ====
} }
