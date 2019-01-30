// This file is part of GQCG-gqcp.
// 
// Copyright (C) 2017-2018  the GQCG developers
// 
// GQCG-gqcp is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// GQCG-gqcp is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with GQCG-gqcp.  If not, see <http://www.gnu.org/licenses/>.
// 
#define BOOST_TEST_MODULE "NumericalDerivator"


#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain

#include "utilities/NumericalDerivator.hpp"



BOOST_AUTO_TEST_CASE ( derive_xcubed ) {

    GQCP::UnaryFunction xcubed = [](double x) { return pow(x, 3);};
    GQCP::NumericalDerivator<4> derivator (xcubed, 0, 0.001);

    BOOST_CHECK(derivator.get_derivative(0) == 0);
    BOOST_CHECK(std::abs(derivator.get_derivative(1) - 1e-06) < 1e-10);
    BOOST_CHECK(std::abs(derivator.get_derivative(2) - 0.006) < 1e-10);
    BOOST_CHECK(std::abs(derivator.get_derivative(3) - 6) < 1e-10);
    BOOST_CHECK(std::abs(derivator.get_derivative(4) - 0) < 1e-10);
}