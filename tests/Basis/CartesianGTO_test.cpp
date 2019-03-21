// This file is part of GQCG-gqcp.
// 
// Copyright (C) 2017-2019  the GQCG developers
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
#define BOOST_TEST_MODULE "CartesianGTO"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise clang++ will complain


#include "Basis/CartesianGTO.hpp"


BOOST_AUTO_TEST_CASE ( constructor ) {

    auto exponents = GQCP::CartesianExponents({0, 0, 0});
    GQCP::Vector<double, 3> center = GQCP::Vector<double, 3>::Zero();

    BOOST_CHECK_THROW(GQCP::CartesianGTO gto (-1.0, exponents, center), std::invalid_argument);  // exponent in the exponential cannot be negative
}


BOOST_AUTO_TEST_CASE ( calculateNormalizationFactor ) {

    auto exponents1 = GQCP::CartesianExponents({1, 0, 1});
    GQCP::Vector<double, 3> center = GQCP::Vector<double, 3>::Zero();
    GQCP::CartesianGTO gto1 (1.0, exponents1, center);

    BOOST_CHECK(std::abs(gto1.calculateNormalizationFactor() - 2.8508218814) < 1.0e-09);  // 'manual' calculation

    auto exponents2 = GQCP::CartesianExponents({1, 2, 3});
    GQCP::CartesianGTO gto2 (2.5, exponents2, center);

    BOOST_CHECK(std::abs(gto2.calculateNormalizationFactor() - 211.2315772257) < 1.0e-09);  // 'manual' calculation
}


BOOST_AUTO_TEST_CASE ( operator_call ) {

    auto exponents1 = GQCP::CartesianExponents({1, 0, 1});
    GQCP::Vector<double, 3> center1;
    center1 << 1.0, 0.0, -0.5;
    GQCP::Vector<double, 3> r1;
    r1 << 0.0, 1.0, 0.0;
    GQCP::CartesianGTO gto1 (1.0, exponents1, center1);

    BOOST_CHECK(std::abs(gto1(r1) - (-0.1502372078)) < 1.0e-09);  // 'manual' calculation


    auto exponents2 = GQCP::CartesianExponents({1, 2, 3});
    GQCP::Vector<double, 3> center2;
    center2 << 0.0, 1.0, 0.0;
    GQCP::Vector<double, 3> r2;
    r2 << -1.0, -1.0, 1.0;
    GQCP::CartesianGTO gto2 (2.5, exponents2, center2);

    BOOST_CHECK(std::abs(gto2(r2) - (-0.0002584649185)) < 1.0e-09);  // 'manual' calculation
}


BOOST_AUTO_TEST_CASE ( calculateDerivative ) {

    auto exponents1 = GQCP::CartesianExponents({1, 0, 1});
    GQCP::Vector<double, 3> center = GQCP::Vector<double, 3>::Zero();
    GQCP::CartesianGTO gto1 (1.0, exponents1, center);

    // GTO1 - x-component
    auto x_derivative1 = gto1.calculateDerivative(GQCP::CartesianDirection::x);
    double ref_coeff1_x1 = -2.0;
    double ref_coeff2_x1 = 1.0;
    auto ref_exp1_x1 = GQCP::CartesianExponents({2, 0, 1});
    auto ref_exp2_x1 = GQCP::CartesianExponents({0, 0, 1});

    BOOST_CHECK(std::abs(x_derivative1.get_coefficients()[0] - ref_coeff1_x1) < 1.0e-12);
    BOOST_CHECK(std::abs(x_derivative1.get_coefficients()[1] - ref_coeff2_x1) < 1.0e-12);
    BOOST_CHECK(x_derivative1.get_functions()[0].get_cartesian_exponents() == ref_exp1_x1);
    BOOST_CHECK(x_derivative1.get_functions()[1].get_cartesian_exponents() == ref_exp2_x1);


    // GTO1 - y-component
    auto y_derivative1 = gto1.calculateDerivative(GQCP::CartesianDirection::y);
    double ref_coeff1_y1 = -2.0;
    auto ref_exp1_y1 = GQCP::CartesianExponents({1, 1, 1});

    BOOST_CHECK(std::abs(y_derivative1.get_coefficients()[0] - ref_coeff1_y1) < 1.0e-12);
    BOOST_CHECK(y_derivative1.get_functions()[0].get_cartesian_exponents() == ref_exp1_y1);


    auto exponents2 = GQCP::CartesianExponents({1, 2, 3});
    GQCP::CartesianGTO gto2 (2.5, exponents2, center);

    // GTO2 - x-component
    auto x_derivative2 = gto2.calculateDerivative(GQCP::CartesianDirection::x);
    double ref_coeff1_x2 = -5.0;
    double ref_coeff2_x2 = 1.0;
    auto ref_exp1_x2 = GQCP::CartesianExponents({2, 2, 3});
    auto ref_exp2_x2 = GQCP::CartesianExponents({0, 2, 3});

    BOOST_CHECK(std::abs(x_derivative2.get_coefficients()[0] - ref_coeff1_x2) < 1.0e-12);
    BOOST_CHECK(std::abs(x_derivative2.get_coefficients()[1] - ref_coeff2_x2) < 1.0e-12);
    BOOST_CHECK(x_derivative2.get_functions()[0].get_cartesian_exponents() == ref_exp1_x2);
    BOOST_CHECK(x_derivative2.get_functions()[1].get_cartesian_exponents() == ref_exp2_x2);

    // GTO2 - z-component
    auto z_derivative2 = gto2.calculateDerivative(GQCP::CartesianDirection::z);
    double ref_coeff1_z2 = -5.0;
    double ref_coeff2_z2 = 3.0;
    auto ref_exp1_z2 = GQCP::CartesianExponents({1, 2, 4});
    auto ref_exp2_z2 = GQCP::CartesianExponents({1, 2, 2});

    BOOST_CHECK(std::abs(z_derivative2.get_coefficients()[0] - ref_coeff1_z2) < 1.0e-12);
    BOOST_CHECK(std::abs(z_derivative2.get_coefficients()[1] - ref_coeff2_z2) < 1.0e-12);
    BOOST_CHECK(z_derivative2.get_functions()[0].get_cartesian_exponents() == ref_exp1_z2);
    BOOST_CHECK(z_derivative2.get_functions()[1].get_cartesian_exponents() == ref_exp2_z2);
}
