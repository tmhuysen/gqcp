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
#define BOOST_TEST_MODULE "ShellSet"

#include <boost/test/unit_test.hpp>

#include "Basis/ScalarBasis/GTOShell.hpp"
#include "Basis/ScalarBasis/GTOBasisSet.hpp"
#include "Basis/ScalarBasis/ShellSet.hpp"


BOOST_AUTO_TEST_CASE ( constructor_basisset ) {

    // Create an STO-3G basisset on (a weird geometry of) H2O
    GQCP::Nucleus h1 (1,  0.0, 0.0, 0.0);
    GQCP::Nucleus o  (8,  0.0, 0.0, 1.0);
    GQCP::Nucleus h2 (1,  0.0, 0.0, 2.0);

    GQCP::ShellSet<GQCP::GTOShell> ref_shellset {
        GQCP::GTOShell(0, h1, {  3.42525091,  0.62391373, 0.16885540}, { 0.15432897, 0.53532814, 0.44463454}, false),
        GQCP::GTOShell(0, o,  {130.7093200,  23.8088610,  6.4436083},  { 0.15432897, 0.53532814, 0.44463454}, false),
        GQCP::GTOShell(0, o,  {  5.0331513,   1.1695961,  0.3803890},  {-0.09996723, 0.39951283, 0.70011547}, false),
        GQCP::GTOShell(1, o,  {  5.0331513,   1.1695961,  0.3803890},  { 0.15591627, 0.60768372, 0.39195739}, false),
        GQCP::GTOShell(0, h2, {  3.42525091,  0.62391373, 0.16885540}, { 0.15432897, 0.53532814, 0.44463454}, false)
    };


    GQCP::Molecule h2o ({h1, o, h2});
    const auto shellset = GQCP::GTOBasisSet("STO-3G").generate(h2o);
    BOOST_CHECK(ref_shellset.asVector() == shellset.asVector());
}
