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
#define BOOST_TEST_MODULE "FFCI"


#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain

#include "HamiltonianBuilder/SelectedCI.hpp"
#include "HamiltonianBuilder/FCI.hpp"
#include "HamiltonianBuilder/FFCI.hpp"
#include "HamiltonianBuilder/DOCI.hpp"

#include "HamiltonianParameters/HamiltonianParameters.hpp"
#include "Molecule.hpp"




BOOST_AUTO_TEST_CASE ( FFCI_constructor ) {

    // Check if a correct constructor works
    GQCP::ProductFockSpace product_fock_space (5, 3, 3);
    GQCP::SelectedFockSpace fock_space (product_fock_space);
    BOOST_CHECK_NO_THROW(GQCP::SelectedCI selected_ci (fock_space));
}


BOOST_AUTO_TEST_CASE ( SelectedCI_vs_FFCI ) {

    // Create H-chain HamiltonianParameters to test results
    size_t K = 4;
    GQCP::Molecule H5 = GQCP::Molecule::HChain(K, 1.1);
    auto random_hamiltonian_parameters = GQCP::HamiltonianParameters::Molecular(H5, "STO-3G");

    // Create compatible Fock spaces
    GQCP::ProductFockSpace product_fock_space (K, 3, 3);
    GQCP::SelectedFockSpace fock_space (product_fock_space, 1);

    // The SelectedFockSpace includes the same configurations as the ProductFockSpace
    // These builder instances should return the same results.
    GQCP::SelectedCI random_sci (fock_space);
    GQCP::FFCI random_fci (product_fock_space, 1);

    Eigen::VectorXd sx = random_sci.calculateDiagonal(random_hamiltonian_parameters);
    Eigen::VectorXd fx = random_fci.calculateDiagonal(random_hamiltonian_parameters);



    Eigen::VectorXd s_mv = random_sci.matrixVectorProduct(random_hamiltonian_parameters, sx, sx);
    Eigen::VectorXd f_mv = random_fci.matrixVectorProduct(random_hamiltonian_parameters, fx, fx);

    Eigen::MatrixXd s_ham = random_sci.constructHamiltonian(random_hamiltonian_parameters);
    Eigen::MatrixXd f_ham = random_fci.constructHamiltonian(random_hamiltonian_parameters);

    std::cout<<std::endl<<s_mv<<std::endl<<"----------------------"<<std::endl;
    std::cout<<std::endl<<f_mv<<std::endl<<"- - - - - - - - - - - - - - - - - - - - - -"<<std::endl;

    std::cout<<std::endl<<s_ham<<std::endl<<"----------------------"<<std::endl;
    std::cout<<std::endl<<f_ham<<std::endl;

    BOOST_CHECK(sx.isApprox(fx));
    BOOST_CHECK(s_mv.isApprox(f_mv));
    BOOST_CHECK(s_ham.isApprox(f_ham));
}
