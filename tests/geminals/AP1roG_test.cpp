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
#define BOOST_TEST_MODULE "AP1roG"

#include <boost/test/unit_test.hpp>

#include "Geminals/AP1roG.hpp"

#include "Geminals/AP1roGLagrangianOptimizer.hpp"
#include "RHF/PlainRHFSCFSolver.hpp"
#include "properties/expectation_values.hpp"



BOOST_AUTO_TEST_CASE ( energy_as_contraction ) {

    auto h2 = GQCP::Molecule::Readxyz("data/h2_olsens.xyz");
    auto ao_mol_ham_par = GQCP::HamiltonianParameters<double>::Molecular(h2, "6-31G**");

    GQCP::PlainRHFSCFSolver plain_scf_solver (ao_mol_ham_par, h2);
    plain_scf_solver.solve();
    auto rhf = plain_scf_solver.get_solution();

    auto mol_ham_par = GQCP::HamiltonianParameters<double>(ao_mol_ham_par, rhf.get_C());


    // Optimize the AP1roG PSE Lagrangian with the initial guess of the geminal coefficients being 0
    GQCP::AP1roGLagrangianOptimizer lagrangian_optimizer (h2, mol_ham_par);
    lagrangian_optimizer.solve();
    double electronic_energy = lagrangian_optimizer.get_electronic_energy();
    auto G = lagrangian_optimizer.get_geminal_coefficients();
    auto multipliers = lagrangian_optimizer.get_multipliers();


    // Calculate the 1- and 2-RDM and check the trace with the one- and two-electron integrals
    auto D = GQCP::calculate1RDM(G, multipliers);
    auto d = GQCP::calculate2RDM(G, multipliers);

    double electronic_energy_by_contraction = GQCP::calculateExpectationValue(mol_ham_par, D, d) - mol_ham_par.get_scalar();  // only the electronic energy
    BOOST_CHECK(std::abs(electronic_energy_by_contraction - electronic_energy) < 1.0e-09);
}
