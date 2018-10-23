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
#ifndef AP1roGPSESolver_hpp
#define AP1roGPSESolver_hpp


#include "HamiltonianParameters/HamiltonianParameters.hpp"
#include "Molecule.hpp"
#include "AP1roG/AP1roGGeminalCoefficients.hpp"
#include "AP1roG/AP1roG.hpp"

namespace GQCP {

/**
 *
 */
class AP1roGPSESolver {
private:
    const size_t K;  // the number of special orbitals
    const size_t N_P;  // the number of electron pairs
    const GQCP::AP1roGGeminalCoefficients initial_geminal_coefficients;
    
    GQCP::HamiltonianParameters ham_par;

    GQCP::AP1roG solution;


public:
    // CONSTRUCTORS
    /**
     *  Constructor based on a given number of electron pairs @param N_P, Hamiltonian parameters @param ham_par and an initial guess for the geminal coefficients @param G
     */
    AP1roGPSESolver(size_t N_P, const GQCP::HamiltonianParameters& ham_par, const GQCP::AP1roGGeminalCoefficients& G);

    /**
     *  Constructor based on a given number of electron pairs @param N_P and Hamiltonian parameters @param ham_par
     *
     *  The initial guess for the geminal coefficients is zero
     */
    AP1roGPSESolver(size_t N_P, const GQCP::HamiltonianParameters& ham_par);

    /**
     *  Constructor based on a given @param molecule, Hamiltonian parameters @param ham_par and an initial guess for the geminal coefficients @param G
     */
    AP1roGPSESolver(const GQCP::Molecule& molecule, const GQCP::HamiltonianParameters& ham_par, const GQCP::AP1roGGeminalCoefficients& G);

    /**
     *  Constructor based on a given @param molecule and Hamiltonian parameters @param ham_par
     *
     *  The initial guess for the geminal coefficients is zero
     */
    AP1roGPSESolver(const GQCP::Molecule& molecule, const GQCP::HamiltonianParameters& ham_par);


    // GETTERS
    GQCP::AP1roG get_solution() const { return this->solution; }


    // PUBLIC METHODS
    /**
     *  Calculate the Jacobian element with compound indices (i,a) and (k,c) at the given geminal coefficients @param G
     *
     *      i and k are subscripts, a and c are superscripts
     */
    double calculateJacobianElement(const AP1roGGeminalCoefficients& G, size_t i, size_t a, size_t k, size_t c) const;

    /**
     *  Calculate and return the Jacobian at the given geminal coefficients @param g
     */
    Eigen::MatrixXd calculateJacobian(const Eigen::VectorXd& g) const;

    /**
     *  Calculate the coordinate function at the given geminal coefficients @param G, with given indices.
     *
     *      i is the subscript and a is the superscript
     */
    double calculateCoordinateFunction(const GQCP::AP1roGGeminalCoefficients& G, size_t i, size_t a) const;

    /**
     *  Calculate the coordinate functions for the pSEs at the given geminal coefficients @param g. @returns a vector F in which every entry is one of the coordinate functions
     */
    Eigen::VectorXd calculateCoordinateFunctions(const Eigen::VectorXd& g) const;

    /**
     *  Set up and solve the projected Schrödinger equations for AP1roG
     */
    void solve();
};


}  // namespace GQCP



#endif /* AP1roGPSESolver_hpp */