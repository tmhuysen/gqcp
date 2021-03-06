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
#pragma once


#include "Mathematical/Representation/Matrix.hpp"

#include <cstddef>
#include <utility>


namespace GQCP {


/**
 *  An enum class for the implemented eigenproblem solver types
 */
enum class SolverType {
    DENSE,
    SPARSE,
    DAVIDSON
};



/**
 *  A base struct to specify eigenproblem solver options, whose derived structs can be used with the eigenproblem solvers
 */
struct BaseSolverOptions {
public:
    // MEMBERS
    size_t number_of_requested_eigenpairs = 1;


    // DESTRUCTOR
    virtual ~BaseSolverOptions() = default;


    // PURE VIRTUAL METHODS
    virtual SolverType get_solver_type() const = 0;
};



/**
 *  A struct to specify dense eigenproblem solver options
 */
struct DenseSolverOptions : public BaseSolverOptions {
public:
    // OVERRIDDEN METHODS
    SolverType get_solver_type () const override { return SolverType::DENSE; };
};



/**
 *  A struct to specify sparse eigenproblem solver options
 */
struct SparseSolverOptions : public BaseSolverOptions {
public:
    // OVERRIDDEN METHODS
    SolverType get_solver_type () const override { return SolverType::SPARSE; };
};



/**
 *  A struct to specify Davidson eigenproblem solver options
 */
struct DavidsonSolverOptions : public BaseSolverOptions {
public:
    // MEMBERS
    double convergence_threshold = 1.0e-08;  // the tolerance on the norm of the residual vector
    double correction_threshold = 1.0e-12;  // the threshold used in solving the (approximated) residue correction equation

    size_t maximum_subspace_dimension = 15;
    size_t collapsed_subspace_dimension = 2;
    size_t maximum_number_of_iterations = 128;

    MatrixX<double> X_0;  // MatrixX<double> of initial guesses, or VectorX<double> of initial guess


    // CONSTRUCTORS
    /**
     *  @param initial_guess        the initial guess(es) for the Davidson algorithm, specified as column(s) of the given vector/matrix
     */
    explicit DavidsonSolverOptions(const MatrixX<double>& initial_guess) :
        X_0 (initial_guess)
    {}


    // OVERRIDDEN METHODS
    SolverType get_solver_type () const override { return SolverType::DAVIDSON; };
};


}  // namespace GQCP
