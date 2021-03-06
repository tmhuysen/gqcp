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


namespace GQCP {


/**
 *  A base class for the implementation of minimizers
 */
class BaseMinimizer {
protected:
    size_t maximum_number_of_iterations;
    double convergence_threshold;

    double is_solved = false;

    VectorX<double> x;  // initial guess, current guess or final solution to the problem


public:
    // CONSTRUCTORS

    /**
     *  @param x0                               an initial guess
     *  @param convergence_threshold            the threshold for convergence on the norm of the gradient
     *  @param maximum_number_of_iterations     the maximum number of iterations in the algorithm
     */
    BaseMinimizer(const VectorX<double>& x0, double convergence_threshold = 1.0e-08, size_t maximum_number_of_iterations = 128);


    // DESTRUCTOR
    virtual ~BaseMinimizer() = default;


    // GETTERS
    const VectorX<double>& get_solution() const;


    // PUBLIC PURE VIRTUAL METHODS

    /**
     *  Solve the problem associated to the numerical minimization method
     *
     *  If successful, it sets
     *      - is_solved to true
     *      - x to the found solution
     */
    virtual void solve() = 0;
};


}  // namespace GQCP
