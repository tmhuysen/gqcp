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


#include "Mathematical/Representation/SquareMatrix.hpp"


namespace GQCP {


/**
 *  A class that represents a Hubbard hopping matrix
 */
class HoppingMatrix : public SquareMatrix<double> {
public:
    // CONSTRUCTORS

    /**
     *  @param H        the Hubbard hopping matrix
     */
    HoppingMatrix(const SquareMatrix<double>& H);

    /**
     *  Generate the Hubbard hopping matrix from an adjacency matrix and parameters U and t
     *
     *  @param A        the Hubbard adjacency matrix, specifying the connectivity of the Hubbard lattice
     *  @param t        the Hubbard parameter t. Note that a positive value for t means a negative neighbour hopping term
     *  @param U        the Hubbard parameter U
     */
    HoppingMatrix(const SquareMatrix<double>& A, double t, double U);


    // NAMED CONSTRUCTORS

    /**
     *  @param upper_triangle       the upper triangle (in column-major ordering) that specifies the Hubbard hopping matrix
     *
     *  @return the hopping matrix that corresponds to the given upper triangle
     */
    static HoppingMatrix FromUpperTriangle(const VectorX<double>& upper_triangle);

    /**
     *  @param K        the number of lattice sites
     *
     *  @return a random hopping matrix with elements distributed uniformly in [-1.0, 1.0]
     */
    static HoppingMatrix Random(size_t K);

    /**
     *  @param csline           a comma-separated line that contains the upper triangle (in column-major ordering) of the Hubbard hopping matrix
     * 
     *  @return the hopping matrix that corresponds to the given comma-separated line
     */
    static HoppingMatrix FromCSLine(const std::string& csline);


    // PUBLIC METHODS

    /**
     *  @return the number of lattice sites corresponding to the Hubbard hopping matrix
     */
    size_t numberOfLatticeSites() const { return this->get_dim(); }
};


}  // namespace GQCP
