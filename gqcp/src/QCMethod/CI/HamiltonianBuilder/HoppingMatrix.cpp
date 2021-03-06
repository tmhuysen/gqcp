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
#include "QCMethod/CI/HamiltonianBuilder/HoppingMatrix.hpp"

#include "Utilities/miscellaneous.hpp"

#include <iostream>


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  @param H        the Hubbard hopping matrix
 */
HoppingMatrix::HoppingMatrix(const SquareMatrix<double>& H) :
    SquareMatrix<double>(H)
{
    if (!H.transpose().isApprox(H)) {
        throw std::invalid_argument("HoppingMatrix::HoppingMatrix(const SquareMatrix<double>&): The given hopping matrix must be symmetric.");
    }
}


/**
 *  Generate the Hubbard hopping matrix from an adjacency matrix and parameters U and t
 *
 *  @param A        the Hubbard adjacency matrix, specifying the connectivity of the Hubbard lattice
 *  @param t        the Hubbard parameter t. Note that a positive value for t means a negative neighbour hopping term
 *  @param U        the Hubbard parameter U
 */
HoppingMatrix::HoppingMatrix(const SquareMatrix<double>& A, double t, double U) :
    HoppingMatrix(U * SquareMatrix<double>::Identity(A.get_dim(), A.get_dim()) - t * A)
{}



/*
 *  NAMED CONSTRUCTORS
 */

/**
 *  @param upper_triangle       the upper triangle (in column-major ordering) that specifies the Hubbard hopping matrix
 *
 *  @return the hopping matrix that corresponds to the given upper triangle
 */
HoppingMatrix HoppingMatrix::FromUpperTriangle(const VectorX<double>& upper_triangle) {

    return SquareMatrix<double>::FullFromTriangle(upper_triangle);
}


/**
 *  @param K        the number of lattice sites
 *
 *  @return a random hopping matrix with elements distributed uniformly in [-1.0, 1.0]
 */
HoppingMatrix HoppingMatrix::Random(size_t K) {

    VectorX<double> v = VectorX<double>::Random(K*(K+1)/2);  // random free variables

    return HoppingMatrix::FromUpperTriangle(v);
}


/**
 *  @param csline           a comma-separated line that contains the upper triangle (in column-major ordering) of the Hubbard hopping matrix
 * 
 *  @return the hopping matrix that corresponds to the given comma-separated line
 */
HoppingMatrix HoppingMatrix::FromCSLine(const std::string& csline) {

    if (csline.empty()) {
        throw std::invalid_argument("HoppingMatrix::FromCSLine(const std::string&): Comma-separated line was empty!");
    }

    // Split the comma-separated line into a std::vector
    std::vector<std::string> splitted_line;
    boost::split(splitted_line, csline, boost::is_any_of(","));

    std::vector<double> triagonal_data;
    for (const auto& x : splitted_line) {
        triagonal_data.push_back(std::stod(x));  // immediately convert string to double
    }

    // Map the std::vector<double> into a VectorX<double> to be used into an other constructor
    GQCP::VectorX<double> triagonal = Eigen::Map<Eigen::VectorXd>(triagonal_data.data(), triagonal_data.size());

    return GQCP::HoppingMatrix::FromUpperTriangle(triagonal);
}


}  // namespace GQCP
