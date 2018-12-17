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
#ifndef GQCP_COMMON_HPP
#define GQCP_COMMON_HPP


#include <cstdlib>
#include <vector>

#include <Eigen/Dense>


namespace GQCP {


using Vectoru = std::vector<size_t>;
using Matrixu = std::vector<Vectoru>;
using VectorXs = Eigen::Matrix<size_t, Eigen::Dynamic, 1>;

using VectorFunction = std::function<Eigen::VectorXd (const Eigen::VectorXd&)>;
using MatrixFunction = std::function<Eigen::MatrixXd (const Eigen::VectorXd&)>;

}  // namespace GQCP


#endif  // GQCP_COMMON_HPP
