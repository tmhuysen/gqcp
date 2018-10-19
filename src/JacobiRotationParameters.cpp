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
#include "JacobiRotationParameters.hpp"

#include <stdexcept>


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  Constructor based on a given @param p, @param q and a @param angle expressed in radians
 */
JacobiRotationParameters::JacobiRotationParameters(size_t p, size_t q, double angle) :
    p (p),
    q (q),
    angle (angle)
{
    // Check if p > q
    if (this->p <= this->q) {
        throw std::invalid_argument("Can't construct a JacobiRotationParameter with p < q.");
    }
}


/*
 *  OPERATORS
 */
/**
 *  Overloading of operator<< for GQCP::JacobiRotationParameters to be used with streams
 */
std::ostream& operator<<(std::ostream& os, const GQCP::JacobiRotationParameters& jacobi_rotation_parameters) {

    os << "p: " << jacobi_rotation_parameters.p << ", q: " << jacobi_rotation_parameters.q << ", angle: " << jacobi_rotation_parameters.angle;
    return os;
}



}  // namespace GQCP
