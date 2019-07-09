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
// 
#ifndef GQCP_LIBCINTTWOELECTRONINTEGRALBUFFER_HPP
#define GQCP_LIBCINTTWOELECTRONINTEGRALBUFFER_HPP


#include "Basis/BaseTwoElectronIntegralBuffer.hpp"



namespace GQCP {


/**
 *  A buffer for storing libcint two-electron integrals
 * 
 *  @tparam _IntegralScalar         the scalar representation of an integral
 *  @tparam _N                      the number of components the operator has
 */
template <typename _IntegralScalar, size_t _N>
class LibcintTwoElectronIntegralBuffer : public BaseTwoElectronIntegralBuffer<_IntegralScalar, _N> {
public:
    using IntegralScalar = _IntegralScalar;  // the scalar representation of an integral
    static constexpr auto N = _N;  // the number of components the operator has
};



}  // namespace GQCP



#endif  // GQCP_LIBCINTTWOELECTRONINTEGRALBUFFER_HPP
