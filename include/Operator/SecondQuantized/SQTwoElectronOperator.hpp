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

#include "Mathematical/ChemicalRankFourTensor.hpp"
#include "Operator/SecondQuantized/SQOneElectronOperator.hpp"
#include "OrbitalOptimization/JacobiRotationParameters.hpp"
#include "Utilities/miscellaneous.hpp"

#include <array>


namespace GQCP {


/**
 *  A class that represents a second-quantized two-electron operator: it holds the matrix representation of its parameters, which are (usually) integrals over first-quantized operators
 *
 *  @tparam _Scalar             the scalar type, i.e. the scalar representation of one of the parameters
 *  @tparam _Components         the number of components of the second-quantized operator
 */
template <typename _Scalar, size_t _Components>
class SQTwoElectronOperator {
public:

    using Scalar = _Scalar;
    static constexpr auto Components = _Components;


private:
    std::array<ChemicalRankFourTensor<Scalar>, Components> G;  // all the matrix representations of the parameters (integrals) of the different components of this second-quantized operator

public:

    /*
     *  CONSTRUCTORS
     */

    /**
     *  @param G            all the matrix representations of the parameters (integrals) of the different components of the second-quantized operator
     */
    SQTwoElectronOperator(const std::array<ChemicalRankFourTensor<Scalar>, Components>& G) : 
        G (G)
    {
        // Check if the given matrix representations have the same dimensions
        const auto dimension = this->F[0].dimension();

        for (size_t i = 1; i < Components; i++) {
            if (dimension != this->F[i].dimension()) {
                throw std::invalid_argument("SQOneElectronOperator(const std::array<ChemicalMatrix<Scalar>, Components>&): The given matrix representations did not have the same dimensions.");
            }
        }
    }


    /**
     *  Construct a two-electron operator with zero parameters
     * 
     *  @param dim          the dimension of the matrix representation of the parameters, i.e. the number of orbitals/sites
     */
    SQTwoElectronOperator(const size_t dim) {
        for (size_t i = 0; i < Components; i++) {
            this->G[i] = ChemicalRankFourTensor<Scalar>(dim);
        }
    }


    /*
     * OPERATORS
     */
    // operator+, operator* OtherScalar




    /*
     *  PUBLIC METHODS
     */

    /**
     *  @return the dimension of the matrix representation of the parameters, i.e. the number of orbitals/sites
     */
    size_t dimension() const {
        return this->G[0].dimension();
    }


    /**
     *  @return all the matrix representations of the parameters (integrals) of the different components of this second-quantized operator
     */
    const std::array<ChemicalRankFourTensor<Scalar>, Components> allParameters() const {
        return this->G;
    }


    /**
     *  @param i            the index of the component
     * 
     *  @return the matrix representation of the parameters (integrals) of one of the the different components of this second-quantized operator
     */
    const ChemicalRankFourTensor<Scalar>& parameters(const size_t i = 0) const {
        return this->G[i];
    }


    /**
     *  In-place rotate the matrix representation of the two-electron operator using a unitary Jacobi rotation matrix constructed from the Jacobi rotation parameters. Note that this function is only available for real (double) matrix representations
     *
     *  @param jacobi_rotation_parameters       the Jacobi rotation parameters (p, q, angle) that are used to specify a Jacobi rotation: we use the (cos, sin, -sin, cos) definition for the Jacobi rotation matrix. See transform() for how the transformation matrix between the two bases should be represented
     */
    template<typename Z = Scalar>
    enable_if_t<std::is_same<Z, double>::value> rotate(const JacobiRotationParameters& jacobi_rotation_parameters) {

        /**
         *  While waiting for an analogous Eigen::Tensor Jacobi module, we implement this rotation by constructing a Jacobi rotation matrix and then doing a rotation with it
         */

        auto dim = static_cast<size_t>(this->dimension(0));  // .dimension() returns a long
        auto J = SquareMatrix<double>::FromJacobi(jacobi_rotation_parameters, dim);  // this is sure to return a unitary matrix

        this->rotate(J);
    }


    /**
     *  @return the one-electron operator that is the difference between a two-electron operator (e_pqrs) and a product of one-electron operators (E_pq E_rs)
     */
    SQOneElectronOperator<Scalar, Components> effectiveOneElectronPartition() const {

        // Initialize a zero operator
        const auto K = this->dimension();  // number of orbitals
        SQOneElectronOperator<Scalar, Components> F (K);


        // Use a formula to set the parameters
        for (size_t i = 0; i < Components; i++) {
            for (size_t p = 0; p < K; p++) {
                for (size_t q = 0; q < K; q++) {
                    for (size_t r = 0; r < K; r++) {
                        F.parameters(i)(p,q) -= 0.5 * this->parameters(i)(p, r, r, q);
                    }
                }
            }
        }

        return F;
    }
};



/*
 *  Convenience aliases
 */
template <typename Scalar>
using ScalarSQTwoElectronOperator = SQTwoElectronOperator<Scalar, 1>;

template <typename Scalar>
using VectorSQTwoElectronOperator = SQTwoElectronOperator<Scalar, 3>;


}  // namespace GQCP
