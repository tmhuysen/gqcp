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
#include "FockSpace/BaseFockSpace.hpp"

#include "FockSpace/FockSpace.hpp"
#include "FockSpace/ProductFockSpace.hpp"
#include "FockSpace/SelectedFockSpace.hpp"


namespace GQCP {


/*
 * PROTECTED CONSTRUCTORS
 */

/**
 *  @param K        the number of orbitals
 *  @param dim      the dimension of the Fock space
 */
BaseFockSpace::BaseFockSpace(size_t K, size_t dim) :
    K (K),
    dim (dim)
{}



/*
 *  NAMED CONSTRUCTORS
 */

/**
 *  Clones a derived BaseFockSpace instance to the heap memory
 *
 *  @param fock_space     reference to a derived BaseFockSpace instance to be cloned.
 *
 *  @return a shared pointer owning the heap-cloned Fock space
 */
std::shared_ptr<BaseFockSpace> BaseFockSpace::CloneToHeap(const BaseFockSpace& fock_space) {

    std::shared_ptr<BaseFockSpace> fock_space_ptr;

    switch (fock_space.get_type()) {

        case FockSpaceType::FockSpace: {
            fock_space_ptr = std::make_shared<FockSpace>(FockSpace(dynamic_cast<const FockSpace&>(fock_space)));
            break;
        }

        case FockSpaceType::ProductFockSpace: {
            fock_space_ptr = std::make_shared<ProductFockSpace>(ProductFockSpace(dynamic_cast<const ProductFockSpace&>(fock_space)));
            break;
        }

        case FockSpaceType::SelectedFockSpace: {
            fock_space_ptr = std::make_shared<SelectedFockSpace>(SelectedFockSpace(dynamic_cast<const SelectedFockSpace&>(fock_space)));
            break;
        }

        case FockSpaceType::FrozenFockSpace: {
            fock_space_ptr = std::make_shared<FrozenFockSpace>(FrozenFockSpace(dynamic_cast<const FrozenFockSpace&>(fock_space)));
            break;
        }

        case FockSpaceType::FrozenProductFockSpace: {
            fock_space_ptr = std::make_shared<FrozenProductFockSpace>(FrozenProductFockSpace(dynamic_cast<const FrozenProductFockSpace&>(fock_space)));
            break;
        }
    }

    return fock_space_ptr;
}



/*
 *  PUBLIC
 */

/**
 *  @return the coefficient vector for the Hartree-Fock wave function (i.e. the 'first' ONV/Slater determinant)
 */
VectorX<double> BaseFockSpace::HartreeFockExpansion() const {
    VectorX<double> expansion = VectorX<double>::Zero(this->dim);
    expansion(0) = 1;  // first configuration is position 0 (conventional ordering of the Fock space)
    return expansion;
}


/**
 *  @return a random normalized coefficient vector, with coefficients uniformly distributed in [-1, 1]
 */
VectorX<double> BaseFockSpace::randomExpansion() const {
    VectorX<double> random = VectorX<double>::Random(this->dim);
    random.normalize();
    return random;
}


/**
 *  @return a constant normalized coefficients vector (i.e. all the coefficients are equal)
 */
VectorX<double> BaseFockSpace::constantExpansion() const {
    VectorX<double> constant = VectorX<double>::Ones(this->dim);
    constant.normalize();
    return constant;
}


}  // namespace GQCP
