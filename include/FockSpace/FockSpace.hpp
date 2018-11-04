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
#ifndef GQCP_FOCKSPACE_HPP
#define GQCP_FOCKSPACE_HPP


#include "FockSpace/BaseFockSpace.hpp"

#include <boost/numeric/conversion/converter.hpp>
#include <boost/math/special_functions.hpp>


namespace GQCP {


/**
 *  The full Fock space for a given set of orbitals and number of electrons
 *  where the ONVs and addresses are linked
 *  through a hashing function calculated with an addressing scheme.
 *  Implementation of the addressing scheme from :
 *      Molecular Electronic-Structure Theory (August 2000) by Trygve Helgaker, Poul Jorgensen, and Jeppe Olsen
 */
class FockSpace: public GQCP::BaseFockSpace {
private:
    const size_t N;  // number of electrons
    Matrixu vertex_weights;  // vertex_weights of the addressing scheme


    // PRIVATE METHODS
    /**
     *  @returns a permutation of the representation, giving the next bitstring permutation in reverse lexical ordering.
     *
     *      Examples:
     *          011 -> 101
     *          101 -> 110
     */
    size_t ulongNextPermutation(size_t representation);


public:
    // CONSTRUCTORS
    /**
     *  Constructor given a @param K (spatial orbitals), N (electrons)
     *  on which the dimensions of the Fock space are based
     */
    FockSpace(size_t K, size_t N);


    // DESTRUCTORS
    ~FockSpace() override = default;


    // GETTERS
    size_t get_vertex_weights(size_t p, size_t m) const { return this->vertex_weights[p][m]; }
    const Matrixu& get_vertex_weights() const { return this->vertex_weights; }
    size_t get_N() const { return this->N; }
    FockSpaceType get_type() const override { return FockSpaceType::FockSpace; }

    // STATIC PUBLIC METHODS
    /**
     *  Given a number of spatial orbitals @param K
     *  and a number of electrons  @param N,
     *  @return the dimension of the Fock space
     */
    static size_t calculateDimension(size_t K, size_t N);



/**
 *  @return the ONV with the corresponding address in the considered space
 */
    ONV get_ONV(size_t address) {
        size_t representation;
        if (this->N == 0) {
            representation = 0;
        }

        else {
            representation = 0;
            size_t m = this->N;  // counts the number of electrons in the spin string up to orbital p

            for (size_t p = this->K; p > 0; p--) {  // p is an orbital index
                size_t weight = get_vertex_weights(p-1, m);

                if (weight <= address) {  // the algorithm can move diagonally, so we found an occupied orbital
                    address -= weight;
                    representation |= ((1) << (p - 1));  // set the (p-1)th bit: see (https://stackoverflow.com/a/47990)

                    m--;  // since we found an occupied orbital, we have one electron less
                    if (m == 0) {
                        break;
                    }
                }
            }
        }
        return ONV(this->K, this->N, representation);
    }


/**
 *  sets @param ONV to the next ONV in the space
 *  performs the ulongNextPermutation() function
 *  and updates the corresponding occupation indices
 *  of the ONV occupation vector
 */
    void setNext(ONV& onv) {
        onv.set_representation(ulongNextPermutation(onv.get_unsigned_representation()));
    }


/**
 *  @return the Fock space address (i.e. the ordering number) of the @param onv in reverse lexical ordering, in the fock space.
 */
    size_t getAddress(const ONV& onv) {
        // An implementation of the formula in Helgaker, starting the addressing count from zero
        size_t address = 0;
        size_t electron_count = 0;  // counts the number of electrons in the spin string up to orbital p
        unsigned long unsigned_onv = onv.get_unsigned_representation();  // copy the unsigned_representation of the onv

        while(unsigned_onv != 0) {  // we will remove the least significant bit each loop, we are finished when no bits are left
            size_t p = __builtin_ctzl(unsigned_onv);  // p is the orbital index counter (starting from 1)
            electron_count++;  // each bit is an electron hence we add it up to the electron count
            address += get_vertex_weights(p , electron_count);
            unsigned_onv ^= unsigned_onv & -unsigned_onv;  // flip the least significant bit
        }
        return address;
    }

};


}  // namespace GQCP


#endif  // GQCP_FOCKSPACE_HPP
