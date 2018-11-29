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
#ifndef GQCP_MOLECULE_HPP
#define GQCP_MOLECULE_HPP


#include <stdlib.h>
#include <string>
#include <vector>

#include <Eigen/Dense>

#include "Atom.hpp"



namespace GQCP {



/**
 *  A class that represents a collection of atoms with a number of electrons
 */
class Molecule {
private:
    const std::vector<GQCP::Atom> atoms;  // coordinates in bohr
    const size_t N;  // number of electrons

    /**
     *  @param xyz_filename     the .xyz-file that contains the molecular coordinates in Angstrom
     *
     *  @return a vector of Atoms that are in the given xyz-file
     */
    static std::vector<GQCP::Atom> parseXYZFile(const std::string& xyz_filename);


public:
    // CONSTRUCTORS
    /**
     *  A constructor that creates a neutral molecule
     *
     *  @param atoms     the atoms that make up the molecule, with coordinates in bohr
     *  @param molecular_charge     +1 -> cation (one electron less than the neutral molecule)
     *                               0 -> neutral molecule
     *                              -1 -> anion (one electron more than the neutral molecule)
     */
    Molecule(const std::vector<GQCP::Atom>& atoms, int molecular_charge);

    /**
     *  A constructor that creates a neutral molecule
     *
     *  @param atoms     the atoms that make up the molecule, with coordinates in bohr
     */
    explicit Molecule(const std::vector<GQCP::Atom>& atoms);

    /**
     *  A constructor from that creates a charged molecule
     *
     *  @param xyz_filename     the .xyz-file that contains the molecular coordinates in Angstrom
     *  @param molecular_charge     +1 -> cation (one electron less than the neutral molecule)
     *                               0 -> neutral molecule
     *                              -1 -> anion (one electron more than the neutral molecule)
     */
    Molecule(const std::string& xyz_filename, int molecular_charge);

    /**
     *  A constructor that creates a neutral molecule
     *
     *  @param xyz_filename     the .xyz-file that contains the molecular coordinates in Angstrom
     */
    explicit Molecule(const std::string& xyz_filename);


    // NAMED CONSTRUCTORS
    /**
     *  @param n            the number of H atoms
     *  @param spacing      the internuclear spacing in bohr
     *  @param charge       the total charge
     *
     *  @return a charged H-chain with equal internuclear spacing
     */
    static Molecule HChain(size_t n, double spacing, int charge);

    /**
     *  @param n            the number of H atoms
     *  @param spacing      the internuclear spacing in bohr
     *
     *  @return a neutral H-chain with equal internuclear spacing
     */
    static Molecule HChain(size_t n, double spacing);

    /**
     *  @param n        the number of H2-molecules
     *  @param a        the internuclear distance in bohr
     *  @param b        the intermolecular distance in bohr
     *  @param charge   the total charge
     *
     *  @return a charged H2-chain
     */
    static Molecule H2Chain(size_t n, double a, double b, int charge);

    /**
     *  @param n        the number of H2-molecules
     *  @param a        the internuclear distance in bohr
     *  @param b        the intermolecular distance in bohr
     *
     *  @return a neutral H2-chain
     */
    static Molecule H2Chain(size_t n, double a, double b);


    // OPERATORS
    /**
     *  @param other        the other molecule
     *
     *  @return if this molecule is equal to the other, within the default GQCP::Atom::tolerance_for_comparison for the coordinates of the atoms
     */
    bool operator==(const GQCP::Molecule& other) const;

    /**
     *  Overloading of operator<< for a GQCP::Molecule to be used with streams
     *
     *  @param os           the output stream which the molecule should be concatenated to
     *  @param molecule     the molecule that should be concatenated to the output stream
     *
     *  @return the updated output stream
     */
    friend std::ostream& operator<<(std::ostream& os, const GQCP::Molecule& molecule);


    // GETTERS
    size_t get_N() const { return this->N; }
    const std::vector<GQCP::Atom>& get_atoms() const { return this->atoms; }
    size_t numberOfAtoms() const { return this->atoms.size(); }


    // PUBLIC METHODS
    /**
     *  @param other        the other molecule
     *  @param tolerance    the tolerance for the coordinates of the atoms
     *
     *  @return if this is equal to the other, within the given tolerance
     */
    bool isEqualTo(const GQCP::Molecule& other, double tolerance=GQCP::Atom::tolerance_for_comparison) const;

    /**
     *  @return the sum of all the charges of the nuclei
     */
    size_t calculateTotalNucleicCharge() const;

    /**
     *  @param index1   the index of the first atom
     *  @param index2   the index of the second atom
     *
     *  @return the distance between the two atoms at index1 and index2
     */
    double calculateInternuclearDistance(size_t index1, size_t index2) const;

    /**
     *  @return the internuclear repulsion energy due to the nuclear framework
     */
    double calculateInternuclearRepulsionEnergy() const;

    /**
     *  @return the electrical dipole moment vector generated by the nuclear framework
     */
    Eigen::Vector3d calculateNuclearDipoleMoment() const;
};



}  // namespace GQCP


#endif  // GQCP_MOLECULE_HPP
