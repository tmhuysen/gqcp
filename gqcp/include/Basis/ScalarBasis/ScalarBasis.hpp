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


#include "Basis/Integrals/Interfaces/LibintInterfacer.hpp"
#include "Basis/Integrals/Interfaces/LibcintInterfacer.hpp"
#include "Basis/Integrals/IntegralCalculator.hpp"
#include "Basis/Integrals/IntegralEngine.hpp"
#include "Basis/ScalarBasis/GTOBasisSet.hpp"
#include "Basis/ScalarBasis/GTOShell.hpp"
#include "Basis/ScalarBasis/ShellSet.hpp"
#include "Mathematical/LinearCombination.hpp"
#include "Mathematical/Representation/QCMatrix.hpp"
#include "Mathematical/Representation/QCRankFourTensor.hpp"
#include "Molecule/Molecule.hpp"
#include "Molecule/NuclearFramework.hpp"
#include "Operator/FirstQuantized/Operator.hpp"

#include <type_traits>


namespace GQCP {


/**
 *  A class that represents a scalar basis: it represents a collection of scalar basis functions. It provides an interface to obtain basis functions and calculate integrals over the shell type
 *
 * @tparam _Shell       the type of shell that this scalar basis contains
 */
template <typename _Shell>
class ScalarBasis {
public:
    using Shell = _Shell;
    using BasisFunction = typename Shell::BasisFunction;


private:
    ShellSet<Shell> shell_set;  // a collection of shells that represents this scalar basis


public:

    /*
     *  CONSTRUCTORS
     */

    /**
     *  @param shell_set        a collection of shells that represents this scalar basis
     */
    ScalarBasis(const ShellSet<Shell>& shell_set): 
        shell_set (shell_set)
    {}


    /**
     *  Construct a scalar basis by placing shells corresponding to the basisset specification on every nucleus of the nuclear framework
     *
     *  @param nuclear_framework        the nuclear framework containing the nuclei on which the shells should be centered
     *  @param basisset_name            the name of the basisset, e.g. "STO-3G"
     *
     *  Note that the normalization factors of the spherical (or axis-aligned Cartesian) GTO primitives are embedded in the contraction coefficients of the underlying shells
     * 
     *  @note This constructor is only available for GTOShells (for the std::enable_if, see https://stackoverflow.com/a/17842695/7930415)
     */
    template <typename Z = GTOShell>
    ScalarBasis(const NuclearFramework& nuclear_framework, const std::string& basisset_name,
                typename std::enable_if<std::is_same<Z, GTOShell>::value>::type* = 0) :
        ScalarBasis(GTOBasisSet(basisset_name).generate(nuclear_framework))
    {
        this->shell_set.embedNormalizationFactorsOfPrimitives();
    }


    /**
     *  Construct a scalar basis by placing shells corresponding to the basisset specification on every nucleus of the molecule
     *
     *  @param molecule             the molecule containing the nuclei on which the shells should be centered
     *  @param basisset_name        the name of the basisset, e.g. "STO-3G"
     *
     *  Note that the normalization factors of the spherical (or axis-aligned Cartesian) GTO primitives are embedded in the contraction coefficients of the underlying shells
     * 
     *  @note This constructor is only available for GTOShells (for the std::enable_if, see https://stackoverflow.com/a/17842695/7930415)
     */
    template <typename Z = GTOShell>
    ScalarBasis(const Molecule& molecule, const std::string& basisset_name,
                typename std::enable_if<std::is_same<Z, GTOShell>::value>::type* = 0) :
        ScalarBasis(molecule.nuclearFramework(), basisset_name)
    {}



    /*
     *  PUBLIC METHODS
     */

    /**
     *  @return the underlying set of shells
     */
    const ShellSet<Shell>& shellSet() const { return this->shell_set; }

    /**
     *  @return the number of basis functions that 'are' in this scalar basis
     */
    size_t numberOfBasisFunctions() const { return this->shell_set.numberOfBasisFunctions(); }

    /**
     *  @return the basis functions that 'are' in this scalar basis
     */
    std::vector<LinearCombination<double, BasisFunction>> basisFunctions() const { return this->shell_set.basisFunctions(); }

    /**
     *  @param i            the index of the requested basis function
     * 
     *  @return the basis function with the given index that 'is' in this scalar basis
     */
    LinearCombination<double, BasisFunction> basisFunction(const size_t i) const { return this->basisFunctions()[i]; }


    /*
     *  PUBLIC METHODS - LIBINT2 INTEGRALS
     */

    /**
     *  @param fq_op            the first-quantized operator
     * 
     *  @return the matrix representation (integrals) of the given first-quantized operator in this scalar basis
     */
    QCMatrix<double> calculateLibintIntegrals(const OverlapOperator& fq_op) const {
        static_assert(std::is_same<Shell, GTOShell>::value, "Can only calculate Libint2 integrals over GTOShells");

        // Construct the libint engine
        const auto max_nprim = this->shell_set.maximumNumberOfPrimitives();
        const auto max_l = this->shell_set.maximumAngularMomentum();
        auto engine = IntegralEngine::Libint(fq_op, max_nprim, max_l);  // cannot be const because libint2::Engine::compute() is not a const method


        // Calculate the integrals using the engine
        const auto integrals = IntegralCalculator::calculate(engine, this->shell_set);
        return integrals[0];
    }


    /**
     *  @param fq_op            the first-quantized operator
     * 
     *  @return the matrix representation (integrals) of the given first-quantized operator in this scalar basis
     */
    QCMatrix<double> calculateLibintIntegrals(const KineticOperator& fq_op) const {
        static_assert(std::is_same<Shell, GTOShell>::value, "Can only calculate Libint2 integrals over GTOShells");

        // Construct the libint engine
        const auto max_nprim = this->shell_set.maximumNumberOfPrimitives();
        const auto max_l = this->shell_set.maximumAngularMomentum();
        auto engine = IntegralEngine::Libint(fq_op, max_nprim, max_l);  // cannot be const because libint2::Engine::compute() is not a const method


        // Calculate the integrals using the engine
        const auto integrals = IntegralCalculator::calculate(engine, this->shell_set);
        return integrals[0];
    }


    /**
     *  @param fq_op            the first-quantized operator
     * 
     *  @return the matrix representation (integrals) of the given first-quantized operator in this scalar basis
     */
    QCMatrix<double> calculateLibintIntegrals(const NuclearAttractionOperator& fq_op) const {
        static_assert(std::is_same<Shell, GTOShell>::value, "Can only calculate Libint2 integrals over GTOShells");

        // Construct the libint engine
        const auto max_nprim = this->shell_set.maximumNumberOfPrimitives();
        const auto max_l = this->shell_set.maximumAngularMomentum();
        auto engine = IntegralEngine::Libint(fq_op, max_nprim, max_l);  // cannot be const because libint2::Engine::compute() is not a const method


        // Calculate the integrals using the engine
        const auto integrals = IntegralCalculator::calculate(engine, this->shell_set);
        return integrals[0];
    }


    /**
     *  @param fq_op            the first-quantized operator
     * 
     *  @return the matrix representation (integrals) of the given first-quantized operator in this scalar basis
     */
    std::array<QCMatrix<double>, 3>  calculateLibintIntegrals(const ElectronicDipoleOperator& fq_op) const {
        static_assert(std::is_same<Shell, GTOShell>::value, "Can only calculate Libint2 integrals over GTOShells");

        // Construct the libint engine
        const auto max_nprim = this->shell_set.maximumNumberOfPrimitives();
        const auto max_l = this->shell_set.maximumAngularMomentum();
        auto engine = IntegralEngine::Libint(fq_op, max_nprim, max_l);  // cannot be const because libint2::Engine::compute() is not a const method


        // Calculate the integrals using the engine
        const auto integrals = IntegralCalculator::calculate(engine, this->shell_set);
        return {integrals[0], integrals[1], integrals[2]};
    }


    /**
     *  @param fq_op            the first-quantized operator
     * 
     *  @return the matrix representation (integrals) of the given first-quantized operator in this scalar basis
     */
    QCRankFourTensor<double> calculateLibintIntegrals(const CoulombRepulsionOperator& fq_op) const {
        static_assert(std::is_same<Shell, GTOShell>::value, "Can only calculate Libint2 integrals over GTOShells");

        // Construct the libint engine
        const auto max_nprim = this->shell_set.maximumNumberOfPrimitives();
        const auto max_l = this->shell_set.maximumAngularMomentum();
        auto engine = IntegralEngine::Libint(fq_op, max_nprim, max_l);  // cannot be const because libint2::Engine::compute() is not a const method


        // Calculate the integrals using the engine
        const auto integrals = IntegralCalculator::calculate(engine, this->shell_set);
        return integrals[0];
    }



    /*
     *  PUBLIC METHODS - LIBCINT INTEGRALS
     *  Note that the Libcint integrals should only be used for Cartesian ShellSets
     */

    /**
     *  Calculate the overlap integrals using Libcint: only use this for all-Cartesian ShellSets
     *
     *  @param fq_op            the first-quantized operator
     * 
     *  @return the matrix representation of the overlap operator in this AO basis, using the libcint integral engine
     */
    QCMatrix<double> calculateLibcintIntegrals(const OverlapOperator& fq_op) const {
        static_assert(std::is_same<Shell, GTOShell>::value, "Can only calculate Libint2 integrals over GTOShells");

        auto engine = IntegralEngine::Libcint(fq_op, this->shell_set);  // cannot be const: Libint2 has a non-const compute() method inside its interface
        const auto integrals = IntegralCalculator::calculate(engine, this->shell_set);
        return integrals[0];
    }


    /**
     *  Calculate the kinetic energy integrals using Libcint: only use this for all-Cartesian ShellSets
     *
     *  @param fq_op            the first-quantized operator
     * 
     *  @return the matrix representation of the kinetic energy operator in this AO basis, using the libcint integral engine
     */
    QCMatrix<double> calculateLibcintIntegrals(const KineticOperator& fq_op) const {
        static_assert(std::is_same<Shell, GTOShell>::value, "Can only calculate Libint2 integrals over GTOShells");

        auto engine = IntegralEngine::Libcint(fq_op, this->shell_set);  // cannot be const: Libint2 has a non-const compute() method inside its interface
        const auto integrals = IntegralCalculator::calculate(engine, this->shell_set);
        return integrals[0];
    }


    /**
     *  Calculate the nuclear attraction energy integrals using Libcint: only use this for all-Cartesian ShellSets
     *
     *  @param fq_op            the first-quantized operator
     * 
     *  @return the matrix representation of the nuclear attraction operator in this AO basis, using the libcint integral engine
     */
    QCMatrix<double> calculateLibcintIntegrals(const NuclearAttractionOperator& fq_op) const {
        static_assert(std::is_same<Shell, GTOShell>::value, "Can only calculate Libint2 integrals over GTOShells");

        auto engine = IntegralEngine::Libcint(fq_op, this->shell_set);  // cannot be const: Libint2 has a non-const compute() method inside its interface
        const auto integrals = IntegralCalculator::calculate(engine, this->shell_set);
        return integrals[0];
    }


    /**
     *  Calculate the electrical dipole integrals using Libcint: only use this for all-Cartesian ShellSets
     *
     *  @param fq_op            the first-quantized operator
     *
     *  @return the matrix representation of the Cartesian components of the electrical dipole operator in this AO basis, using the libcint integral engine
     */
    std::array<QCMatrix<double>, 3> calculateLibcintIntegrals(const ElectronicDipoleOperator& fq_op) const {
        static_assert(std::is_same<Shell, GTOShell>::value, "Can only calculate Libint2 integrals over GTOShells");

        auto engine = IntegralEngine::Libcint(fq_op, this->shell_set);  // cannot be const: Libint2 has a non-const compute() method inside its interface
        const auto integrals = IntegralCalculator::calculate(engine, this->shell_set);
        return {integrals[0], integrals[1], integrals[2]};
    }

    /**
     *  Calculate the Coulomb repulsion energy integrals using Libcint: only use this for all-Cartesian ShellSets
     *
     *  @param fq_op            the first-quantized operator
     * 
     *  @return the matrix representation of the Coulomb repulsion operator in this AO basis, using the libcint integral engine
     */
    QCRankFourTensor<double> calculateLibcintIntegrals(const CoulombRepulsionOperator& fq_op) const {
        static_assert(std::is_same<Shell, GTOShell>::value, "Can only calculate Libint2 integrals over GTOShells");

        auto engine = IntegralEngine::Libcint(fq_op, this->shell_set);  // cannot be const: Libint2 has a non-const compute() method inside its interface
        const auto integrals = IntegralCalculator::calculate(engine, this->shell_set);
        return integrals[0];
    }
};


}  // namespace GQCP
