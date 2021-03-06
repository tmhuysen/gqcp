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
#include <pybind11/pybind11.h>

namespace py = pybind11;


/**
 *  As stated in the FAQ (https://pybind11.readthedocs.io/en/stable/faq.html#how-can-i-reduce-the-build-time), it is good practice to split the binding code over multiple files
 */
namespace gqcpy {

void bindQCMethodDOCINewtonOrbitalOptimizer(py::module& module);
void bindQCMethodDOCIRHF(py::module& module);
void bindQCMethodHubbard(py::module& module);
void bindQCMethodFCI(py::module& module);
void bindQCMethodFukuiDysonAnalysis(py::module& module);
void bindMullikenConstrainedFCI(py::module& module);
void bindVersion(py::module& module);

void bindMolecule(py::module& module);
void bindNucleus(py::module& module);
void bindSpinorBasis(py::module& module);
void bindSQOneElectronOperator(py::module& module);
void bindSQTwoElectronOperator(py::module& module);

}  // namespace gqcpy




/**
 *  The actual Python binding into the gqcpy Python module
 */
PYBIND11_MODULE (gqcpy, module) {

    gqcpy::bindVersion(module);

    // Bind basic functionality
    gqcpy::bindMolecule(module);
    gqcpy::bindNucleus(module);
    gqcpy::bindSpinorBasis(module);
    gqcpy::bindSQOneElectronOperator(module);
    gqcpy::bindSQTwoElectronOperator(module);


    // Bind quantum chemical methods
    gqcpy::bindQCMethodDOCINewtonOrbitalOptimizer(module);
    gqcpy::bindQCMethodDOCIRHF(module);
    gqcpy::bindQCMethodHubbard(module);
    gqcpy::bindQCMethodFCI(module);
    gqcpy::bindQCMethodFukuiDysonAnalysis(module);
    gqcpy::bindMullikenConstrainedFCI(module);
}
