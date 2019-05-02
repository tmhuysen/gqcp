#include "Grid/EvaluatedScalarFunctionOnCubicGrid.hpp"



namespace GQCP {



/*
 *  CONSTRUCTORS
 */

/**
 *  Evaluate a scalar function on a cubic grid
 *
 *  @param grid         the cubic grid
 *  @param function     the scalar function
 */
EvaluatedScalarFunctionOnCubicGrid::EvaluatedScalarFunctionOnCubicGrid(const CubicGrid& grid, const ScalarFunction<double, double, 3>& function) :
    grid (grid)
{

    const auto& steps = this->grid.get_steps();
    this->values = Tensor<double, 3>(steps[0], steps[1], steps[2]);
    this->values.setZero();

    for (size_t i = 0; i < steps[0]; i++) {
        for (size_t j = 0; j < steps[1]; j++) {
            for (size_t k = 0; k < steps[2]; k++) {
                const auto& r = grid.position(i, j, k);
                this->values(i,j,k) = function(r);
            }
        }
    }
}


/*
 *  PUBLIC METHODS
 */

/**
 *  Write the evaluated scalar values to a cube file
 *
 *  @param filename     the name of the cubefile that has to be generated
 *  @param molecule     the molecule that should be placed in the cubefile
 */
void EvaluatedScalarFunctionOnCubicGrid::toCubeFile(const std::string& filename, const Molecule& molecule) const {

    std::ofstream cubefile;
    cubefile.open(filename, std::fstream::out);

    const auto& steps = this->grid.get_steps();
    const auto& origin = this->grid.get_origin();
    const auto& step_sizes = this->grid.get_step_sizes();
    cubefile << "Fun times" << std::endl;
    cubefile << "Fun times" << std::endl;
    const auto& atoms = molecule.get_atoms();
    cubefile << std::scientific;
    cubefile << atoms.size() << " " << origin(0) << " " << origin(1) << " " << origin(2) << std::endl;
    cubefile << steps[0] << " " << step_sizes[0] << " " << 0.0 << " " << 0.0 << std::endl;
    cubefile << steps[1] << " " << 0.0 << " " << step_sizes[1] << " " << 0.0 << std::endl;
    cubefile << steps[2] << " " << 0.0 << " " << 0.0 << " " << step_sizes[0] << std::endl;
    for (const auto& atom : atoms) {
        cubefile << atom.atomic_number << " " << 0.0 << " " << atom.position(0) << " " << atom.position(1) << " " << atom.position(2) << std::endl;
    }

    for (size_t i = 0; i < steps[0]; i++) {
        for (size_t j = 0; j < steps[1]; j++) {
            for (size_t k = 0; k < steps[2]; k++) {
                cubefile << this->values(i,j,k) << " ";
                if (k % 6 == 5) { cubefile << std::endl; }  // there can only be 5 values on one line (hence % 6)
            }
        cubefile << std::endl;
        }
    }
}



}  // namespace GQCP
