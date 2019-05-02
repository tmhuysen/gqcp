#include "EvaluatedScalarFunctionOnCubicGrid.hpp"



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
    // Evaluate the 
    auto steps = this->grid.get_steps();
    this->values = Tensor<double, 3>(steps[0], steps[1], steps[2]);
    this->values.setZero();

    for (size_t i = 0; i < steps[0]; i++) {
        for (size_t j = 0; j < steps[1]; j++) {
            for (size_t k = 0; k < steps[2]; k++) {
                auto r = grid.position(i, j, k);
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


    

}


}  // namespace GQCP
