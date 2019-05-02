#include "CubicGrid.hpp"


namespace GQCP {



/*
 *  CONSTRUCTORS
 */






/*
 *  PUBLIC METHODS
 */

/**
 *  Evaluate the given function on every point of the grid and add the result
 *
 *  @param function     the scalar function
 */

void CubicGrid::evaluate(const ScalarFunction<double, double, 3>& function) {

    for (size_t i = 0; i < this->steps[0]; i++) {
        for (size_t j = 0; j < this->steps[1]; j++) {
            for (size_t k = 0; k < this->steps[2]; k++) {

                this->values(i,j,k) += function(r);
            }
        }
    }
}



}  // namespace GQCP
