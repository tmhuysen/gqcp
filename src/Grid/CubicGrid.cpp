#include "Grid/CubicGrid.hpp"


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  @param origin       the origin of the grid
 *  @param steps        the number of steps in the x, y, z-directions
 *  @param step_sizes   the step sizes in the x, y, z-directions
 */
CubicGrid::CubicGrid(Vector<double, 3> origin, std::array<size_t, 3> steps, std::array<double, 3> step_sizes) :
    origin (origin),
    steps (steps),
    step_sizes (step_sizes)
{}



/*
 *  PUBLIC METHODS
 */

/**
 *  @param i        step number in the x-direction
 *  @param j        step number in the y-direction
 *  @param k        step number in the z-direction
 *
 *  @return the position vector associated to the given indices
 */
Vector<double, 3> CubicGrid::position(size_t i, size_t j, size_t k) const {

    double x = this->origin(0) + i * this->step_sizes[0];
    double y = this->origin(1) + j * this->step_sizes[1];
    double z = this->origin(2) + k * this->step_sizes[2];

    return Vector<double, 3>(x, y, z);
}

}  // namespace GQCP 
