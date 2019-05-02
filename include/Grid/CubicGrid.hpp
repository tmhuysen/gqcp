#ifndef CubicGrid_h
#define CubicGrid_h


#include "math/Matrix.hpp"
#include "math/Tensor.hpp"


namespace GQCP {


class CubicGrid {
private:
    Vector<double, 3> origin;  // origin of the grid
    std::array<size_t, 3> steps;  // steps in x, y, z-directions
    std::array<double, 3> step_sizes;  // step sizes in x, y, z-directions

public:
    // CONSTRUCTORS
    /**
     *  @param origin       the origin of the grid
     *  @param steps        the number of steps in the x, y, z-directions
     *  @param step_sizes   the step sizes in the x, y, z-directions
     */
    CubicGrid(Vector<double, 3> origin, std::array<size_t, 3> steps, std::array<double, 3> step_sizes);


    // GETTERS
    const Vector<double, 3>& get_origin() const { return this->origin; }
    const std::array<size_t, 3>& get_steps() const { return this->steps; }
    const std::array<double, 3>& get_step_sizes() const { return this->step_sizes; }


    // PUBLIC METHODS
    /**
     *  @param i        step number in the x-direction
     *  @param j        step number in the y-direction
     *  @param k        step number in the z-direction
     *
     *  @return the position vector associated to the given indices
     */
    Vector<double, 3> position(size_t i, size_t j, size_t k) const;
};



}  // namespace GQCP




#endif  /* CubicGrid_h */
