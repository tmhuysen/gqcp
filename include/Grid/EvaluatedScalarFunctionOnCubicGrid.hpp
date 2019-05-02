#ifndef EvaluatedScalarFunctionOnCubicGrid_hpp
#define EvaluatedScalarFunctionOnCubicGrid_hpp



#include "math/ScalarFunction.hpp"
#include "Grid/CubicGrid.hpp"
#include "Molecule.hpp"




namespace GQCP {


/**
 *  A class that can represent 'property' densities on a cubic grid
 */
class EvaluatedScalarFunctionOnCubicGrid {
private:
    CubicGrid grid;
    Tensor<double, 3> values;


public:
    // CONSTRUCTORS
    /**
     *  Evaluate a scalar function on a cubic grid
     *
     *  @param grid         the cubic grid
     *  @param function     the scalar function
     */
    EvaluatedScalarFunctionOnCubicGrid(const CubicGrid& grid, const ScalarFunction<double, double, 3>& function);


    // PUBLIC METHODS
    /**
     *  Write the evaluated scalar values to a cube file
     *
     *  @param filename     the name of the cubefile that has to be generated
     *  @param molecule     the molecule that should be placed in the cubefile
     */
    void toCubeFile(const std::string& filename, const Molecule& molecule) const;
};


}  // namespace GQCP



#endif  /* EvaluatedScalarFunctionOnCubicGrid_hpp */
