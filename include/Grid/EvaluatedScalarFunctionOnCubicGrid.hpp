#ifndef CubicGrid_hpp
#define CubicGrid_hpp



#include "math/ScalarFunction.hpp"



namespace GQCP {


/**
 *  A class that can represent 'property' densities on a cubic grid
 */
class EvaluatedScalarFunctionOnCubicGrid {
private:
    Tensor<double, 3> values;


public:
    // CONSTRUCTORS
    EvaluatedScalarFunctionOnCubicGrid(


    /**
     *  Evaluate the given function on every point of the grid and add the result
     *
     *  @param function     the scalar function
     */
    void evaluate(const ScalarFunction<double, double, 3>& function);
};



}  // namespace GQCP


#endif /* CubicGrid_hpp */
