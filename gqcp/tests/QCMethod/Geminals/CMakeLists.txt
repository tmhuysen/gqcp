list(APPEND test_target_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/AP1roG_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/AP1roGGeminalCoefficients_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/AP1roGJacobiOrbitalOptimizer_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/AP1roGLagrangianNewtonOrbitalOptimizer_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/AP1roGPSESolver_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/APIGGeminalCoefficients_test.cpp
)

set(test_target_sources ${test_target_sources} PARENT_SCOPE)
