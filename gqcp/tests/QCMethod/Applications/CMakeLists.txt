list(APPEND test_target_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/AtomicDecompositionParameters_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/CISolver_DOCI_constrained_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/constrained_RHF_test.cpp
)

set(test_target_sources ${test_target_sources} PARENT_SCOPE)
