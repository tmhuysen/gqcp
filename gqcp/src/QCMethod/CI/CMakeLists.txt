target_sources(gqcp
    PRIVATE
        CISolver.cpp
        DOCINewtonOrbitalOptimizer.cpp
)

add_subdirectory(HamiltonianBuilder)
