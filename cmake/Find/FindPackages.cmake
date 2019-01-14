# In this CMake file, we will find all required packages


# Find the Boost package
find_package(Boost REQUIRED COMPONENTS program_options)

# Find Eigen3
find_package(Eigen3 3.3.5 REQUIRED)

# Find Libint
find_package(Libint2 REQUIRED)

# Find Spectra
find_package(Spectra REQUIRED)

# Find doxygen
if(BUILD_DOCS)
    find_package(Doxygen REQUIRED dot)
endif(BUILD_DOCS)

# Find MKL
if(USE_MKL)
    find_package(MKL)
endif(USE_MKL)

# Find google benchmarks
if(BUILD_BENCHMARKS)
    find_package(benchmark)
endif(BUILD_BENCHMARKS)