/**
 *  A benchmark executable for the DOCI matvec
 */

#include "benchmark/benchmark.h"

#include "HamiltonianParameters/HamiltonianParameters_constructors.hpp"
#include "HamiltonianBuilder/DOCI.hpp"


static void matvec(benchmark::State& state) {
    // Prepare parameters
    size_t K = state.range(0);
    size_t N = state.range(1);
    GQCP::FockSpace fock_space (K, N);
    GQCP::DOCI doci (fock_space);

    GQCP::HamiltonianParameters ham_par = GQCP::constructRandomHamiltonianParameters(K);
    Eigen::VectorXd diagonal = doci.calculateDiagonal(ham_par);
    Eigen::VectorXd x = fock_space.randomExpansion();

    // Code inside this loop is measured repeatedly
    for (auto _ : state) {
        Eigen::VectorXd matvec = doci.matrixVectorProduct(ham_par, x, diagonal);
        // Make sure the variable is not optimized away by compiler
        benchmark::DoNotOptimize(matvec);
    }

    state.counters["Orbitals"] = K;
    state.counters["Electron pairs"] = N;
    state.counters["Dimension"] = fock_space.get_dimension();
}


static void CustomArguments(benchmark::internal::Benchmark* b) {
    for (int i = 5; i < 9; ++i){
        // b-Args({Orbitals, Electrons})
        b->Args({28,i});
    }
}

// Perform the benchmarks
BENCHMARK(matvec)->Unit(benchmark::kMillisecond)->Apply(CustomArguments);


BENCHMARK_MAIN();
