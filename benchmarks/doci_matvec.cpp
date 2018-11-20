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

    // Random ham_par
    GQCP::HamiltonianParameters ham_par = GQCP::constructRandomHamiltonianParameters(K);
    Eigen::VectorXd diagonal = doci.calculateDiagonal(ham_par);
    Eigen::VectorXd random = Eigen::VectorXd::Random(diagonal.rows());

    // Code inside this loop is measured repeatedly
    for (auto _ : state) {
        Eigen::VectorXd matvec = doci.matrixVectorProduct(ham_par, random, diagonal);
        // Make sure the variable is not optimized away by compiler
        benchmark::DoNotOptimize(matvec);
    }

    state.counters["Orbitals"] = K;
    state.counters["Electron pairs"] = N;
}


static void CustomArguments(benchmark::internal::Benchmark* b) {
    b->Args({10, 7});
    b->Args({11, 7});
    b->Args({12, 7});
    b->Args({28, 7});
    b->Args({28, 8});
    b->Args({28, 9});
    b->Args({28, 10});
    b->Args({28, 11});
    b->Args({28, 12});
    b->Args({28, 13});
    b->Args({28, 14});
}

// Perform the benchmarks
BENCHMARK(matvec)->Unit(benchmark::kMillisecond)->Apply(CustomArguments);


BENCHMARK_MAIN();
// --benchmark_counters_tabular=true --benchmark_out=<docimatvec_benchmark>
// --benchmark_out=<docimatvec_benchmark> --benchmark_counters_tabular=true
