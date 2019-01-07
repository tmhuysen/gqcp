/**
 *  A benchmark executable for the construction of the FCI Hamiltonian
 */

#include <benchmark/benchmark.h>

#include "HamiltonianParameters/HamiltonianParameters_constructors.hpp"
#include "HamiltonianBuilder/FCI.hpp"


static void initializeFCI(benchmark::State& state) {
    // Prepare parameters
    size_t K = state.range(0);
    size_t N = state.range(1);
    GQCP::ProductFockSpace fock_space (K, N, N);
    // Code inside this loop is measured repeatedly
    for (auto _ : state) {
        GQCP::FCI fci (fock_space);
        //benchmark::DoNotOptimize(fci);  // make sure the variable is not optimized away by compiler
    }

    state.counters["Orbitals"] = K;
    state.counters["Electron pairs"] = N;
    state.counters["Dimension"] = fock_space.get_dimension();
}


static void CustomArguments(benchmark::internal::Benchmark* b) {

    b->Args({12,6});
}


// Perform the benchmarks
BENCHMARK(initializeFCI)->Unit(benchmark::kMillisecond)->Apply(CustomArguments);
BENCHMARK_MAIN();