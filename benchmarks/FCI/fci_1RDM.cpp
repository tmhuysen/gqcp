/**
 *  A benchmark executable for the FCI 1RDMs
 */

#include <benchmark/benchmark.h>

#include "HamiltonianParameters/HamiltonianParameters.hpp"
#include "RDM/FCIRDMBuilder.hpp"


static void matvec(benchmark::State& state) {
    // Prepare parameters
    size_t K = state.range(0);
    size_t N = state.range(1);
    GQCP::ProductFockSpace fock_space (K, N, N);
    GQCP::FCIRDMBuilder rdm_builder (fock_space);
    GQCP::VectorX<double> x = fock_space.randomExpansion();

    // Code inside this loop is measured repeatedly
    for (auto _ : state) {
        GQCP::OneRDMs<double> one_rdms = rdm_builder.calculate1RDMs(x);

        benchmark::DoNotOptimize(one_rdms);  // make sure the variable is not optimized away by compiler
    }

    state.counters["Orbitals"] = K;
    state.counters["Electron pairs"] = N;
    state.counters["Dimension"] = fock_space.get_dimension();
}


static void CustomArguments(benchmark::internal::Benchmark* b) {
    for (int i = 2; i < 7; ++i) {  // need int instead of size_t
        b->Args({12, i});  // orbitals, electron pairs
    }
}


// Perform the benchmarks
BENCHMARK(matvec)->Unit(benchmark::kMillisecond)->Apply(CustomArguments);
BENCHMARK_MAIN();
