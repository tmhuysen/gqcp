/**
 *  A benchmark executable for the FCI 1RDM
 */

#include <benchmark/benchmark.h>

#include "HamiltonianParameters/HamiltonianParameters_constructors.hpp"
#include "HamiltonianBuilder/FCI.hpp"
#include "RDM/FCIRDMBuilder.hpp"


static void rdm(benchmark::State& state) {
    // Prepare parameters
    size_t K = state.range(0);
    size_t N = state.range(1);
    GQCP::ProductFockSpace fock_space (K, N, N);
    GQCP::FCI fci (fock_space);

    Eigen::VectorXd x = fock_space.randomExpansion();
    GQCP::FCIRDMBuilder fcirdm (fock_space);
    // Code inside this loop is measured repeatedly
    for (auto _ : state) {
        GQCP::OneRDMs var = fcirdm.calculate1RDMs(x);

        benchmark::DoNotOptimize(var);  // make sure the variable is not optimized away by compiler
    }

    state.counters["Orbitals"] = K;
    state.counters["Electron pairs"] = N;
    state.counters["Dimension"] = fock_space.get_dimension();
}


static void CustomArguments(benchmark::internal::Benchmark* b) {
    for (int i = 2; i < 6; ++i) {  // need int instead of size_t
        b->Args({10, i});  // orbitals, electron pairs
    }
    //b->Args({14, 7});
}


// Perform the benchmarks
BENCHMARK(rdm)->Unit(benchmark::kMillisecond)->Apply(CustomArguments);
BENCHMARK_MAIN();
