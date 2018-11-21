/**
 *  A benchmark executable for the DOCI constructHamiltonian
 */

#include "benchmark/benchmark.h"

#include "HamiltonianParameters/HamiltonianParameters_constructors.hpp"
#include "HamiltonianBuilder/DOCI.hpp"


static void constructHamiltonian(benchmark::State& state) {
    // Prepare parameters
    size_t K = state.range(0);
    size_t N = state.range(1);
    GQCP::FockSpace fock_space (K, N);
    GQCP::DOCI doci (fock_space);
    GQCP::HamiltonianParameters ham_par = GQCP::constructRandomHamiltonianParameters(K);

    // Code inside this loop is measured repeatedly
    for (auto _ : state) {
        Eigen::MatrixXd hamiltonian = doci.constructHamiltonian(ham_par);
        // Make sure the variable is not optimized away by compiler
        benchmark::DoNotOptimize(hamiltonian);
    }

    state.counters["Orbitals"] = K;
    state.counters["Electron pairs"] = N;
    state.counters["Dimension"] = fock_space.get_dimension();
}


static void CustomArguments(benchmark::internal::Benchmark* b) {
    for (int i = 5; i < 9; ++i){
        // b-Args({Orbitals, Electrons})
        b->Args({16,i});
    }
}

// Perform the benchmarks
BENCHMARK(constructHamiltonian)->Unit(benchmark::kMillisecond)->Apply(CustomArguments);


BENCHMARK_MAIN();
