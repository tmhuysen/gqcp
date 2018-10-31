/**
 *  A benchmark executable for the DOCI matvec
 */

#include "benchmark/benchmark.h"

#include "HamiltonianParameters/HamiltonianParameters_constructors.hpp"
#include "HamiltonianBuilder/DOCI.hpp"





static void matvec(benchmark::State& state) {
    // Code inside this loop is measured repeatedly
    size_t K = state.range(0);
    size_t N = state.range(1);
    GQCP::FockSpace fock_space (K, N);
    GQCP::DOCI doci (fock_space);
    GQCP::HamiltonianParameters ham_par = GQCP::constructRandomHamiltonianParameters(K);
    Eigen::VectorXd diagonal = doci.calculateDiagonal(ham_par);
    Eigen::VectorXd random = Eigen::VectorXd::Random(diagonal.rows());
    for (auto _ : state) {
        Eigen::VectorXd matvec = doci.matrixVectorProduct(ham_par, random, diagonal);
        // Make sure the variable is not optimized away by compiler
        benchmark::DoNotOptimize(matvec);
    }
}

static void matvec2(benchmark::State& state) {
    // Code inside this loop is measured repeatedly
    size_t K = state.range(0);
    size_t N = state.range(1);
    GQCP::FockSpace fock_space (K, N);
    GQCP::DOCI doci (fock_space);
    GQCP::HamiltonianParameters ham_par = GQCP::constructRandomHamiltonianParameters(K);
    Eigen::VectorXd diagonal = doci.calculateDiagonal(ham_par);
    Eigen::VectorXd random = Eigen::VectorXd::Random(diagonal.rows());
    for (auto _ : state) {
        Eigen::VectorXd matvec = doci.matrixVectorProduct2(ham_par, random, diagonal);
        // Make sure the variable is not optimized away by compiler
        benchmark::DoNotOptimize(matvec);
    }
}

static void CustomArguments(benchmark::internal::Benchmark* b) {
    for (int i = 10; i <= 28; ++i){
        b->Args({i, i/2});
    }



}

BENCHMARK(matvec)->Apply(CustomArguments);
BENCHMARK(matvec2)->Apply(CustomArguments)->UseRealTime();


BENCHMARK_MAIN();