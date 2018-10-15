#include "RDM/TwoRDM.hpp"


namespace GQCG {


/*
 *  CONSTRUCTORS
 */

TwoRDM::TwoRDM(Eigen::Tensor<double, 4> two_rdm) :
    BaseRDM (two_rdm.dimensions()[0]),
    two_rdm (two_rdm)
{}

/**
 * Constructor where two_rdm = @param two_rdm_aaaa + @param two_rdm_bbbb + @param two_rdm_aabb  + @param two_rdm_bbaa
 */
TwoRDM::TwoRDM(Eigen::Tensor<double, 4> two_rdm_aaaa, Eigen::Tensor<double, 4> two_rdm_bbbb, Eigen::Tensor<double, 4> two_rdm_aabb, Eigen::Tensor<double, 4> two_rdm_bbaa) :
    BaseRDM (two_rdm_aaaa.dimensions()[0]),
    two_rdm (two_rdm_aaaa + two_rdm_bbbb + two_rdm_aabb + two_rdm_bbaa),
    two_rdm_aaaa (two_rdm_aaaa),
    two_rdm_bbbb (two_rdm_bbbb),
    two_rdm_aabb (two_rdm_aabb),
    two_rdm_bbaa (two_rdm_bbaa)
{}

/*
 *  PUBLIC METHODS
 */

/**
 *  @return the trace of this->two_rdm
 */
double TwoRDM::trace() {
    // TODO: when Eigen3 releases tensor.trace(), use it to implement the reduction

    auto K = static_cast<size_t>(this->two_rdm.dimension(1));

    double trace = 0.0;
    for (size_t p = 0; p < K; p++) {
        for (size_t q = 0; q < K; q++) {
            trace += this->two_rdm(p,p,q,q);
        }
    }

    return trace;
}

/**
 *  @return Eigen::MatrixXd D, the reduced-over 2-RDM : D(p,q) = this->two_rdm(p,q,r,r)
 */
Eigen::MatrixXd TwoRDM::reduce_2RDM() {
    // TODO: when Eigen3 releases tensor.trace(), use it to implement the reduction

    auto K = static_cast<size_t>(this->two_rdm.dimension(1));

    Eigen::MatrixXd D = Eigen::MatrixXd::Zero(K, K);
    for (size_t p = 0; p < K; p++) {
        for (size_t q = 0; q < K; q++) {
            for (size_t r = 0; r < K; r++) {
                D(p,q) += this->two_rdm(p,q,r,r);
            }
        }
    }

    return D;
}

}  // namespace GQCG

