#ifndef GQCP_DOCIRDMBUILDER_HPP
#define GQCP_DOCIRDMBUILDER_HPP


#include "FockSpace/FockSpace.hpp"
#include "RDM/RDMBuilder.hpp"
#include "RDM/RDMs.hpp"


namespace GQCP {


/**
 *  DOCIRDMBuilder is a class for the calculation of a density matrix from a given wave function
 *  or coefficient expansion in a doubly occupied or single Fock space
 */
class DOCIRDMBuilder : public RDMBuilder {
    FockSpace fock_space;  // both the alpha and beta Fock space


public:
    // CONSTRUCTOR
    explicit DOCIRDMBuilder(const FockSpace& fock_space);


    // DESTRUCTOR
    ~DOCIRDMBuilder() = default;


    // OVERRIDDEN PUBLIC METHODS
    /**
     *  @return all 1-RDMs from a coefficient vector @param x
     */
    OneRDMs calculate1RDMs(const Eigen::VectorXd& x) override;

    /**
     *  @return all 2-RDMs from a coefficient vector @param x
     */
    TwoRDMs calculate2RDMs(const Eigen::VectorXd& x) override;

    /**
     *  @return the Fock space of the RDMBuilder
     */
    BaseFockSpace* get_fock_space() override { return &fock_space; }
};


}  // namespace GQCP


#endif  // GQCP_DOCIRDMBUILDER_HPP
