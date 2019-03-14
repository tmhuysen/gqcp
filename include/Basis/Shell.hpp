#ifndef Shell_hpp
#define Shell_hpp


#include "Atom.hpp"
#include "Basis/Contraction.hpp"


namespace GQCP {


class Shell {
private:
    Atom atom;
    std::vector<double> exponents;  // shared for every contraction
    std::vector<Contraction> contractions;


public:
    // CONSTRUCTORS
    Shell(const Atom& atom, const std::vector<double>& exponents, const std::vector<Contraction>& contractions);

    // GETTERS
    const Atom& get_atom() const { return this->atom; }

    const std::vector<double>& get_exponents() const {
        return this->exponents;
    }

    const std::vector<Contraction>& get_contractions() const { return this->contractions; }
};


}  // namespace GQCP


#endif  /* Shell_hpp */
