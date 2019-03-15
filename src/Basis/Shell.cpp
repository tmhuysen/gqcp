//
//  Shell.cpp
//  gqcp
//
//  Created by Laurent Lemmens on 14/03/2019.
//  Copyright © 2019 Laurent Lemmens. All rights reserved.
//

#include "Basis/Shell.hpp"


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

Shell::Shell(const Atom& atom, const std::vector<double>& exponents, const std::vector<Contraction>& contractions) :
    atom (atom),
    exponents (exponents),
    contractions (contractions)
{
    for (const auto& contraction : this->contractions) {
        if (contraction.coefficients.size() != this->exponents.size()) {
            throw std::invalid_argument("Shell(Atom, std::vector<double>, std::vector<Contraction>): the exponents and contractions must match in size.");
        }
    }
}



/*
 *  PUBLIC METHODS
 */

size_t Shell::numberOfContractions() const {
    return this->contractions.size();
}


}  // namespace GQCP
