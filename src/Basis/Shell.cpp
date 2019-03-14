//
//  Shell.cpp
//  gqcp
//
//  Created by Laurent Lemmens on 14/03/2019.
//  Copyright © 2019 Laurent Lemmens. All rights reserved.
//

#include "Basis/Shell.hpp"


namespace GQCP {


Shell::Shell(const Atom& atom, const std::vector<double>& exponents, const std::vector<Contraction>& contractions) :
    atom (atom),
    exponents (exponents),
    contractions (contractions)
{
    if (exponents.size() != contractions.size()) {
        throw std::invalid_argument("Shell(Atom, std::vector<double>, std::vector<Contraction>): the exponents and contractions must match in size.");
    }
}



}  // namespace GQCP
