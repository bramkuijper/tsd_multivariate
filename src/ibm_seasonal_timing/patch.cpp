#include "patch.hpp"
#include "individual.hpp"

Patch::Patch(Parameters const &par) :
    females(par.n/2, Individual(par, true)),
    males(par.n/2, Individual(par, false))
{}
