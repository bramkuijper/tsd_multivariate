#include "patch.hpp"
#include "individual.hpp"

Patch::Patch(Parameters const &par) :
    females(par.n[female], Individual(par, true)),
    males(par.n[male], Individual(par, false))
{}
