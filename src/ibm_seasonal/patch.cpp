#include "patch.hpp"
#include "individual.hpp"

Patch::Patch(Parameters const &par) :
    females(par.n[female], Individual(true, 
            par.init_z, par.init_timing)),
    males(par.n[male], Individual(false,
            par.init_z, par.init_timing))

{}
