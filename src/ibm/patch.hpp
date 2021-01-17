#ifndef PATCH_HPP_
#define PATCH_HPP_

#include <vector>
#include "individual.hpp"


// an individual patch of the metapopulation

class Patch {

    public:
        std::vector<Individual> breedersF;
        std::vector<Individual> breedersM;

        std::vector<Individual> phil_juvsF;
        std::vector<Individual> phil_juvsM;

        Patch();

        Patch(Patch const &other);

        void operator=(Patch const &other);
}; // end class Patch


#endif
