#ifndef _PATCH_HPP_
#define _PATCH_HPP_

#include <vector>
#include "individual.hpp"
#include "parameters.hpp"

class Patch
{

    public: 
        std::vector <Individual> females{};
        std::vector <Individual> female_survivors{};
        std::vector <Individual> males{};
        std::vector <Individual> male_survivors{};
        std::vector <Individual> juveniles{};
        Patch(Parameters const &par);


};

#endif
