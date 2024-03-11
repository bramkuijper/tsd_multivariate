#ifndef _PATCH_HPP_
#define _PATCH_HPP_

#include <vector>
#include "individual.hpp"

class Patch
{
    std::vector <Individual> females{};
    std::vector <Individual> males{};
    std::vector <Individual> juveniles{};

    public: 
        Patch();


};

#endif
