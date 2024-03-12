#ifndef _INDIVIDUAL_HPP_
#define _INDIVIDUAL_HPP_

#include <random>
#include "parameters.hpp"

class Individual
{

    public:
        bool is_female{true};
        double z{0.0};
        int t{1};

        Individual(bool const is_female, double const z, int const t);
        Individual(Individual const &other);

        Individual(Individual const &mom, 
                Individual const &dad,
                Parameters const &par,
                std::mt19937 &rng_r);

        void operator=(Individual const &other);
};


#endif 
