#ifndef _INDIVIDUAL_HPP_
#define _INDIVIDUAL_HPP_

#include <random>
#include "parameters.hpp"

class Individual
{

    public:
        bool is_female{true}; // whether individual is female or not
        double z{0.0}; // the sex allocation temperature threshold: above the threshold male, below female
        int t{1}; // the time step trait at which to reproduce. t = 0 reproduce never. t = 1 reproduce every time step, t = 2 reproduce every 2nd time step etc
       
        // the amount of reproductive effort invested per time step
        double effort_per_timestep{0.0};

        // reproductive resources available
        double resources{0.0};

        Individual(Parameters & const pars);
        Individual(Individual const &other);

        Individual(Individual const &mom, 
                Individual const &dad,
                Parameters const &par,
                std::mt19937 &rng_r);

        void operator=(Individual const &other);
};


#endif 
