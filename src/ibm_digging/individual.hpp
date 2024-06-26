#ifndef _INDIVIDUAL_HPP_
#define _INDIVIDUAL_HPP_

#include <random>
#include "parameters.hpp"

class Individual
{

    public:
        bool is_female{true}; // whether individual is female or not
        double a{0.0}; // the sex allocation temperature intercept
        double b{0.0}; // the sex allocation temperature slope
        int t{1}; // the time step in the season at which the individual reproduces
        double depth{0.0}; // the sex allocation temperature slope
       
        Individual(Parameters const &pars, bool const is_female);
        Individual(Individual const &other);

        Individual(Individual const &mom, 
                Individual const &dad,
                Parameters const &par,
                std::mt19937 &rng_r);

        void operator=(Individual const &other);

        double determine_sex(double const temperature);
};


#endif 
