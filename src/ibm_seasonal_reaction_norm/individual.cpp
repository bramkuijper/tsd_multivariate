#include <algorithm>
#include <random>
#include <cmath>
#include "individual.hpp"
#include "parameters.hpp"

Individual::Individual(Parameters const &params, bool const is_female) :
    is_female{is_female},
    a{params.init_a}, // pivotal temp 
    b{params.init_b}, // temp slope (sigmoidal)
    t{params.init_t}, // timing intercept
    tb{params.init_tb} // timing slope on perceived temp
{}

// copy constructor
Individual::Individual(Individual const &other) :
    is_female{other.is_female},
    a{other.a},
    b{other.b},
    t{other.t},
    tb{other.tb},
    attempted_to_mate{other.attempted_to_mate}
{
}

// birth constructor
Individual::Individual(Individual const &mom,
        Individual const &dad,
        Parameters const &par,
        std::mt19937 &rng_r) :
    is_female{true} // by default individual is female until sex is determined
{
    std::uniform_real_distribution uniform{0.0,1.0};
    std::normal_distribution normal{0.0, par.sdmu};

    // inherit haploid traits
    a = uniform(rng_r) < 0.5 ? dad.a : mom.a; // pivotal temperature
    b = uniform(rng_r) < 0.5 ? dad.b : mom.b; // sex alloc slope on temperature
    t = uniform(rng_r) < 0.5 ? dad.t : mom.t; // timing of reproduction, intercept
    tb = uniform(rng_r) < 0.5 ? dad.tb : mom.tb; // timing of reproduction
                                              //
    // mutate the intercept
    if (uniform(rng_r) < par.mu_a)
    {
        a += normal(rng_r);

        a = std::clamp(a, par.ab_range[0], par.ab_range[1]);
    }

    // mutate the threshold above which individuals develop as male
    if (uniform(rng_r) < par.mu_b)
    {
        b += normal(rng_r);
        b = std::clamp(b, par.ab_range[0], par.ab_range[1]);
    }

    if (uniform(rng_r) < par.mu_t)
    {
        t+=normal(rng_r);
    }
    
    if (uniform(rng_r) < par.mu_tb)
    {
        tb += normal(rng_r);
    }
} // end birth constructor


// assignment operator
void Individual::operator=(Individual const &other)
{
    is_female = other.is_female;
    a = other.a;
    b = other.b;
    t = other.t;
    tb = other.tb;

    attempted_to_mate = other.attempted_to_mate;
}

// determine sex of this individual
double Individual::determine_sex(double const temperature)
{
    double p_female = 1.0 / (1.0 + std::exp(-b *(temperature - a)));
    
    return(p_female);
} // end determine_sex

