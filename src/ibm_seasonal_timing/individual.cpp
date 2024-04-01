#include <algorithm>
#include <random>
#include <cmath>
#include "individual.hpp"
#include "parameters.hpp"

Individual::Individual(Parameters const &params, bool const is_female) :
    is_female{is_female},
    a{params.init_a},
    b{params.init_b},
    t{params.init_t}
{
} // end Individual()

// copy constructor
Individual::Individual(Individual const &other) :
    is_female{other.is_female},
    a{other.a},
    b{other.b},
    t{other.t}
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
    a = uniform(rng_r) < 0.5 ? dad.a : mom.a; // sex alloc threshold
    b = uniform(rng_r) < 0.5 ? dad.b : mom.b; // sex alloc threshold
    t = uniform(rng_r) < 0.5 ? dad.t : mom.t; // timing of reproduction
                                              //
    // mutate the intercept
    if (uniform(rng_r) < par.mu_a)
    {
        a += normal(rng_r);
    }

    // mutate the threshold above which individuals develop as male
    if (uniform(rng_r) < par.mu_b)
    {
        b += normal(rng_r);
    }

    // mutate the timestep at which an individual will reproduce
    double delta_t_double{0.0};
    double delta_t_int{0.0};

    if (uniform(rng_r) < par.mu_t)
    {
        // either increment or decrement with values from a certain range
        delta_t_double = uniform(rng_r) * par.unif_range_sdmu_t;
    
        delta_t_int = std::lround(delta_t_double);

        t += uniform(rng_r) < 0.5 ? -delta_t_int : delta_t_int;
    
        t = std::clamp(t,0,par.max_t);
    }

} // end birth constructor

// assignment operator
void Individual::operator=(Individual const &other)
{
    is_female = other.is_female;
    a = other.a;
    b = other.b;
    t = other.t;
}

// determine sex of this individual
double Individual::determine_sex(double const temperature)
{
    double p_female = 1.0 / (1.0 + std::exp(-(a + b * temperature)));

    return(p_female);
} // end determine_sex

