#include <algorithm>
#include <random>
#include "individual.hpp"
#include "parameters.hpp"

Individual::Individual(Parameters const &params, is_female) :
    is_female{is_female},
    z{params.init_z},
    t{params.init_t},
    effort_per_timestep{params.init_effort_per_timestep},
    resources{params.init_resources}
{
} // end Individual()

// copy constructor
Individual::Individual(Individual const &other) :
    is_female{other.is_female},
    z{other.z},
    t{other.t},
    effort_per_timestep{params.init_effort_per_timestep},
    resources{other.resources}
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

    // inherit haploid traits
    z = uniform(rng_r) < 0.5 ? dad.z : mom.z; // sex alloc threshold
    t = uniform(rng_r) < 0.5 ? dad.t : mom.t; // timing of reproduction
                                              //
    effort_per_timestep = uniform(rng_r) < 0.5 ? 
        dad.effort_per_timestep 
        : 
        mom.effort_per_timestep;

    // mutate the threshold above which individuals develop as male
    if (uniform(rng_r) < par.mu_z)
    {
        std::normal_distribution normal{0.0, par.sdmu_z};

        z += normal(rng_r);
    }

    // mutate the timestep at which an individual will reproduce
    double delta_t_double{0.0};
    double delta_t_int{0.0};

    if (uniform(rng_r) < par.mu_t)
    {
        // either increment or decrement with values from a certain range
        delta_t_double = uniform(rng_r) * par.unif_range_sdmu_t;
    
        delta_t_int = lround(delta_t_double);

        t += uniform(rng_r) < 0.5 ? -delta_t_int : delta_t_int;
    
        t = std::clamp(t,0,par.max_threshold);
    }


    // mutate the amount of reproductive resources an individual will invest
    // per timestep
    if (uniform(rng_r) < par.mu_effort_per_timestep)
    {
        effort_per_timestep += uniform(rng_r) * par.unif_range_sdmu_effort_timestep;

        effort_per_timestep = std::clamp(effort_per_timestep, 0.0, par.init_resources);
    }

    resources = params.init_resources;
} // end birth constructor

// assignment operator
void Individual::operator=(Individual const &other)
{
    is_female = other.is_female;
    z = other.z;
    t = other.t;
    effort_per_timestep = other.effort_per_timestep;
    resources = other.resources;
}
