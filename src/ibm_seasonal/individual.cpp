#include <algorithm>
#include <random>
#include "individual.hpp"
#include "parameters.hpp"

Individual::Individual(double const init_z, int const init_t) :
    z{init_z},
    t{init_t}
{
} // end Individual()

// copy constructor
Individual::Individual(Individual const &other) :
    z{other.z},
    t{other.t}
{

}

// birth constructor
Individual::Individual(Individual const &mom,
        Individual const &dad,
        Parameters const &par,
        std::mt19937 &rng_r)
{
    std::uniform_real_distribution uniform{0.0,1.0};

    // inherit haploid traits
    z = uniform(rng_r) < 0.5 ? dad.z : mom.z;
    t = uniform(rng_r) < 0.5 ? dad.t : mom.t;

    // mutate the threshold above which individuals develop as male
    if (uniform(rng_r) < par.mu_z)
    {
        std::normal_distribution normal{0.0, par.sdmu_z};

        z += normal(rng_r);
    }

    double delta_t_double{0.0};
    double delta_t_int{0.0};

    if (uniform(rng_r) < par.mu_t)
    {
        // either increment or decrement with values from a certain range
        delta_t_double = uniform(rng_r) * par.unif_range_sdmu_t;
    
        delta_t_int = lround(delta_t_double);

        t += uniform(rng_r) < 0.5 ? -delta_t_int : delta_t_int;
    
        t = std::clamp(t,0,max_t);
    }
} // end birth constructor

// assignment operator
void Individual::operator=(Individual const &other)
{
    z = other.z;
    t = other.t;
}
