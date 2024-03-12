#include "tsd_seasonal.hpp"
#include "parameters.hpp"

// constructor
TSDSeasonal::TSDSeasonal(Parameters const &par) :
    par{par},
    data_file{par.file_name},
    uniform{0.0,1.0},
    patch_sampler{0, static_cast<int>(par.npatches) - 1},
    female_sampler{0, static_cast<int>(par.n[female]) - 1},
    male_sampler{0, static_cast<int>(par.n[male]) - 1},
    rd{},
    seed{rd()},
    rng_r{seed},
    metapopulation(500,Patch(par))
{
    // now run the simulation
    for (time_step = 0; 
            time_step < par.max_time; ++time_step)
    {
        // mortality 
        survive();
        reproduce();
        replace();
    }
} // end constructor

void TSDSeasonal::survive()
{
    // dist to sample number of female survivors
    // this could be in constructor but putting it here now 
    // as things could potentially depend on local envt
    std::binomial_distribution<unsigned> 
        female_survivor_sampler(
                par.n[female], 
                par.survival_prob[female]);
    
    std::binomial_distribution<unsigned> 
        male_survivor_sampler(
                par.n[male], 
                par.survival_prob[male]);

    for (std::vector<Patch>::iterator patch_iter = metapopulation.begin();
            patch_iter != metapopulation.end();
            ++patch_iter)
    {
        patch_iter->female_survivors.clear();
        patch_iter->male_survivors.clear();

        // get number of female survivors
        int n_female_survivors = female_survivor_sampler(rng_r);
        int n_male_survivors = male_survivor_sampler(rng_r);

        // sample surviving f and put them into female survivors
        std::sample(patch_iter->females.begin(), 
                patch_iter->females.end(),
                std::back_inserter(patch_iter->female_survivors),
                n_female_survivors,
                rng_r);
        
        // sample surviving m and put them into male survivors
        std::sample(patch_iter->males.begin(), 
                patch_iter->males.end(),
                std::back_inserter(patch_iter->male_survivors),
                n_male_survivors,
                rng_r);
    }
} // end survive()

void TSDSeasonal::reproduce()
{
    for (std::vector<Patch>::iterator patch_iter = metapopulation.begin();
            patch_iter != metapopulation.end();
            ++patch_iter)
    {
    }
} // end reproduce()

void TSDSeasonal::replace()
{
}
