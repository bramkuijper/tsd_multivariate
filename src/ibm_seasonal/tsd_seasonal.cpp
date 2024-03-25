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
        update_environment();
    }
} // end constructor

void TSDSeasonal::update_environment()
{
    temperature = params.temperature_amplitude * sin(time_step * params.frequency)
}

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

        assert(patch_iter->females.size() == par.n[female];
        assert(patch_iter->males.size() == par.n[male];

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


// calculate fecundity for a given female
int TSDSeasonal::calculate_fecundity(Individual &mother)
{
    double resources_to_current_effort = mother.effort_per_timestep;

    // if too few resources left to make a full clutch, use a fraction of resources
    if (resources_to_current_effort > mother.resources)
    {
        resources_to_current_effort = mother.resources;

        mother.resources = 0.0;
    }
    else
    {
        mother.resources -= resources_to_current_effort;
        assert(mother.resources >= 0.0);
        assert(mother.resources < params.init_resources);
    }

    n_eggs = std::floor(resources_to_current_effort);

    // stochastic rounding, if effort_per_timestep = 5.3
    // then we have 5 eggs + 0.3 prob of having another
    if (uniform(rng_r) < resources_to_current_effort - n_eggs)
    {
        ++n_eggs;
    }

    return(n_eggs);
} // end calculate_fecundity

void TSDSeasonal::reproduce()
{
    // go through all survivors and assess whether they are breeding
    for (int patch_idx = 0;
            patch_idx < metapopulation.size();
            ++patch_idx)
    {
        for (int female_survivor_idx = 0; 
                metapopulation[patch_idx].female_survivors.size();
                ++female_survivor_idx)
        {
            assert(metapopulation[patch_idx].
                    female_survivors[female_survivor_idx].is_female);
            // ok 
            if (metapopulation[patch_idx].
                    female_survivors[female_survivor_idx].t % timestep == 0)
            {
                // fecundity
                n_eggs = calculate_fecundity(female_survivors[female_survivor_idx]);


            }
        } // end for female_idx
    } // end for patch_idx
} // end reproduce()

void TSDSeasonal::replace()
{
}
