#include <cassert>
#include <vector>
#include <cmath>
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
    temperature = par.temperature_amplitude * std::sin(time_step * par.frequency);
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

        assert(patch_iter->females.size() == par.n[female]);
        assert(patch_iter->males.size() == par.n[male]);

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
    assert(mother.resources <= par.init_resources);

    double resources_to_current_effort = mother.effort_per_timestep * mother.resources;

    mother.resources -= resources_to_current_effort;

    // may dip to n've coz rounding errors, hence make bound
    if (mother.resources < 0.0)
    {
        mother.resources = 0.0;

    }

    int n_eggs = std::floor(resources_to_current_effort);

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
    int n_eggs;

    // empty vector with available local males
    std::vector <unsigned int> available_local_males{};

    // go through all survivors and assess whether they are breeding
    for (unsigned int patch_idx = 0;
            patch_idx < metapopulation.size();
            ++patch_idx)
    {
        available_local_males.clear();
        juveniles.clear();
        
        for (unsigned int male_survivor_idx = 0; 
                metapopulation[patch_idx].male_survivors.size();
                ++male_survivor_idx)
        {
            if (metapopulation[patch_idx].
                    male_survivors[male_survivor_idx].t % time_step == 0)
            {
                available_local_males.push_back(male_survivor_idx);
            }
        }

        // no male available: no offspring produced.
        if (available_local_males.size() < 1)
        {
            continue;
        }

        // mix the list of available males
        std::shuffle(available_local_males.begin(),
                available_local_males.end(),
                rng_r);


        for (unsigned int female_survivor_idx = 0; 
                metapopulation[patch_idx].female_survivors.size();
                ++female_survivor_idx)
        {
            assert(metapopulation[patch_idx].
                    female_survivors[female_survivor_idx].is_female);
            // ok 
            if (metapopulation[patch_idx].
                    female_survivors[female_survivor_idx].t % time_step == 0)
            {
                // get fecundity
                n_eggs = calculate_fecundity(
                        metapopulation[patch_idx].
                            female_survivors[female_survivor_idx]);

                for (int egg_idx = 0; egg_idx < n_eggs; ++n_eggs)
                {
                    // obtain father_idx
                    unsigned int father_idx = available_local_males[
                        egg_idx % available_local_males.size()];

                    assert(father_idx < metapopulation[patch_idx].
                            male_survivors.size());

                    Individual Kid(metapopulation[patch_idx].
                            female_survivors[female_survivor_idx],
                                metapopulation[patch_idx].
                                female_survivors[female_survivor_idx],
                                par,
                                rng_r);

                    p_female = Kid.determine_sex(environment);

                    Kid.is_female = uniform(rng_r) < p_female;

                    if (uniform(rng_r) < 
                            calculate_survival(Kid.is_female ? female : male))
                    {
                        juveniles.push_back(Kid);
                    }
                } // end for()
            } // end if
        } // end for female_idx
    } // end for patch_idx
} // end reproduce()

// calculate survival
bool TSDSeasonal::calculate_survival(Sex const the_sex)
{
    return(std::exp(-0.5 * (par.t_opt[the_sex] - temperature) * 
                (par.t_opt[the_sex] - temperature) / par.omega_t[the_sex]));
} // end calculate_survival

void TSDSeasonal::replace()
{
    // go through all survivors and assess whether they are breeding
    for (unsigned int patch_idx = 0;
            patch_idx < metapopulation.size();
            ++patch_idx)
    {
        n_vacant_f = metapopulation[patch_idx].female_survivors();
        n_vacant_m = metapopulation[patch_idx].male_survivors();

        for (unsigned int juvenile_idx = 0; 
                metapopulation[patch_idx].juveniles.size();
                ++juvenile_idx)
        {
            if (metapopulation[patch_idx].juveniles[juvenile_idx].is_female
                    && )
            {
                metapopulation[patch_idx].female_survivors.append(
            }
        }
    }
} // end replace()
