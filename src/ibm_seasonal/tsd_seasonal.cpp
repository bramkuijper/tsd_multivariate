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
    temperature = par.amplitude * std::sin(time_step * par.frequency);
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

    global_productivity = 0;

    // go through all survivors and assess whether they are breeding
    for (unsigned int patch_idx = 0;
            patch_idx < metapopulation.size();
            ++patch_idx)
    {
        available_local_males.clear();

        // clear male and female juvs
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

        // now go through all females and make offspring if 
        // season is alright
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

                    // calculate individual SR
                    double p_female = Kid.determine_sex(environment);

                    // realise sex determination
                    Kid.is_female = uniform(rng_r) < p_female;

                    // survive and add to pool of juveniles
                    if (Kid.is_female)
                    {
                        if (uniform(rng_r) < 
                                calculate_survival(female))
                        {
                            juvenile_females.push_back(Kid);
                        }
                    }
                    else
                    {
                        if (uniform(rng_r) < 
                                calculate_survival(male))
                        {
                            juvenile_males.push_back(Kid);
                        }
                    }
                } // end for()
            } // end if
        } // end for female_idx

        // shuffle the vectors
        std::shuffle(
                metapopulation[patch_idx].juveniles.begin(),
                    metapopulation[patch_idx].juveniles.end(),
                    rng_r
                    );

        // add this to calculation of global productivity
        // which is needed to calculate dispersal
        global_productivity += 
            metapopulation[patch_idx].juveniles.size();

    } // end for patch_idx
} // end reproduce()

// calculate survival
bool TSDSeasonal::calculate_survival(Sex const the_sex)
{
    return(std::exp(-0.5 * (par.t_opt[the_sex] - temperature) * 
                (par.t_opt[the_sex] - temperature) / par.omega_t[the_sex]));
} // end calculate_survival


// replace vacancies with newborn juveniles
void TSDSeasonal::replace()
{
    int n_vacancies;

    std::vector<double> probabilities{};

    // go through all survivors and assess whether they are breeding
    for (unsigned int patch_idx = 0;
            patch_idx < metapopulation.size();
            ++patch_idx)
    {
        // calculate number of vacancies
        n_vacancies = par.n - etapopulation[patch_idx].
            juveniles_f.size() - metapopulation[patch_idx].juveniles_m.size();

        // event 1: sample local male
        probabilities.push_back(
                (1.0 - d[male]) * metapopulation[patch_idx].
                juvenile_males.size());

        // event 2: sample dispersing male
        probabilities.push_back(
                d[male] * productivity[male] / par.npatches);

        // event 3: sample philopatric female
        probabilities.push_back((1.0 - d[female]) * metapopulation[patch_idx].
                juvenile_females.size());

        // event 4: sample dispersing female
        probabilities.push_back(d[female] * productivity[female] / par.npatches);

        std::discrete_distribution(probabilities.begin(), probabilities.end());

        for (int replace_idx = 0; 
                replace_idx < n_vacancies;
                ++replace_idx)
        {

            if (metapopulation[patch_idx].juveniles.size() == 0)
            {
                break;
            }

            if (metapopulation[patch_idx].juveniles.back().is_female)
            {
                metapopulation[patch_idx].female_survivors.push_back(
                        metapopulation[patch_idx].juveniles.back());
            }
            else
            {
                metapopulation[patch_idx].male_survivors.push_back(
                        metapopulation[patch_idx].juveniles.back()
                        );
            }

            metapopulation[patch_idx].juveniles.pop_back();
        } // end replace_idx
    } // end for patch_idx
} // end replace()

void TSDSeasonal::write_parameters()
{
    data_file << std::endl << std::endl
        << "seed;" << seed << std::endl
        << "npatches;" << par.npatches << std::endl
        << "n;" << par.n << std::endl
        << "df;" << par.d[female] << std::endl
        << "dm;" << par.d[male] << std::endl
        << "sf;" << par.survival_prob[female] << std::endl
        << "sm;" << par.survival_prob[male] << std::endl
        << "toptf;" << par.t_opt[female] << std::endl
        << "toptm;" << par.t_opt[male] << std::endl
        << "omegaf;" << par.omega_t[female] << std::endl
        << "omegam;" << par.omega_t[male] << std::endl
        << "frequency;" << par.frequency << std::endl
        << "amplitude;" << par.temperature_amplitude << std::endl
        << "init_t;" << par.init_t << std::endl
        << "init_resources;" << par.init_resources << std::endl
        << "init_effort_per_timestep;" << par.init_effort_per_timestep << std::endl
        << "init_a;" << par.init_a << std::endl
        << "init_b;" << par.init_b << std::endl
        << "mu_a;" << par.mu_a << std::endl
        << "mu_b;" << par.mu_b << std::endl
        << "mu_effort_per_timestep;" << par.mu_effort_per_timestep << std::endl
        << "mu_t;" << par.mu_t << std::endl
        << "unif_range_sdmu_t;" << par.unif_range_sdmu_t << std::endl
        << "sdmu;" << par.sdmu << std::endl; 
} // end write_parameters
