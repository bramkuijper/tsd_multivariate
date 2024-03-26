#include <cassert>
#include <stdexcept>
#include <vector>
#include <cmath>
#include "tsd_seasonal.hpp"
#include "parameters.hpp"

// constructor
TSDSeasonal::TSDSeasonal(Parameters const &par) :
    par{par},
    data_file{par.file_name},
    uniform{0.0,1.0},
    rd{},
    seed{rd()},
    rng_r{seed},
    metapopulation(par.npatches,Patch(par))
{
    // now run the simulation
    for (time_step = 0; 
            time_step < par.max_simulation_time; ++time_step)
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
    for (std::vector<Patch>::iterator patch_iter = metapopulation.begin();
            patch_iter != metapopulation.end();
            ++patch_iter)
    {
        patch_iter->female_survivors.clear();
        patch_iter->male_survivors.clear();

        std::binomial_distribution<unsigned> 
            female_survivor_sampler(
                    patch_iter->females.size(),
                    par.survival_prob[female]);
        
        std::binomial_distribution<unsigned> 
            male_survivor_sampler(
                    patch_iter->males.size(),
                    par.survival_prob[male]);


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

        metapopulation[patch_idx].male_juveniles.clear();
        metapopulation[patch_idx].female_juveniles.clear();
        
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
                    double p_female = Kid.determine_sex(temperature);

                    // realise sex determination
                    Kid.is_female = uniform(rng_r) < p_female;

                    // survive and add to pool of juveniles
                    if (Kid.is_female)
                    {
                        if (uniform(rng_r) < 
                                calculate_survival(female))
                        {
                            metapopulation[patch_idx].female_juveniles.push_back(Kid);
                        }
                    }
                    else
                    {
                        if (uniform(rng_r) < 
                                calculate_survival(male))
                        {
                            metapopulation[patch_idx].male_juveniles.push_back(Kid);
                        }
                    }
                } // end for()
            } // end if
        } // end for female_idx

    } // end for patch_idx
} // end reproduce()


void TSDSeasonal::calculate_patch_productivities()
{
    // vectors to store all local productivities
    // of males and females, allowing us to make 
    // a productivity distribution
    std::vector<int> local_productivity_male{};
    std::vector<int> local_productivity_female{};
     
    // reset sum of all productivities to 0
    global_productivity[male] = 0;
    global_productivity[female] = 0;

    int njuv_males,njuv_females;

    // go through all survivors and assess whether they are breeding
    for (unsigned int patch_idx = 0;
            patch_idx < metapopulation.size();
            ++patch_idx)
    {
        njuv_males = metapopulation[patch_idx].male_juveniles.size();
        njuv_females = metapopulation[patch_idx].female_juveniles.size();

        local_productivity_male.push_back(njuv_males);
        local_productivity_female.push_back(njuv_females);
        
        // add this to calculation of global productivity
        // which is needed to calculate dispersal
        global_productivity[male] += njuv_males;

        global_productivity[female] += njuv_females;
    }// end for patch_idx


    // now update the distribution of patch productivities
    // which is used for remote sampling
    std::discrete_distribution<unsigned>::param_type 
        distribution_param_juvs_m(
                local_productivity_male.begin(),
                local_productivity_male.end());

    std::discrete_distribution<unsigned>::param_type 
        distribution_param_juvs_f(
                local_productivity_female.begin(),
                local_productivity_female.end()
                );

    productivity_distribution_f.param(distribution_param_juvs_f);
    productivity_distribution_m.param(distribution_param_juvs_m);


} // end calculate_patch_productivities


// calculate survival
bool TSDSeasonal::calculate_survival(Sex const the_sex)
{
    return(std::exp(-0.5 * (par.t_opt[the_sex] - temperature) * 
                (par.t_opt[the_sex] - temperature) / par.omega_t[the_sex]));
} // end calculate_survival

// add juvenile to stack of survivors
void TSDSeasonal::add_juv_to_survivors(std::vector<Individual> &from,
        std::vector<Individual> &to)
{
    assert(from.size() > 0);

    // sample from the local juveniles
    std::sample(from.begin(), 
            from.end(),
            std::back_inserter(to),
            1,
            rng_r);
} // end add_juv_to_survivors()

// replace vacancies with newborn juveniles
void TSDSeasonal::replace()
{
    unsigned n_vacancies;

    unsigned event_idx, patch_origin;
    
    std::vector<double> probabilities{};

    // go through all survivors and assess whether they are breeding
    for (unsigned int patch_idx = 0;
            patch_idx < metapopulation.size();
            ++patch_idx)
    {

        // event 0: sample local male
        probabilities.push_back(
                (1.0 - par.d[male]) * metapopulation[patch_idx].
                male_juveniles.size());

        // event 1: sample dispersing male
        probabilities.push_back(
                par.d[male] * 
                    static_cast<double>(
                        global_productivity[male]) / par.npatches);

        // event 2: sample philopatric female
        probabilities.push_back(
                (1.0 - par.d[female]) * metapopulation[patch_idx].
                    female_juveniles.size());

        // event 3: sample dispersing female
        probabilities.push_back(
                par.d[female] * 
                    static_cast<double>(
                        global_productivity[female]) / par.npatches);

        // make a discrete distribution of the different 
        // events that can take place
        std::discrete_distribution<unsigned> recruitment_event(
                probabilities.begin(), probabilities.end());

        // calculate number of vacancies
        n_vacancies = par.n - metapopulation[patch_idx].
            female_juveniles.size() - metapopulation[patch_idx].male_juveniles.size();

        assert(n_vacancies <= par.n);
        assert(n_vacancies >= 0);

        // replace juveniles
        for (unsigned int replace_idx = 0; 
                replace_idx < n_vacancies;
                ++replace_idx)
        {
            event_idx = recruitment_event(rng_r);

            switch(event_idx)
            {
                case 0: // recruit local male
                    {
                        patch_origin = patch_idx;

                        unsigned int sizet0 = metapopulation[patch_idx].male_survivors.size();

                        add_juv_to_survivors(
                            metapopulation[patch_origin].male_juveniles,
                            metapopulation[patch_idx].male_survivors);

                        assert(metapopulation[patch_idx].male_survivors.size() == sizet0 + 1);
                    }
                    break;
                    
                case 1:  // recruit remote male
                    {
                        patch_origin = productivity_distribution_m(rng_r); 
                        
                        unsigned int sizet0 = metapopulation[patch_idx].male_survivors.size();

                        add_juv_to_survivors(
                            metapopulation[patch_origin].male_juveniles,
                            metapopulation[patch_idx].male_survivors);

                        assert(metapopulation[patch_idx].male_survivors.size() == sizet0 + 1);
                    }
                    break;
                
                case 2: // recruit local female
                    {
                        patch_origin = patch_idx;

                        unsigned int sizet0 = metapopulation[patch_idx].female_survivors.size();

                        add_juv_to_survivors(
                            metapopulation[patch_origin].female_juveniles,
                            metapopulation[patch_idx].female_survivors);

                        assert(metapopulation[patch_idx].female_survivors.size() == sizet0 + 1);
                    }
                    break;

                case 3: // recruit remote female
                    {
                        patch_origin = productivity_distribution_f(rng_r); 
                        
                        unsigned int sizet0 = metapopulation[patch_idx].female_survivors.size();

                        add_juv_to_survivors(
                            metapopulation[patch_origin].female_juveniles,
                            metapopulation[patch_idx].female_survivors);

                        assert(metapopulation[patch_idx].female_survivors.size() == sizet0 + 1);
                    }
                    break;

                default:
                    throw std::runtime_error("juvenile sampler event out of bounds.");
                    break;
            } // end switch;

        } // end replace_idx

        metapopulation[patch_idx].females = metapopulation[patch_idx].female_survivors;
        metapopulation[patch_idx].males = metapopulation[patch_idx].male_survivors;

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
        << "amplitude;" << par.amplitude << std::endl
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
