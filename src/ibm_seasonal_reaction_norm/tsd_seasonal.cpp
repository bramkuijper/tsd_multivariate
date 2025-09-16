#include <cassert>
#include <stdexcept>
#include <algorithm>
#include <vector>
#include <cmath>
#include "tsd_seasonal.hpp"
#include "parameters.hpp"

// the evolution of TSD in a seasonal environment
// when individuals can simply evolve a time at which to reproduce
// rather than relying on a cue indicative of seasonality

// constructor
TSDSeasonal::TSDSeasonal(Parameters const &par) :
    par{par},
    data_file{par.file_name},
    uniform{0.0,1.0},
    rd{},
    seed{rd()},
    rng_r{seed},
    metapopulation(par.npatches,Patch(par)) // initialize meta population
{
    // start simulation by writing headers to output file
    write_headers();

    // now run the simulation
    for (time_step = 0;  
            time_step <= par.max_simulation_time; ++time_step)
    {

        // reproduction happens every time step to the 
        // subset of the individuals willing to reproduce
        reproduce();
        update_environment();

        // once per season replace adults
        if (time_step > 0 && time_step % par.max_t_season == 0)
        {
            adult_survival();
            fill_vacancies();
            replace();

            // clear the stack of juveniles and be ready for the
            // next season
            clear_juveniles();
            reset_adult_breeding_status();
        }
    
        if (time_step % par.skip_output == 0 || 
                (time_step > par.max_simulation_time / 2 - par.interval / 2 && 
                 time_step < par.max_simulation_time / 2 + par.interval / 2))
        {
            write_data();
        }
    }

    write_parameters();
} // end constructor

// once an individual is selected to breed at a
// particular timestep (based on its seasonal reaction norm)
// it cannot breed another time during the season.
// Only once the season starts again, will the breeding status
// of this individual be reset
void TSDSeasonal::reset_adult_breeding_status()
{
    n_already_attempted[male] = 0;
    n_already_attempted[female] = 0;

    // update the environment in each patch wi
    for (unsigned int patch_idx = 0; 
            patch_idx < metapopulation.size();
            ++patch_idx)
    {
        // reset males 
        for (auto male_iterator = metapopulation[patch_idx].males.begin(); 
                male_iterator != metapopulation[patch_idx].males.end();
                ++male_iterator)
        {
            male_iterator->attempted_to_mate = false;
        }

        // reset females
        for (auto female_iterator = metapopulation[patch_idx].females.begin(); 
                female_iterator != metapopulation[patch_idx].females.end();
                ++female_iterator)
        {
            female_iterator->attempted_to_mate = false;
        }
    } // patch_idx
} // end reset_adult_breeding_status()

// envt is sinusoidal varying over time from -1 to 1 + error
void TSDSeasonal::update_environment()
{
    // update the intercept to a different one once we are halfway
    // the simulation time
    double intercept = time_step < par.max_simulation_time / 2 ? 
        par.temperature_intercept 
        : 
        par.temperature_intercept_change;

    // update the environment in each patch wi
    for (unsigned int patch_idx = 0; 
            patch_idx < metapopulation.size();
            ++patch_idx)
    {
        metapopulation[patch_idx].temperature = 
            intercept + 
            par.amplitude * std::sin(time_step * 2 * M_PI / par.max_t_season) +
            standard_normal(rng_r) * par.temp_error_sd;
    }
}// end update_environment()

void TSDSeasonal::adult_survival()
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
                    static_cast<unsigned>(patch_iter->females.size()),
                    par.survival_prob[female]);
        
        std::binomial_distribution<unsigned> 
            male_survivor_sampler(
                    static_cast<unsigned>(patch_iter->males.size()),
                    par.survival_prob[male]);


        // get number of female survivors
        unsigned n_female_survivors = female_survivor_sampler(rng_r);
        unsigned n_male_survivors = male_survivor_sampler(rng_r);

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


void TSDSeasonal::clear_juveniles()
{
    // go through all survivors and assess whether they are breeding
    for (unsigned int patch_idx = 0;
            patch_idx < metapopulation.size();
            ++patch_idx)
    {
        metapopulation[patch_idx].male_juveniles.clear();
        metapopulation[patch_idx].female_juveniles.clear();
    }
} // end clear_juveniles()

void TSDSeasonal::reproduce()
{
    // reset fecundity variable
    fecundity = 0;

    // aux variable for timing of reproduction
    double prob_reproduce, cue;

    // empty vector with available local males
    std::vector <unsigned int> available_local_males{};
    std::vector <unsigned int> already_attempted_to_mate{};

    // reset stats
    n_available_adults[male] = 0;
    n_available_adults[female] = 0;
    
    // go through all survivors and assess whether they are breeding
    for (unsigned int patch_idx = 0;
            patch_idx < metapopulation.size();
            ++patch_idx)
    {
        available_local_males.clear();

        for (unsigned int male_idx = 0; 
                male_idx < metapopulation[patch_idx].males.size();
                ++male_idx)
        {
            if (!metapopulation[patch_idx].males[male_idx].attempted_to_mate)
            {
                cue = metapopulation[patch_idx].temperature +
                          standard_normal(rng_r) * par.cue_error;

                prob_reproduce = metapopulation[patch_idx].males[male_idx].t +
                    metapopulation[patch_idx].males[male_idx].tb * cue;

                if (uniform(rng_r) < prob_reproduce)
                {
                    available_local_males.push_back(male_idx);
                    metapopulation[patch_idx].males[male_idx].attempted_to_mate = true;
                }
            }
        } // end male_idx

        // update counter for the stats function on available number of 
        // individuals at this time step
        n_available_adults[male] += static_cast<unsigned>(available_local_males.size());
        n_already_attempted[male] += static_cast<unsigned>(available_local_males.size());

        // mix the list of available males
        std::shuffle(available_local_males.begin(),
                available_local_males.end(),
                rng_r);

        // now go through all females and make offspring if 
        // season is alright
        for (unsigned int female_idx = 0; 
                female_idx < metapopulation[patch_idx].females.size();
                ++female_idx)
        {
            assert(metapopulation[patch_idx].females[female_idx].is_female);

            if (metapopulation[patch_idx].females[female_idx].attempted_to_mate)
            {
                continue;
            }
            
            cue = metapopulation[patch_idx].temperature +
                      standard_normal(rng_r) * par.cue_error;

            prob_reproduce = metapopulation[patch_idx].females[female_idx].t +
                metapopulation[patch_idx].females[female_idx].tb * cue;

            // check whether a female indeed reproduces at this time step
            if (uniform(rng_r) < prob_reproduce)
            {
                metapopulation[patch_idx].females[female_idx].attempted_to_mate = true;
                
                // update counter for the stats function on available number of 
                // individuals at this time step
                ++n_available_adults[female];
                ++n_already_attempted[female];

                if (available_local_males.size() > 0)
                {
                    for (unsigned egg_idx = 0; egg_idx < par.fecundity; ++egg_idx)
                    {
                        // obtain father_idx
                        unsigned int father_idx = available_local_males[
                            egg_idx % static_cast<unsigned>(available_local_males.size())];

                        assert(father_idx >= 0);
                        assert(father_idx < metapopulation[patch_idx].males.size());

                        Individual Kid(metapopulation[patch_idx].females[female_idx],
                                    metapopulation[patch_idx].males[father_idx],
                                    par,
                                    rng_r);

                        // calculate individual SR
                        double p_female = Kid.determine_sex(
                                metapopulation[patch_idx].temperature);

                        // realise sex determination
                        Kid.is_female = uniform(rng_r) < p_female;

                        // survive and add to pool of juveniles
                        if (Kid.is_female)
                        {
                            if (uniform(rng_r) < 
                                    calculate_survival(patch_idx, female))
                            {
                                metapopulation[patch_idx].female_juveniles.push_back(Kid);
                            }
                        }
                        else
                        {
                            if (uniform(rng_r) < 
                                    calculate_survival(patch_idx, male))
                            {
                                metapopulation[patch_idx].male_juveniles.push_back(Kid);
                            }
                        }
                    } // end for()
                } // if local_available_males.size() > 0
            } // end if uniform() <  prob_reproduce
        } // end for female_idx
    } // end for patch_idx
} // end reproduce()


void TSDSeasonal::calculate_patch_productivities()
{
    // vectors to store all local productivities
    // of males and females, allowing us to make 
    // a productivity distribution
    std::vector<unsigned> local_productivity_male{};
    std::vector<unsigned> local_productivity_female{};
     
    // reset sum of all productivities to 0
    global_productivity[male] = 0;
    global_productivity[female] = 0;

    unsigned njuv_males,njuv_females;

    // go through all survivors and assess whether they are breeding
    for (unsigned int patch_idx = 0;
            patch_idx < metapopulation.size();
            ++patch_idx)
    {
        njuv_males = static_cast<unsigned>(
                metapopulation[patch_idx].male_juveniles.size()
                );

        njuv_females = static_cast<unsigned>(
                metapopulation[patch_idx].female_juveniles.size()
                );

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
double TSDSeasonal::calculate_survival(unsigned int const patch_idx, Sex const the_sex)
{
    return(std::exp(-0.5 * 
                (par.t_opt[the_sex] - metapopulation[patch_idx].temperature) * 
                (par.t_opt[the_sex] - metapopulation[patch_idx].temperature) / 
                    par.omega_t[the_sex]));

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
// this should only happen at the end of each season
// after adult survival
void TSDSeasonal::fill_vacancies()
{
    // first calculate all the productivities
    calculate_patch_productivities();

    unsigned n_vacancies;

    unsigned event_idx, patch_origin;
    
    std::vector<double> probabilities{};

    double sumprob{0.0};


    // go through all survivors and assess whether they are breeding
    for (unsigned int patch_idx = 0;
            patch_idx < metapopulation.size();
            ++patch_idx)
    {
        probabilities.clear();
        sumprob = 0.0;

        double plocal_male = 
                (1.0 - par.d[male]) * metapopulation[patch_idx].
                male_juveniles.size();

        // event 0: sample local male
        probabilities.push_back(plocal_male);
        sumprob += plocal_male;

        // event 1: sample dispersing male
        double premote_male = 
                par.d[male] * 
                    static_cast<double>(
                        global_productivity[male]) / par.npatches;

        probabilities.push_back(premote_male);
        sumprob += premote_male;

        // event 2: sample philopatric female
        double plocal_female = (1.0 - par.d[female]) * metapopulation[patch_idx].
                    female_juveniles.size();

        probabilities.push_back(plocal_female);
        sumprob += plocal_female;

        // event 3: sample dispersing female
        double premote_female = par.d[female] * 
                    static_cast<double>(
                        global_productivity[female]) / par.npatches;

        probabilities.push_back(premote_female);
        sumprob += premote_female;
        
        if (sumprob == 0.0)
        {
            continue;
        }


        // make a discrete distribution of the different 
        // events that can take place
        std::discrete_distribution<unsigned> recruitment_event(
                probabilities.begin(), probabilities.end());

        // calculate number of vacancies
        n_vacancies = par.n - 
            static_cast<unsigned>(
                    metapopulation[patch_idx].female_survivors.size()
                    ) - 
            static_cast<unsigned>(
                    metapopulation[patch_idx].male_survivors.size()
                    );

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

                        unsigned int sizet0 = static_cast<unsigned>(
                                metapopulation[patch_idx].male_survivors.size()
                                );

                        add_juv_to_survivors(
                            metapopulation[patch_origin].male_juveniles,
                            metapopulation[patch_idx].male_survivors
                            );

                        assert(metapopulation[patch_idx].
                                male_survivors.size() == sizet0 + 1);
                    }
                    break;
                    
                case 1:  // recruit remote male
                    {
                        patch_origin = productivity_distribution_m(rng_r); 
                        
                        unsigned int sizet0 = static_cast<unsigned>(
                                metapopulation[patch_idx].male_survivors.size()
                                );

                        add_juv_to_survivors(
                            metapopulation[patch_origin].male_juveniles,
                            metapopulation[patch_idx].male_survivors);

                        assert(metapopulation[patch_idx].male_survivors.size() == sizet0 + 1);
                    }
                    break;
                
                case 2: // recruit local female
                    {
                        patch_origin = patch_idx;

                        unsigned int sizet0 = static_cast<unsigned>(
                                metapopulation[patch_idx].female_survivors.size()
                                );

                        add_juv_to_survivors(
                            metapopulation[patch_origin].female_juveniles,
                            metapopulation[patch_idx].female_survivors);

                        assert(metapopulation[patch_idx].female_survivors.size() == sizet0 + 1);
                    }
                    break;

                case 3: // recruit remote female
                    {
                        patch_origin = productivity_distribution_f(rng_r); 

                        assert(patch_origin >= 0);
                        assert(patch_origin < par.npatches);
                        
                        unsigned int sizet0 = static_cast<unsigned>(
                                metapopulation[patch_idx].female_survivors.size()
                                );

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

    } // end for patch_idx
    

} // end fill_vacancies()


// now replace breeders with surviving breeders (which includes filled vacancies)
void TSDSeasonal::replace()
{
    survivors[male] = 0;
    survivors[female] = 0;

    for (auto patch_iter = metapopulation.begin();
            patch_iter != metapopulation.end();
            ++patch_iter)
    {
        survivors[male] += patch_iter->male_survivors.size();
        survivors[female] += patch_iter->female_survivors.size();

        patch_iter->females = patch_iter->female_survivors;
        patch_iter->males = patch_iter->male_survivors;
    }
    // population extinct, we are done
    if (survivors[male] == 0 || survivors[female]==0)
    {
        std::cout << "extinct :-( in generation " << time_step << std::endl;
        write_data();
        write_parameters();
        exit(1);
    }
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
        << "ab_min;" << par.ab_range[0] << std::endl
        << "ab_max;" << par.ab_range[1] << std::endl
        << "amplitude;" << par.amplitude << std::endl
        << "intercept;" << par.temperature_intercept << std::endl
        << "intercept_change;" << par.temperature_intercept_change << std::endl
        << "init_t;" << par.init_t << std::endl
        << "max_t_season;" << par.max_t_season << std::endl
        << "init_a;" << par.init_a << std::endl
        << "init_b;" << par.init_b << std::endl
        << "mu_a;" << par.mu_a << std::endl
        << "mu_b;" << par.mu_b << std::endl
        << "mu_t;" << par.mu_t << std::endl
        << "mu_tb;" << par.mu_tb << std::endl
        << "sdmu;" << par.sdmu << std::endl
        << "temp_error_sd;" << par.temp_error_sd << std::endl
        << "cue_error;" << par.cue_error << std::endl; 
} // end write_parameters

void TSDSeasonal::write_data()
{
    // first calculate all the productivities
    // again to make sure numbers are ok
    calculate_patch_productivities();

    double meana{0.0};
    double ssa{0.0};

    double meanb{0.0};
    double ssb{0.0};
    
    double meant{0.0};
    double sst{0.0};
    
    double meantb{0.0};
    double sstb{0.0};

    double x;

    // number of adult females and males
    unsigned nf{0};
    unsigned nm{0};
    // number of adult females and males
    // which have survived
    unsigned nf_surv{0};
    unsigned nm_surv{0};

    double adult_sr{0.0};
    double adult_sr_surv{0.0};

    double mean_temperature{0.0};
    double ss_temperature{0.0};

    

    // go through all survivors and assess whether they are breeding
    for (auto patch_iter = metapopulation.begin();
            patch_iter != metapopulation.end();
            ++patch_iter)
    {
        nf += patch_iter->females.size();
        nm += patch_iter->males.size();
        nf_surv += patch_iter->female_survivors.size();
        nm_surv += patch_iter->male_survivors.size();

        mean_temperature += patch_iter->temperature;
        ss_temperature += patch_iter->temperature * patch_iter->temperature;

        for (auto female_iter = patch_iter->females.begin();
                female_iter != patch_iter->females.end();
                ++female_iter)
        {
            x = female_iter->a;
            meana += x;
            ssa += x * x;

            x = female_iter->b;
            meanb += x;
            ssb += x * x;

            x = female_iter->t;
            meant += x;
            sst += x * x;
            
            x = female_iter->tb;
            meantb += x;
            sstb += x * x;
        }
        
        for (auto male_iter = patch_iter->males.begin();
                male_iter != patch_iter->males.end();
                ++male_iter)
        {
            x = male_iter->a;
            meana += x;
            ssa += x * x;

            x = male_iter->b;
            meanb += x;
            ssb += x * x;

            x = male_iter->t;
            meant += x;
            sst += x * x;
            
            x = male_iter->tb;
            meantb += x;
            sstb += x * x;
        }
    }

    // adult sex ratio
    adult_sr = (nf + nm) == 0 ? 0 : static_cast<double>(nm) / (nf + nm);
    adult_sr_surv = (nf_surv + nm_surv) == 0 ? 0 : 
        static_cast<double>(nm_surv) / (nf_surv + nm_surv);

    meana /= nf + nm;
    double vara = ssa / (nf + nm) - meana * meana;

    meanb /= nf + nm;
    double varb = ssb / (nf + nm) - meanb * meanb;
    
    meant /= nf + nm;
    double vart = sst / (nf + nm) - meant * meant;
    
    meantb /= nf + nm;
    double vartb = sstb / (nf + nm) - meantb * meantb;

    // calculate sex ratio among surviving juveniles
    double global_juv_sr_after_survival = global_productivity[male] + global_productivity[female] == 0 ? 0 :
        static_cast<double>(global_productivity[male]) / 
            (global_productivity[male] + global_productivity[female]);

    mean_temperature /= par.npatches;
    double var_temperature = ss_temperature / par.npatches - 
        mean_temperature * mean_temperature;

    data_file << time_step << ";" << 
        meana << ";" <<
        vara << ";" <<

        meanb << ";" <<
        varb << ";" <<
        
        meant << ";" <<
        vart << ";" <<

        meantb << ";" <<
        vartb << ";" <<

        // surviving male juvs
        static_cast<double>(global_productivity[male])/par.npatches << ";" << 
        
        // surviving female juvs
        static_cast<double>(global_productivity[female])/par.npatches << ";" << 

        // surviving female adults
        static_cast<double>(nf_surv)/par.npatches << ";" << 
        
        // surviving male adults
        static_cast<double>(nm_surv)/par.npatches << ";" << 

        // available males per time step
        static_cast<double>(n_available_adults[male])/par.npatches << ";" <<
        
        // available females per time step
        static_cast<double>(n_available_adults[female])/par.npatches << ";" <<

        // prior attempts male
        static_cast<double>(n_already_attempted[male])/par.npatches << ";" <<
        
        // prior attempts male
        static_cast<double>(n_already_attempted[female])/par.npatches << ";" <<

        // juvenile sr
        global_juv_sr_after_survival << ";" << 

        //  adult sr
        adult_sr << ";" << 

        //  adult sr
        adult_sr_surv << ";" << 
        static_cast<double>(nf)/par.npatches  << ";" << 
        static_cast<double>(nm)/par.npatches  << ";" << 
        mean_temperature << ";" <<
        var_temperature << ";" <<
        std::endl;
} // end write_data()

void TSDSeasonal::write_headers()
{
    data_file << "time;a;var_a;b;var_b;t;var_t;tb;var_tb;surviving_male_juvs;surviving_female_juvs;surviving_male_adults;surviving_female_adults;available_males;available_females;already_attempted_males;already_attempted_females;surviving_juv_sr;adult_sr;adult_surv_sr;nf;nm;mean_environment;var_environment;" 
        << std::endl;
} // write_headers()
