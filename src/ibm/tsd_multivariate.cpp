//

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <cassert>
#include <random>
#include <unistd.h>
#include "individual.hpp"
#include "patch.hpp"

// set up the random number generator using a good way of getting random seed
std::random_device rd;
unsigned seed = rd();
std::mt19937 rng_r(seed);

// create random number sampler for a uniform [0,1] distribution
std::uniform_real_distribution<> uniform(0.0,1.0);

// create random number distribution for discrete 0 or 1 distribution
std::bernoulli_distribution discrete01(0.5);

// parameters -- just the initial values here, 
// we set them later using argc, argv (command line)

int clutch_max = 50; // maximum clutch size
int n_patches = 50; // maximum clutch size
int max_generations = 100;
int skip = 10;

int n[2] = {10,10}; // nf, nm: females and males per patch

double init_d[2] = {0.0,0.0}; // initialize 
double init_b = 0.0; // initial value 
double init_sr[2] = {0.0,0.0};

// juvenile survival of males, 
// females in their respective envts
double v[2][2] = {{0.0,0.0},{0,0}}; 

// mutation rates
double mu_sr = 0.0; 
double mu_b = 0.0; 
double mu_d[2] = {0.0,0.0}; 
double sdmu = 0.0; 

// rates of environmental change
double s[2] = {{0.1,0.5}};

// solely for stats purposes
// frequency of envt 2
double p2 = 0.0;

// base name used for output files
std::string file_basename;

// (empty) vector of dispersing female juveniles
std::vector <Individual> disp_juvsF;
// (empty) vector of dispersing male juveniles
std::vector <Individual> disp_juvsM;

// the total meta population 
std::vector <Patch> meta_population;

// initializes parameters from the command line
void init_pars_from_cmd(int arc, char **argv)
{
    n[Female] = atoi(argv[1]);
    n[Male] = atoi(argv[2]);

    init_d[Female] = atof(argv[3]);
    init_d[Male] = atof(argv[3]);
    init_b = atof(argv[5]);
    init_sr = atof(argv[5]);

    d[Female] = atof(argv[7]); 
    d[Male] = atof(argv[8]); 
   
    
    v[Female][0] = atof(argv[9]);
    v[Female][1] = atof(argv[10]);
    v[Male][0] = atof(argv[11]);
    v[Male][1] = atof(argv[12]);

    s[0] = atof(argv[13]);
    s[1] = atof(argv[14]);

    p2 = s[0] / (s[0] + s[1]);

    mu_sr = atof(argv[20]);
    mu_b = atof(argv[21]);
    mu_d = atof(argv[22]);
    sdmu = atof(argv[23]);
    file_basename = argv[26];
} // end init_pars_from_cmd()



// initialize the population by creating individuals
// and assigning them traits
void initialize_population()
{
    // loop through all the patches
    for (int patch_idx = 0; patch_idx < n_patches; ++patch_idx)
    {
        // initialize an empty patch
        Patch deme_i{};

        // loop through all the females in the patch
        for (int female_idx = 0; female_idx < nf; ++female_idx)
        {
            // initialize female
            Individual female_i{};

            // now assign initial values to each allele
            for (int allele_idx = 0; allele_idx < 2; ++allele_idx)
            {
                female_i.sr[0][allele_idx] = init_sr[0];
                female_i.sr[1][allele_idx] = init_sr[1];

                female_i.d[Female][allele_idx] = init_d[Female];
                female_i.d[Male][allele_idx] = init_d[Male];

                female_i.b[allele_idx] = init_b;
            }

            deme_i.breedersF.push_back(female_i);
        } // end for female_idx
        
        // loop through all the males in the patch
        for (int male_idx = 0; male_idx < nm; ++male_idx)
        {
            // initialize male
            Individual male_i{};

            for (int allele_idx = 0; allele_idx < 2; ++allele_idx)
            {
                male_i.sr[0][allele_idx] = init_sr[0];
                male_i.sr[1][allele_idx] = init_sr[1];

                male_i.d[Female][allele_idx] = init_d[Female];
                male_i.d[Male][allele_idx] = init_d[Male];

                male_i.b[allele_idx] = init_b;
            }

            male_i.envt_quality_high = 
                uniform(rng_r) < p_high;


            deme_i.breedersM.push_back(male_i);
        } 

        meta_population.push_back(deme_i);
    } // end for patch_idx

    assert(meta_population[n_patches - 1].breedersF.size() == nf);
    assert(meta_population[n_patches - 1].breedersM.size() == nm);
    assert(meta_population.size() > 0);
    assert(meta_population.size() == n_patches);
} // end initialize_population()

// function to write parameters 
// as a two column table
// to file data_file
void write_parameters(std::ofstream &data_file)
{
    data_file << 
        std::endl << 
        std::endl <<
        "p2;" << p2 << std::endl <<
        "init_df;" << init_d[Female] << std::endl <<
        "init_dm;" << init_d[Male] << std::endl <<
        "init_b;" << init_b << std::endl <<
        "init_sr;" << init_sr << std::endl <<
        "seed;" << seed << std::endl <<
        "nm;" << n[Males] << std::endl <<
        "nf;" << n[Females] << std::endl <<
        "s1;" << s[0] << std::endl <<
        "s2;" << s[1] << std::endl <<
        "vf1;" << v[Female][0] << std::endl <<
        "vf2;" << v[Female][1] << std::endl <<
        "vm1;" << v[Male][0] << std::endl <<
        "vm2;" << v[Male][1] << std::endl <<
        "mu_sr;" << mu_sr << std::endl <<
        "mu_b;" << mu_p << std::endl <<
        "mu_d;" << mu_d << std::endl <<
        "sdmu;" << sdmu_p << std::endl <<
        "basename;" << file_basename << std::endl <<
        std::endl;
}

// TODO



void write_stats_headers(std::ofstream &data_file)
{
    data_file << "generation;t;p;tprime;var_t;var_p;var_tprime" << std::endl;
}

// write statistics to file data_file
void write_stats_per_timestep(int time_step, std::ofstream &data_file)
{
    double mean_t = 0.0;
    double mean_tprime = 0.0;
    double mean_p = 0.0;
    double ss_t = 0.0;
    double ss_tprime = 0.0;
    double ss_p = 0.0;

    double z;

    for (int patch_idx = 0; patch_idx < n_patches; ++patch_idx)
    {
        assert(meta_population[patch_idx].breedersF.size() == nf);
        assert(meta_population[patch_idx].breedersM.size() == nm);

        for (int female_idx = 0; female_idx < nf; ++female_idx)
        {
            for (int allele_idx = 0; allele_idx < 2; ++allele_idx)
            {
                z = meta_population[patch_idx].breedersF[female_idx].t[allele_idx];
                mean_t += z;
                ss_t += z * z;

                z = meta_population[patch_idx].breedersF[female_idx].p[allele_idx];
                mean_p += z;
                ss_p += z * z;

                z = meta_population[patch_idx].breedersF[female_idx].tprime[allele_idx];
                mean_tprime += z;
                ss_tprime += z * z;
            }
        }
            
        
        for (int male_idx = 0; male_idx < nm; ++male_idx)
        {

            for (int allele_idx = 0; allele_idx < 2; ++allele_idx)
            {
                z = meta_population[patch_idx].breedersM[male_idx].t[allele_idx];
                mean_t += z;
                ss_t += z * z;

                z = meta_population[patch_idx].breedersM[male_idx].p[allele_idx];
                mean_p += z;
                ss_p += z * z;

                z = meta_population[patch_idx].breedersM[male_idx].tprime[allele_idx];
                mean_tprime += z;
                ss_tprime += z * z;
            }
        }
    } // end for (int patch_idx

    mean_t /= (nm + nf) * n_patches;
    mean_p /= (nm + nf) * n_patches;
    mean_tprime /= (nm + nf) * n_patches;

    double var_t = ss_t / ((nm + nf) * n_patches) - mean_t * mean_t;
    double var_p = ss_p / ((nm + nf) * n_patches) - mean_p * mean_p;
    double var_tprime = ss_tprime / ((nm + nf) * n_patches) - mean_tprime * mean_tprime;

    data_file << time_step << ";"
        << mean_t <<  ";"
        << mean_p <<  ";"
        << mean_tprime <<  ";"
        << var_t <<  ";"
        << var_p <<  ";"
        << var_tprime <<  std::endl;
} // end void write_stats_per_timestep()

double mutate(double const val, double const mu, double const sdmu)
{
    // no mutation
    if (uniform(rng_r) > mu)
    {
        return(val);
    }

    std::normal_distribution<double> mutation_values(0.0, sdmu);

    // otherwise return mutated value
    return(val + mutation_values(rng_r));
}

void create_offspring(
        Individual &offspring,
        Individual const &mother,
        Individual const &father
        )
{
    for (int allele_idx = 0; allele_idx < 2; ++allele_idx)
    {
        // inherit alleles
        offspring.p[0] = mutate(mother.p[discrete01(rng_r)], mu_p, sdmu_p);
        offspring.p[1] = mutate(father.p[discrete01(rng_r)], mu_p, sdmu_p);
        offspring.t[0] = mutate(mother.t[discrete01(rng_r)], mu_t, sdmu_t);
        offspring.t[1] = mutate(father.t[discrete01(rng_r)], mu_t, sdmu_t);
        offspring.tprime[0] = mutate(mother.tprime[discrete01(rng_r)], mu_tprime, sdmu_tprime);
        offspring.tprime[1] = mutate(father.tprime[discrete01(rng_r)], mu_tprime, sdmu_tprime);
    }

    // set environment for this offspring
    //
    offspring.envt_quality_high = uniform(rng_r) < p_high;

    offspring.p_phen = 0.0;

    // express loci dependent on who is in control
    switch(control_p) {

        case maternal:
            offspring.p_phen = mother.p[0] + mother.p[1];
            break;
        case paternal:
            offspring.p_phen = father.p[0] + father.p[1];
            break;
        case madumnal:
            offspring.p_phen = offspring.p[0];
            break;
        case padumnal:
            offspring.p_phen = offspring.p[1];
            break;
        case offspr:
            offspring.p_phen = offspring.p[0] + offspring.p[1];
            break;
        default:
            std::cout << "p_phen: should not have happened.." << std::endl;
            break;
    }

    // express baseline ornament
    double t_phen = 0.0;

    switch(control_t) {

        case maternal:
            t_phen = mother.t[0] + mother.t[1];
            break;
        case paternal:
            t_phen = father.t[0] + father.t[1];
            break;
        case madumnal:
            t_phen = offspring.t[0];
            break;
        case padumnal:
            t_phen = offspring.t[1];
            break;
        case offspr:
            t_phen = offspring.t[0] + offspring.t[1];
            break;
        default:
            std::cout << "t_phen: should not have happened.." << std::endl;
            break;
    }
    
    // express ornament plasticity
    double tprime_phen = 0.0;

    switch(control_tprime) {

        case maternal:
            tprime_phen = mother.tprime[0] + mother.tprime[1];
            break;
        case paternal:
            tprime_phen = father.tprime[0] + father.tprime[1];
            break;
        case madumnal:
            tprime_phen = offspring.tprime[0];
            break;
        case padumnal:
            tprime_phen = offspring.tprime[1];
            break;
        case offspr:
            tprime_phen = offspring.tprime[0] + offspring.tprime[1];
            break;
        default:
            std::cout << "tprime_phen: should not have happened.." << std::endl;
            break;
    }

    offspring.s = t_phen + tprime_phen * offspring.envt_quality_high;
} // end create offspring

void mate_produce_offspring()
{
    // clear existing stacks of dispersing juveniles
    disp_juvsM.clear();
    disp_juvsF.clear();

    // vector to store the odds of choosing this mate according 
    // to open-ended preferences
    // i.e., exp(a p s) 
    std::vector <double> choice_odds(nm,0.0);

    // auxiliary variable storing current female preference
    double p, s;

    // auxiliary variable to check whether kid is male or female
    bool kid_is_female;

    int father_idx;

    // loop through all patches of the metapopulation
    for (int patch_idx = 0; patch_idx < n_patches; ++patch_idx)
    {
        assert(meta_population[patch_idx].breedersM.size() > 0);
        assert(meta_population[patch_idx].breedersM.size() == nm);
        
        assert(meta_population[patch_idx].breedersF.size() > 0);
        assert(meta_population[patch_idx].breedersF.size() == nf);


        // loop through all females in a particular site
        for (int female_idx = 0; female_idx < nf; ++female_idx)
        {
            // reset local counts of juveniles in this stack
            // as we will produce them now
            meta_population[patch_idx].phil_juvsF.clear();
            meta_population[patch_idx].phil_juvsM.clear();

            // express preference
            p = meta_population[patch_idx].breedersF[female_idx].p_phen; 

            // loop through all males in the site
            // and express ornaments and calculate mate choice odds
            for (int male_idx = 0; male_idx < nm; ++male_idx)
            {
                s = meta_population[patch_idx].breedersM[male_idx].s; 

                choice_odds[male_idx] = exp(a * s * p);
            }

            // make distribution of males to choose from
            std::discrete_distribution <int> mate_choice_distribution(
                    choice_odds.begin()
                    ,choice_odds.end());

            for (int egg_i = 0; egg_i < clutch_max; ++egg_i)
            {
                // sample male to sire egg_i
                // according to mate choice distribution
                // more attractive males are more likely to be chosen
                // (at least when p != 0)
                father_idx = mate_choice_distribution(rng_r);

                // do bounds checking on the father
                // to see whether there are no bugs and
                // resulting buffer overflows
                assert(father_idx >= 0);
                assert(father_idx < nm);

                // check if mate is high/low quality and 
                // have offspring survive accordingly
                if (meta_population[patch_idx].breedersM[father_idx].envt_quality_high)
                {
                    // offspring does not survive
                    // continue to the next egg
                    if (uniform(rng_r) > fh)
                    {
                        continue;
                    }
                }
                else
                {
                    if (uniform(rng_r) > fl)
                    {
                        continue;
                    }
                }

                // offspring has survived, now allocate it and decide whether 
                // it disperses or not
                Individual Kid{};

                // make the offspring
                create_offspring(Kid, 
                        meta_population[patch_idx].breedersF[female_idx],
                        meta_population[patch_idx].breedersM[father_idx]);

                // decide sex randomly
                kid_is_female = discrete01(rng_r);

//                std::cout << kid_is_female << std::endl;

                if (uniform(rng_r) < d[kid_is_female]) // kid disperses
                {
                    if (kid_is_female)
                    {
                        disp_juvsF.push_back(Kid);
                    }
                    else
                    {
                        disp_juvsM.push_back(Kid);
                    }
                }
                else // stays in natal site
                {
                    if (kid_is_female)
                    {
                        // assign kid to stack of philopatric juvenile females
                        meta_population[patch_idx].
                            phil_juvsF.push_back(Kid);
                    }
                    else
                    {
                        // assign kid to stack of philopatric juvenile males
                        meta_population[patch_idx].
                            phil_juvsM.push_back(Kid);

                    }
                    
                    assert(meta_population[patch_idx].phil_juvsM.size() < nf * clutch_max);
                    assert(meta_population[patch_idx].phil_juvsF.size() < nf * clutch_max);
                }

            } // end for for (int egg_i = 0; egg_i < clutch_max; ++egg_i)

        } // end for for (int female_idx = 0; female_idx < nf; ++female_idx)

    } // end for (int patch_idx = 0; patch_idx < n_patches; ++patch_idx)

    // randomly shuffle dispersers
    shuffle(disp_juvsM.begin(), disp_juvsM.end(), rng_r);
    shuffle(disp_juvsF.begin(), disp_juvsF.end(), rng_r);

} // end mate_produce_offspring()


void adult_mortality_replacement()
{
    // auxiliary variable storing 
    // preferences, ornaments
    // male quality
    double p,s;

    bool envt_high;

    double mortality_prob_males;

    // auxiliary variables for the number 
    // of juvenile male and female immigrants
    // available in the local patch
    int nm_imm_local, nf_imm_local;

    // auxiliary variables for the number 
    // of juvenile male and female immigrants
    // globally available
    int nm_imm_total = disp_juvsM.size();
    int nf_imm_total = disp_juvsF.size();

    // loop through all patches of the metapopulation
    for (int patch_idx = 0; patch_idx < n_patches; ++patch_idx)
    {
        nm_imm_local = nm_imm_total / n_patches;
        nf_imm_local = nf_imm_total / n_patches;
                    
        // loop through all females in a particular site
        for (int female_idx = 0; female_idx < nf; ++female_idx)
        {
            p = meta_population[patch_idx].breedersF[female_idx].p_phen;

            // female survival affected by cost of preference
            if (uniform(rng_r) < 
                    base_mort + (1.0 - base_mort) * (1.0 - (base_surv + (1.0 - base_surv) * exp(-cp * p * p))))
            {
                // female dies.
//                std::cout << meta_population[patch_idx].njuvsF << std::endl;

                // choose immigrant rathern than local female
                if (uniform(rng_r) < (double) nf_imm_local  / 
                        (nf_imm_local + meta_population[patch_idx].phil_juvsF.size()))
                {
                    assert(disp_juvsF.size() > 0);

                    meta_population[patch_idx].breedersF[female_idx] = 
                        disp_juvsF.back();

                    disp_juvsF.pop_back();

                    // reduce number of 
                    // immigrant juvenile females for this local patch
                    --nf_imm_local;
                }
                else
                {
                    assert(meta_population[patch_idx].phil_juvsF.size() > 0);

                    std::uniform_int_distribution<int> 
                        local_female_sampler(0, meta_population[patch_idx].phil_juvsF.size() - 1);

                    int sampled_juvF = local_female_sampler(rng_r);

                    assert(sampled_juvF >= 0);
                    assert(sampled_juvF < meta_population[patch_idx].phil_juvsF.size());

                    // assign one of the local female juveniles to this spot
                    meta_population[patch_idx].breedersF[female_idx] = 
                        meta_population[patch_idx].phil_juvsF[sampled_juvF];

                    // delete individual by overwriting it with the last
                    // individual from the stack and then reducing
                    // the count
                    meta_population[patch_idx].phil_juvsF[sampled_juvF] = 
                        meta_population[patch_idx].phil_juvsF.back();

                    meta_population[patch_idx].phil_juvsF.pop_back();
                }
            } // end if (uniform(rng_r)
        } // end female_idx

        for (int male_idx = 0; male_idx < nm; ++male_idx)
        {
            // save ornament size
            s = meta_population[patch_idx].breedersM[male_idx].s;   

            // save environmentally-induced quality of the male, v
            envt_high = meta_population[patch_idx].breedersM[male_idx].envt_quality_high;

            // for survival function, see 2nd term of eq. (2a)
            // in Iwasa & Pomiankowski 1999 JTB
            //
            // obviously the mortality probability 
            // c s^2/(1+kv) will be larger than 1
            // when s has a large magnitude, but beyond  
            // changing the curvature of the function somewhat
            // the effect is the same: the individual will die
            mortality_prob_males = 1.0 - std::clamp(cs * s * s/(1.0 + k * envt_high),0.0,1.0);


            if (uniform(rng_r) < 
                    base_mort + (1.0 - base_mort) * mortality_prob_males)
            {
                // male dies
                assert(nm_imm_local >= 0);

                // choose immigrant rathern than local female
                if (uniform(rng_r) < (double) nm_imm_local  / 
                        (nm_imm_local + meta_population[patch_idx].phil_juvsM.size()))
                {
                    assert(disp_juvsM.size() > 0);

                    meta_population[patch_idx].breedersM[male_idx] = 
                        disp_juvsM.back();

                    disp_juvsM.pop_back();

                    // reduce number of 
                    // immigrant juvenile males for this local patch
                    --nm_imm_local;
                }
                else
                {
                    assert(meta_population[patch_idx].phil_juvsM.size() > 0);

                    assert(meta_population[patch_idx].phil_juvsM.size() < nf * clutch_max);

                    std::uniform_int_distribution<int> 
                        local_male_sampler(0, meta_population[patch_idx].phil_juvsM.size() - 1);

                    int sampled_juvM = local_male_sampler(rng_r);

                    assert(sampled_juvM >= 0);
                    assert(sampled_juvM < meta_population[patch_idx].phil_juvsM.size());

                    // assign one of the local female juveniles to this spot
                    meta_population[patch_idx].breedersM[male_idx] = 
                        meta_population[patch_idx].phil_juvsM[sampled_juvM];

                    // delete individual by overwriting it with the last
                    // individual from the stack and then reducing
                    // the count
                    meta_population[patch_idx].phil_juvsM[sampled_juvM] = 
                        meta_population[patch_idx].phil_juvsM.back();

                    // reduce the number of local female juveniles
                    meta_population[patch_idx].phil_juvsM.pop_back();
                }
            } // end if if (uniform(rng_r) < base_surv + ...
        }
    }
} // end adult_mortality_replacement()

// the main part of the code
int main(int argc, char **argv)
{
    // variable storing the name of the yaml parameter file
    std::string yaml_file{};
 
    // first
    int opt;

    std::stringstream usage;
    
    usage << "Usage of this programme: " << argv[0] << " -f location_of_yaml_parameter_file.yaml" << std::endl;

    while ((opt = getopt(argc, argv, "f:")) != -1) 
    {
        if (opt == 'f')
        {
            yaml_file = optarg;
        }
        else
        {
            std::cerr << usage.str();

            exit(EXIT_FAILURE);
        }
    } // end while

    if (argc < 2)
    {
        std::cerr << usage.str();
        exit(EXIT_FAILURE);
    }

    init_pars_from_yaml(yaml_file);

//    init_pars_from_cmd(argc, argv);
    // ok parameters initialized, now lets do some work

    // initialize output files
    std::ofstream data_file(file_basename + ".csv");

    // write headers to the data file
    write_stats_headers(data_file);

    initialize_population();

    for (int generation_idx = 0; generation_idx < max_generations; ++generation_idx)
    {
        mate_produce_offspring();
        adult_mortality_replacement();

        if (generation_idx % skip == 0)
        {
            write_stats_per_timestep(generation_idx, data_file);
        }
    }

    write_parameters(data_file);
}
