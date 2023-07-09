// Temperature-dependent sex determination
// coevolving with phylopatry, burrowing, phenology

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
unsigned int clutch_max = 50; // maximum clutch size
unsigned int n_patches = 50; // maximum clutch size
unsigned int max_generations = 100;
unsigned int generation_perturb = 50;
unsigned int skip = 1;

unsigned int n[2] = {10,10}; // n[Female], nm: females and males per patch

double init_d[2] = {0.0,0.0}; // initialize dispersal
double init_b = 0.0; // initial value of the burrowing trait
double init_sr[2] = {0.0,0.0};

bool d_cue_is_envt = false;

// juvenile survival of males, 
// females in their respective envts
// first index is sex, 2nd index is envt
double v[2][2] = {{0.0,0.0},{0,0}}; 

// juvenile survival of males, 
// females in their respective envts before perturbation
double v_orig[2][2] = {{0.0,0.0},{0,0}}; 

// juvenile survival of males, 
// females post perturbation
double v_pert[2][2] = {{0.0,0.0},{0,0}}; 

// does borrowing depth affect survival, y/n
bool burrow_mod_survival = 0;

// mutation rates
double mu_sr = 0.0; 
double mu_b = 0.0; 
double mu_d[2] = {0.0,0.0}; 
double sdmu = 0.0; 

// rates of environmental change
double s[2] = {0.1,0.5};

// rates of environmental change pre perturbation
double s_orig[2] = {0.1,0.5};

// rates of environmental change post perturbation
double s_pert[2] = {0.1,0.5};


// whether environmental variation is spatial or not
bool spatial = false;

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

double wbar = 0.0;
double var_wbar = 0.0;

// initializes parameters from the command line
void init_pars_from_cmd(int argc, char **argv)
{
    n[Female] = atoi(argv[1]); // number of males and females per patch
    n[Male] = atoi(argv[2]);
    n_patches = atoi(argv[3]); // number of patches
    clutch_max = atoi(argv[4]); // max clutch size

    init_d[Female] = atof(argv[5]); // initial dispersal rates
    init_d[Male] = atof(argv[6]); 
    init_b = atof(argv[7]); // initial value of the burrowing trait
    init_sr[0] = atof(argv[8]); // initial sex ratio
    init_sr[1] = atof(argv[9]);

    v[Female][0] = v_orig[Female][0] = atof(argv[10]); // survival probabilities per sex pre perturbation
    v[Female][1] = v_orig[Female][1] = atof(argv[11]);
    v[Male][0] = v_orig[Male][0] = atof(argv[12]);
    v[Male][1] = v_orig[Male][1] = atof(argv[13]);

    v_pert[Female][0] = atof(argv[14]); // survival probabilities per sex
    v_pert[Female][1] = atof(argv[15]);
    v_pert[Male][0] = atof(argv[16]);
    v_pert[Male][1] = atof(argv[17]);

    burrow_mod_survival = atoi(argv[18]); // does burrowing affect survival?

    s_orig[0] = s[0] = atof(argv[19]); // environmental switch rates pre perturbation
    s_orig[1] = s[1] = atof(argv[20]);

    s_pert[0] = atof(argv[21]);
    s_pert[1] = atof(argv[22]);

    spatial = atoi(argv[23]);

    d_cue_is_envt = atoi(argv[24]);

    p2 = s_orig[0] / (s_orig[0] + s_orig[1]); // fraction of environments at a certain frequency

    mu_sr = atof(argv[25]); // mutation rates
    mu_b = atof(argv[26]);
    mu_d[Female] = atof(argv[27]);
    mu_d[Male] = atof(argv[28]);
    sdmu = atof(argv[29]);
    max_generations = atoi(argv[30]);
    generation_perturb = atoi(argv[31]);
    file_basename = argv[32];
} // end init_pars_from_cmd()



// initialize the population by creating individuals
// and assigning them traits
void initialize_population()
{
    meta_population.reserve(n_patches);

    // loop through all the patches
    for (unsigned int patch_idx = 0; patch_idx < n_patches; ++patch_idx)
    {
        // initialize an empty patch
        Patch deme_i{};

        deme_i.breedersM.reserve(n[Male]);
        deme_i.breedersF.reserve(n[Female]);

        // loop through all the females in the patch
        for (unsigned int female_idx = 0; female_idx < n[Female]; ++female_idx)
        {
            // initialize female
            Individual female_i{};

            // now assign initial values to each allele
            for (unsigned int allele_idx = 0; allele_idx < 2; ++allele_idx)
            {
                female_i.sr[0][allele_idx] = 0.5 * init_sr[0];
                female_i.sr[1][allele_idx] = 0.5 * init_sr[1];

                female_i.d[0][allele_idx] = 0.5 * init_d[0];
                female_i.d[1][allele_idx] = 0.5 * init_d[1];

                female_i.b[allele_idx] = 0.5 * init_b;
            }

            deme_i.breedersF.push_back(female_i);
        } // end for female_idx
        
        // loop through all the males in the patch
        for (unsigned int male_idx = 0; male_idx < n[Male]; ++male_idx)
        {
            // initialize male
            Individual male_i{};

            for (unsigned int allele_idx = 0; allele_idx < 2; ++allele_idx)
            {
                male_i.sr[0][allele_idx] = 0.5 * init_sr[0];
                male_i.sr[1][allele_idx] = 0.5 * init_sr[1];

                male_i.d[0][allele_idx] = 0.5 * init_d[0];
                male_i.d[1][allele_idx] = 0.5 * init_d[1];

                male_i.b[allele_idx] = 0.5 * init_b;
            }

            deme_i.breedersM.push_back(male_i);
        } 

        meta_population.push_back(deme_i);
    } // end for patch_idx

    assert(meta_population[n_patches - 1].breedersF.size() == n[Female]);
    assert(meta_population[n_patches - 1].breedersM.size() == n[Male]);
    assert(meta_population.size() > 0);
    assert(meta_population.size() == n_patches);

    disp_juvsM.reserve(n_patches * n[Female] * clutch_max);
    disp_juvsF.reserve(n_patches * n[Female] * clutch_max);
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
        "init_sr1;" << init_sr[0] << std::endl <<
        "init_sr2;" << init_sr[1] << std::endl <<
        "seed;" << seed << std::endl <<
        "nm;" << n[Male] << std::endl <<
        "nf;" << n[Female] << std::endl <<
        "npatches;" << n_patches << std::endl <<
        "clutch_max;" << clutch_max << std::endl <<
        "spatial;" << spatial << std::endl <<
        "d_cue_is_envt;" << d_cue_is_envt << std::endl <<
        "s1;" << s_orig[0] << std::endl <<
        "s2;" << s_orig[1] << std::endl <<
        "s_pert1;" << s_pert[0] << std::endl <<
        "s_pert2;" << s_pert[1] << std::endl <<
        "vf1;" << v_orig[Female][0] << std::endl <<
        "vf2;" << v_orig[Female][1] << std::endl <<
        "vm1;" << v_orig[Male][0] << std::endl <<
        "vm2;" << v_orig[Male][1] << std::endl <<
        "v_pertf1;" << v_pert[Female][0] << std::endl <<
        "v_pertf2;" << v_pert[Female][1] << std::endl <<
        "v_pertm1;" << v_pert[Male][0] << std::endl <<
        "v_pertm2;" << v_pert[Male][1] << std::endl <<
        "burrow_mod_survival;" << burrow_mod_survival << std::endl <<
        "mu_sr;" << mu_sr << std::endl <<
        "mu_b;" << mu_b << std::endl <<
        "mu_df;" << mu_d[Female] << std::endl <<
        "mu_dm;" << mu_d[Male] << std::endl <<
        "sdmu;" << sdmu << std::endl <<
        "generation_perturb;" << generation_perturb << std::endl <<
        "basename;" << file_basename << std::endl <<
        std::endl;
} // end write_parameters

// write the headers for the file with the statistics
void write_stats_headers(std::ofstream &data_file)
{
    data_file << "generation;sr1;sr2;varsr1;varsr2;b;varb;df;dm;vardf;vardm;freq1;varfreq1;wbar;var_wbar;" << std::endl;
} // end write_stats_headers

// write statistics to file data_file
void write_stats_per_timestep(unsigned int const time_step, std::ofstream &data_file)
{
    // mean values
    double mean_sr[2] = {0.0,0.0};
    double mean_d[2] = {0.0,0.0};
    double mean_b = 0.0;

    // sum of squares
    double ss_sr[2] = {0.0,0.0};
    double ss_d[2] = {0.0,0.0};
    double ss_b = 0.0;

    double z;
    
    double mean_freq1 = 0.0;
    double ss_freq1 = 0.0;

    for (unsigned int patch_idx = 0; patch_idx < n_patches; ++patch_idx)

    {
        assert(meta_population[patch_idx].breedersF.size() == n[Female]);
        assert(meta_population[patch_idx].breedersM.size() == n[Male]);

        z = meta_population[patch_idx].envt_hi;
        mean_freq1 += z;
        ss_freq1 += z*z;

        for (unsigned int female_idx = 0; female_idx < n[Female]; ++female_idx)
        {

            z = meta_population[patch_idx].breedersF[female_idx].sr[0][0]
                    + meta_population[patch_idx].breedersF[female_idx].sr[0][1];
                
            mean_sr[0] += z;
            ss_sr[0] += z * z;

            z = meta_population[patch_idx].breedersF[female_idx].sr[1][0]
                    + meta_population[patch_idx].breedersF[female_idx].sr[1][1];

            mean_sr[1] += z;
            ss_sr[1] += z * z;

            z = meta_population[patch_idx].breedersF[female_idx].d[0][0]
                    + meta_population[patch_idx].breedersF[female_idx].d[0][1];

            mean_d[Female] += z;
            ss_d[Female] += z * z;
            
            z = meta_population[patch_idx].breedersF[female_idx].d[1][0]
                    + meta_population[patch_idx].breedersF[female_idx].d[1][1];

            mean_d[Male] += z;
            ss_d[Male] += z * z;
            
            z = meta_population[patch_idx].breedersF[female_idx].b[0]
                    + meta_population[patch_idx].breedersF[female_idx].b[1];
            mean_b += z;
            ss_b += z * z;
        }
        
        for (unsigned int male_idx = 0; male_idx < n[Male]; ++male_idx)
        {
            z = meta_population[patch_idx].breedersM[male_idx].sr[0][0]
                    + meta_population[patch_idx].breedersM[male_idx].sr[0][1];
                
            mean_sr[0] += z;
            ss_sr[0] += z * z;

            z = meta_population[patch_idx].breedersM[male_idx].sr[1][0]
                    + meta_population[patch_idx].breedersM[male_idx].sr[1][1];

            mean_sr[1] += z;
            ss_sr[1] += z * z;

            z = meta_population[patch_idx].breedersM[male_idx].d[0][0]
                    + meta_population[patch_idx].breedersM[male_idx].d[0][1];

            mean_d[Female] += z;
            ss_d[Female] += z * z;
            
            z = meta_population[patch_idx].breedersM[male_idx].d[1][0]
                    + meta_population[patch_idx].breedersM[male_idx].d[1][1];

            mean_d[Male] += z;
            ss_d[Male] += z * z;
            
            z = meta_population[patch_idx].breedersM[male_idx].b[0]
                    + meta_population[patch_idx].breedersM[male_idx].b[1];
            mean_b += z;
            ss_b += z * z;
        } // end for (int male_idx = 0; male_idx < n[Male]; ++male_idx)

    } // end for (int patch_idx

    mean_sr[0] /= (n[Male] + n[Female]) * n_patches;
    mean_sr[1] /= (n[Male] + n[Female]) * n_patches;
    mean_d[Female] /= (n[Male] + n[Female]) * n_patches;
    mean_d[Male] /= (n[Male] + n[Female]) * n_patches;
    mean_b /= (n[Male] + n[Female]) * n_patches;

    mean_freq1 /= n_patches;

    double var_freq1 = ss_freq1 / n_patches - mean_freq1 * mean_freq1;

    double var_sr[2];

    var_sr[0] = ss_sr[0] / ((n[Male] + n[Female]) * n_patches) 
        - mean_sr[0] * mean_sr[0];

    var_sr[1] = ss_sr[1] / ((n[Male] + n[Female]) * n_patches) 
        - mean_sr[1] * mean_sr[1];

    double var_d[2];

    var_d[Female] = ss_d[Female] / ((n[Male] + n[Female]) * n_patches) 
        - mean_d[Female] * mean_d[Female];

    var_d[Male] = ss_d[Male] / ((n[Male] + n[Female]) * n_patches) 
        - mean_d[Male] * mean_d[Male];

    double var_b = ss_b / ((n[Male] + n[Female]) * n_patches) 
        - mean_b * mean_b;

    data_file << time_step << ";"
        
        << mean_sr[0] <<  ";"
        << mean_sr[1] <<  ";"
        << var_sr[0] <<  ";"
        << var_sr[1] <<  ";"

        << mean_b << ";"
        << var_b << ";"

        << mean_d[Female] <<  ";"
        << mean_d[Male] <<  ";"
        << var_d[Female] <<  ";"
        << var_d[Male] << ";"
        << mean_freq1 <<  ";"
        << var_freq1  << ";"
        << wbar << ";"
        << var_wbar << ";"
        << std::endl;
} // end void write_stats_per_timestep()

// mutational probability
double mutate(double const val, double const mu, std::normal_distribution<double> &mutation_values)
{
    // no mutation
    if (uniform(rng_r) > mu)
    {
        return(val);
    }

    // otherwise return mutated value
    return(val + mutation_values(rng_r));
}

// create a single offspring
void create_offspring(
        Individual &offspring,
        Individual const &mother,
        Individual const &father,
        bool const natal_envt // environment at birth in absence of burrowing
        )
{
    std::normal_distribution<double> mutational_distribution(0.0, sdmu);

    for (unsigned int allele_idx = 0; allele_idx < 2; ++allele_idx)
    {
        // inherit alleles
        offspring.sr[0][0] = std::clamp(
                mutate(
                    mother.sr[0][discrete01(rng_r)]
                    ,mu_sr
                    ,mutational_distribution)
                ,0.0,0.5);

        offspring.sr[0][1] = std::clamp(
                mutate(
                    father.sr[0][discrete01(rng_r)]
                    ,mu_sr
                    ,mutational_distribution)
                ,0.0,0.5);
        
        offspring.sr[1][0] = std::clamp(
                mutate(
                    mother.sr[1][discrete01(rng_r)]
                    ,mu_sr
                    ,mutational_distribution)
                ,0.0,0.5);

        offspring.sr[1][1] = std::clamp(
                mutate(
                    father.sr[1][discrete01(rng_r)]
                    ,mu_sr
                    ,mutational_distribution)
                ,0.0,0.5);

        offspring.d[0][0] = std::clamp(
                mutate(
                    mother.d[0][discrete01(rng_r)]
                    ,mu_d[0]
                    ,mutational_distribution)
                ,0.0,0.5);

        offspring.d[0][1] = std::clamp(
                mutate(
                    father.d[0][discrete01(rng_r)]
                    ,mu_d[0]
                    ,mutational_distribution)
                ,0.0,0.5);
        
        offspring.d[1][0] = std::clamp(
                mutate(
                    mother.d[1][discrete01(rng_r)]
                    ,mu_d[1]
                    ,mutational_distribution)
                ,0.0,0.5);

        offspring.d[1][1] = std::clamp(
                mutate(
                    father.d[1][discrete01(rng_r)]
                    ,mu_d[1]
                    ,mutational_distribution)
                ,0.0,0.5);
        
        offspring.b[0] = std::clamp(
                mutate(
                    mother.b[discrete01(rng_r)]
                    ,mu_b
                    ,mutational_distribution)
                ,0.0,0.5);

        offspring.b[1] = std::clamp(
                mutate(
                    father.b[discrete01(rng_r)]
                    ,mu_b
                    ,mutational_distribution)
                ,0.0,0.5);

    } // end for (int allele_idx = 0; allele_idx < 2; ++allele_idx)


    bool current_envt = natal_envt;

    // if in envt 2, burrowing matters
    // so any nonzero borrowing depth 
    // may switch environment back to 1
    if (natal_envt && uniform(rng_r) < mother.b[0] + mother.b[1])
    {
        current_envt = !natal_envt;
    }
    
    // sex determination - note that sr measures proportions sons 
    // as in most sex allocation papers - so that if any uniform number
    // is larger than the sr expression loci, it must be a daughter
    offspring.sex = Female;

    // sex allocation decision making 
    if (uniform(rng_r) < 
            offspring.sr[current_envt][0] + offspring.sr[current_envt][1])
    {
        offspring.sex = Male;
    }

    offspring.viability = v[offspring.sex][
        burrow_mod_survival == 0 ? natal_envt : current_envt];

} // end create offspring

void mate_produce_offspring()
{
    // clear existing stacks of dispersing juveniles
    disp_juvsM.clear();
    disp_juvsF.clear();

    // sample a random local male on the patch
    std::uniform_int_distribution<unsigned int> male_sampler(0,n[Male] - 1);

    // aux variable to store index of father
    unsigned int father_idx;

    bool local_envt_hi, cue;

    wbar = 0.0;
    var_wbar = 0.0;

    // loop through all patches of the metapopulation
    for (unsigned int patch_idx = 0; patch_idx < n_patches; ++patch_idx)
    {
        assert(meta_population[patch_idx].breedersM.size() > 0);
        assert(meta_population[patch_idx].breedersM.size() == n[Male]);
        
        assert(meta_population[patch_idx].breedersF.size() > 0);
        assert(meta_population[patch_idx].breedersF.size() == n[Female]);

        // reset local counts of juveniles in this stack
        // as we will produce them now
        meta_population[patch_idx].phil_juvsF.clear();
        meta_population[patch_idx].phil_juvsM.clear();

        local_envt_hi = meta_population[patch_idx].envt_hi;

        // loop through all females in a particular site
        for (unsigned int female_idx = 0; female_idx < n[Female]; ++female_idx)
        {
            // loop through a female's clutch
            for (unsigned int egg_i = 0; egg_i < clutch_max; ++egg_i)
            {
                // sample male to sire egg_i
                // according to mate choice distribution
                // more attractive males are more likely to be chosen
                // (at least when p != 0)
                father_idx = male_sampler(rng_r);

                // do bounds checking on the father
                // to see whether there are no bugs and
                // resulting buffer overflows
                assert(father_idx >= 0);
                assert(father_idx < n[Male]);

                // offspring has survived, now allocate it and decide whether 
                // it disperses or not
                Individual Kid{};

                // make the offspring
                create_offspring(Kid, 
                        meta_population[patch_idx].breedersF[female_idx]
                        ,meta_population[patch_idx].breedersM[father_idx]
                        ,local_envt_hi
                        );

                // offspring does not survive, create the next
                if (uniform(rng_r) > Kid.viability)
                {
                    continue;
                }

                ++wbar;
                ++var_wbar;

                cue = d_cue_is_envt ? local_envt_hi : Kid.sex;

                if (uniform(rng_r) < Kid.d[cue][0] + Kid.d[cue][1]) // kid disperses
                {
                    if (Kid.sex == Female)
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
                    if (Kid.sex == Female)
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
                    
                    assert(meta_population[patch_idx].phil_juvsM.size() < 
                            n[Female] * clutch_max);
                    assert(meta_population[patch_idx].phil_juvsF.size() < 
                            n[Female] * clutch_max);
                }

            } // end for for (int egg_i = 0; egg_i < clutch_max; ++egg_i)

        } // end for for (int female_idx = 0; female_idx < n[Female]; ++female_idx)

    } // end for (int patch_idx = 0; patch_idx < n_patches; ++patch_idx)

    wbar /= n_patches * n[Female];
    var_wbar = var_wbar / (n_patches * n[Female]) - wbar * wbar;

    // randomly shuffle dispersers
    shuffle(disp_juvsM.begin(), disp_juvsM.end(), rng_r);
    shuffle(disp_juvsF.begin(), disp_juvsF.end(), rng_r);
} // end mate_produce_offspring()


void adult_mortality_replacement()
{
    // auxiliary variables for the number 
    // of juvenile male and female immigrants
    // available in the local patch
    unsigned int nm_imm_local, nf_imm_local;

    // auxiliary variables for the number 
    // of juvenile male and female immigrants
    // globally available
    unsigned int nm_imm_total = disp_juvsM.size();
    unsigned int nf_imm_total = disp_juvsF.size();

    // if only temporal variation (i.e., all patches
    // switch between envt'al states at the same moment
    bool envt_change = !spatial && uniform(rng_r) < s[meta_population[0].envt_hi];

    // loop through all patches of the metapopulation
    for (unsigned int patch_idx = 0; patch_idx < n_patches; ++patch_idx)
    {
        nm_imm_local = nm_imm_total / n_patches;
        nf_imm_local = nf_imm_total / n_patches;
                    
        // loop through all females in a particular site
        // and replace
        for (unsigned int female_idx = 0; female_idx < n[Female]; ++female_idx)
        {
            assert(nf_imm_local > 0 || 
                    meta_population[patch_idx].phil_juvsF.size() > 0);

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

                std::uniform_int_distribution<unsigned int> 
                    local_female_sampler(0, meta_population[patch_idx].phil_juvsF.size() - 1);

                unsigned int sampled_juvF = local_female_sampler(rng_r);

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
        } // end female_idx

        // loop through all males and replace
        for (unsigned int male_idx = 0; male_idx < n[Male]; ++male_idx)
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

                assert(meta_population[patch_idx].phil_juvsM.size() < n[Female] * clutch_max);

                std::uniform_int_distribution<unsigned int> 
                    local_male_sampler(0, meta_population[patch_idx].phil_juvsM.size() - 1);

                unsigned int sampled_juvM = local_male_sampler(rng_r);

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
            } // end if else local vs global recruitment
        }  // end for (int male_idx = 0; male_idx < n[Male]; ++male_idx)

        // environmental change in case of spatial environments
        if ((spatial && uniform(rng_r) < s[meta_population[patch_idx].envt_hi])
                || envt_change) 
        {
                meta_population[patch_idx].envt_hi = 
                    !meta_population[patch_idx].envt_hi;
        }
    } // end for (int patch_idx = 0; patch_idx < n_patches; ++patch_idx)
} // end adult_mortality_replacement()

void environmental_perturbation()
{
    for (int envt_idx = 0; envt_idx < 2; ++envt_idx)
    {
        s[envt_idx] = s_pert[envt_idx];

        for (int sex_idx = 0; sex_idx < 2; ++sex_idx)
        {
            v[sex_idx][envt_idx]  = v_pert[sex_idx][envt_idx];
        }
    }
} // end environmental_perturbation()

// the main part of the code
int main(int argc, char **argv)
{
    init_pars_from_cmd(argc, argv);
    // ok parameters initialized, now lets do some work

    // initialize output files
    std::ofstream data_file(file_basename + ".csv");

    // write headers to the data file
    write_stats_headers(data_file);

    initialize_population();

    for (unsigned int generation_idx = 0; 
            generation_idx < max_generations; ++generation_idx)
    {
        mate_produce_offspring();

        adult_mortality_replacement();

        if (generation_idx % skip == 0)
        {
            write_stats_per_timestep(generation_idx, data_file);
        }

        if (generation_idx == generation_perturb)
        {
            environmental_perturbation();
        }
    }

    write_parameters(data_file);
}



