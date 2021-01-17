// sexual selection in a spatially structured population
// assessing the scope for intragenomic conflict
//

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <cassert>
#include <random>
#include <unistd.h>
#include <yaml-cpp/yaml.h>
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

// parameters
int clutch_max = 50; // maximum clutch size
int n_patches = 50; // maximum clutch size
int max_generations = 100;
int skip = 10;

int nm = 5; // actual number of males/patch
int nf = 5; // actual number of females/patch
double d[2] = {0.0,0.0}; //dm, df: male and female dispersal probabilities
double init_t = 0.0; // initial value of baseline ornament
double init_p = 0.0; // initial value of preference
double init_tprime = 0.0; // initial value of ornament cond-dep 
double p_high = 0.0; // prob male is high quality
double l = 0.0; // probability female mates locally
double base_surv = 0.0; // baseline mortality
double base_mort = 0.0; // baseline mortality
double a = 0.0; // efficacy of mate choice
double cp = 0.0; // cost of female preference
double cs = 0.0; // cost of male ornament
double k = 1.0; // rate at which ornamentation costs decrease with increasing male quality
double fl = 0.0; // juvenile survival of an offspring born from a low quality male
double fh = 0.0; // juvenile survival of an offspring born from a hi quality male
double mu_t = 0.0; 
double mu_tprime = 0.0; 
double mu_p = 0.0; 
double sdmu_t = 0.0; 
double sdmu_tprime = 0.0; 
double sdmu_p = 0.0; 

// base name used for output files
std::string file_basename;


// variable that characterizes which individual/allele
// is controlling expression of a phenotype
enum Expression {

    offspr = 0,
    maternal = 1,
    paternal = 2,
    madumnal = 3,
    padumnal = 4
};

// who controls expression of p
// again, these are parameters that are
// set once the simulation reads in a yaml file 
Expression control_p = offspr;
// who controls expression of t
Expression control_t = offspr;
// who controls expression of t'
Expression control_tprime = offspr;

// (empty) vector of dispersing female juveniles
std::vector <Individual> disp_juvsF;
// (empty) vector of dispersing male juveniles
std::vector <Individual> disp_juvsM;

std::vector <Patch> meta_population;

// initializes parameters from the command line
void init_pars_from_cmd(int arc, char **argv)
{
    nm = atoi(argv[1]);
    nf = atoi(argv[2]);
    init_t = atof(argv[3]);
    init_p = atof(argv[4]);
    init_tprime = atof(argv[5]);
    p_high = atof(argv[6]);
    d[0] = atof(argv[7]); // dm
    d[1] = atof(argv[8]); // df
    base_surv = atof(argv[9]);
    base_mort = atof(argv[10]);
    a = atof(argv[11]);
    cp = atof(argv[12]);
    cs = atof(argv[13]);
    k = atof(argv[14]);
    fl = atof(argv[15]);
    fh = atof(argv[16]);
    control_p = static_cast<Expression>(atoi(argv[17]));
    control_t = static_cast<Expression>(atoi(argv[18]));
    control_tprime = static_cast<Expression>(atoi(argv[19]));
    mu_t = atof(argv[20]);
    mu_p = atof(argv[21]);
    mu_tprime = atof(argv[22]);
    sdmu_t = atof(argv[23]);
    sdmu_p = atof(argv[24]);
    sdmu_tprime = atof(argv[25]);
    file_basename = argv[26];
}

// preferred method: initialize parameter
// from yaml file with parameter
void init_pars_from_yaml(std::string const &yaml_file_name)
{
    // make the yaml file into an object so that we
    // can extract things
    YAML::Node param_file = YAML::LoadFile(yaml_file_name);

    // variables to set if reading the YAML
    // file fails
    bool fail = false;
    std::string fail_param{};

    // check whether the nm parameter exists
    if (param_file["nm"]) {

        // if yes, assign it to the variable
        // nm as an int
        nm = param_file["nm"].as<int>();
    }
    else
    {
        // if not set failure, but continue
        fail = true;
        fail_param = "nm";
    }

    if (param_file["nf"]) {
        nf = param_file["nf"].as<int>();
    }
    else
    {
        fail = true;
        fail_param = "nf";
    }
    
    if (param_file["init_t"]) {
        init_t = param_file["init_t"].as<double>();
    }
    else
    {
        fail = true;
        fail_param = "init_t";
    }
    
    if (param_file["init_p"]) {
        init_p = param_file["init_p"].as<double>();
    }
    else
    {
        fail = true;
        fail_param = "init_p";
    }
    
    if (param_file["init_tprime"]) {
        init_tprime = param_file["init_tprime"].as<double>();
    }
    else
    {
        fail = true;
        fail_param = "init_tprime";
    }
    
    if (param_file["p_high"]) {
        p_high = param_file["p_high"].as<double>();
    }
    else
    {
        fail = true;
        fail_param = "p_high";
    }


    if (param_file["dm"]) {
        d[0] = param_file["dm"].as<double>();
    }
    else
    {
        fail = true;
        fail_param = "dm";
    }
    
    if (param_file["df"]) {
        d[1] = param_file["df"].as<double>();
    }
    else
    {
        fail = true;
        fail_param = "df";
    }
    
    if (param_file["base_surv"]) {
        base_surv = param_file["base_surv"].as<double>();
    }
    else
    {
        fail = true;
        fail_param = "base_surv";
    }
    
    if (param_file["base_mort"]) {
        base_mort = param_file["base_mort"].as<double>();
    }
    else
    {
        fail = true;
        fail_param = "base_mort";
    }

    if (param_file["a"]) {
        a = param_file["a"].as<double>();
    }
    else
    {
        fail = true;
        fail_param = "a";
    }
    
    if (param_file["cp"]) {
        cp = param_file["cp"].as<double>();
    }
    else
    {
        fail = true;
        fail_param = "cp";
    }
    
    if (param_file["cs"]) {
        cs = param_file["cs"].as<double>();
    }
    else
    {
        fail = true;
        fail_param = "cs";
    }
    
    if (param_file["k"]) {
        k = param_file["k"].as<double>();
    }
    else
    {
        fail = true;
        fail_param = "k";
    }
    
    if (param_file["fl"]) {
        fl = param_file["fl"].as<double>();
    }
    else
    {
        fail = true;
        fail_param = "fl";
    }
    
    if (param_file["fh"]) {
        fh = param_file["fh"].as<double>();
    }
    else
    {
        fail = true;
        fail_param = "fh";
    }
    
    if (param_file["control_p"]) {
        control_p = static_cast<Expression>(param_file["control_p"].as<int>());
    }
    else
    {
        fail = true;
        fail_param = "control_p";
    }
    
    if (param_file["control_t"]) {
        control_t = static_cast<Expression>(param_file["control_t"].as<int>());
    }
    else
    {
        fail = true;
        fail_param = "control_t";
    }

    if (param_file["control_tprime"]) {
        control_tprime = static_cast<Expression>(param_file["control_tprime"].as<int>());
    }
    else
    {
        fail = true;
        fail_param = "control_tprime";
    }

    if (param_file["mu_t"]) {
        mu_t = param_file["mu_t"].as<double>();
    }
    else
    {
        fail = true;
        fail_param = "mu_t";
    }
    
    if (param_file["mu_p"]) {
        mu_p = param_file["mu_p"].as<double>();
    }
    else
    {
        fail = true;
        fail_param = "mu_p";
    }
    
    if (param_file["mu_tprime"]) {
        mu_tprime = param_file["mu_tprime"].as<double>();
    }
    else
    {
        fail = true;
        fail_param = "mu_tprime";
    }
    
    if (param_file["sdmu_t"]) {
        sdmu_t = param_file["sdmu_t"].as<double>();
    }
    else
    {
        fail = true;
        fail_param = "sdmu_t";
    }
    
    if (param_file["sdmu_p"]) {
        sdmu_p = param_file["sdmu_p"].as<double>();
    }
    else
    {
        fail = true;
        fail_param = "sdmu_p";
    }
    
    if (param_file["sdmu_tprime"]) {
        sdmu_tprime = param_file["sdmu_tprime"].as<double>();
    }
    else
    {
        fail = true;
        fail_param = "sdmu_tprime";
    }
    
    if (param_file["file_basename"]) {
        file_basename = param_file["file_basename"].as<std::string>();
    }
    else
    {
        fail = true;
        fail_param = "file_basename";
    }
    
    if (param_file["max_generations"]) {
        max_generations = param_file["max_generations"].as<int>();
    }
    else
    {
        fail = true;
        fail_param = "max_generations";
    }
    
    if (param_file["n_patches"]) {
        n_patches = param_file["n_patches"].as<int>();
    }
    else
    {
        fail = true;
        fail_param = "n_patches";
    }
    
    if (param_file["clutch_max"]) {
        clutch_max = param_file["clutch_max"].as<int>();
    }
    else
    {
        fail = true;
        fail_param = "clutch_max";
    }

    if (fail)
    {
        std::cout << "Cannot find '" << fail_param <<  "' in " << yaml_file_name << std::endl;
        throw "Cannot find '" + fail_param + "' in " + yaml_file_name;
    }

} // end init_pars_from_yaml()

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
                female_i.t[allele_idx] = init_t;
                female_i.p[allele_idx] = init_p;
                female_i.tprime[allele_idx] = init_tprime;
            }

            deme_i.breedersF.push_back(female_i);
        }
        
        for (int male_idx = 0; male_idx < nm; ++male_idx)
        {
            // initialize male
            Individual male_i{};

            for (int allele_idx = 0; allele_idx < 2; ++allele_idx)
            {
                male_i.t[allele_idx] = init_t;
                male_i.p[allele_idx] = init_p;
                male_i.tprime[allele_idx] = init_tprime;
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
        "p_high;" << p_high << std::endl <<
        "init_tprime;" << init_tprime << std::endl <<
        "init_t;" << init_t << std::endl <<
        "init_p;" << init_p << std::endl <<
        "seed;" << seed << std::endl <<
        "nm;" << nm << std::endl <<
        "nf;" << nf << std::endl <<
        "dm;" << d[0] << std::endl <<
        "df;" << d[1] << std::endl <<
        "base_surv;" << base_surv << std::endl <<
        "base_mort;" << base_mort << std::endl <<
        "a;" << a << std::endl <<
        "cp;" << cp << std::endl <<
        "cs;" << cs << std::endl <<
        "k;" << k << std::endl <<
        "fl;" << fl << std::endl <<
        "fh;" << fh << std::endl <<
        "control_p;" << control_p << std::endl <<
        "control_t;" << control_t << std::endl <<
        "control_tprime;" << control_tprime << std::endl <<
        "mu_t;" << mu_t << std::endl <<
        "mu_p;" << mu_p << std::endl <<
        "mu_tprime;" << mu_tprime << std::endl <<
        "sdmu_t;" << sdmu_t << std::endl <<
        "sdmu_p;" << sdmu_p << std::endl <<
        "sdmu_tprime;" << sdmu_tprime << std::endl <<
        std::endl;
}


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
