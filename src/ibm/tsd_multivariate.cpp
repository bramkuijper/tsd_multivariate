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

int n[2] = {10,10}; // n[Female], nm: females and males per patch

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
        for (int female_idx = 0; female_idx < n[Female]; ++female_idx)
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

    assert(meta_population[n_patches - 1].breedersF.size() == n[Female]);
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
        "n[Female];" << n[Females] << std::endl <<
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
} // end write_parameters

// write the headers for the file with the statistics
void write_stats_headers(std::ofstream &data_file)
{
    data_file << "generation;sr1;sr2;varsr1;varsr2;b;varb;df;dm;vardf;vardm" << std::endl;
} // end write_stats_headers

// write statistics to file data_file
void write_stats_per_timestep(int time_step, std::ofstream &data_file)
{
    // mean values
    double mean_sr[2] = {0.0,0.0};
    double mean_d[2] = {0.0,0.0};
    double mean_b = 0.0;

    // sum of squares
    double ss_sr[2] = {0.0,0.0};
    double ss_d[2] = {0.0,0.0};
    double ss_b = 0.0;

    for (int patch_idx = 0; patch_idx < n_patches; ++patch_idx)
    {
        assert(meta_population[patch_idx].breedersF.size() == n[Female]);
        assert(meta_population[patch_idx].breedersM.size() == nm);

        for (int female_idx = 0; female_idx < n[Female]; ++female_idx)
        {
            for (int allele_idx = 0; allele_idx < 2; ++allele_idx)
            {
                z = meta_population[patch_idx].breedersF[female_idx].sr[0][allele_idx]
                mean_sr[0] += z;
                ss_sr[0] += z * z;

                z = meta_population[patch_idx].breedersF[female_idx].sr[1][allele_idx]
                mean_sr[1] += z;
                ss_sr[1] += z * z;

                z = meta_population[patch_idx].breedersF[female_idx].d[Female][allele_idx]
                mean_d[Female] += z;
                ss_d[Female] += z * z;
                
                z = meta_population[patch_idx].breedersF[female_idx].d[Male][allele_idx]
                mean_d[Male] += z;
                ss_d[Male] += z * z;
                
                z = meta_population[patch_idx].breedersF[female_idx].b[allele_idx]
                mean_b += z;
                ss_b += z * z;
            }
        } // end for female_idx
            
        
        for (int male_idx = 0; male_idx < n[Male]; ++male_idx)
        {
            for (int allele_idx = 0; allele_idx < 2; ++allele_idx)
            {
                z = meta_population[patch_idx].breedersM[male_idx].sr[0][allele_idx]
                mean_sr[0] += z;
                ss_sr[0] += z * z;

                z = meta_population[patch_idx].breedersM[male_idx].sr[1][allele_idx]
                mean_sr[1] += z;
                ss_sr[1] += z * z;

                z = meta_population[patch_idx].breedersM[male_idx].d[Female][allele_idx]
                mean_d[Female] += z;
                ss_d[Female] += z * z;
                
                z = meta_population[patch_idx].breedersM[male_idx].d[Male][allele_idx]
                mean_d[Male] += z;
                ss_d[Male] += z * z;
                
                z = meta_population[patch_idx].breedersM[male_idx].b[allele_idx]
                mean_b += z;
                ss_b += z * z;
            }
        } // end for (int male_idx = 0; male_idx < n[Male]; ++male_idx)
    } // end for (int patch_idx

    mean_sr[0] /= (n[Male] + n[Female]) * n_patches;
    mean_sr[1] /= (n[Male] + n[Female]) * n_patches;
    mean_df[Female] /= (n[Male] + n[Female]) * n_patches;
    mean_df[Male] /= (n[Male] + n[Female]) * n_patches;
    mean_b /= (n[Male] + n[Female]) * n_patches;

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
        << var_d[Male] << std::endl;
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
        bool const natal_envt // current environment
        )
{
    std::normal_distribution<double> mutational_distribution(0.0, sdmu);

    for (int allele_idx = 0; allele_idx < 2; ++allele_idx)
    {
        // inherit alleles
        offspring.sr[0][0] = mutate(
                mother.sr[0][discrete01(rng_r)]
                ,mu_sr
                ,mutational_distribution);

        offspring.sr[0][1] = mutate(
                father.sr[0][discrete01(rng_r)]
                ,mu_sr
                ,mutational_distribution);
        
        offspring.sr[1][0] = mutate(
                mother.sr[1][discrete01(rng_r)]
                ,mu_sr
                ,mutational_distribution);

        offspring.sr[1][1] = mutate(
                father.sr[1][discrete01(rng_r)]
                ,mu_sr
                ,mutational_distribution);

        offspring.d[Female][0] = mutate(
                mother.d[Female][discrete01(rng_r)]
                ,mu_d[Female]
                ,mutational_distribution);

        offspring.d[Female][1] = mutate(
                father.d[Female][discrete01(rng_r)]
                ,mu_d[Female]
                ,mutational_distribution);
        
        offspring.d[Male][0] = mutate(
                mother.d[Male][discrete01(rng_r)]
                ,mu_d[Male]
                ,mutational_distribution);

        offspring.d[Male][1] = mutate(
                father.d[Male][discrete01(rng_r)]
                ,mu_d[Male]
                ,mutational_distribution);
        
        offspring.b[0] = mutate(
                mother.b[discrete01(rng_r)]
                ,mu_b
                ,mutational_distribution);

        offspring.b[1] = mutate(
                father.b[discrete01(rng_r)]
                ,mu_b
                ,mutational_distribution);
    } // end for (int allele_idx = 0; allele_idx < 2; ++allele_idx)

    // sex determination - note that sr measures proportions sons 
    // as in most sex allocation papers - so that if any uniform number
    // is larger than the sr expression loci, it must be a daughter
    offspring.sex = Male;

    bool current_envt = natal_envt;

    if (natal_envt && uniform(rng_r) < mother.b)
    {
        current_envt = !natal_envt;
    }

    // sex allocation decision making 
    if (uniform(rng_r) > 
            offspring.sr[current_envt][0] + offspring.sr[current_envt][1])
    {
        offspring.sex = Female;

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
    std::uniform_int_distribution<int> male_sampler(0,n[Male] - 1);

    // aux variable to store local patch environment
    bool local_envt_hi; 

    // loop through all patches of the metapopulation
    for (int patch_idx = 0; patch_idx < n_patches; ++patch_idx)
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
        for (int female_idx = 0; female_idx < n[Female]; ++female_idx)
        {
            // loop through a female's clutch
            for (int egg_i = 0; egg_i < clutch_max; ++egg_i)
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
                assert(father_idx < n[Males]);

                // offspring has survived, now allocate it and decide whether 
                // it disperses or not
                Individual Kid{};

                // make the offspring
                create_offspring(Kid, 
                        meta_population[patch_idx].breedersF[female_idx]
                        ,meta_population[patch_idx].breedersM[father_idx]
                        ,meta_population[patch_idx].envt_hi
                        );

                if (uniform(rng_r) > offspring.viability)
                {
                    continue;
                }

                if (uniform(rng_r) < d[Kid.sex]) // kid disperses
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
                    
                    assert(meta_population[patch_idx].phil_juvsM.size() < n[Female] * clutch_max);
                    assert(meta_population[patch_idx].phil_juvsF.size() < n[Female] * clutch_max);
                }

            } // end for for (int egg_i = 0; egg_i < clutch_max; ++egg_i)

        } // end for for (int female_idx = 0; female_idx < n[Female]; ++female_idx)

    } // end for (int patch_idx = 0; patch_idx < n_patches; ++patch_idx)

    // randomly shuffle dispersers
    shuffle(disp_juvsM.begin(), disp_juvsM.end(), rng_r);
    shuffle(disp_juvsF.begin(), disp_juvsF.end(), rng_r);

} // end mate_produce_offspring()


void adult_mortality_replacement()
{
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
        for (int female_idx = 0; female_idx < n[Female]; ++female_idx)
        {
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
        } // end female_idx

        for (int male_idx = 0; male_idx < n[Male]; ++male_idx)
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
