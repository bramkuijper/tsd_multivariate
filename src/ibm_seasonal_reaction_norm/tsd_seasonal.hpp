#ifndef _TSD_SEASONAL_HPP_
#define _TSD_SEASONAL_HPP_

#include <vector>
#include <iostream>
#include <fstream>
#include <random>
#include <string>
#include "parameters.hpp"
#include "patch.hpp"

class TSDSeasonal
{
    private:
        Parameters par;
        std::ofstream data_file;

        std::uniform_real_distribution<double> uniform;
        std::normal_distribution<double> standard_normal{0.0,1.0};
        
        // productivity distributions for females and males, 
        // initialized with some dummy values
        std::discrete_distribution<unsigned> productivity_distribution_f{1,1};
        std::discrete_distribution<unsigned> productivity_distribution_m{1,1};

        long unsigned time_step{0};

        std::random_device rd;
        unsigned int seed;
        std::mt19937 rng_r;

        // all the patches of the metapopoulation
        std::vector<Patch> metapopulation;

        // total productivity for males and females
        int global_productivity[2]{0,0};
        unsigned int survivors[2]{0,0};
        unsigned int fecundity{0};

        void adult_survival();
        void reproduce();
        void replace();
        void fill_vacancies();

        void write_headers();
        void write_data();
        void write_parameters();

        // reset an individual's breeding statuses
        // at the start of each season
        void reset_adult_breeding_status();

        // update the state of the environment
        void update_environment();

        // calculate fecundity of an individual
        // and update resource levels
        int calculate_fecundity(Individual &mother);

        double calculate_survival(unsigned int const patch_idx, Sex const the_sex);

        void calculate_patch_productivities();

        void add_juv_to_survivors(std::vector<Individual> &from,
                std::vector<Individual> &to);

        // clear the stack of juveniles and be ready for the
        // next season
        void clear_juveniles();

    public:
        TSDSeasonal(Parameters const &par);


};


#endif
