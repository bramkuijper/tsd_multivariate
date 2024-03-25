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
        std::uniform_int_distribution<int> patch_sampler;
        std::uniform_int_distribution<int> female_sampler;
        std::uniform_int_distribution<int> male_sampler;

        long unsigned time_step{0};

        double temperature{0.0};

        std::random_device rd;
        unsigned int seed;
        std::mt19937 rng_r;

        // all the patches of the metapopoulation
        std::vector<Patch> metapopulation;

        // total productivity for males and females
        int global_productivity[2]{0,0};

        void survive();
        void reproduce();
        void replace();

        void write_data_headers();
        void write_data();
        void write_parameters();

        // update the state of the environment
        void update_environment();

        // calculate fecundity of an individual
        // and update resource levels
        int calculate_fecundity(Individual &mother);

        bool calculate_survival(Sex const the_sex);

    public:
        TSDSeasonal(Parameters const &par);


};


#endif
