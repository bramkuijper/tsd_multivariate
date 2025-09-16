#ifndef _PARAMETERS_HPP_
#define _PARAMETERS_HPP_

#include <string>

enum Sex
{
    female = 0,
    male = 1
};

class Parameters
{
    public:
        unsigned npatches{250};
        unsigned n{10};

        double d[2]{.5,.5};
        double survival_prob[2]{0.1,0.1};

        double t_opt[2]{-0.5,0.5};
        double omega_t[2]{0.1,0.1};

        int fecundity{10};

        // temperature sinusoidal
        double amplitude{1.0};
        double temperature_intercept{0.0};
        double temperature_intercept_change{0.8};

        // max number of timesteps in one season
        long unsigned max_t_season{50};

        // timing value, where t is time in the season of reproduction
        int init_t{3};

        // maximum duration of the simulation
        long unsigned max_simulation_time{100};

        // interval surrounding the perturbation
        // during which all data is plotted at every time step
        long unsigned interval{300};

        // initial sex allocation threshold 
        double init_a{0.0};
        double init_b{0.0};

        double mu_a{0.02};
        double mu_b{0.02};
        double ab_range[2]{-10,10};
        
        double mu_t{0.05};

        double unif_range_sdmu_t{2.0};
        double sdmu{0.05};

        unsigned int skip_output{1};

        double temp_error_sd{0.1};
        double cue_error{0.1};

        std::string file_name{"sim_tsd_seasonal"};

};

#endif 

