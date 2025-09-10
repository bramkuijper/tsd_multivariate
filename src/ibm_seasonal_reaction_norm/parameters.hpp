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

        // the optimal seasonal temperature to reproduce
        // (temperature varies as a sinusoidal, so starts
        // at 0 at the start of the season, then increases
        // then decreases)
        double t_opt[2]{-0.5,0.5};

        // the width of the seasonal function
        double omega_t[2]{0.1,0.1};

        int fecundity{30};

        // amplitude of the temperature sinusoidal
        double amplitude{1.0};
        double temperature_intercept{0.0};
        double temperature_intercept_change{0.8};

        // max number of timesteps in one season
        long unsigned max_t_season{50};

        // timing value, where t is time in the season of reproduction
        double init_t{2.0/50};
        double init_tb{0.0}; // initial slope of the reaction norm

        // the duration of the simulation
        long unsigned max_simulation_time{100};

        // initial sex allocation threshold 
        double init_a{0.0};
        double init_b{0.0};

        double mu_a{0.02};
        double mu_b{0.02};
        double ab_range[2]{-10,10};
        
        double mu_t{0.05};
        double mu_tb{0.05};

        double sdmu{0.05};

        unsigned int skip_output{1};

        double temp_error_sd{0.1};
        double cue_error{0.01};

        std::string file_name{"sim_tsd_seasonal"};

};

#endif 

