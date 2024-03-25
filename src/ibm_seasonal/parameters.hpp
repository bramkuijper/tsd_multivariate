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
        unsigned npatches{100};
        unsigned n[2]{5,5};

        double d[2]{.5,.5};
        double survival_prob[2]{.5,.5};

        double t_opt[2]{-0.1,0.0};
        double omega_t[2]{0.1,0.1};

        double freq{0.05};

        // timing value where t % timing == 0 determines breeding
        int init_t{1};

        double init_resources{10};

        double init_effort_per_timestep{2.0};

        // max value of the sex allocation threshold
        int max_threshold{20};

        unsigned max_time{100};

        // initial sex allocation threshold 
        double init_z{0.5};

        double mu_z{0.02};
        double mu_effort_per_timestep{0.02};
        double mu_t{0.05};
        double unif_range_sdmu_t{2.0};
        double sdmu{0.02};

        // temperature sinusoidal
        double temperature_amplitude{1.0};
        double frequency{1.0};

        std::string file_name{"sim_tsd_seasonal"};
};

#endif 

