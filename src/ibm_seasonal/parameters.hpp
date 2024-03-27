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

        double t_opt[2]{-0.1,0.0};
        double omega_t[2]{0.1,0.1};

        // temperature sinusoidal
        double amplitude{1.0};
        double frequency{1.0};

        // timing value where t % timing == 0 determines breeding
        int init_t{1};
        int max_t{20};

        double init_resources{10};

        double init_effort_per_timestep{0.1};

        long unsigned max_simulation_time{100};

        // initial sex allocation threshold 
        double init_a{0.0};
        double init_b{0.0};

        double mu_a{0.02};
        double mu_b{0.02};

        double mu_effort_per_timestep{0.02};
        double mu_t{0.05};

        double unif_range_sdmu_t{2.0};
        double sdmu{0.02};

        unsigned int skip_output = 100;

        std::string file_name{"sim_tsd_seasonal"};
};

#endif 

