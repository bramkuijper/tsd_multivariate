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


        double freq{0.05};

        // reproductive resources allocated when reproducing
        double resources{30};
        
        // timing value where t % timing == 0 determines breeding
        int init_timing{1};

        // max threshold
        int max_threshold{20};
        unsigned max_time{100};

        // initial threshold 
        double init_z{0.5};

        double mu_z{0.02};
        double mu_t{0.05};
        double unif_range_sdmu_t{2.0};
        double sdmu_z{0.01};

        std::string file_name{"sim_tsd_seasonal"};
};

#endif 

