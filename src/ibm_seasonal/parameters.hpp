#ifndef _PARAMETERS_HPP_
#define _PARAMETERS_HPP_

class Parameters
{
    public:
        unsigned npatches{100};
        unsigned n[2]{5,5};

        double d[2]{.5,.5};

        double freq{0.05};

        // reproductive resources allocated when reproducing
        double resources{30};
        
        // timing value where t % timing == 0 determines breeding
        int init_timing{1};

        int max_t = 20;

        // initial threshold 
        double init_z{0.5};

        double mu_z{0.02};
        double mu_t{0.05};
        double unif_range_sdmu_t{2.0};
        double sdmu_z{0.01};
};

#endif 

