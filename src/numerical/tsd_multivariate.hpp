#ifndef TSD_MULTIVARIATE_HPP
#include <string>

enum Sex {
    Male = 0,
    Female = 1
};

// struct with parameters
struct parstruct{

    double surv[2][2];
    double d[2];
    double b;
    double s[2];
    double sigma[2][2];
    double p[2];
    double v[2][2];
    double u[2][2];
    std::string base;
};

class TSD_Multivariate
{
    private:
        // sex specific survival probabilities
        double surv[2][2];

        // sex-specific dispersal
        double d[2];

        // burrow depth
        double b;

        // proportions sons in either envt
        double s[2];

        // environmental change
        double sigma[2][2];

        // frequency of either envt
        double p[2];

        // eigenvalue
        double lambda;

        // reproductive values
        double v[2][2];

        // stable class frequencies
        double u[2][2];

        std::string base;

        static constexpr long int max_time = 1e08;
        static constexpr double vanish_bound = 1e-08;

        // private functions

        // calculate eigenvectors
        void eigenvectors(bool const output);

        // total number of competing offspring in the local site
        double C(Sex const sex_t1, bool const envt_t1);

        // total fecundity * survival of offsprign born in envt_t1 and
        // of sex sex_t1
        double fecundity_survival(bool const envt_t1, Sex const sex_t1);

    public: 
        
        TSD_Multivariate();
        
        TSD_Multivariate(parstruct const &parstruct);

        // run the simulation
        void run();

        // initialize using command line arguments
        void init_arguments(int argc, char **argv);
        
        // element a((envt_t,sex_t) -> (envt_t1,sex_t1)) of the resident
        // transition matrix
        double A_resident(
                bool const envt_t
                ,Sex const sex_t
                ,bool const envt_t1
                ,Sex const sex_t1);
};

#endif
