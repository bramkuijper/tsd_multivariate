#ifndef TSD_MULTIVARIATE_HPP
#include <string>
#include <fstream>

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
    double n[2];
    double eul_d;
    double eul_sr;
    double eul_b;
    bool delta_surv;
    std::string base;
};

class TSD_Multivariate
{
    private:
        // sex specific survival probabilities
        double surv_t0[2][2]; // initial survival probabilities
        double surv_tend[2][2]; // initial survival probabilities
        double surv[2][2]; // current survival probabilities
        

        // population sizes per patch
        double n[2];

        // sex-specific dispersal
        double d[2];

        // burrow depth
        double b;

        // proportions sons in either envt
        double s[2];

        // environmental change
        double sigma_t0[2][2]; // initial sigma
        double sigma_tend[2][2]; // new sigma
        double sigma[2][2]; // current sigma

        bool delta_surv;

        // frequency of either envt
        double p[2];

        // eigenvalue
        double lambda;

        // euler's constant
        double eul_d;
        double eul_sr;
        double eul_b;

        // reproductive values
        // first index sex
        // second index environment
        double v[2][2];

        // stable class frequencies
        // first index sex
        // second index environment
        double u[2][2];

        // relatedness coefficients between focal adult
        // and random offspring in local patch
        // first index is sex of focal adult
        // second index is environment of focal adult
        double rj[2][2];

        // coefficients of consanguinity between two randomly
        // sampled focal adults
        // index is environment
        double Qff[2];
        double Qfm[2];
        double Qmm[2];

        // the basename of the output file
        std::string base;

        static constexpr long int max_time = 1e08;
        static constexpr double vanish_bound = 1e-07;

        // private functions

        // calculate eigenvectors
        void eigenvectors(bool const output);

        // total number of competing offspring in the local site
        double C(Sex const sex_t1, bool const envt_t1);

        // derivative of the local competition function wrt to 
        // the mutant dispersal trait d[sex_d]
        double dCddlocal(Sex const sex_t1, 
                bool const envt_t1, 
                Sex const sex_d);
        
        // derivative of the local competition function wrt to 
        // the mutant sex ratio trait s[envt_d]
        double dCdsrlocal(Sex const sex_t1, 
                bool const envt_t1, 
                bool const envt_d);
        
        // derivative of the local competition function wrt to 
        // the burrowing trait b
        double dCdblocal(Sex const sex_t1, 
                bool const envt_t1);

        // total fecundity * survival of offsprign born in envt_t1 and
        // of sex sex_t1
        double fecundity_survival(bool const envt_t1, Sex const sex_t1);

        // derivative of fecundity functin wrt to trait s[envt_dsr]
        double dfecundity_survival_ds(
                bool const envt_t1
                ,Sex const sex_t1
                ,bool const envt_dsr
                );

        // derivative of fecundity functin wrt to trait b
        double dfecundity_survival_db(bool const envt_t1,Sex const sex_t1);

        void write_data(std::ofstream &output_file, int const generation);

        void write_parameters(std::ofstream &output_file);

        void write_data_headers(std::ofstream &output_file);

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

        double dB_dsr(
                bool const envt_t
                ,Sex const sex_t
                ,bool const envt_t1
                ,Sex const sex_t1
                ,bool const envt_dsr);

        double dB_db(
                bool const envt_t
                ,Sex const sex_t
                ,bool const envt_t1
                ,Sex const sex_t1);
        
        double dB_dd(
                bool const envt_t
                ,Sex const sex_t
                ,bool const envt_t1
                ,Sex const sex_t1
                ,Sex const sex_d);

        bool update_traits();

        void iterate_relatedness();
};

#endif
