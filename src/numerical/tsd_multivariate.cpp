#include <cstdlib>
#include <cmath>
#include <fstream>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>


#include "tsd_multivariate.hpp"

TSD_Multivariate::TSD_Multivariate() :
    surv{{0.0,0.0},{0.0,0.0}}
    ,d{0.0,0.0}
    ,b{0.0}
    ,s{0.0,0.0}
    ,sigma{{0.0,0.0},{0.0,0.0}}
    ,p{0.0,0.0}
    ,lambda{0.0}
    ,v{{0.0,0.0},{0.0,0.0}}
    ,u{{0.0,0.0},{0.0,0.0}}
    ,delta_surv{false}
    ,base{}
{
} // end void TSD_Multivariate::TSD_Multivariate

// constructor with arguments
TSD_Multivariate::TSD_Multivariate(
        parstruct const &pstruct) :
            surv{{pstruct.surv[0][0],
                    pstruct.surv[0][1]},
                {pstruct.surv[1][0],
                    pstruct.surv[1][1]}}
            ,d{pstruct.d[0],d[1]}
            ,b{pstruct.b}
            ,s{pstruct.s[0],pstruct.s[1]}
            ,sigma{
                {pstruct.sigma[0][0]
                    ,pstruct.sigma[0][1]},
                {pstruct.sigma[1][0]
                    ,pstruct.sigma[1][1]}}
        ,p{pstruct.p[0],pstruct.p[1]}
        ,lambda{1.0}
        ,v{
            {pstruct.v[0][0],pstruct.v[0][1]}
            ,{pstruct.v[1][0],pstruct.v[1][1]}}
        ,u{
            {pstruct.u[0][0],pstruct.u[0][1]}
            ,{pstruct.u[1][0],pstruct.u[1][1]}}
        ,delta_surv{pstruct.delta_surv}
        ,base{pstruct.base}
{
}

// initialize arguments
void TSD_Multivariate::init_arguments(int argc, char **argv)
{
    surv[Male][0] = std::atof(argv[1]);
    surv[Male][1] = std::atof(argv[2]);
    surv[Female][0] = std::atof(argv[3]);
    surv[Female][1] = std::atof(argv[4]);
    d[Male] = std::atof(argv[5]);
    d[Female] = std::atof(argv[5]);
    b = std::atof(argv[5]);
    s[0] = std::atof(argv[5]);
    s[1] = std::atof(argv[5]);
    sigma[0][1] = std::atof(argv[5]);
    sigma[0][0] = 1.0 - sigma[0][1];
    sigma[1][0] = std::atof(argv[5]);
    sigma[1][0] = 1.0 - sigma[1][0];
    lambda = std::atof(argv[6]);
    delta_surv = std::atoi(argv[6]);
    base = argv[7];
    p[0] = sigma[1][0] / (sigma[1][0] + sigma[0][1]);
    p[1] = 1.0 - p[0];
} // end init_arguments()

// run the iteration
void TSD_Multivariate::run()
{
    double b_tplus1 = 0.0;
    double s_tplus1[2] = {0.0,0.0};

    base = base + ".csv";

    // initialize the file only when you run the thing
    std::ofstream output_file{base};

    write_parameters(output_file);


    for (long int time_step = 0; time_step < max_time; ++time_step)
    {
        // update the eigenvectors during every time step
        eigenvectors(true);
    }
} // end run()

void TSD_Multivariate::write_data_headers(std::ofstream &output_file)
{
    output_file << "generation;b;lambda;";

    for (int iter_i = 0; iter_i < 2; ++iter_i)
    {
        output_file << "d" << (static_cast<Sex>(iter_i) == Male ? "m" : "f") << ";";
        output_file << "sr" << iter_i + 1 << ";";

        for (int iter_j = 0; iter_j < 2; ++iter_j)
        {
            output_file << "v" << (static_cast<Sex>(iter_i) == Male ? "m" : "f") 
                << iter_j + 1 << ";"; 

            output_file << "u" << (static_cast<Sex>(iter_i) == Male ? "m" : "f") << iter_j + 1
                << iter_j + 1 << ";";
        }
    }
}

void TSD_Multivariate::write_data(std::ofstream &output_file, int const generation)
{
    output_file << generation << ";" << b << ";" << lambda << ";";

    for (int iter_i = 0; iter_i < 2; ++iter_i)
    {
        output_file << d[iter_i] << ";";
        output_file << s[iter_i] << ";";

        for (int iter_j = 0; iter_j < 2; ++iter_j)
        {
            output_file << v[iter_i][iter_j] << ";" << u[iter_i][iter_j] << ";";
        }
    }
}

void TSD_Multivariate::write_parameters(std::ofstream &output_file)
{
    output_file << std::endl
        << std::endl
        << "p1;" << p[0] << std::endl;

    for (int iter_i = 0; iter_i < 2; ++iter_i)
    {
        for (int iter_j = 0; iter_j < 2; ++iter_j)
        {
            output_file << "surv" 
                << (static_cast<Sex>(iter_i) == Male ? "m" : "f") 
                << (iter_j + 1) 
                << ";"
                << surv[iter_i][iter_j] 
                << std::endl;
            
            output_file << "sigma" 
                << (iter_i + 1) 
                << (iter_j + 1) 
                << ";"
                << sigma[iter_i][iter_j]
                << std::endl;
        }
    }
}

// total fecundity * survival of offsprign born in envt_t1 and
// of sex sex_t1
double TSD_Multivariate::fecundity_survival(bool const envt_t1, Sex const sex_t1) 
{
    if (sex_t1 == Male) 
    {
        return(envt_t1 == 0 ? 
                s[envt_t1] * surv[sex_t1][envt_t1]
                :
                (1.0 - b) * s[envt_t1] * surv[sex_t1][envt_t1] + 
                    b * s[!envt_t1] * (delta_surv * surv[sex_t1][!envt_t1] 
                                        + (1 - delta_surv) * surv[sex_t1][envt_t1]));
    }

    return(envt_t1 == 0 ?
            (1.0 - s[envt_t1]) * surv[sex_t1][envt_t1]
            :
            (1.0 - b) * (1.0 - s[envt_t1]) * surv[sex_t1][envt_t1] + 
                    b * (1.0 - s[!envt_t1]) * (delta_surv * surv[sex_t1][!envt_t1] 
                                        + (1 - delta_surv) * surv[sex_t1][envt_t1]));
            
}

// total number of competing offspring
double TSD_Multivariate::C(Sex const sex_t1, bool const envt_t1)
{
    return((1.0 - d[sex_t1]) * fecundity_survival(envt_t1, sex_t1) 
            + d[sex_t1] * 
                    (p[0] * fecundity_survival(0, sex_t1)
                     + (1.0 - p[1]) * fecundity_survival(1, sex_t1)));
}

// transition matrix for a rare sex allocation mutant
double TSD_Multivariate::Bsr(
        bool const envt_t
        ,Sex const sex_t
        ,bool const envt_t1
        ,Sex const sex_t1)
{
    double bij = sigma[envt_t][envt_t1] * (1.0 - d[sex_t]) *
        dfecundity_sur
}


// get entries of the resident transition matrix
double TSD_Multivariate::A_resident(
        bool const envt_t
        ,Sex const sex_t
        ,bool const envt_t1
        ,Sex const sex_t1)
{
    // philopatric fitness component
    double aij = sigma[envt_t][envt_t1] * (1.0 - d[sex_t]) * 
                                fecundity_survival(envt_t1, sex_t1) / C(sex_t1, envt_t1)
               // remote fitness component 
                    + (
                        p[envt_t1] * (1.0 - sigma[envt_t1][!envt_t1]) 
                        + 
                        p[!envt_t1] * sigma[!envt_t1][envt_t1]
                    ) * d[sex_t] *
                                fecundity_survival(envt_t1, sex_t1) / C(sex_t1, envt_t1);

    return(aij);                                        
}

// calculate left and right eigenvectors
void TSD_Multivariate::eigenvectors(bool const output)
{
    size_t ndim = 4;

    // allocate GSL things to work on the eigenvectors
    //
    // allocate a vector to calculate eigenvalues
    gsl_vector_complex *eval = gsl_vector_complex_alloc(ndim);

    // allocate a matrix to calculate eigenvectors
    gsl_matrix_complex *evec = gsl_matrix_complex_alloc(ndim, ndim);
     
    // allocate a vector to calculate eigenvalues for the tranpose
    gsl_vector_complex *evalT = gsl_vector_complex_alloc(ndim);
    
    // allocate a matrix to calculate eigenvectors for the tranpose
    gsl_matrix_complex *evecT = gsl_matrix_complex_alloc(ndim, ndim);

    gsl_matrix *m = gsl_matrix_alloc(ndim,ndim);
    gsl_matrix *mT = gsl_matrix_alloc(ndim,ndim);
    gsl_matrix *m_stats = gsl_matrix_alloc(ndim,ndim);
     
    gsl_eigen_nonsymmv_workspace * workspace_w = gsl_eigen_nonsymmv_alloc(ndim);
    gsl_eigen_nonsymmv_workspace * workspace_wT = gsl_eigen_nonsymmv_alloc(ndim);

    // auxiliary variable to see whether stuff has been converged
    bool converged;

    // auxiliary variables to store the lhs
    double v_i_tplus1, u_i_tplus1, ev_tplus1;

    double Atmp;

    size_t row_i, col_j;
        
    // now calculate solutions of the eigenvectors
    // until convergence
    for (long int iter_time_step = 0; 
            iter_time_step < max_time;
            ++iter_time_step)
    {
        row_i = 0;
        col_j = 0;

        for (size_t envt_i = 0; envt_i < 2; ++envt_i)
        {
            for (size_t sex_i = 0; sex_i < 2; ++sex_i)
            {
                for (size_t envt_j = 0; envt_j < 2; ++envt_j)
                {
                    for (size_t sex_j = 0; sex_j < 2; ++sex_j)
                    {
                        Atmp = A_resident(
                                envt_i
                                ,static_cast<Sex>(sex_i)
                                ,envt_j
                                ,static_cast<Sex>(sex_j));

                        gsl_matrix_set(m, row_i, col_j, Atmp);
                        gsl_matrix_set(mT, col_j, row_i, Atmp);
                        gsl_matrix_set(m_stats, row_i, col_j, Atmp);

                        // increment column after 
                        // each (sex_i, envt_i) combination
                        ++col_j;

                        // if maximum column number has been reached
                        // reset to 0
                        if (col_j > 3)
                        {
                            col_j = 0;
                        }
                        
                    }// end for sex_j
                } // end for envt_j

                // increment row after each (sex_i, envt_i) combination
                ++row_i;

                // if max row number has been reached
                // reset to 0
                if (row_i > 3)
                {
                    row_i = 0;
                }
            } // end for envt_i
        } //end for sex_i
        
        // now calculate evs
        // set up the workspaces
        gsl_eigen_nonsymmv(m, eval, evec, workspace_w);
        gsl_eigen_nonsymmv(mT, evalT, evecT, workspace_wT);
                
        

        gsl_eigen_nonsymmv_sort (eval, evec, 
                                GSL_EIGEN_SORT_ABS_DESC);
            
        // get the first column from the eigenvector matrix
        // giving us the right eigenvector
        gsl_vector_complex_view evec_i 
                   = gsl_matrix_complex_column(evec, 0);

        // get the first column from the eigenvector matrix
        // of the transpose of the left eigenvector
        // giving us the left eigenvector
        gsl_eigen_nonsymmv_sort (evalT, evecT, 
                                GSL_EIGEN_SORT_ABS_DESC);

        // get the last column from the eigenvector
        gsl_vector_complex_view evecT_i 
                   = gsl_matrix_complex_column(evecT, 0);

        double product_vu = 0;
        double total_u = 0;

        // make variables that normalize both eigenvectors
        for (size_t iter1 = 0; iter1 < ndim; ++iter1)
        {
            // norm for right ev
            total_u += fabs(GSL_REAL(
                   gsl_vector_complex_get(&evec_i.vector,iter1)));

            // norm for left ev
            product_vu += fabs(GSL_REAL(
                   gsl_vector_complex_get(&evec_i.vector,iter1)))
                   * fabs(GSL_REAL(
                   gsl_vector_complex_get(&evecT_i.vector,iter1)));
        }


        product_vu /= total_u;
       
        // by default we assume convergence until proved otherwise
        converged = true;

        bool envt_idx;
        Sex sex_idx;

        // get the results from the eigenvector
        for (size_t dim_i = 0; dim_i < ndim; ++dim_i)
        {
            // get the right eigenvalue out of the eigenvector of A
            u_i_tplus1 = fabs(GSL_REAL(
                        gsl_vector_complex_get(&evec_i.vector,dim_i)))/
                total_u;
            
            // get the left eigenvalue out of the eigenvector of A^T
            v_i_tplus1 = fabs(GSL_REAL(
                        gsl_vector_complex_get(&evecT_i.vector,dim_i)))/
                product_vu;

            // translate dimensions into
            // envt sex, with
            // 0: envt = 0, sex = 0 
            sex_idx = static_cast<Sex>(dim_i % 2);
            envt_idx = static_cast<Sex>(std::floor(dim_i/2));

            if (fabs(u[sex_idx][envt_idx] - u_i_tplus1) > vanish_bound)
            {
                converged = false;
            }

            if (fabs(v[sex_idx][envt_idx] - v_i_tplus1) > vanish_bound)
            {
                converged = false;
            }

            v[sex_idx][envt_idx] = v_i_tplus1;
            u[sex_idx][envt_idx] = u_i_tplus1;

            assert(std::isnan(v[sex_idx][envt_idx]) == 0);
            assert(std::isnan(u[sex_idx][envt_idx]) == 0);
        } // for dim_i
        
        // now get the eigenvalue
        ev_tplus1 = GSL_REAL(gsl_vector_complex_get(eval, 0));

        if (fabs(lambda - ev_tplus1) > vanish_bound)
        {
            converged = false;
        }

        lambda = ev_tplus1;

        if (converged)
        {
            break;
        }
    } // end iteration


    // free all the stuff we've used
    gsl_eigen_nonsymmv_free (workspace_w);
    gsl_eigen_nonsymmv_free (workspace_wT);
    gsl_vector_complex_free(evalT);
    gsl_matrix_complex_free(evecT);
    gsl_matrix_free(mT);
    gsl_matrix_free(m);
    gsl_matrix_free(m_stats);

}

