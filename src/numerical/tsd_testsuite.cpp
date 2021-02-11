#define BOOST_TEST_MODULE testsuite
#include "tsd_multivariate.hpp"
#include <boost/test/unit_test.hpp>
#include <string>

struct TestStructure
{
    parstruct parameters;

    TestStructure() 
    {
        parameters.surv[0][0] = 0.5;
        parameters.surv[1][0] = 0.5;
        parameters.surv[0][1] = 0.5;
        parameters.surv[1][1] = 0.5;
        parameters.d[0] = 0.2;
        parameters.d[1] = 0.8;
        parameters.b = 0.25;
        parameters.s[0] = 0.25;
        parameters.s[1] = 0.25;
        parameters.p[0] = 0.5;
        parameters.p[1] = 0.5;
        parameters.sigma[0][0] = 0.5;
        parameters.sigma[1][0] = 0.5;
        parameters.sigma[0][1] = 0.5;
        parameters.sigma[1][1] = 0.5;
        parameters.v[0][0] = 1.0;
        parameters.v[1][0] = 1.0;
        parameters.v[0][1] = 1.0;
        parameters.v[1][1] = 1.0;
        parameters.u[0][0] = 0.25;
        parameters.u[1][0] = 0.25;
        parameters.u[0][1] = 0.25;
        parameters.u[1][1] = 0.25;
//        double burying_depth = 0.3;
//
//        double s01 = 0.3;
//        double s10 = 0.8;
//
//        double s1 = 0.8;
//        double s2 = 0.3;
//
//        double p0 = s10 / (s10 + s01);
//
//        double d[2] = {dm,df};
//        double s[2] = {s1,s2};
//        double sigma[2][2] = {{1.0 - s01,s10},{s01,1.0 - s10}};
//        double p[2] = {p0,1.0 - p0};
//        double v[2][2] = {{1.1,0.9},{0.8,1.2}};
//        double u[2][2] = {{0.3,0.2},{0.25,0.25}};

        parameters.base = "test_output";
    } // end TestStructure()

    void TestARes()
    {
        TSD_Multivariate testsim(parameters);

        std::cout << testsim.A_resident(0,static_cast<Sex>(0),0,static_cast<Sex>(0)) << std::endl;
    }

    // destruct - remove files etc
    ~TestStructure() 
    {
    } // end ~TestStructure
};

BOOST_FIXTURE_TEST_SUITE(test_case1, TestStructure)

BOOST_AUTO_TEST_CASE(test_case2)
{
    TestARes();
}

BOOST_AUTO_TEST_SUITE_END()
