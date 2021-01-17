#include "tsd_turtle.hpp"
#include <boost/test/unit_test.hpp>
#include <string>

struct TestStructure
{
    TestStructure() {
        double surv_female_0 = 0.5;
        double surv_female_1 = 1.0;

        double df = 0.1;
        double dm = 0.8;

        double burying_depth = 0.3;

        double s01 = 0.3;
        double s10 = 0.8;

        double s1 = 0.8;
        double s2 = 0.3;

        double p0 = s10 / (s10 + s01);

        double surv[2][2] = {{1.0,1.0},{surv_female_0,surv_female_1}};
        double d[2] = {dm,df};
        double s[2] = {s1,s2};
        double sigma[2][2] = {{1.0 - s01,s10},{s01,1.0 - s10}};
        double p[2] = {p0,1.0 - p0};
        double v[2][2] = {{1.1,0.9},{0.8,1.2}};
        double u[2][2] = {{0.3,0.2},{0.25,0.25}};

        std::string base_str{"test_output"};
    } // end TestStructure()


    // destructure - remove files etc
    ~TestStructure() {
    } // end ~TestStructure


};

BOOST_FIXTURE_TEST_SUITE(
    TSDTurtle testsim{
            surv
            ,d
            ,burying_depth
            ,s
            ,sigma
            ,p
            ,v
            ,u
            ,base_str};
}
