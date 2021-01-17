#ifndef INDIVIDUAL_HPP_
#define INDIVIDUAL_HPP_

class Individual
{
    public:
        // baseline ornament investment
        double t[2];

        // quality dependent ornament investment
        double tprime[2];

        // preference locus, p = 0 means random mating
        double p[2];

        // preference phenotype
        double p_phen;

        // total size of the ornament
        double s;

        // whether this male lives in a 
        // high quality environment
        bool envt_quality_high;

        // constructor (i.e., build this individual)
        Individual();

        // copy constructor (i.e., copy this individual from another one)
        Individual(Individual const &other);

        // change the assignment operator so that also individuals
        // can be assigned, i..e., Individual1 = Individual2;
        void operator=(Individual const &other);

}; // end class individual





#endif
