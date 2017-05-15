#ifndef RANDOM_H
#define RANDOM_H

#include <boost/random/mersenne_twister.hpp>

class Random {
    private:

        boost::mt19937 generator;

    public:
        // constructor
        Random();

        void seed();//

        void seed(unsigned int seed);

        // continuous uniform random number in (a,b)
        double continuousUniform(double a, double b);

        // exponential random number
        double exponential(double rate);

        // poisson random number
        double poisson(double mean);
    
};//end class Random

#endif
