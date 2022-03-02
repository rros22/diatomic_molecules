#include "rand.hpp"
#include <iostream>

#include <random>
#include <ctime>

double rand_range(double min, double max){

    static std::default_random_engine generator(unsigned(time(nullptr)));
    std::uniform_real_distribution<double> distribution(min, max);

    return distribution(generator);

}
