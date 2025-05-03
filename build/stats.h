#ifndef STATS_H
#define STATS_H

#include <iostream>
using namespace std;

class stats {
public:
    static double hmgt;           // how many generated tracks (generated TR2)
    static double hmgthTR1;       // how many generated tracks hitted TR1
    static double hmgthL2;        // how many generated tracks hitted L2
    static double hmgthL1;        // how many generated tracks hitted L1
    static double hmgthL0;        // how many generated tracks hitted L0
    static double hmrt;           // how many reco tracks
    static double hmgthL012;      // how many generated tracks hitted all 3 layers

    stats(); // Default constructor

    friend std::ostream &operator<<(std::ostream &output, const stats &s);
};

#endif