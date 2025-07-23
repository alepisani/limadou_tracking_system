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
    static int hmrtar;         // how many reco tracks are real
    static double hmgthL012;      // how many generated tracks hitted all 3 layers
    static bool hitL0;
    static bool hitL1;
    static bool hitL2;
    static double hmgth1L;        //how many generated tarcks hitted only 1 layer
    static double hmgth2L;        //how many generated tarcks hitted only 2 layer
    static double hmgth0L;        //how many generated tarcks hitted only 0 layer
    static double hmgthdcL2;      //how many generated tracks hitted dead chip on layer 2
    static double hmgthdcL1;      //how many generated tracks hitted dead chip on layer 1
    static double hmgthdcL0;      //how many generated tracks hitted dead chip on layer 0
    static double fakehit;        //how many fake hit
    static double hmthL2;         //how many tracks (beam test) hitted layer 2  
    static double hmthL1;         //how many tracks (beam test) hitted layer 1
    static double hmthL0;         //how many tracks (beam test) hitted layer 0  
    static double hmbh3L;         //how many beam hitted all 3 layer
    static double hmbh2L;         //how many beam hitted only 2 layer
    static double hmbh1L;         //how many beam hitted only 1 layer
    static double hmbh0L;         //how many beam hitted only 0 layer
    

    stats(); // Default constructor
    void reset();

    friend std::ostream &operator<<(std::ostream &output, const stats &s);
};

#endif