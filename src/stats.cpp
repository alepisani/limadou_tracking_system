#include <string>
#include <array> 
#include <map>
#include <iostream>
#include <cmath>
#include <TApplication.h>
#include <TCanvas.h>
#include <TView.h>
#include <TList.h>
#include <TPolyLine3D.h>
#include "TH1F.h"
#include "TLine.h"
#include "TRandom3.h"
#include "TMarker3DBox.h"
#include "../include/stats.h"
using namespace std;

double stats::hmgt = 0;
double stats::hmgthTR1 = 0;
double stats::hmgthL2 = 0;
double stats::hmgthL1 = 0;
double stats::hmgthL0 = 0;
double stats::hmgthL012 = 0;
double stats::hmrt = 0;
int stats::hmrtar = 0;
bool stats::hitL0;
bool stats::hitL1;
bool stats::hitL2;
double stats::hmgth1L = 0;
double stats::hmgth2L = 0;
double stats::hmgth0L = 0;
double stats::hmgthdcL2 = 0;
double stats::hmgthdcL1 = 0;
double stats::hmgthdcL0 = 0;
double stats::fakehit = 0;
double stats::hmthL2 = 0;
double stats::hmthL1 = 0;
double stats::hmthL0 = 0;
double stats::hmbh3L = 0;
double stats::hmbh2L = 0;
double stats::hmbh1L = 0;
double stats::hmbh0L = 0;



// Constructor
stats::stats(){
}


    void stats::reset(){
        stats::hmgt = 0;
        stats::hmgthTR1 = 0;
        stats::hmgthL2 = 0;
        stats::hmgthL1 = 0;
        stats::hmgthL0 = 0;
        stats::hmgthL012 = 0;
        stats::hmrt = 0;
        stats::hmrtar = 0;
        stats::hitL0;
        stats::hitL1;
        stats::hitL2;
        stats::hmgth1L = 0;
        stats::hmgth2L = 0;
        stats::hmgth0L = 0;
        stats::hmgthdcL2 = 0;
        stats::hmgthdcL1 = 0;
        stats::hmgthdcL0 = 0;
        stats::fakehit = 0;
        stats::hmthL2 = 0;
        stats::hmthL1 = 0;
        stats::hmthL0 = 0;
        stats::hmbh3L = 0;
        stats::hmbh2L = 0;
        stats::hmbh1L = 0;
        stats::hmbh0L = 0;
    }

// Overload operator<<
std::ostream &operator<<(std::ostream &output, const stats &s) {
    output << "how many generated tracks hitted TR2: " << stats::hmgt << endl;
    output << "how many generated tracks hitted TR1: " << stats::hmgthTR1 << endl;
    output << "how many generated tracks hitted L2: " << stats::hmgthL2-stats::hmgthdcL2 << endl;
    output << "how many generated tracks hitted L1: " << stats::hmgthL1-stats::hmgthdcL1 << endl;
    output << "how many generated tracks hitted L0: " << stats::hmgthL0-stats::hmgthdcL0 << endl;
    output << "how many generated tracks hitted all  3 layer: " << stats::hmgthL012 << endl;
    output << "how many generated tracks hitted only 2 layer: " << stats::hmgth2L << endl;
    output << "how many generated tracks hitted only 1 layer: " << stats::hmgth1L << endl;
    output << "how many generated tracks hitted      0 layer: " << stats::hmgth0L << endl;
    output << "how many generated tracks hitted dead chip on layer 2: " << stats::hmgthdcL2 << endl;
    output << "how many generated tracks hitted dead chip on layer 1: " << stats::hmgthdcL1 << endl;
    output << "how many generated tracks hitted dead chip on layer 0: " << stats::hmgthdcL0 << endl;
    output << "how many fake hit: " << stats::fakehit << endl;
    output << "how many reco tracks: " << stats::hmrt << endl;
    output << "how many reco tracks are real: " << stats::hmrtar << endl;
    output << "how many cluster (beam test) on layer 2: " << stats::hmthL2 << endl;
    output << "how many cluster (beam test) on layer 1: " << stats::hmthL1 << endl;
    output << "how many cluster (beam test) on layer 0: " << stats::hmthL0 << endl;
    output << "how many beam hitted all  3 layer: " << stats::hmbh3L << endl;
    output << "how many beam hitted only 2 layer: " << stats::hmbh2L << endl;
    output << "how many beam hitted only 1 layer: " << stats::hmbh1L << endl;
    output << "how many beam hitted      0 layer: " << stats::hmbh0L << endl;
    
    return output;
}




