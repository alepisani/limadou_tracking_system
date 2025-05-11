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
bool stats::hitL0;
bool stats::hitL1;
bool stats::hitL2;
double stats::hmgth1L = 0;
double stats::hmgth2L = 0;
double stats::hmgth0L = 0;
double stats::hmgthdcL2 = 0;
double stats::hmgthdcL1 = 0;
double stats::hmgthdcL0 = 0;



// Constructor
stats::stats(int events){
    hmgt = events;
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
    output << "how many generated tracks hitted only 0 layer: " << stats::hmgth0L << endl;
    output << "how many generated tracks hitted dead chip on layer 2: " << stats::hmgthdcL2 << endl;
    output << "how many generated tracks hitted dead chip on layer 1: " << stats::hmgthdcL1 << endl;
    output << "how many generated tracks hitted dead chip on layer 0: " << stats::hmgthdcL0 << endl;
    output << "how many reco tracks: " << stats::hmrt << endl;
    return output;
}


