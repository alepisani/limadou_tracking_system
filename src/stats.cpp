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
#include "../build/stats.h"
using namespace std;


// Constructor
stats::stats(){}

double stats::hmgt = 0;
double stats::hmgthTR1 = 0;
double stats::hmgthL2 = 0;
double stats::hmgthL1 = 0;
double stats::hmgthL0 = 0;
double stats::hmgthL012 = 0;
double stats::hmrt = 0;

// Overload operator<<
std::ostream &operator<<(std::ostream &output, const stats &s) {
    output << "how many generated tracks: " << stats::hmgt << endl;
    output << "how many generated tracks hitted TR1: " << stats::hmgthTR1 << endl;
    output << "how many generated tracks hitted L2: " << stats::hmgthL2 << endl;
    output << "how many generated tracks hitted L1: " << stats::hmgthL1 << endl;
    output << "how many generated tracks hitted L0: " << stats::hmgthL0 << endl;
    output << "how many generated tracks hitted each layer: " << stats::hmgthL012 << endl;
    output << "how many reco tracks: " << stats::hmrt << endl;
    return output;
}


