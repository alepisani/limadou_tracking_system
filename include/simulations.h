#ifndef SIMULATIONS_H
#define SIMULATIONS_H

#include <string>
#include <array> 
#include <map>
#include <vector>
#include <iostream>
#include <cmath>
#include <TCanvas.h> // Include ROOT's TCanvas header
#include <TView.h>   // Include ROOT's TView header
#include "TMath.h"
#include "./LTrackerTrack.h"
#include "../include/eventdata.h"
#include "../include/chip.h"
#include "../include/display.h"
#include "../include/stats.h"
#include "../include/LTrackerCluster.h"
using namespace std;

class simulations{
public:
    //classe per stipare i dati dalle simulazioni

    std::vector<int> gen_tracks;
    std::vector<double> real_time;
    std::vector<double> cpu_time;
    std::vector<int> reco_trk;
    std::vector<int> hmgth3l;
    std::vector<double> reco_hmgth3l;
    std::vector<int> tracklet;

    simulations();
    static double mean(const vector<double>& v);
    void sim(int);
    friend std::ostream& operator<<(std::ostream& os, const simulations& sim);

};




#endif