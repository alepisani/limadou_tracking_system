#ifndef SIMULATIONS_H
#define SIMULATIONS_H

#include <string>
#include <array> 
#include <map>
#include <vector>
#include <iostream>
#include <cmath>
#include <chrono>
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
    //double radius;
    std::vector<double> real_time;
    std::vector<double> cpu_time;
    std::vector<double> reco_trk;
    std::vector<double> reco_trk_real;
    std::vector<double> hmgth3l;
    std::vector<double> reco_hmgth3l;
    std::vector<double> tracklet;
    std::vector<double> radius;
    std::vector<double> chi2cut;
    std::vector<double> delta_gentrk_recoreal;
    std::vector<double> delta_gentrk_reco;
    std::vector<double> efficiency;

    simulations();
    static double mean(const vector<double>& v);
    void sim_only_trk_3L(int);
    void sim_old_algo(int);
    void sim_trk_32L(int);
    void printProgressBarWithETA(int, int, std::chrono::steady_clock::time_point, int);
      

};




#endif