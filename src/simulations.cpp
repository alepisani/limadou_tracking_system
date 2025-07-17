#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <iomanip>
#include <TApplication.h>
#include <TCanvas.h>
#include <TView.h>
#include <TList.h>
#include <TPolyLine3D.h>
#include "TMarker3DBox.h"
#include "TH1F.h"
#include <numeric>
#include <TStopwatch.h>
#include "../include/eventdata.h"
#include "../include/chip.h"
#include "../include/display.h"
#include "../include/stats.h"
#include "../include/LTrackerCluster.h"
#include "../include/LTrackerTrack.h"
#include "../include/simulations.h"
#include "simulations.h"


simulations::simulations(){
    simulations::gen_tracks = {2,3,4,5,6,7,8,9,10,12,14,16,18,20,25,30,35,40,50};
    //simulations::gen_tracks = {2,3,4,5,6,7,8,9,10};
    
}


double simulations::mean(const vector<double> &v) {
    if (v.empty()) return 0.0;
    return accumulate(v.begin(), v.end(), 0.0) / v.size();
}


void simulations::sim(int iteration_per_event){
/* 
    fai una barra di caricamento
voglio che faccia X iterationi dello stesso processo
es. X iterazioni di 2/3/4/5/6/7/8/9/10/12/14/16/18/20/25/30/40/50 generated gen_tracks
    per ogni i-esima iterazione mi deve estrapolare dei dati e metterli in un vector
    alla fine di questa iterazione fai la media, cancella i singoli valori e conserva la media nei vector */

    display simu;
    LTrackerTrack ltt;
    stats stats;
    std::vector<double> r_time;
    std::vector<double> c_time;
    std::vector<double> reco;
    std::vector<double> gen_tr3L;
    std::vector<double> trkl;

    simu.take_angle_distribution();
    


    for(int i=0; i < gen_tracks.size(); ++i){
        cout << "~";
        r_time.clear();
        c_time.clear();
        reco.clear();
        gen_tr3L.clear();
        trkl.clear();
        stats.reset();
        ltt.Reset();
        
        for(int j=0; j<iteration_per_event; ++j){

            stats.reset();
            ltt.Reset();
            simu.tracks_no_print_hist(gen_tracks[i], ltt);
            //simu.tracks(gen_tracks[i], ltt, real_tracks);
            TStopwatch t;
            t.Start();
            ltt.computeTracklets();
            ltt.new_computing();
            //ltt.printRecoTracks_new_alg(real_tracks, gen_tracks[i]);
            t.Stop();

            r_time.push_back(t.RealTime());
            c_time.push_back(t.CpuTime());
            reco.push_back(stats::hmrt);
            gen_tr3L.push_back(stats::hmgthL012);
            trkl.push_back(ltt.tracklet_lay02.size());
            cout << "+-";

        }
        real_time.push_back(mean(r_time));
        cpu_time.push_back(mean(c_time));
        reco_trk.push_back(mean(reco));
        hmgth3l.push_back(mean(gen_tr3L));
        tracklet.push_back(mean(trkl));



    }

    cout << *this << endl;
    //cout << real_time.size() << endl;
}



std::ostream& operator<<(std::ostream& os, const simulations& sim) {
    size_t n = sim.real_time.size(); // supponiamo tutti i vector abbiano stessa size
    os << "======== Simulation Results ========\n";
    for (size_t i = 0; i < n; ++i) {
        os << "Tracks:  "       << sim.gen_tracks[i]                     << " | "
           << "Real Time: "     << std::fixed << std::setprecision(5) << sim.real_time[i]    << " s | "
           << "CPU Time: "      << std::fixed << std::setprecision(5) << sim.cpu_time[i]     << " s | "
           << "Reco Trk: "      << sim.reco_trk[i]                     << " | "
           << "Hmgth3l: "       << sim.hmgth3l[i]                      << " | "
           << "Tracklet: "      << sim.tracklet[i]
           << "\n";
    }

    return os;
}
