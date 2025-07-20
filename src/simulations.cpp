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
#include <chrono>
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
    //simulations::gen_tracks = {2,3,4,5,6,7,8,9,10,12,14,16,18,20,25,30};
    //simulations::gen_tracks = {2,3,4,5};
    //simulations::gen_tracks = {8};
    radius = 7.;
}


double simulations::mean(const vector<double> &v) {
    if (v.empty()) return 0.0;
    return accumulate(v.begin(), v.end(), 0.0) / v.size();
}

void simulations::printProgressBarWithETA(int current, int total, std::chrono::steady_clock::time_point start_time, int barWidth = 30) {
    using namespace std::chrono;

    float progress = static_cast<float>(current) / total;
    int pos = static_cast<int>(barWidth * progress);

    auto now = steady_clock::now();
    auto elapsed = duration_cast<seconds>(now - start_time).count();
    int eta = (progress > 0.0) ? static_cast<int>(elapsed / progress - elapsed) : 0;

    std::ostringstream oss;
    oss << "\r[";

    for (int i = 0; i < barWidth; ++i) {
        if (i < pos) oss << "█";
        else oss << " ";
    }

    oss << "] ";
    oss << std::setw(3) << int(progress * 100.0) << "% ";
    oss << "(i=" << current << "/" << total << ") ";
    oss << "⏱ " << elapsed << "s ";
    oss << "ETA: " << eta << "s";

    // Aggiungi spazi per sovrascrivere eventuali residui
    std::string output = oss.str();
    size_t terminal_width = 80;
    if (output.size() < terminal_width)
        output += std::string(terminal_width - output.size(), ' ');

    std::cout << output << std::flush;
}


void simulations::sim_only_trk_3L(int iteration_per_event){

    auto start_time = std::chrono::steady_clock::now();
    display simu;
    LTrackerTrack ltt;
    stats stats;
    chips cc;
    std::vector<double> r_time;
    std::vector<double> c_time;
    std::vector<double> reco;
    std::vector<double> gen_tr3L;
    std::vector<double> trkl;

    /* TCanvas* real_tracks = new TCanvas("MC_tracks", "3D View_mc", 800, 600);
    TView* rt = TView::CreateView(1);
    rt->SetRange(-100, -100, 0, 100, 100, 70);
    rt->ShowAxis();
    simu.draw_TR12(real_tracks);
    cc.print_all_chips(cc, real_tracks); */




    simu.take_distributions();
    
    for(int i=0; i < gen_tracks.size(); ++i){
        r_time.clear();
        c_time.clear();
        reco.clear();
        gen_tr3L.clear();
        trkl.clear();
        stats.reset();
        ltt.Reset();
        
        printProgressBarWithETA(i, gen_tracks.size(), start_time);
        
        for(int j=0; j<iteration_per_event; ++j){

            //cout << "i=" << i << ", j=" << j << std::endl << std::flush;

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

        }
        real_time.push_back(mean(r_time));
        cpu_time.push_back(mean(c_time));
        reco_trk.push_back(mean(reco));
        hmgth3l.push_back(mean(gen_tr3L));
        tracklet.push_back(mean(trkl));

    }
    cout << endl;
    cout << *this << endl;
    //cout << real_time.size() << endl;
}



std::ostream& operator<<(std::ostream& os, const simulations& sim) {
    // Header
    os << "\n=== Simulation Results ===\n";
    os << std::setw(10) << "Gen Trk"
       << std::setw(12) << "RealTime"
       << std::setw(12) << "CPUTime"
       << std::setw(12) << "RecoTrk"
       << std::setw(12) << "Gen3L"
       << std::setw(12) << "Tracklet"
       << "\n";

    // Separator
    os << std::string(70, '-') << "\n";

    // Corpo della tabella
    for (size_t i = 0; i < sim.real_time.size(); ++i) {
        os << std::setw(10) << sim.gen_tracks[i]
           << std::setw(12) << std::fixed << std::setprecision(6) << sim.real_time[i]
           << std::setw(12) << std::fixed << std::setprecision(6) << sim.cpu_time[i]
           << std::setw(12) << std::fixed << std::setprecision(3) << sim.reco_trk[i]
           << std::setw(12) << std::fixed << std::setprecision(3) << sim.hmgth3l[i]
           << std::setw(12) << std::fixed << std::setprecision(3) << sim.tracklet[i]
           << "\n";
    }

    os << std::string(70, '-') << "\n";

    return os;
}
