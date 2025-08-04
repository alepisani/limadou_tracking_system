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
#include <fstream>
#include <filesystem>
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

simulations::simulations()
{
    // simulations::gen_tracks = {2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 25, 30};
    // radius = {6, 4, 2, 1.8, 1.6, 1.4, 1.2, 1., 0.8, 0.6, 0.4, 0.2, 0.1, 0.05};
    simulations::gen_tracks = {2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 25, 30, 40, 50};
    radius = {0.5, 0.45, 0.4, 0.35, 0.3, 0.25, 0.23, 0.21, 0.19, 0.17, 0.15, 0.13, 0.11, 0.09, 0.07, 0.05, 0.03, 0.01};

    // limite massimo dimensione del chip ~13.7
    // limite minimo dimensione singolo pixel ~0.029
}

double simulations::mean(const vector<double> &v)
{
    if (v.empty())
        return 0.0;
    return accumulate(v.begin(), v.end(), 0.0) / v.size();
}

void simulations::printProgressBarWithETA(int current, int total, std::chrono::steady_clock::time_point start_time, int barWidth = 30)
{
    using namespace std::chrono;

    float progress = static_cast<float>(current) / total;
    int pos = static_cast<int>(barWidth * progress);

    auto now = steady_clock::now();
    auto elapsed = duration_cast<seconds>(now - start_time).count();
    int eta = (progress > 0.0) ? static_cast<int>(elapsed / progress - elapsed) : 0;

    std::ostringstream oss;
    oss << "\r[";

    for (int i = 0; i < barWidth; ++i)
    {
        if (i < pos)
            oss << "█";
        else
            oss << " ";
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

void simulations::sim_only_trk_3L(int iteration_per_event)
{

    auto start_time = std::chrono::steady_clock::now();
    display simu;
    LTrackerTrack ltt;
    stats stats;
    chips cc;

    simu.take_distributions();

    // std::string path = "/mnt/c/Users/user/Desktop/limadou_sim.csv";
    // std::ofstream file(path, std::ios::out);
    // file << "GenTrk,Raggio,Efficiency,Gen3L,RecoTrk,RecoReal,Tracklet,RealTime,CPUTime\n";

    std::string path = "/mnt/c/Users/user/Desktop/limadou_sim.csv";
    bool file_exists = std::filesystem::exists(path);
    std::ofstream file(path, std::ios::app);

    // Scrivi intestazione solo se il file non esisteva prima
    if (!file_exists)
    {
        file << "GenTrk,Raggio,Efficiency,Gen3L,RecoTrk,RecoReal,Tracklet,RealTime,CPUTime\n";
    }

    for (int i = 0; i < gen_tracks.size(); ++i)
    {
        std::cout << "==== GenTracks = " << gen_tracks[i] << " ====" << std::endl;

        for (int m = 0; m < radius.size(); ++m)
        {
            cout << endl;
            std::cout << " raggio = " << radius[m] << " mm" << std::endl;

            std::vector<double> r_time, c_time, reco, reco_real, gen_tr3L, trkl;
            stats.reset();
            ltt.Reset();

            for (int j = 0; j < iteration_per_event; ++j)
            {
                stats.reset();
                ltt.Reset();
                simu.tracks_no_print_hist(gen_tracks[i], ltt);
                TStopwatch t;
                t.Start();
                ltt.computeTracklets();
                ltt.new_algo(radius[m]);
                t.Stop();

                r_time.push_back(t.RealTime());
                c_time.push_back(t.CpuTime());
                reco.push_back(stats::hmrt);
                reco_real.push_back(stats::hmrtar);
                gen_tr3L.push_back(stats::hmgthL012);
                trkl.push_back(ltt.tracklet_lay02.size());

                printProgressBarWithETA(j + 1, iteration_per_event, start_time);
            }

            // Scrivi una riga per ogni combinazione GenTrack-Raggio
            file << gen_tracks[i] << ","
                 << radius[m] << ","
                 << std::fixed << std::setprecision(3) << mean(reco) / mean(gen_tr3L) << ","
                 << std::fixed << std::setprecision(3) << mean(gen_tr3L) << ","
                 << std::fixed << std::setprecision(3) << mean(reco) << ","
                 << std::fixed << std::setprecision(3) << mean(reco_real) << ","
                 << std::fixed << std::setprecision(3) << mean(trkl) << ","
                 << std::fixed << std::setprecision(6) << mean(r_time) << ","
                 << std::fixed << std::setprecision(6) << mean(c_time) << "\n";
        }

        std::cout << std::endl;
        std::cout << *this << std::endl;
    }

    file.close();
    std::cout << "File CSV scritto correttamente in: " << path << "\n";
}

void simulations::sim_old_algo(int iteration_per_event)
{

    auto start_time = std::chrono::steady_clock::now();
    display simu;
    LTrackerTrack ltt;
    stats stats;
    chips cc;

    simu.take_distributions();

    std::string path = "/mnt/c/Users/user/Desktop/limadou_sim.csv";
    std::ofstream file(path, std::ios::out);
    file << "GenTrk,Efficiency,Gen3L,RecoTrk,RecoReal,Tracklet,RealTime,CPUTime\n";

    for (int i = 0; i < gen_tracks.size(); ++i)
    {
        std::cout << "==== GenTracks = " << gen_tracks[i] << " ====" << std::endl;

        std::vector<double> r_time, c_time, reco, reco_real, gen_tr3L, trkl;
        stats.reset();
        ltt.Reset();

        for (int j = 0; j < iteration_per_event; ++j)
        {
            stats.reset();
            ltt.Reset();
            simu.tracks_no_print_hist(gen_tracks[i], ltt);
            TStopwatch t;
            t.Start();
            ltt.computeTracklets();
            ltt.computeTrackCandidates();
            t.Stop();

            r_time.push_back(t.RealTime());
            c_time.push_back(t.CpuTime());
            reco.push_back(stats::hmrt);
            reco_real.push_back(stats::hmrtar);
            gen_tr3L.push_back(stats::hmgthL012);
            trkl.push_back(ltt.tracklet_lay02.size() + ltt.tracklet_lay01.size() + ltt.tracklet_lay12.size());

            printProgressBarWithETA(j + 1, iteration_per_event, start_time);
        }

        // Scrivi una riga per ogni combinazione GenTrack-Raggio
        file << gen_tracks[i] << ","
             << std::fixed << std::setprecision(3) << mean(reco_real) / mean(gen_tr3L) << ","
             << std::fixed << std::setprecision(3) << mean(gen_tr3L) << ","
             << std::fixed << std::setprecision(3) << mean(reco) << ","
             << std::fixed << std::setprecision(3) << mean(reco_real) << ","
             << std::fixed << std::setprecision(3) << mean(trkl) << ","
             << std::fixed << std::setprecision(6) << mean(r_time) << ","
             << std::fixed << std::setprecision(6) << mean(c_time) << "\n";

        std::cout << std::endl;
        std::cout << *this << std::endl;
    }

    file.close();
    std::cout << "File CSV scritto correttamente in: " << path << "\n";
}

void simulations::sim_trk_32L(int iteration_per_event)
{

    auto start_time = std::chrono::steady_clock::now();
    display simu;
    LTrackerTrack ltt;
    stats stats;
    chips cc;

    simu.take_distributions();

    // std::string path = "/mnt/c/Users/user/Desktop/limadou_sim.csv";
    // std::ofstream file(path, std::ios::out);
    // file << "GenTrk,Raggio,Efficiency,Gen3L,RecoTrk,RecoReal,Tracklet,RealTime,CPUTime\n";

    std::string path = "/mnt/c/Users/user/Desktop/limadou_sim.csv";
    bool file_exists = std::filesystem::exists(path);
    std::ofstream file(path, std::ios::app);

    // Scrivi intestazione solo se il file non esisteva prima
    if (!file_exists)
    {
        file << "GenTrk,Raggio,Eff,Eff_real,fake_reco_trk,GenTrk,RecoTrk,RecoReal,Tracklet,RealTime,CPUTime\n";
    }

    for (int i = 0; i < gen_tracks.size(); ++i)
    {
        std::cout << "==== GenTracks = " << gen_tracks[i] << " ====" << std::endl;

        for (int m = 0; m < radius.size(); ++m)
        {
            cout << endl;
            std::cout << " raggio = " << radius[m] << " mm" << std::endl;

            std::vector<double> r_time, c_time, reco, reco_real, gen_trk, trkl;
            stats.reset();
            ltt.Reset();

            for (int j = 0; j < iteration_per_event; ++j)
            {
                stats.reset();
                ltt.Reset();
                simu.tracks_no_print_hist(gen_tracks[i], ltt);
                TStopwatch t;
                t.Start();
                ltt.computeTracklets();
                ltt.new_algo(radius[m]);
                t.Stop();

                r_time.push_back(t.RealTime());
                c_time.push_back(t.CpuTime());
                reco.push_back(stats::hmrt);
                reco_real.push_back(stats::hmrtar);
                gen_trk.push_back(stats::hmgthL012 + stats::hmgth2L);
                trkl.push_back(ltt.tracklet_lay02.size());

                printProgressBarWithETA(j + 1, iteration_per_event, start_time);
            }

            // Scrivi una riga per ogni combinazione GenTrack-Raggio
            file << gen_tracks[i] << ","
                 << radius[m] << ","
                 << std::fixed << std::setprecision(3) << mean(reco) / mean(gen_trk) << ","
                 << std::fixed << std::setprecision(3) << mean(reco_real) / mean(gen_trk) << ","
                 << std::fixed << std::setprecision(3) << (mean(reco) - mean(reco_real)) / mean(gen_trk) << ","
                 << std::fixed << std::setprecision(3) << mean(gen_trk) << ","
                 << std::fixed << std::setprecision(3) << mean(reco) << ","
                 << std::fixed << std::setprecision(3) << mean(reco_real) << ","
                 << std::fixed << std::setprecision(3) << mean(trkl) << ","
                 << std::fixed << std::setprecision(6) << mean(r_time) << ","
                 << std::fixed << std::setprecision(6) << mean(c_time) << "\n";
        }

        std::cout << std::endl;
        std::cout << *this << std::endl;
    }

    file.close();
    std::cout << "File CSV scritto correttamente in: " << path << "\n";
}

std::ostream &operator<<(std::ostream &os, const simulations &sim)
{
    // Header
    os << "\n=== Simulation Results === \n";
    os << std::setw(10) << "Gen Trk"
       << std::setw(12) << "Gen Trk - reco"
       //<< std::setw(12) << "Gen Trk - recoreal"
       << std::setw(12) << "Gen3L"
       << std::setw(12) << "RecoTrk"
       << std::setw(12) << "RecoReal"
       << std::setw(12) << "Tracklet"
       << std::setw(12) << "RealTime"
       << std::setw(12) << "CPUTime"
       << "\n";

    // Separator
    os << std::string(82, '-') << "\n";

    // Corpo della tabella
    for (size_t i = 0; i < sim.real_time.size(); ++i)
    {
        os << std::setw(10) << sim.gen_tracks[i]
           << std::setw(12) << std::fixed << std::setprecision(3) << sim.delta_gentrk_reco[i]
           //<< std::setw(12) << std::fixed << std::setprecision(3) << sim.delta_gentrk_recoreal[i]
           << std::setw(12) << std::fixed << std::setprecision(3) << sim.hmgth3l[i]
           << std::setw(12) << std::fixed << std::setprecision(3) << sim.reco_trk[i]
           << std::setw(12) << std::fixed << std::setprecision(3) << sim.reco_trk_real[i]
           << std::setw(12) << std::fixed << std::setprecision(3) << sim.tracklet[i]
           << std::setw(12) << std::fixed << std::setprecision(6) << sim.real_time[i]
           << std::setw(12) << std::fixed << std::setprecision(6) << sim.cpu_time[i]
           << "\n";
    }

    os << std::string(70, '-') << "\n";

    return os;
}
