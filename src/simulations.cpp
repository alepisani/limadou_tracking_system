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
#include <TTree.h>
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
#include <cmath>
#include <cstdint>
#include <string>






simulations::simulations()
{
    simulations::gen_tracks = {10};
    // radius = {6, 4, 2, 1.8, 1.6, 1.4, 1.2, 1., 0.8, 0.6, 0.4, 0.2, 0.1, 0.05};
    // simulations::gen_tracks = {2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 25, 30, 40, 50};
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

    std::string path = "../data/limadou_sim3L.csv";
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
        // std::cout << *this << std::endl;
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

    std::string path = "../data/limadou_sim_oldalgo.csv";
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
        // std::cout << *this << std::endl;
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

    std::vector<double> r_time, c_time, reco, reco_real, gen_trk;
    r_time.reserve(iteration_per_event);
    c_time.reserve(iteration_per_event);
    reco.reserve(iteration_per_event);
    reco_real.reserve(iteration_per_event);
    gen_trk.reserve(iteration_per_event);

    // TH1F *htheta = new TH1F("htheta", "theta;#theta;counts", (iteration_per_event * radius.size() * gen_tracks.size()) / 10, -370, 370);
    // TH1F *hphi = new TH1F("hphi", "phi;#phi;counts", (iteration_per_event * radius.size() * gen_tracks.size()) / 10, -370, 370);

    simu.take_distributions();

    std::string path = "../data/limadou_sim32L.csv";
    bool file_exists = std::filesystem::exists(path);
    std::ofstream file(path, std::ios::app);

    // Scrivi intestazione solo se il file non esisteva prima
    if (!file_exists)
    {
        file << "GenTrk, Raggio, Eff, err_e, Eff_real, err_er, fake_reco_trk, err_frt, GenTrk, RecoTrk, RecoReal, RealTime, CPUTime\n";
    }

    for (int i = 0; i < gen_tracks.size(); ++i)
    {
        std::cout << "==== GenTracks = " << gen_tracks[i] << " ====" << std::endl;
        r_time.clear();
        c_time.clear();
        reco.clear();
        gen_trk.clear();
        reco_real.clear();

        for (int m = 0; m < radius.size(); ++m)
        {
            cout << endl;
            std::cout << " raggio = " << radius[m] << " mm" << std::endl;
            r_time.clear();
            c_time.clear();
            reco.clear();
            gen_trk.clear();
            reco_real.clear();
            stats.reset();
            ltt.Reset();

            for (int j = 0; j < iteration_per_event; ++j)
            {
                stats.reset();
                ltt.Reset();

                //LTrackerTrack ltt;
                //stats stats;

                simu.tracks_no_print_hist(gen_tracks[i], ltt);
                TStopwatch t;
                t.Start();
                ltt.computeTracklets();
                ltt.new_algo(radius[m]);
                t.Stop();

                //for (int i = 0; i < ltt.tracks.size(); ++i)
                //{
                //    float &remap_phi = ltt.tracks[i].phi;
                //    if (ltt.tracks[i].phi < -180.)
                //    {
                //        remap_phi = ltt.tracks[i].phi + 180.;
                //    }
                //    if (ltt.tracks[i].phi > 180.)
                //    {
                //        remap_phi = ltt.tracks[i].phi - 180.;
                //    }
                //    hphi->Fill(remap_phi);
                //    htheta->Fill(ltt.tracks[i].theta);
                //}

                r_time.push_back(t.RealTime());
                c_time.push_back(t.CpuTime());
                reco.push_back(stats::hmrt);
                reco_real.push_back(stats::hmrtar);
                gen_trk.push_back(stats::hmgthL012 + stats::hmgth2L);

                printProgressBarWithETA(j + 1, iteration_per_event, start_time);
            }

            printf("\ndimensione dei vector: reco_real %lu bytes\n", sizeof(reco_real));
            printf("dimesione display %lu bytes\n", sizeof(simu));
            printf("dimesione ltt %lu bytes\n", sizeof(ltt));
            printf("dimesione stats %lu bytes\n", sizeof(stats  ));


            // errori statistici
            double eff = mean(reco) / mean(gen_trk);
            double eff_real = mean(reco_real) / mean(gen_trk);
            double ineff_fake = eff - eff_real;
            double err_reco = TMath::RMS(reco.begin(), reco.end()) / TMath::Sqrt(iteration_per_event);
            double err_reco_real = TMath::RMS(reco_real.begin(), reco_real.end()) / TMath::Sqrt(iteration_per_event);
            double err_gen_trk = TMath::RMS(gen_trk.begin(), gen_trk.end()) / TMath::Sqrt(iteration_per_event);
            // propagazione errori
            double err_eff = TMath::Sqrt(pow(err_reco / mean(gen_trk), 2) + pow((mean(reco) * err_gen_trk) / (pow(mean(gen_trk), 2)), 2));
            double err_effreal = TMath::Sqrt(pow(err_reco_real / mean(gen_trk), 2) + pow((mean(reco_real) * err_gen_trk) / (pow(mean(gen_trk), 2)), 2));
            double err_ineffake = TMath::Sqrt(pow(err_eff, 2) + pow(err_effreal, 2));

            // Scrivi una riga per ogni combinazione GenTrack-Raggio
            file << gen_tracks[i] << ","
                 << radius[m] << ","
                 << std::fixed << std::setprecision(3) << eff << ","
                 << std::fixed << std::setprecision(3) << err_eff << ","
                 << std::fixed << std::setprecision(3) << eff_real << ","
                 << std::fixed << std::setprecision(3) << err_effreal << ","
                 << std::fixed << std::setprecision(3) << ineff_fake << ","
                 << std::fixed << std::setprecision(3) << err_ineffake << ","
                 << std::fixed << std::setprecision(3) << mean(gen_trk) << ","
                 << std::fixed << std::setprecision(3) << mean(reco) << ","
                 << std::fixed << std::setprecision(3) << mean(reco_real) << ","
                 << std::fixed << std::setprecision(6) << mean(r_time) << ","
                 << std::fixed << std::setprecision(6) << mean(c_time) << "\n";
        }

        std::cout << std::endl;
        // std::cout << *this << std::endl;
    }

    file.close();
    std::cout << "File CSV scritto correttamente in: " << path << "\n";

    // char fil[200];
    // sprintf(fil, "../data/stats_new_algo_reco.root");
    // TFile f(fil, "RECREATE");
    // htheta->Write();
    // hphi->Write();
    // f.Close();
}



























/* struct RunningStat {
    std::uint64_t n = 0;
    double mean = 0.0;
    double M2 = 0.0;

    void add(double x) {
        ++n;
        double delta = x - mean;
        mean += delta / double(n);
        double delta2 = x - mean;
        M2 += delta * delta2;
    }

    double get_mean() const { return mean; }

    // varianza campionaria (n>1), altrimenti 0
    double get_variance() const { return (n > 1) ? (M2 / double(n - 1)) : 0.0; }

    double get_stddev() const { return std::sqrt(get_variance()); }

    // standard error: stddev / sqrt(N); se N==0 ritorna 0
    double get_std_error() const { return (n > 0) ? (get_stddev() / std::sqrt((double)n)) : 0.0; }

    void reset() { n = 0; mean = 0.0; M2 = 0.0; }
};

void simulations::sim_trk_32L(int iteration_per_event)
{
    auto start_time = std::chrono::steady_clock::now();
    display simu;
    LTrackerTrack ltt;
    stats stats;

    // non accumuliamo più tutti i singoli valori: usiamo statistiche online
    simu.take_distributions();

    std::string path = "../data/limadou_sim32L.csv";
    bool file_exists = std::filesystem::exists(path);
    std::ofstream file(path, std::ios::app);

    if (!file_exists) {
        file << "GenTrk, Raggio, Eff, err_e, Eff_real, err_er, fake_reco_trk, err_frt, GenTrk, RecoTrk, RecoReal, RealTime, CPUTime\n";
    }

    for (size_t ig = 0; ig < gen_tracks.size(); ++ig) {
        std::cout << "==== GenTracks = " << gen_tracks[ig] << " ====" << std::endl;

        for (size_t m = 0; m < radius.size(); ++m) {
            std::cout << std::endl << " raggio = " << radius[m] << " mm" << std::endl;

            // running stats per combinazione GenTracks-Radius
            RunningStat rs_rtime, rs_ctime, rs_reco, rs_recoreal, rs_gentrk;

            // assicurati che Reset rilasci memoria interna se necessario (vedi suggerimento sotto)
            stats.reset();
            ltt.Reset();

            for (int j = 0; j < iteration_per_event; ++j) {
                stats.reset();
                ltt.Reset();

                simu.tracks_no_print_hist(gen_tracks[ig], ltt);

                TStopwatch t;
                t.Start();
                ltt.computeTracklets();
                ltt.new_algo(radius[m]);
                t.Stop();

                // leggi i valori come nel tuo codice originale (presumo membri statici o globali)
                rs_rtime.add(t.RealTime());
                rs_ctime.add(t.CpuTime());
                rs_reco.add(stats::hmrt);
                rs_recoreal.add(stats::hmrtar);
                rs_gentrk.add(stats::hmgthL012 + stats::hmgth2L);

                printProgressBarWithETA(j + 1, iteration_per_event, start_time);

                // (opzionale) se ltt tiene grandi buffer, ogni X iterazioni
                // puoi forzare una compattezza interna (implementa compact_internal_vectors() in LTrackerTrack)
                // if ((j+1) % 10000 == 0) { ltt.compact_internal_vectors(); }
            }

            // estrai medie
            double mean_reco = rs_reco.get_mean();
            double mean_recoreal = rs_recoreal.get_mean();
            double mean_gentrk = rs_gentrk.get_mean();

            // evita divisione per zero
            if (mean_gentrk == 0.0) {
                std::cerr << "Warning: mean_gentrk == 0 for gen_tracks=" << gen_tracks[ig] << " radius=" << radius[m] << "\n";
                mean_gentrk = 1.0; // per evitare NaN, ma segnala logicamente il caso
            }

            double eff = mean_reco / mean_gentrk;
            double eff_real = mean_recoreal / mean_gentrk;
            double ineff_fake = eff - eff_real;

            // errori (std error dalle running stats)
            double err_reco = rs_reco.get_std_error();
            double err_reco_real = rs_recoreal.get_std_error();
            double err_gen_trk = rs_gentrk.get_std_error();

            // propagazione errori (stessa formula usata in origine)
            double err_eff = std::sqrt(std::pow(err_reco / mean_gentrk, 2) +
                                      std::pow((mean_reco * err_gen_trk) / (mean_gentrk * mean_gentrk), 2));
            double err_effreal = std::sqrt(std::pow(err_reco_real / mean_gentrk, 2) +
                                          std::pow((mean_recoreal * err_gen_trk) / (mean_gentrk * mean_gentrk), 2));
            double err_ineffake = std::sqrt(err_eff * err_eff + err_effreal * err_effreal);

            file << gen_tracks[ig] << ","
                 << radius[m] << ","
                 << std::fixed << std::setprecision(3) << eff << ","
                 << std::fixed << std::setprecision(3) << err_eff << ","
                 << std::fixed << std::setprecision(3) << eff_real << ","
                 << std::fixed << std::setprecision(3) << err_effreal << ","
                 << std::fixed << std::setprecision(3) << ineff_fake << ","
                 << std::fixed << std::setprecision(3) << err_ineffake << ","
                 << std::fixed << std::setprecision(3) << mean_gentrk << ","
                 << std::fixed << std::setprecision(3) << mean_reco << ","
                 << std::fixed << std::setprecision(3) << mean_recoreal << ","
                 << std::fixed << std::setprecision(6) << rs_rtime.get_mean() << ","
                 << std::fixed << std::setprecision(6) << rs_ctime.get_mean() << "\n";
        } // radius
        std::cout << std::endl;
    } // gen_tracks

    file.close();
    std::cout << "File CSV scritto correttamente in: " << path << "\n";
} */
