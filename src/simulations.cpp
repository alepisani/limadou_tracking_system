#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <iomanip>
#include <TApplication.h>
#include <TCanvas.h>
#include <TView.h>
#include <TList.h>
#include "TLegend.h"
#include <TPolyLine3D.h>
#include "TMarker3DBox.h"
#include "TH1F.h"
#include "TH2.h"
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
#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include <algorithm>
#include <malloc.h>
#include "THStack.h"
#include "TH1.h"        // for TH1::AddDirectory(false)
#include "TROOT.h"      // declares gROOT
#include "TDirectory.h" // declares gDirectory

simulations::simulations()
{
    // radius in mm
    // simulations::gen_tracks = {1};
    radius = {0.4};
    // simulations::gen_tracks = {2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 25, 30, 40, 50};
    simulations::gen_tracks = {2, 3, 4, 5, 6, 7, 8, 9, 10};
    // radius = {2, 1.5, 1, 0.9, 0.8, 0.7, 0.65, 0.6, 0.55, 0.5, 0.47, 0.45, 0.43, 0.4, 0.37, 0.35, 0.33, 0.3, 0.27, 0.25, 0.23, 0.2, 0.15, 0.1, 0.05, 0.01, 0}; //mm

    // simulations::gen_tracks = {100};
    // radius = {0.3};

    // limite massimo dimensione del chip ~13.7 mm --> raggio massimo ~ 6mm = 6000 microm
    // limite minimo dimensione singolo pixel ~0.029 mm --> raggio minimo ~ 0.015 = 15 microm
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
            gen_tr3L.push_back(stats::hmgthL012 + stats::hmgth2L);
            // gen_tr3L.push_back(stats::hmgthL012);
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
    float pi = TMath::Pi();
    float degtorad = TMath::DegToRad();
    float radtodeg = TMath::RadToDeg();
    std::vector<float> alltheta;
    std::vector<float> allphi;
    float nbins = (iteration_per_event * radius.size() * gen_tracks.size());
    TH1F *htheta_real = new TH1F("htheta_real", "#theta;#theta (deg);counts", 45, 0, 90);
    TH1F *htheta_reco = new TH1F("htheta_reco", "#theta;#theta (deg);counts", 45, 0, 90);
    TH1F *h_theta_diff = new TH1F("h_theta_diff", "#theta;(#theta_{real} - #theta_{reco}) (deg);counts", 100, 0, 0.03);
    TH1F *h_phi_diff = new TH1F("h_phi_diff", "#phi;(#phi_{real} - #phi_{reco}) (deg);counts", 1000, -190, 190);
    TH1F *hphi_real = new TH1F("hphi_real", "#phi;#phi;counts", 50, -190, 190);
    TH1F *hphi_reco = new TH1F("hphi_reco", "#phi;#phi;counts", 50, -190, 190);
    TH2D *h_real = new TH2D("h_theta_vs_phi_real", "#theta vs #phi;#phi (deg);#theta (deg)", 90, -185, 185, 20, 0, 90);
    TH2D *h_reco = new TH2D("h_theta_vs_phi_reco", "#theta vs #phi;#phi (deg);#theta (deg)", 90, -185, 185, 20, 0, 90);
    TH1F *hx_real = new TH1F("hx_real", "x;x;counts", 50, -90, +90);
    TH1F *hx_reco = new TH1F("hx_reco", "x;x;counts", 50, -90, +90);
    TH2D *hx_reco_theta = new TH2D("hx_reco_theta", "#theta vs x (RECO); x (mm); #theta (deg)", 90, -77, 77, 20, 0, 90);
    TH2D *hx_real_theta = new TH2D("hx_real_theta", "#theta vs x (REAL); x (mm); #theta (deg)", 90, -77, 77, 20, 0, 90);
    TH2D *hy_reco_theta = new TH2D("hy_reco_theta", "#theta vs y (RECO); y (mm); #theta (deg)", 90, -85, 85, 20, 0, 90);
    TH2D *hy_real_theta = new TH2D("hy_real_theta", "#theta vs y (REAL); y (mm); #theta (deg)", 90, -85, 85, 20, 0, 90);
    TH2D *hx_reco_phi = new TH2D("hx_reco_phi", "#phi vs x (RECO); x (mm); #phi (deg)", 90, -77, 77, 50, -185, 185);
    TH2D *hx_real_phi = new TH2D("hx_real_phi", "#phi vs x (REAL); x (mm); #phi (deg)", 90, -77, 77, 50, -185, 185);
    TH2D *hy_reco_phi = new TH2D("hy_reco_phi", "#phi vs y (RECO); y (mm); #phi (deg)", 90, -85, 85, 50, -185, 185);
    TH2D *hy_real_phi = new TH2D("hy_real_phi", "#phi vs y (REAL); y (mm); #phi (deg)", 90, -85, 85, 50, -185, 185);

    TH2D *h_dtheta_theta = new TH2D("h_dtheta_theta", "#Delta_{#theta} vs #theta; #theta (deg);   #Delta_{#theta} (deg)", 180, 0, 90, 100, -12, 12);
    TH2D *h_dtheta_phi = new TH2D("h_dtheta_phi", "#Delta_{#theta} vs #phi;   #phi (deg);     #Delta_{#theta} (deg)", 150, -200, 200, 100, -12, 12);
    TH2D *h_dtheta_x = new TH2D("h_dtheta_x", "#Delta_{#theta} vs x;      x (mm);         #Delta_{#theta} (deg)", 90, -77, 77, 100, -12, 12);
    TH2D *h_dtheta_y = new TH2D("h_dtheta_y", "#Delta_{#theta} vs y;      y (mm);         #Delta_{#theta} (deg)", 90, -85, 85, 100, -12, 12);

    TH2D *h_dphi_theta = new TH2D("h_dphi_theta", "#Delta_{#phi} vs #theta;   #theta (deg);   #Delta_{#phi} (deg)", 360, 0, 90, 150, -300, 300);
    TH2D *h_dphi_phi = new TH2D("h_dphi_phi", "#Delta_{#phi} vs #phi;     #phi (deg);     #Delta_{#phi} (deg)", 150, -200, 200, 150, -300, 300);
    TH2D *h_dphi_x = new TH2D("h_dphi_x", "#Delta_{#phi} vs x;        x (mm);         #Delta_{#phi} (deg)", 90, -77, 77, 150, -300, 300);
    TH2D *h_dphi_y = new TH2D("h_dphi_y", "#Delta_{#phi} vs y;        y (mm);         #Delta_{#phi} (deg)", 90, -85, 85, 150, -300, 300);

    TH2D *h_dx_theta = new TH2D("h_dx_theta", "#Delta_x vs #theta;        #theta (deg);   #Delta_x (mm)", 180, 0, 90, 100, -5, 5);
    TH2D *h_dx_phi = new TH2D("h_dx_phi", "#Delta_x vs #phi;          #phi (deg);     #Delta_x (mm)", 150, -200, 200, 100, -5, 5);
    TH2D *h_dx_x = new TH2D("h_dx_x", "#Delta_x vs x;             x (mm);         #Delta_x (mm)", 90, -77, 77, 100, -5, 5);
    TH2D *h_dx_y = new TH2D("h_dx_y", "#Delta_x vs y;             y (mm);         #Delta_x (mm)", 90, -85, 85, 100, -5, 5);

    TH2D *h_dy_theta = new TH2D("h_dy_theta", "#Delta_x vs #theta;        #theta (deg);   #Delta_x (mm)", 180, 0, 90, 100, -5, 5);
    TH2D *h_dy_phi = new TH2D("h_dy_phi", "#Delta_x vs #phi;          #phi (deg);     #Delta_x (mm)", 150, -200, 200, 100, -5, 5);
    TH2D *h_dy_x = new TH2D("h_dy_x", "#Delta_x vs x;             x (mm);         #Delta_x (mm)", 90, -77, 77, 100, -5, 5);
    TH2D *h_dy_y = new TH2D("h_dy_y", "#Delta_x vs y;             y (mm);         #Delta_x (mm)", 90, -85, 85, 100, -5, 5);

    auto start_time = std::chrono::steady_clock::now();
    display simu;
    LTrackerTrack ltt;
    stats stats;

    simu.take_distributions();

    std::vector<double> r_time, c_time, reco, reco_real, gen_trk, reco_real3, reco_real2, fake3, fake2;
    r_time.reserve(iteration_per_event);
    c_time.reserve(iteration_per_event);
    reco.reserve(iteration_per_event);
    reco_real.reserve(iteration_per_event);
    reco_real3.reserve(iteration_per_event);
    reco_real2.reserve(iteration_per_event);
    gen_trk.reserve(iteration_per_event);
    fake3.reserve(iteration_per_event);
    fake2.reserve(iteration_per_event);

    std::string path = "../data/limadou_sim32L.csv";
    bool file_exists = std::filesystem::exists(path);
    std::ofstream file(path, std::ios::app);

    // Scrivi intestazione solo se il file non esisteva prima
    if (!file_exists)
    {
        file << "GenTrk, Raggio, Eff, err_e, Eff_real, err_er, fake_reco_trk, err_frt, GenTrk, RecoTrk, RecoReal, RealTime, CPUTime, eff3hit, eff2hit, fake3, fake2\n";
    }

    for (int i = 0; i < gen_tracks.size(); ++i)
    {
        std::cout << "==== GenTracks = " << gen_tracks[i] << " ====" << std::endl;
        r_time.clear();
        c_time.clear();
        reco.clear();
        gen_trk.clear();
        reco_real.clear();
        reco_real3.clear();
        reco_real2.clear();
        fake2.clear();
        fake3.clear();

        for (int m = 0; m < radius.size(); ++m)
        {
            cout << endl;
            std::cout << " raggio = " << radius[m] << " mm" << std::endl;
            r_time.clear();
            c_time.clear();
            gen_trk.clear();
            reco.clear();
            reco_real.clear();
            reco_real3.clear();
            reco_real2.clear();
            fake2.clear();
            fake3.clear();
            stats.reset();
            ltt.Reset();

            for (int j = 0; j < iteration_per_event; ++j)
            {
                stats.reset();
                ltt.Reset();

                double theta_reco, theta_real, phi_reco, phi_real, x_reco, x_real, y_reco, y_real;
                simu.tracks_no_print_hist(gen_tracks[i], ltt);
                for (int k = 0; k < simu.generated_tracks.size(); ++k)
                {
                    htheta_real->Fill(simu.generated_tracks[k].theta * radtodeg);
                    hphi_real->Fill(simu.generated_tracks[k].phi * radtodeg);
                    h_real->Fill(simu.generated_tracks[k].phi * radtodeg, simu.generated_tracks[k].theta * radtodeg);
                    hx_real->Fill(simu.generated_tracks[k].x0);
                    hx_real_theta->Fill(simu.generated_tracks[k].x0, simu.generated_tracks[k].theta * radtodeg);
                    hy_real_theta->Fill(simu.generated_tracks[k].y0, simu.generated_tracks[k].theta * radtodeg);
                    hx_real_phi->Fill(simu.generated_tracks[k].x0, simu.generated_tracks[k].phi * radtodeg);
                    hy_real_phi->Fill(simu.generated_tracks[k].y0, simu.generated_tracks[k].phi * radtodeg);
                    theta_real = simu.generated_tracks[k].theta * radtodeg;
                    phi_real = simu.generated_tracks[k].phi * radtodeg;
                    x_real = simu.generated_tracks[k].x0;
                    y_real = simu.generated_tracks[k].y0;
                }
                TStopwatch t;
                t.Start();
                ltt.computeTracklets();
                ltt.new_algo(radius[m]);
                t.Stop();

                for (int g = 0; g < ltt.tracks.size(); ++g)
                {
                    htheta_reco->Fill(ltt.tracks[g].theta * radtodeg);
                    hphi_reco->Fill(ltt.tracks[g].phi * radtodeg);
                    h_reco->Fill(ltt.tracks[g].phi * radtodeg, ltt.tracks[g].theta * radtodeg);
                    hx_reco->Fill(ltt.tracks[g].x0);
                    hx_reco_theta->Fill(ltt.tracks[g].x0, ltt.tracks[g].theta * radtodeg);
                    hy_reco_theta->Fill(ltt.tracks[g].y0, ltt.tracks[g].theta * radtodeg);
                    hx_reco_phi->Fill(ltt.tracks[g].x0, ltt.tracks[g].phi * radtodeg);
                    hy_reco_phi->Fill(ltt.tracks[g].y0, ltt.tracks[g].phi * radtodeg);
                    theta_reco = ltt.tracks[g].theta * radtodeg;
                    phi_reco = ltt.tracks[g].phi * radtodeg;
                    x_reco = ltt.tracks[g].x0;
                    y_reco = ltt.tracks[g].y0;
                }

                simu.generated_tracks.clear();
                r_time.push_back(t.RealTime());
                c_time.push_back(t.CpuTime());
                reco.push_back(stats::hmrt);
                reco_real.push_back(stats::hmrtar);
                reco_real3.push_back(stats::hmrtar3);
                reco_real2.push_back(stats::hmrtar2);
                // gen_trk.push_back(stats::hmgthL012 + stats::hmgth2L);
                gen_trk.push_back(stats::hmgthL012);
                fake3.push_back(stats::hmrtaf3);
                fake2.push_back(stats::hmrtaf2);

                double dtheta = theta_real - theta_reco;
                double dphi = phi_real - phi_reco;
                double dx = x_real - x_reco;
                double dy = y_real - y_reco;

                if (ltt.tracks.size())
                {
                    h_theta_diff->Fill(dtheta);
                    h_phi_diff->Fill(dphi);

                    h_dtheta_theta->Fill(theta_real, dtheta);
                    h_dtheta_phi->Fill(phi_real, dtheta);
                    h_dtheta_x->Fill(x_real, dtheta);
                    h_dtheta_y->Fill(y_real, dtheta);

                    h_dphi_theta->Fill(theta_real, dphi);
                    h_dphi_phi->Fill(phi_real, dphi);
                    h_dphi_x->Fill(x_real, dphi);
                    h_dphi_y->Fill(y_real, dphi);

                    h_dx_theta->Fill(theta_real, dx);
                    h_dx_phi->Fill(phi_real, dx);
                    h_dx_x->Fill(x_real, dx);
                    h_dx_y->Fill(y_real, dx);

                    h_dy_theta->Fill(theta_real, dy);
                    h_dy_phi->Fill(phi_real, dy);
                    h_dy_x->Fill(x_real, dy);
                    h_dy_y->Fill(y_real, dy);
                }

                printProgressBarWithETA(j + 1, iteration_per_event, start_time);
            }

            // errori statistici
            double eff = mean(reco) / mean(gen_trk);
            double eff_real = mean(reco_real) / mean(gen_trk);
            double ineff_fake = eff - eff_real;
            double eff3hit = mean(reco_real3) / mean(gen_trk);
            double eff2hit = mean(reco_real2) / mean(gen_trk);
            double ineff_fake3 = mean(fake3) / mean(gen_trk);
            double ineff_fake2 = mean(fake2) / mean(gen_trk);
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
                 << std::fixed << std::setprecision(6) << eff << ","
                 << std::fixed << std::setprecision(6) << err_eff << ","
                 << std::fixed << std::setprecision(6) << eff_real << ","
                 << std::fixed << std::setprecision(6) << err_effreal << ","
                 << std::fixed << std::setprecision(6) << ineff_fake << ","
                 << std::fixed << std::setprecision(6) << err_ineffake << ","
                 << std::fixed << std::setprecision(6) << mean(gen_trk) << ","
                 << std::fixed << std::setprecision(6) << mean(reco) << ","
                 << std::fixed << std::setprecision(6) << mean(reco_real) << ","
                 << std::fixed << std::setprecision(6) << mean(r_time) << ","
                 << std::fixed << std::setprecision(6) << mean(c_time) << ","
                 << std::fixed << std::setprecision(6) << eff3hit << ","
                 << std::fixed << std::setprecision(6) << eff2hit << ","
                 << std::fixed << std::setprecision(6) << ineff_fake3 << ","
                 << std::fixed << std::setprecision(6) << ineff_fake2 << "\n";
        }

        std::cout << std::endl;
    }

    file.close();
    std::cout << "File CSV scritto correttamente in: " << path << "\n";

    char fil[200];
    sprintf(fil, "../data/simulations_angle_reco.root");
    TFile *f = TFile::Open(fil, "RECREATE");

    // --- write your histograms individually ---
    htheta_reco->Write("htheta_reco", TObject::kOverwrite);
    htheta_real->Write("htheta_real", TObject::kOverwrite);
    hphi_reco->Write("hphi_reco", TObject::kOverwrite);
    hphi_real->Write("hphi_real", TObject::kOverwrite);
    h_reco->Write("h_reco", TObject::kOverwrite);
    h_real->Write("h_real", TObject::kOverwrite);
    hx_reco->Write("hx_reco", TObject::kOverwrite);
    hx_real->Write("hx_real", TObject::kOverwrite);
    hx_reco_theta->Write("hx_reco_theta", TObject::kOverwrite);
    hx_real_theta->Write("hx_real_theta", TObject::kOverwrite);
    hy_reco_theta->Write("hy_reco_theta", TObject::kOverwrite);
    hy_real_theta->Write("hy_real_theta", TObject::kOverwrite);
    hx_reco_phi->Write("hx_reco_phi", TObject::kOverwrite);
    hx_real_phi->Write("hx_real_phi", TObject::kOverwrite);
    hy_reco_phi->Write("hy_reco_phi", TObject::kOverwrite);
    hy_real_phi->Write("hy_real_phi", TObject::kOverwrite);

    h_dtheta_theta->SetStats(0);
    h_dtheta_theta->Write("h_dtheta_theta", TObject::kOverwrite);
    h_dtheta_phi->SetStats(0);
    h_dtheta_phi->Write("h_dtheta_phi", TObject::kOverwrite);
    h_dtheta_x->SetStats(0);
    h_dtheta_x->Write("h_dtheta_x", TObject::kOverwrite);
    h_dtheta_y->SetStats(0);
    h_dtheta_y->Write("h_dtheta_y", TObject::kOverwrite);

    h_dphi_theta->SetStats(0);
    h_dphi_theta->Write("h_dphi_theta", TObject::kOverwrite);
    h_dphi_phi->SetStats(0);
    h_dphi_phi->Write("h_dphi_phi", TObject::kOverwrite);
    h_dphi_x->SetStats(0);
    h_dphi_x->Write("h_dphi_x", TObject::kOverwrite);
    h_dphi_y->SetStats(0);
    h_dphi_y->Write("h_dphi_y", TObject::kOverwrite);

    h_dx_theta->SetStats(0);
    h_dx_theta->Write("h_dx_theta", TObject::kOverwrite);
    h_dx_phi->SetStats(0);
    h_dx_phi->Write("h_dx_phi", TObject::kOverwrite);
    h_dx_x->SetStats(0);
    h_dx_x->Write("h_dx_x", TObject::kOverwrite);
    h_dx_y->SetStats(0);
    h_dx_y->Write("h_dx_y", TObject::kOverwrite);

    h_dy_theta->SetStats(0);
    h_dy_theta->Write("h_dy_theta", TObject::kOverwrite);
    h_dy_phi->SetStats(0);
    h_dy_phi->Write("h_dy_phi", TObject::kOverwrite);
    h_dy_x->SetStats(0);
    h_dy_x->Write("h_dy_x", TObject::kOverwrite);
    h_dy_y->SetStats(0);
    h_dy_y->Write("h_dy_y", TObject::kOverwrite);

    // normalization
    if (htheta_reco->Integral() > 0)
        htheta_reco->Scale(1.0 / htheta_reco->Integral("width"));
    if (htheta_real->Integral() > 0)
        htheta_real->Scale(1.0 / htheta_real->Integral("width"));
    if (hphi_reco->Integral() > 0)
        hphi_reco->Scale(1.0 / hphi_reco->Integral("width"));
    if (hphi_real->Integral() > 0)
        hphi_real->Scale(1.0 / hphi_real->Integral("width"));
    if (h_reco->Integral() > 0)
        h_reco->Scale(1.0 / h_reco->Integral("width"));
    if (h_real->Integral() > 0)
        h_real->Scale(1.0 / h_real->Integral("width"));

    // re-write after scaling
    h_reco->Write("h_reco_normalized", TObject::kOverwrite);
    h_real->Write("h_real_normalized", TObject::kOverwrite);
    h_theta_diff->Write("h_theta_diff", TObject::kOverwrite);
    h_phi_diff->Write("h_phi_diff", TObject::kOverwrite);

    // --- make overlay canvas for theta ---
    TCanvas *c_theta = new TCanvas("c_theta_overlay", "theta reco vs real", 800, 600);
    htheta_reco->SetLineColor(kBlue);
    htheta_real->SetLineColor(kRed);
    htheta_reco->Draw("HIST");
    htheta_real->Draw("HISTSAME");

    TLegend *leg1 = new TLegend(0.6, 0.7, 0.8, 0.8);
    leg1->AddEntry(htheta_reco, "Reco", "l");
    leg1->AddEntry(htheta_real, "Real", "l");
    leg1->Draw();

    c_theta->Write(); // saves the overlay canvas

    // --- make overlay canvas for phi ---
    TCanvas *c_phi = new TCanvas("c_phi_overlay", "phi reco vs real", 800, 600);
    hphi_reco->SetLineColor(kBlue);
    hphi_real->SetLineColor(kRed);
    hphi_reco->Draw("HIST");
    hphi_real->Draw("HISTSAME");

    TLegend *leg2 = new TLegend(0.6, 0.7, 0.8, 0.8);
    leg2->AddEntry(hphi_reco, "Reco", "l");
    leg2->AddEntry(hphi_real, "Real", "l");
    leg2->Draw();

    c_phi->Write(); // saves the overlay canvas

    TCanvas *alldistros = new TCanvas("alldistros", "16 histograms", 1600, 1600);
    alldistros->Divide(4, 4); // 4 columns, 4 rows

    // First row: dtheta
    alldistros->cd(1);
    h_dtheta_theta->Draw("colz");
    alldistros->cd(2);
    h_dtheta_phi->Draw("colz");
    alldistros->cd(3);
    h_dtheta_x->Draw("colz");
    alldistros->cd(4);
    h_dtheta_y->Draw("colz");
    // Second row: dphi
    alldistros->cd(5);
    h_dphi_theta->Draw("colz");
    alldistros->cd(6);
    h_dphi_phi->Draw("colz");
    alldistros->cd(7);
    h_dphi_x->Draw("colz");
    alldistros->cd(8);
    h_dphi_y->Draw("colz");
    // Third row: dx
    alldistros->cd(9);
    h_dx_theta->Draw("colz");
    alldistros->cd(10);
    h_dx_phi->Draw("colz");
    alldistros->cd(11);
    h_dx_x->Draw("colz");
    alldistros->cd(12);
    h_dx_y->Draw("colz");
    // Fourth row: dy
    alldistros->cd(13);
    h_dy_theta->Draw("colz");
    alldistros->cd(14);
    h_dy_phi->Draw("colz");
    alldistros->cd(15);
    h_dy_x->Draw("colz");
    alldistros->cd(16);
    h_dy_y->Draw("colz");

    alldistros->Write();
    alldistros->SaveAs("../data/alldistros.png");

    // finish
    f->Flush();
    f->ls();
    f->Close();

}
