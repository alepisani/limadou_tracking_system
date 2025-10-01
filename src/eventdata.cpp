#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <TApplication.h>
#include <TCanvas.h>
#include <TView.h>
#include <TList.h>
#include <TPolyLine3D.h>
#include "TMarker3DBox.h"
#include <chrono>
#include "TH1F.h"
#include "TH2.h"
#include "TLegend.h"
#include "../include/eventdata.h"
#include "../include/chip.h"
#include "../include/display.h"
#include "../include/stats.h"
#include "../include/LTrackerCluster.h"
#include "../include/LTrackerTrack.h"
#include "../include/simulations.h"
#include "eventdata.h"
using namespace std;

// std::string input_filename = "../data/HEPD02-FM_m-Exp-20250907-071441-Events-00351_01656-p01_L2.root";
// std::string input_filename = "../data/HEPD02-FM_m-Exp-20250907-000001-Events-00351_01437-p01_L2.root";

// std::string input_filename = "../data/HEPD02-FM_m-Exp-20250907-045249-Events-00351_01585-p01_L2.root";      // this file has weird peaks
// std::string input_filename = "../data/HEPD02-FM_m-Exp-20250907-002417-Events-00351_01449-p01_L2.root";
std::string input_filename = "../../data_beam_test/TEST_MUONS_m_MAIN_1000.0MeV_-999.0deg_-0.05V_boot207_run510_L2.root";

eventdata::eventdata() {}

std::unordered_map<int, eventdata> alldata;

std::ostream &operator<<(std::ostream &output, const eventdata &ev)
{
    for (size_t i = 0; i < ev.cls_mean_x.size(); i++)
    {
        // output << "turret_idx: " << (int)ev.DIR_turret_idx[i] << " ";
        // output << "stave_idx: " << (int)ev.DIR_stave_idx[i] << " ";
        // output << "chip_idx0: " << (int)ev.DIR_chip_idx0[i] << " ";
        // output << "chip_idx1: " << (int)ev.DIR_chip_idx1[i] << " ";
        // output << "chip_id:   " << ev.DIR_chip_id[i] << " ";
        // output << "xpos:      " << ev.DIR_xpos[i] << " || ";
        // output << "ypos:      " << ev.DIR_ypos[i] << " || ";
        // output << "zpos:      " << ev.DIR_zpos[i] << " || ";
        // output << "cls_idx:   " << ev.DIR_cls_idx[i];
        // output << "total events: " << alldata.size();
        // output << "hmthL2: " << ;
        // output << "hmthL1: " << ;
        // output << "hmthL0: " << ;
        output << endl;
    }
    return output;
}

void eventdata::takedata()
{

    TFile *file = TFile::Open(input_filename.c_str()); // open for reading

    if (!file || file->IsZombie())
    {
        std::cerr << "Errore nell'aprire il file ROOT\n";
        return;
    }
    TTree *tree = (TTree *)file->Get("L2;1");
    if (!tree)
    {
        std::cerr << "TTree non trovato nel file\n";
        return;
    }

    // Oggetti temporanei per collegare i branch
    std::vector<unsigned int> *cls_size = nullptr;
    std::vector<float> *cls_mean_x = nullptr;
    std::vector<float> *cls_mean_y = nullptr;
    std::vector<float> *cls_mean_z = nullptr;
    std::vector<float> *cls_mean_x_err = nullptr;
    std::vector<float> *cls_mean_y_err = nullptr;

    tree->SetBranchAddress("cls_size", &cls_size); // in pixel
    tree->SetBranchAddress("cls_mean_x", &cls_mean_x);
    tree->SetBranchAddress("cls_mean_y", &cls_mean_y);
    tree->SetBranchAddress("cls_mean_z", &cls_mean_z);
    tree->SetBranchAddress("cls_mean_x_err", &cls_mean_x_err);
    tree->SetBranchAddress("cls_mean_y_err", &cls_mean_y_err);

    Long64_t nEntries = tree->GetEntries();

    for (Long64_t i = 0; i < nEntries; ++i)
    {
        tree->GetEntry(i);
        eventdata ev;
        ev.cls_size = *cls_size;
        ev.cls_mean_x = *cls_mean_x;
        ev.cls_mean_y = *cls_mean_y;
        ev.cls_mean_z = *cls_mean_z;
        ev.cls_mean_x_err = *cls_mean_x_err;
        ev.cls_mean_y_err = *cls_mean_y_err;
        alldata[i] = ev;
    }

    file->Close();
}

void eventdata::analize_data()
{
    /** canvas = false  -->  used for the all dataset, store all the stats in the ttree
     *  canvas = true   -->  used for a batch, just a way to fix the algo
     */
    bool print_canvas = false;

    double radtodeg = TMath::RadToDeg();
    takedata();
    int cls_max = 0;
    int index = 0;
    auto start_time = std::chrono::steady_clock::now();
    float pi = TMath::Pi();
    int nbins = alldata.size() / 200;
    int theta_bins = 50;
    int phi_bins = 50;
    TH1F *h_dx0 = new TH1F("hdx0", "#dx0;#dx0;counts", 1000, -5, 5);
    TH1F *h_dx1 = new TH1F("hdx1", "#dx1;#dx1;counts", 1000, -5, 5);
    TH1F *h_dx2 = new TH1F("hdx2", "#dx2;#dx2;counts", 1000, -5, 5);
    TH1F *h_dy0 = new TH1F("hdy0", "#dy0;#dy0;counts", 1000, -5, 5);
    TH1F *h_dy1 = new TH1F("hdy1", "#dy1;#dy1;counts", 1000, -5, 5);
    TH1F *h_dy2 = new TH1F("hdy2", "#dy2;#dy2;counts", 1000, -5, 5);

    TH1F *h_errx0 = new TH1F("herrx0", "errx;errx;counts", 100, 0, 0.1);
    TH1F *h_erry0 = new TH1F("herry0", "erry;erry;counts", 100, 0, 0.1);

    TH1F *htheta = new TH1F("htheta", "#theta;#theta;counts", theta_bins, -5, 90);
    TH1F *hphi = new TH1F("hphi", "#phi;#phi;counts", phi_bins, -185, 185);
    TH2D *h = new TH2D("h_theta_vs_phi", "#theta vs #phi;#phi (deg);#theta (deg)", nbins, -185, 185, nbins, 0, 90);
    TH1F *hchi2 = new TH1F("hchi2_newalgo", "#chi2;#chi2;counts", 500, 0, 50000);
    TH2D *h_chi2_theta = new TH2D("h_chi2_theta", "#chi^2 vs #theta;        #theta (deg);   #chi^2", 30, 0, 90, 50, -2, 500);
    h_chi2_theta->SetStats(0);
    TH2D *h_chi2_phi = new TH2D("h_chi2_phi", "#chi^2 vs #phi;          #phi (deg);     #chi^2", 30, -200, 200, 50, -2, 500);
    h_chi2_phi->SetStats(0);

    TH2D *h_chi2_dx0 = new TH2D("h_chi2_dx0", "#chi^{2} vs dx0;        dx0;   #chi^{2}", 10, -1, 1, 502, -2, 500);
    TH2D *h_chi2_dx1 = new TH2D("h_chi2_dx1", "#chi^{2} vs dx1;        dx1;   #chi^{2}", 10, -1, 1, 502, -2, 500);
    TH2D *h_chi2_dx2 = new TH2D("h_chi2_dx2", "#chi^{2} vs dx2;        dx2;   #chi^{2}", 10, -1, 1, 502, -2, 500);
    TH2D *h_chi2_dy0 = new TH2D("h_chi2_dy0", "#chi^{2} vs dy0;        dy0;   #chi^{2}", 10, -1, 1, 502, -2, 500);
    TH2D *h_chi2_dy1 = new TH2D("h_chi2_dy1", "#chi^{2} vs dy1;        dy1;   #chi^{2}", 10, -1, 1, 502, -2, 500);
    TH2D *h_chi2_dy2 = new TH2D("h_chi2_dy2", "#chi^{2} vs dy2;        dy2;   #chi^{2}", 10, -1, 1, 502, -2, 500);

    TH2D *h_chi2_errx0 = new TH2D("h_chi2_errx0", "#chi^{2} vs errx0;        errx0;   #chi^{2}", 100, 0, 0.05, 100, -2, 500);
    TH2D *h_chi2_erry0 = new TH2D("h_chi2_erry0", "#chi^{2} vs erry0;        erry0;   #chi^{2}", 100, 0, 0.05, 100, -2, 500);

    TH2D *h_err_dx = new TH2D("herrdx", "err vs dx;        err;   dx", 100, 0, 0.05, 100, -1, 1);

    TCanvas *canvas = new TCanvas("MC_tracks", "3D View_mc", 800, 600);
    TView *rt = TView::CreateView(1);
    rt->SetRange(-100, -100, 0, 100, 100, 70);
    rt->ShowAxis();

    display d;
    chips cc;
    stats s;
    simulations sim;
    if (print_canvas)
    {
        cc.print_all_chips(cc, canvas);
        d.draw_TR12(canvas);
    }

    // selecting with the index the event you want to make the reco
    int n;
    if (!print_canvas)
        n = alldata.size();
    if (print_canvas)
        n = 1000;
    for (int i = 0; i < n; ++i)
    {
        LTrackerTrack ltt;
        eventdata ev;
        LTrackerCluster cl;
        ev = alldata[i];

        ev.from_ev_to_cluster(cl, ev);
        for (int j = 0; j < cl.cls_mean_z.size(); ++j)
        {
            LCluster c;
            c.x = cl.cls_mean_x[j];
            c.y = cl.cls_mean_y[j];
            c.z = cl.cls_mean_z[j];
            c.errx = cl.cls_mean_err_x[j];
            c.erry = cl.cls_mean_err_y[j];
            c.id = i;
            // cout << "\n cluster: \n" << c;
            if (c.z < 36. && c.z > 30.)
            {
                ltt.tidy_clusters_lay2.try_emplace(j, c);
                if (print_canvas)
                {
                    TMarker3DBox *p = new TMarker3DBox(c.x, c.y, c.z, 2, 2, 0, 0, 0);
                    p->SetLineWidth(1.4);
                    p->Draw();
                }
            }
            if (c.z < 29. && c.z > 22.)
            {
                ltt.tidy_clusters_lay1.try_emplace(j, c);
                if (print_canvas)
                {
                    TMarker3DBox *p = new TMarker3DBox(c.x, c.y, c.z, 2, 2, 0, 0, 0);
                    p->SetLineWidth(1.4);
                    p->Draw();
                }
            }
            if (c.z < 20. && c.z > 15.)
            {
                ltt.tidy_clusters_lay0.try_emplace(j, c);
                if (print_canvas)
                {
                    TMarker3DBox *p = new TMarker3DBox(c.x, c.y, c.z, 2, 2, 0, 0, 0);
                    p->SetLineWidth(1.4);
                    p->Draw();
                }
            }
        }

        ltt.computeTracklets();
        ltt.new_algo(0.4);

        if (!print_canvas)
        {
            for (int m = 0; m < ltt.tracks.size(); ++m)
            {
                htheta->Fill(ltt.tracks[m].theta * radtodeg);
                hphi->Fill(ltt.tracks[m].phi * radtodeg);
                hchi2->Fill(ltt.tracks[m].chi2);
                h->Fill(ltt.tracks[m].phi * radtodeg, ltt.tracks[m].theta * radtodeg);
                h_chi2_theta->Fill(ltt.tracks[m].theta * radtodeg, ltt.tracks[m].chi2);
                h_chi2_phi->Fill(ltt.tracks[m].phi * radtodeg, ltt.tracks[m].chi2);

                // if(ltt.tracks[m].n_clus == 3){}
                h_chi2_dx0->Fill(ltt.tracks[m].dx0, ltt.tracks[m].chi2);
                h_chi2_dx1->Fill(ltt.tracks[m].dx1, ltt.tracks[m].chi2);
                h_chi2_dx2->Fill(ltt.tracks[m].dx2, ltt.tracks[m].chi2);
                h_chi2_dy0->Fill(ltt.tracks[m].dy0, ltt.tracks[m].chi2);
                h_chi2_dy1->Fill(ltt.tracks[m].dy1, ltt.tracks[m].chi2);
                h_chi2_dy2->Fill(ltt.tracks[m].dy2, ltt.tracks[m].chi2);

                h_dx0->Fill(ltt.tracks[m].dx0);
                h_dx1->Fill(ltt.tracks[m].dx1);
                h_dx2->Fill(ltt.tracks[m].dx2);
                h_dy0->Fill(ltt.tracks[m].dy0);
                h_dy1->Fill(ltt.tracks[m].dy1);
                h_dy2->Fill(ltt.tracks[m].dy2);

                h_errx0->Fill(ltt.tracks[m].err_x0);
                h_erry0->Fill(ltt.tracks[m].err_y0);

                h_chi2_errx0->Fill(ltt.tracks[m].err_x0, ltt.tracks[m].chi2);
                h_chi2_erry0->Fill(ltt.tracks[m].err_y0, ltt.tracks[m].chi2);

                h_err_dx->Fill(ltt.tracks[m].err_x0, ltt.tracks[m].dx0);
            }
        }

        if (print_canvas)
        {
            ltt.printRecoTracks_new_alg(canvas);
            // ltt.print_all_tracklet(ltt);
        }
        sim.printProgressBarWithETA(i + 1, n, start_time, 30);

        if (print_canvas)
        {
            for (int m = 0; m < ltt.tracks.size(); ++m)
            {
                if (ltt.tracks[m].chi2 > 400.)
                {
                    // cout << ltt.tracks[m] << endl;
                }
            }
        }
    }

    if (!print_canvas)
    {
        TH1F *h_theta = new TH1F("h_theta", "h_theta", theta_bins, -5, 90);
        h_theta->SetStats(0);
        TH1F *h_theta_m2 = new TH1F("h_theta_m2", "h_theta_m2", theta_bins, -5, 90);
        h_theta_m2->SetStats(0);
        TH1F *h_phi = new TH1F("h_phi", "h_phi", phi_bins, -185, 185);
        h_phi->SetStats(0);
        TH1F *h_phi_m2 = new TH1F("h_phi_m2", "h_phi_m2", phi_bins, -185, 185);
        h_phi_m2->SetStats(0);
        TH1F *h_chi2 = new TH1F("h_chi2", "#chi2;#chi2;counts", 500, 0, 500);
        h_chi2->SetStats(0);
        TH1F *h_chi2_m2 = new TH1F("h_chi2_m2", "#chi2;#chi2;counts", 500, 0, 500);
        h_chi2_m2->SetStats(0);

        TFile *fIn = new TFile(input_filename.c_str());
        TTree *oldTree = (TTree *)fIn->Get("L2");

        TFile *fOut = new TFile("../data/overlay_plots.root", "RECREATE");

        // --- Input branch (std::vector<float>) ---
        std::vector<float> *theta = nullptr;
        std::vector<float> *phi = nullptr;
        std::vector<float> *chi2 = nullptr;
        std::vector<float> *theta_m2 = nullptr;
        std::vector<float> *phi_m2 = nullptr;
        std::vector<float> *chi2_m2 = nullptr;
        oldTree->SetBranchAddress("theta", &theta);
        oldTree->SetBranchAddress("phi", &phi);
        oldTree->SetBranchAddress("chi2", &chi2);
        oldTree->SetBranchAddress("theta_m2", &theta_m2);
        oldTree->SetBranchAddress("phi_m2", &phi_m2);
        oldTree->SetBranchAddress("chi2_m2", &chi2_m2);

        Long64_t nentries = oldTree->GetEntries();
        for (Long64_t i = 0; i < nentries; i++)
        {
            oldTree->GetEntry(i);

            // Add size checks before accessing vectors
            if (theta && phi && theta_m2 && phi_m2 && chi2 && chi2_m2)
            {
                for (size_t j = 0; j < theta->size(); ++j)
                {
                    // mask for passing throgh the triggers
                    if (theta->at(j) < 75.3)
                    {
                        h_theta->Fill(theta->at(j));
                        h_phi->Fill(phi->at(j));
                        h_chi2->Fill(chi2->at(j));
                    }
                }
                for (size_t j = 0; j < theta_m2->size(); ++j)
                {
                    if (theta_m2->at(j) < 75.3)
                    {
                        h_theta_m2->Fill(theta_m2->at(j));
                        h_phi_m2->Fill(phi_m2->at(j));
                        h_chi2_m2->Fill(chi2_m2->at(j));
                    }
                }
            }
        }
        h_theta->Scale(1.0 / h_theta->Integral("width"));
        h_theta_m2->Scale(1.0 / h_theta_m2->Integral("width"));
        htheta->Scale(1.0 / htheta->Integral("width"));
        h_phi->Scale(1.0 / h_phi->Integral("width"));
        h_phi_m2->Scale(1.0 / h_phi_m2->Integral("width"));
        hphi->Scale(1.0 / hphi->Integral("width"));

        h_theta->Write();
        h_theta_m2->Write();
        htheta->Write();
        h_phi->Write();
        h_phi_m2->Write();
        hphi->Write();
        h_chi2->Write();
        h_chi2_m2->Write();
        hchi2->Write();
        h_chi2_theta->Write();
        h_chi2_phi->Write();

        h_chi2_dx0->Write();
        h_chi2_dx1->Write();
        h_chi2_dx2->Write();
        h_chi2_dy0->Write();
        h_chi2_dy1->Write();
        h_chi2_dy2->Write();

        h_errx0->Write();
        h_erry0->Write();

        h_chi2_errx0->Write();
        h_chi2_erry0->Write();

        h_err_dx->Write();

        // Find the maximum y value among all histograms
        double maxY1 = h_theta->GetMaximum();
        double maxY2 = h_theta_m2->GetMaximum();
        double maxY3 = htheta->GetMaximum();
        double maxY = maxY1;
        if (maxY2 > maxY)
            maxY = maxY2;
        if (maxY3 > maxY)
            maxY = maxY3;

        // Set the y-axis range with 10% padding
        h_theta->GetYaxis()->SetRangeUser(0, maxY * 1.1);
        TCanvas *c_theta = new TCanvas("c_theta_overlay", "compare_theta", 800, 600);
        h_theta->SetTitle("#theta comparison");
        h_theta->GetXaxis()->SetTitle("#theta (deg)");
        h_theta->GetYaxis()->SetTitle("Counts");
        h_theta->SetLineColor(kBlue);
        h_theta_m2->SetLineColor(kRed);
        htheta->SetLineColor(kBlack);
        h_theta->Draw("HIST");
        h_theta_m2->Draw("HISTSAME");
        htheta->Draw("HISTSAME");
        h->GetXaxis()->SetTitle("Theta [rad]");
        h->GetYaxis()->SetTitle("Counts");

        TLegend *leg1 = new TLegend(0.6, 0.7, 0.9, 0.9);
        leg1->AddEntry(h_theta, Form("hough transform (N=%.0f)", h_theta->GetEntries()), "l");
        leg1->AddEntry(h_theta_m2, Form("old_algo (N=%.0f)", h_theta_m2->GetEntries()), "l");
        leg1->AddEntry(htheta, Form("new_algo (N=%.0f)", htheta->GetEntries()), "l");
        leg1->Draw();
        c_theta->Write();

        // Find the maximum y value among all histograms
        double maxY1phi = h_phi->GetMaximum();
        double maxY2phi = h_phi_m2->GetMaximum();
        double maxY3phi = hphi->GetMaximum();
        double maxYphi = maxY1phi;
        if (maxY2phi > maxYphi)
            maxYphi = maxY2phi;
        if (maxY3phi > maxYphi)
            maxYphi = maxY3phi;

        // Set the y-axis range with 10% padding
        h_phi->GetYaxis()->SetRangeUser(0, maxYphi * 1.1);

        TCanvas *c_phi = new TCanvas("c_phi_overlay", "compare_phi", 800, 600);
        h_phi->SetTitle("#phi comparison");
        h_phi->GetXaxis()->SetTitle("#phi (deg)");
        h_phi->GetYaxis()->SetTitle("Counts");
        h_phi->SetLineColor(kBlue);
        h_phi_m2->SetLineColor(kRed);
        hphi->SetLineColor(kBlack);
        h_phi->Draw("HIST");
        h_phi_m2->Draw("HISTSAME");
        hphi->Draw("HISTSAME");
        TLegend *leg2 = new TLegend(0.6, 0.7, 0.9, 0.9);
        leg2->AddEntry(h_phi, Form("hough transform (N=%.0f)", h_phi->GetEntries()), "l");
        leg2->AddEntry(h_phi_m2, Form("old_algo (N=%.0f)", h_phi_m2->GetEntries()), "l");
        leg2->AddEntry(hphi, Form("new_algo (N=%.0f)", hphi->GetEntries()), "l");
        leg2->Draw();
        c_phi->Write();

        h_dx0->Write();
        h_dx1->Write();
        h_dx2->Write();
        h_dy0->Write();
        h_dy1->Write();
        h_dy2->Write();

        fOut->Close();
        fIn->Close();
    }
}

void eventdata::print_data_on_canvas(TCanvas *can)
{

    int negative_clus;
    can->cd();
    display d;
    d.draw_TR12(can);
    chips c;
    c.print_all_chips(c, can);

    // int events = alldata.size();
    int events = 1;
    TH1F *hx = new TH1F("x", "x;x;counts", events, -display::TR2Size[0] * 2.5, display::TR2Size[0] * 2.5);
    TH1F *hy = new TH1F("y", "y;y;counts", events, -100, 100);
    TH1F *hz = new TH1F("z", "z;z;counts", events, 10, 40);

    eventdata ev;
    std::vector<LTrackerCluster> clusters;

    for (int i = 0; i < events; i++)
    {
        LTrackerCluster cluster;
        ev = alldata[i];
        // cluster.CalculateClusterPosition(ev);
        ev.from_ev_to_cluster(cluster, ev);

        clusters.push_back(cluster);
    }

    cout << "events: " << events << endl;

    for (int i = 0; i < clusters.size(); i++)
    {
        // cout << "~~~~~~~~~~~~~~~~~~~~~~" << endl;
        // cout << clusters[i] << endl;
        bool hit_L2 = false;
        bool hit_L1 = false;
        bool hit_L0 = false;
        if (clusters[i].cls_mean_x.empty())
        {
            stats::hmbh0L++;
        }
        for (int j = 0; j < clusters[i].cls_mean_x.size(); j++)
        {

            TMarker3DBox *p = new TMarker3DBox(clusters[i].cls_mean_x[j], clusters[i].cls_mean_y[j], clusters[i].cls_mean_z[j], 2, 2, 0, 0, 0);
            p->SetLineColor(kRed);
            p->SetLineWidth(3);
            p->Draw();
            hx->Fill(clusters[i].cls_mean_x[j]);
            hy->Fill(clusters[i].cls_mean_y[j]);
            hz->Fill(clusters[i].cls_mean_z[j]);

            // stats

            if (clusters[i].cls_mean_z[j] < 36. && clusters[i].cls_mean_z[j] > 30.)
            {
                stats::hmthL2++;
                hit_L2 = true;
            }
            if (clusters[i].cls_mean_z[j] < 29. && clusters[i].cls_mean_z[j] > 22.)
            {
                stats::hmthL1++;
                hit_L1 = true;
            }
            if (clusters[i].cls_mean_z[j] < 20. && clusters[i].cls_mean_z[j] > 15.)
            {
                stats::hmthL0++;
                hit_L0 = true;
            }
            if (clusters[i].cls_mean_z[j] < 0.)
            {
                negative_clus++;
            }
        }
        if (hit_L2 && hit_L1 && hit_L0)
        {
            stats::hmbh3L++;
        }
        if ((hit_L2 && hit_L1 && !hit_L0) || (hit_L2 && !hit_L1 && hit_L0) || (!hit_L2 && hit_L1 && hit_L0))
        {
            stats::hmbh2L++;
        }
        if ((hit_L2 && !hit_L1 && !hit_L0) || (!hit_L2 && hit_L1 && !hit_L0) || (!hit_L2 && !hit_L1 && hit_L0))
        {
            stats::hmbh1L++;
        }
    }

    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   negative_clus  " << negative_clus << endl;
    char file[200];
    sprintf(file, "../data_beam_test/data.root");
    TFile f(file, "RECREATE");
    hx->Write();
    hy->Write();
    hz->Write();

    f.Close();
}

void eventdata::from_ev_to_cluster(LTrackerCluster &cluster, eventdata &ev)
{

    cluster.cls_mean_x = ev.cls_mean_x;
    cluster.cls_mean_y = ev.cls_mean_y;
    cluster.cls_mean_z = ev.cls_mean_z;
    cluster.cls_mean_err_x = ev.cls_mean_x_err;
    cluster.cls_mean_err_y = ev.cls_mean_y_err;
    cluster.cls_idx = {-1, -1, -1};

    // values to be discuss and maybe change
}
