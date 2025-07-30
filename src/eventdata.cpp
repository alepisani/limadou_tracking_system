#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <TApplication.h>
#include <TCanvas.h>
#include <TView.h>
#include <TList.h>
#include <TPolyLine3D.h>
#include "TMarker3DBox.h"
#include "TH1F.h"
#include "../include/eventdata.h"
#include "../include/chip.h"
#include "../include/display.h"
#include "../include/stats.h"
#include "../include/LTrackerCluster.h"
#include "../include/LTrackerTrack.h"
#include "eventdata.h"
using namespace std;

eventdata::eventdata(){}

std::unordered_map<int, eventdata> alldata;


std::ostream &operator<<(std::ostream &output, const eventdata &ev) {
    for (size_t i = 0; i < ev.cls_mean_x.size(); i++) {
        // output << "turret_idx: " << (int)ev.DIR_turret_idx[i] << " ";
        // output << "stave_idx: " << (int)ev.DIR_stave_idx[i] << " ";
        // output << "chip_idx0: " << (int)ev.DIR_chip_idx0[i] << " ";
        // output << "chip_idx1: " << (int)ev.DIR_chip_idx1[i] << " ";
        // output << "chip_id:   " << ev.DIR_chip_id[i] << " ";
        //output << "xpos:      " << ev.DIR_xpos[i] << " || ";
        //output << "ypos:      " << ev.DIR_ypos[i] << " || ";
        //output << "zpos:      " << ev.DIR_zpos[i] << " || ";
        // output << "cls_idx:   " << ev.DIR_cls_idx[i];
        // output << "total events: " << alldata.size();
        //output << "hmthL2: " << ;
        //output << "hmthL1: " << ;
        //output << "hmthL0: " << ;
        output << endl;
    }
    return output;
}



void eventdata::takedata(){

    //cambia nome file in base a cosa vuoi
    //TFile *file = TFile::Open("../data_beam_test/");
    
    TFile *file = TFile::Open("../data_beam_test/TEST_MUONS_m_MAIN_1000.0MeV_-999.0deg_-0.05V_boot207_run510_L2.root");
    if (!file || file->IsZombie()) {
        std::cerr << "Errore nell'aprire il file ROOT\n";
        return;
    }
    TTree *tree = (TTree*)file->Get("L2;1");
    if (!tree) {
        std::cerr << "TTree non trovato nel file\n";
        return;
    }

    // Oggetti temporanei per collegare i branch
    /* std::vector<unsigned char> *DIR_turret_idx = nullptr;
    std::vector<unsigned char> *DIR_stave_idx = nullptr;
    std::vector<unsigned char> *DIR_chip_idx0 = nullptr;
    std::vector<unsigned char> *DIR_chip_idx1 = nullptr;
    std::vector<unsigned int>  *DIR_chip_id = nullptr;
    std::vector<float> *DIR_xpos = nullptr;
    std::vector<float> *DIR_ypos = nullptr;
    std::vector<float> *DIR_zpos = nullptr;
    std::vector<int> *DIR_cls_idx = nullptr; */
    std::vector<unsigned int> *cls_size = nullptr;
    std::vector<float> *cls_mean_x = nullptr;
    std::vector<float> *cls_mean_y = nullptr;
    std::vector<float> *cls_mean_z = nullptr;
    std::vector<float> *cls_mean_x_err = nullptr;
    std::vector<float> *cls_mean_y_err = nullptr;
    

    /* tree->SetBranchAddress("DIR_turret_idx", &DIR_turret_idx);
    tree->SetBranchAddress("DIR_stave_idx",  &DIR_stave_idx);
    tree->SetBranchAddress("DIR_chip_idx0",  &DIR_chip_idx0);
    tree->SetBranchAddress("DIR_chip_idx1",  &DIR_chip_idx1);
    tree->SetBranchAddress("DIR_chip_id",    &DIR_chip_id);
    tree->SetBranchAddress("DIR_xpos",       &DIR_xpos);
    tree->SetBranchAddress("DIR_ypos",       &DIR_ypos);
    tree->SetBranchAddress("DIR_zpos",       &DIR_zpos);
    tree->SetBranchAddress("DIR_cls_idx",    &DIR_cls_idx); */
    tree->SetBranchAddress("cls_size",  &cls_size);    //in pixel
    tree->SetBranchAddress("cls_mean_x",  &cls_mean_x);
    tree->SetBranchAddress("cls_mean_y",  &cls_mean_y);
    tree->SetBranchAddress("cls_mean_z",  &cls_mean_z);
    tree->SetBranchAddress("cls_mean_x_err",  &cls_mean_x_err);
    tree->SetBranchAddress("cls_mean_y_err",  &cls_mean_y_err);

    


    Long64_t nEntries = tree->GetEntries();

    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);
        eventdata ev;
        /* ev.DIR_turret_idx = *DIR_turret_idx;
        ev.DIR_stave_idx  = *DIR_stave_idx;
        ev.DIR_chip_idx0  = *DIR_chip_idx0;
        ev.DIR_chip_idx1  = *DIR_chip_idx1;
        ev.DIR_chip_id    = *DIR_chip_id;
        ev.DIR_xpos       = *DIR_xpos;
        ev.DIR_ypos       = *DIR_ypos;
        ev.DIR_zpos       = *DIR_zpos;
        ev.DIR_cls_idx    = *DIR_cls_idx; */
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

void eventdata::analize_data(){

    
    TCanvas* canvas = new TCanvas("MC_tracks", "3D View_mc", 800, 600);
    TView* rt = TView::CreateView(1);
    rt->SetRange(-100, -100, 0, 100, 100, 70);
    rt->ShowAxis();

    takedata();
    display d;
    LTrackerTrack ltt;
    /* for(int i=0; i < alldata.size(); ++i){
        eventdata ev;
        LTrackerCluster cl;
        LCluster c;
        ev = alldata[i];
        ev.from_ev_to_cluster(cl, ev, c);
        if(cl.cls_mean_z < 36. && cl.cls_mean_z > 30.){
            ltt.tidy_clusters_lay2.push_back(c);
        }
        if(cl.cls_mean_z < 29. && cl.cls_mean_z > 22.){
            ltt.tidy_clusters_lay1.push_back(c);
        }
        if(cl.cls_mean_z < 20. && cl.cls_mean_z > 15.){
            ltt.tidy_clusters_lay0.push_back(c);
        }
        
        
        ltt.computeTracklets();
        ltt.new_computing(2.);
    } */
    eventdata ev;
    LTrackerCluster cl;
    LCluster c;
    chips cc;
    stats s;
    ev = alldata[0];
    cc.print_all_chips(cc,canvas);
    d.draw_TR12(canvas);
    ev.from_ev_to_cluster(cl, ev);
    for (int j = 0; j < cl.cls_mean_z.size(); ++j) {
        LCluster c;
        c.x = cl.cls_mean_x[j];
        c.y = cl.cls_mean_y[j];
        c.z = cl.cls_mean_z[j];
        c.errx = cl.cls_mean_err_x[j];
        c.erry = cl.cls_mean_err_y[j];
        c.id = j;
        if (c.z < 36. && c.z > 30.) {
            ltt.tidy_clusters_lay2.try_emplace(j, c);
        }
        if (c.z < 29. && c.z > 22.) {
            ltt.tidy_clusters_lay1.try_emplace(j, c);
        }
        if (c.z < 20. && c.z > 15.) {
            ltt.tidy_clusters_lay0.try_emplace(j, c);
        }
    }

    ltt.computeTracklets();
    ltt.new_computing(2.);
    ltt.printRecoTracks_new_alg(canvas);

}

void eventdata::print_data_on_canvas(TCanvas* can){

    int negative_clus;
    can->cd();
    display d;
    d.draw_TR12(can);
    chips c;
    c.print_all_chips(c, can);

    //int events = alldata.size();
    int events = 1;
    TH1F* hx = new TH1F("x", "x;x;counts", events, -display::TR2Size[0]*2.5, display::TR2Size[0]*2.5);
    TH1F* hy = new TH1F("y", "y;y;counts", events, -100, 100);
    TH1F* hz = new TH1F("z", "z;z;counts", events, 10, 40);

    eventdata ev;
    std::vector<LTrackerCluster> clusters;

    for (int i = 0; i < events; i++){
        LTrackerCluster cluster;
        ev = alldata[i];
        //cluster.CalculateClusterPosition(ev);
        ev.from_ev_to_cluster(cluster, ev);

        clusters.push_back(cluster);
    }    

    cout << "events: " << events << endl;

    for(int i=0; i<clusters.size(); i++){
        //cout << "~~~~~~~~~~~~~~~~~~~~~~" << endl;
        //cout << clusters[i] << endl;
        bool hit_L2 = false;
        bool hit_L1 = false;
        bool hit_L0 = false;
        if(clusters[i].cls_mean_x.empty()){
        stats::hmbh0L++;
        }
        for(int j=0; j<clusters[i].cls_mean_x.size(); j++){

            TMarker3DBox *p = new TMarker3DBox(clusters[i].cls_mean_x[j],clusters[i].cls_mean_y[j],clusters[i].cls_mean_z[j],2,2,0,0,0);
            p->SetLineColor(kRed);
            p->SetLineWidth(3);
            p->Draw();
            hx->Fill(clusters[i].cls_mean_x[j]);
            hy->Fill(clusters[i].cls_mean_y[j]);
            hz->Fill(clusters[i].cls_mean_z[j]);    
                 
            //stats

            if(clusters[i].cls_mean_z[j] < 36. && clusters[i].cls_mean_z[j] > 30.){
                stats::hmthL2++;
                hit_L2 = true;
            }
            if(clusters[i].cls_mean_z[j] < 29. && clusters[i].cls_mean_z[j] > 22.){
                stats::hmthL1++;
                hit_L1 = true;
            }
            if(clusters[i].cls_mean_z[j] < 20. && clusters[i].cls_mean_z[j] > 15.){
                stats::hmthL0++;
                hit_L0 = true;
            }
            if(clusters[i].cls_mean_z[j] < 0.){
                negative_clus++;
            }
        }   
        if(hit_L2 && hit_L1 && hit_L0){
            stats::hmbh3L++;}
        if((hit_L2 && hit_L1 && !hit_L0) || (hit_L2 && !hit_L1 && hit_L0) || (!hit_L2 && hit_L1 && hit_L0)){
            stats::hmbh2L++;}
        if((hit_L2 && !hit_L1 && !hit_L0) || (!hit_L2 && hit_L1 && !hit_L0) || (!hit_L2 && !hit_L1 && hit_L0)){
            stats::hmbh1L++;}
    }
    


    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   negative_clus  " << negative_clus << endl;
    char file[200];
    sprintf(file,"../data_beam_test/data.root");
    TFile f(file,"RECREATE");
    hx->Write();
    hy->Write();
    hz->Write();
    
    f.Close();
}

void eventdata::from_ev_to_cluster(LTrackerCluster& cluster, eventdata& ev){

    cluster.cls_mean_x = ev.cls_mean_x;
    cluster.cls_mean_y = ev.cls_mean_y;
    cluster.cls_mean_z = ev.cls_mean_z;
    cluster.cls_mean_err_x = ev.cls_mean_x_err;
    cluster.cls_mean_err_y = ev.cls_mean_y_err;
    cluster.cls_idx = {-1,-1,-1};


    //values to be discuss and maybe change

}

