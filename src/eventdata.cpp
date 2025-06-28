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
using namespace std;

eventdata::eventdata(){}

std::unordered_map<int, eventdata> alldata;


std::ostream &operator<<(std::ostream &output, const eventdata &ev) {
    for (size_t i = 0; i < ev.DIR_xpos.size(); i++) {
        // output << "turret_idx: " << (int)ev.DIR_turret_idx[i] << " ";
        // output << "stave_idx: " << (int)ev.DIR_stave_idx[i] << " ";
        // output << "chip_idx0: " << (int)ev.DIR_chip_idx0[i] << " ";
        // output << "chip_idx1: " << (int)ev.DIR_chip_idx1[i] << " ";
        // output << "chip_id:   " << ev.DIR_chip_id[i] << " ";
        output << "xpos:      " << ev.DIR_xpos[i] << " || ";
        output << "ypos:      " << ev.DIR_ypos[i] << " || ";
        output << "zpos:      " << ev.DIR_zpos[i] << " || ";
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

    // Carica il TTree (sostituisci con il nome corretto se diverso)
    TTree *tree = (TTree*)file->Get("L2;1");
    if (!tree) {
        std::cerr << "TTree non trovato nel file\n";
        return;
    }

    // Oggetti temporanei per collegare i branch
    std::vector<unsigned char> *DIR_turret_idx = nullptr;
    std::vector<unsigned char> *DIR_stave_idx = nullptr;
    std::vector<unsigned char> *DIR_chip_idx0 = nullptr;
    std::vector<unsigned char> *DIR_chip_idx1 = nullptr;
    std::vector<unsigned int>  *DIR_chip_id = nullptr;
    std::vector<float> *DIR_xpos = nullptr;
    std::vector<float> *DIR_ypos = nullptr;
    std::vector<float> *DIR_zpos = nullptr;
    std::vector<int> *DIR_cls_idx = nullptr;

    tree->SetBranchAddress("DIR_turret_idx", &DIR_turret_idx);
    tree->SetBranchAddress("DIR_stave_idx",  &DIR_stave_idx);
    tree->SetBranchAddress("DIR_chip_idx0",  &DIR_chip_idx0);
    tree->SetBranchAddress("DIR_chip_idx1",  &DIR_chip_idx1);
    tree->SetBranchAddress("DIR_chip_id",    &DIR_chip_id);
    tree->SetBranchAddress("DIR_xpos",       &DIR_xpos);
    tree->SetBranchAddress("DIR_ypos",       &DIR_ypos);
    tree->SetBranchAddress("DIR_zpos",       &DIR_zpos);
    tree->SetBranchAddress("DIR_cls_idx",    &DIR_cls_idx);

    Long64_t nEntries = tree->GetEntries();

    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);
        eventdata ev;
        ev.DIR_turret_idx = *DIR_turret_idx;
        ev.DIR_stave_idx  = *DIR_stave_idx;
        ev.DIR_chip_idx0  = *DIR_chip_idx0;
        ev.DIR_chip_idx1  = *DIR_chip_idx1;
        ev.DIR_chip_id    = *DIR_chip_id;
        ev.DIR_xpos       = *DIR_xpos;
        ev.DIR_ypos       = *DIR_ypos;
        ev.DIR_zpos       = *DIR_zpos;
        ev.DIR_cls_idx    = *DIR_cls_idx;
        alldata[i] = ev;

        
    }
    file->Close();
}


void eventdata::print_data_on_canvas(TCanvas* can){

    can->cd();
    display d;
    d.draw_TR12(can);
    chips c;
    c.print_all_chips(c, can);

    int events = alldata.size();
    //int events = 20;
    TH1F* hx = new TH1F("x", "x;x;counts", events, -display::TR2Size[0]*2.5, display::TR2Size[0]*2.5);
    TH1F* hy = new TH1F("y", "y;y;counts", events, -100, 100);
    TH1F* hz = new TH1F("z", "z;z;counts", events, 10, 40);

    eventdata ev;
    LTrackerCluster cluster;
    for (int i = 0; i < events; i++){
        ev = alldata[i];
        cluster.CalculateClusterPosition(ev);
    }    

    cout << "quanti cluster: " << cluster.cls_mean_x.size() << endl;
    cout << "events: " << events << endl;
    for(int i=0; i<cluster.cls_mean_x.size(); i++){
        //cout << "p: " << cluster.cls_mean_x[i] << " , y: " << cluster.cls_mean_y[i] << "  , z: " << cluster.cls_mean_z[i] << endl;
        TMarker3DBox *p = new TMarker3DBox(cluster.cls_mean_x[i],cluster.cls_mean_y[i],cluster.cls_mean_z[i],2,2,0,0,0);
        p->SetLineColor(kRed);
        p->SetLineWidth(3);
        p->Draw();
        hx->Fill(cluster.cls_mean_x[i]);
        hy->Fill(cluster.cls_mean_y[i]);
        hz->Fill(cluster.cls_mean_z[i]);

        //stats
        if(cluster.cls_mean_z[i] < 36 && cluster.cls_mean_z[i] > 30){stats::hmthL2++;}
        if(cluster.cls_mean_z[i] < 29 && cluster.cls_mean_z[i] > 22){stats::hmthL1++;}
        if(cluster.cls_mean_z[i] < 20 && cluster.cls_mean_z[i] > 15){stats::hmthL0++;}
    }
    stats::hmbh0L = events - stats::hmbh3L - stats::hmbh2L - stats::hmbh1L; 




    char file[200];
    sprintf(file,"../data_beam_test/data.root");
    TFile f(file,"RECREATE");
    hx->Write();
    hy->Write();
    hz->Write();
    
    f.Close();
}


