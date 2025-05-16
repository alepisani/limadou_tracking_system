#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include "../include/eventdata.h"
using namespace std;

eventdata::eventdata(){}

std::unordered_map<int, eventdata> alldata;

void eventdata::takedata(){

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



    //printiamo un evento a caso
    eventdata evento;
    evento = alldata[300];
    for(int i = 0; i < evento.DIR_xpos.size(); i++){
        cout << "z_pos:     " << evento.DIR_zpos[i] << endl;
    }
    cout << "alldata size: " << alldata.size() << endl;
    

}




