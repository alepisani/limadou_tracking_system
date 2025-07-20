#include "../include/display.h"
#include "../include/chip.h"
#include "../include/stats.h"
#include "../include/LTrackerTrack.h"
#include "../include/LTrackerCluster.h"
#include "../include/eventdata.h"
#include "../include/simulations.h"
#include <string>
#include <array> 
#include <map>
#include <iostream>
#include <cmath>
#include <TTree.h>
#include "TFile.h"
#include "TSystem.h" // Include ROOT's TSystem header
#include "TApplication.h"
#include <thread>
#include <TStopwatch.h>
using namespace std;


void run(int events){

    TCanvas* real_tracks = new TCanvas("MC_tracks", "3D View_mc", 800, 600);
    TView* rt = TView::CreateView(1);
    rt->SetRange(-100, -100, 0, 100, 100, 70);
    rt->ShowAxis();

    
    stats s;
    display generated_tracks;
    chips cc;
    LTrackerTrack tracker;
    LTrackCandidate tc;
    
    //interfaccia rivelatore
    cc.print_all_chips(cc, real_tracks);
    generated_tracks.draw_TR12(real_tracks);

    //funzione MC
    generated_tracks.take_distributions();
    generated_tracks.tracks(events, tracker, real_tracks);


    //algoritmo di ricostruzione tracce
    TStopwatch t;
    t.Start();
    tracker.computeTracklets();
    //tracker.computeTrackCandidates();
    tracker.new_computing();
    t.Stop();

    //tracker.printRecoTracks_old_alg(real_tracks, events);
    tracker.printRecoTracks_new_alg(real_tracks, events);

    cout << "-+-+-+-+-+-+-+-+-+-+-+-+-+-+" << endl;
    cout << "Real time: " << t.RealTime() << " s\n";
    cout << "CPU time:  " << t.CpuTime()  << " s\n";
    cout << "-+-+-+-+-+-+-+-+-+-+-+-+-+-+" << endl;


    

    cout << "stats \n" << s << endl;


}


int main(int argc, char** argv) {
    TApplication app("ROOT Application", &argc, argv);

    //take data from beam test
    
    /* TCanvas* can = new TCanvas("can", "3D View", 800, 600);
    TView* rt = TView::CreateView(1);
    rt->SetRange(-100, -100, 0, 100, 100, 70);
    rt->ShowAxis();
    eventdata e;
    e.takedata();
    e.print_data_on_canvas(can);
    stats s;
    cout << "stats \n" << s << endl; */
    

    

    //track simulation
    //int events = 10;
    //run(events);

    simulations sim;
    sim.sim_only_trk_3L(100);


     
    //prove varie
    //TFile *file = TFile::Open("../data_beam_test/TEST_MUONS_m_MAIN_1000.0MeV_-999.0deg_-0.05V_boot207_run510_L2.root");
    /* 
    
    TTree *tree = (TTree*)file->Get("L2;1");
    TTree *tree = (TTree*)file->Get("Tmd;1");
    tree->Print();

    
    Long64_t nentries = tree->GetEntries();
    Bool_t trig_conf_flag[6];
    tree->SetBranchAddress("trig_conf_flag[6]", &trig_conf_flag);
    for (Long64_t i=0; i<nentries; ++i) {
    //for (Long64_t i=10000; i<10500; ++i) {
        tree->GetEntry(i);
        cout << "Evento " << i << ": ";
        //cout << trigger_mask;

        for (int j=0; j<6; ++j) {
            cout << "flag[" << j << "]=" << trig_conf_flag[j] << " ";
        }


        cout << endl;
    } 

    cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;

      
    UInt_t trigger_mask;
    TTree *tr = (TTree*)file->Get("Tmd;1");
    Long64_t nentr = tr->GetEntries();
    tr->SetBranchAddress("trigger_mask", &trigger_mask);
    for (Long64_t i=0; i<nentr; ++i) {
        tr->GetEntry(i);
        cout << "Evento " << i << ": ";
        cout << trigger_mask;
        cout << endl;
    }
     */
    






    // Thread secondario per leggere l'input senza bloccare l'interfaccia grafica
    std::thread inputThread([]() {
        std::cin.get(); // Attende INVIO
        gSystem->Exit(0); // Termina il loop di app.Run()
    });

    // Avvia il loop degli eventi di ROOT (canvas interattivo)
    app.Run(); // Bloccante ma necessario per l'interazione

    inputThread.join(); // Attende la fine del thread di input

    return 0;
}