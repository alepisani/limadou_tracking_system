#include "../include/display.h"
#include "../include/stats.h"
#include "../include/LTrackerTrack.h"
#include "../include/LTrackerCluster.h"
#include <string>
#include <array> 
#include <map>
#include <iostream>
#include <cmath>
#include "TSystem.h" // Include ROOT's TSystem header
#include "TApplication.h"
#include <thread>
using namespace std;


void printKeysDirectly(const std::unordered_map<int, LCluster>& map) {
    std::cout << "Keys in the unordered_map: ";
    for (const auto& pair : map) {
        std::cout << pair.first << " "; // Print the key directly
    }
    std::cout << std::endl;
}

void run(){
    /*
    gSystem->Load("libHist");
    gSystem->Load("libGraf");
    gSystem->Load("libGraf3d");
    gSystem->Load("libMathCore");
    gSystem->Load("libRIO");

    gSystem->AddIncludePath("-Ibuild");
    gSystem->CompileMacro("./src/stats.cpp", "kg", "", "build");
    gSystem->CompileMacro("./src/LTrackerCluster.cpp", "kg", "", "build");
    gSystem->CompileMacro("./src/LTrackerTrack.cpp", "kg", "", "build");
    gSystem->CompileMacro("./src/display.cpp", "kg", "", "build");

    gSystem->Load("build/LTrackerTrack_cpp.so");
    */

    TCanvas* real_tracks = new TCanvas("MC_tracks", "3D View_mc", 800, 600);
    TView* rt = TView::CreateView(1);
    rt->SetRange(-100, -100, 0, 100, 100, 70);
    rt->ShowAxis();

    TCanvas* recon_tracks = new TCanvas("reco_tracks", "3D View_recot", 800, 600);
    TView* recot = TView::CreateView(1);
    recot->SetRange(-100, -100, 0, 100, 100, 70);
    recot->ShowAxis();

    int events = 50;
    stats s(events);
    
    display generated_tracks;
    LTrackerTrack tracker;
    generated_tracks.draw_TR12(real_tracks);
    generated_tracks.layers(real_tracks);
    generated_tracks.tracks(events, false, tracker, real_tracks);

    display reco_tracks;
    reco_tracks.draw_TR12(recon_tracks);
    reco_tracks.layers(recon_tracks);    
    tracker.computeTracklets();
    tracker.computeTrackCandidates(recon_tracks);

    real_tracks->Draw();
    real_tracks->Update();
    recon_tracks->Draw();
    recon_tracks->Update();

    cout << "stats \n" << s << endl;
}


/*
    for (int i = 0; i < events; i++) {
        if (tracker.tidy_clusters_lay2.find(i) != tracker.tidy_clusters_lay2.end()) {
            std::cout << "Cluster: " << tracker.tidy_clusters_lay2.at(i) << std::endl;
        } else {
            std::cout << "No cluster found for key: " << i << std::endl;
        }
    }
    */



int main(int argc, char** argv) {
    TApplication app("ROOT Application", &argc, argv);

    run(); 

    std::cout << "press ENTER to close everything" << std::endl;

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