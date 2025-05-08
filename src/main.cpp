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


void run(){

    TCanvas* real_tracks = new TCanvas("MC_tracks", "3D View_mc", 800, 600);
    TView* rt = TView::CreateView(1);
    rt->SetRange(-100, -100, 0, 100, 100, 70);
    rt->ShowAxis();

    //TCanvas* recon_tracks = new TCanvas("reco_tracks", "3D View_recot", 800, 600);
    //TView* recot = TView::CreateView(1);
    //recot->SetRange(-100, -100, 0, 100, 100, 70);
    //recot->ShowAxis();

    int events = 10;
    stats s(events);
    
    display generated_tracks;
    LTrackerTrack tracker;
    generated_tracks.draw_TR12(real_tracks);
    generated_tracks.layers(real_tracks);
    generated_tracks.tracks(events, tracker, real_tracks);

    cout << "stats \n" << s << endl;
}


int main(int argc, char** argv) {
    TApplication app("ROOT Application", &argc, argv);

    run(); 

    std::cout << "press ENTER to close" << std::endl;

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