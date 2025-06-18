#include "../include/display.h"
#include "../include/chip.h"
#include "../include/stats.h"
#include "../include/LTrackerTrack.h"
#include "../include/LTrackerCluster.h"
#include "../include/eventdata.h"
#include <string>
#include <array> 
#include <map>
#include <iostream>
#include <cmath>
#include "TSystem.h" // Include ROOT's TSystem header
#include "TApplication.h"
#include <thread>
#include <TStopwatch.h>
using namespace std;


void run(){

    TCanvas* real_tracks = new TCanvas("MC_tracks", "3D View_mc", 800, 600);
    TView* rt = TView::CreateView(1);
    rt->SetRange(-100, -100, 0, 100, 100, 70);
    rt->ShowAxis();

    int events = 2;
    stats s;
    display generated_tracks;
    chips cc;
    LTrackerTrack tracker;
    LTrackCandidate tc;
    
    //interfaccia rivelatore
    cc.print_all_chips(cc, real_tracks);
    generated_tracks.draw_TR12(real_tracks);

    //funzione MC
    generated_tracks.tracks(events, tracker, real_tracks);


    //algoritmo di ricostruzione tracce
    tracker.computeTracklets();
    tracker.computeTrackCandidates(real_tracks);
    tracker.printRecoTracks(real_tracks);
    

    

    

    cout << "stats \n" << s << endl;


}


int main(int argc, char** argv) {
    TApplication app("ROOT Application", &argc, argv);
    TStopwatch timer;       // crea il cronometro
    timer.Start();          // avvia il timer

    //take data from beam test
    /*
    TCanvas* can = new TCanvas("can", "3D View", 800, 600);
    TView* rt = TView::CreateView(1);
    rt->SetRange(-100, -100, 0, 100, 100, 70);
    rt->ShowAxis();
    eventdata e;
    e.takedata();
    e.print_data_on_canvas(can);
    stats s;
    cout << "stats \n" << s << endl;
    */

    

    //track simulation
    run();
    


    


    timer.Stop();           // ferma il timer

    std::cout << "Real time: " << timer.RealTime() << " s\n";
    std::cout << "CPU time:  " << timer.CpuTime()  << " s\n";

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