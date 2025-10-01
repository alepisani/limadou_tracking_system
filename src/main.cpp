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
#include "TH1.h"
#include "TROOT.h"
using namespace std;

void run(int *events)
{

    TCanvas *real_tracks = new TCanvas("MC_tracks", "3D View_mc", 800, 600);
    TView *rt = TView::CreateView(1);
    rt->SetRange(-100, -100, 0, 100, 100, 70);
    rt->ShowAxis();

    stats s;
    display gen_tracks;
    chips cc;
    LTrackerTrack tracker;
    LTrackCandidate tc;

    // interfaccia rivelatore
    cc.print_all_chips(cc, real_tracks);
    gen_tracks.draw_TR12(real_tracks);

    // funzione MC
    gen_tracks.take_distributions();
    gen_tracks.tracks(events, tracker, real_tracks);
    // for(int i = 0; i < gen_tracks.generated_tracks.size(); ++i){
    //     double theta = gen_tracks.generated_tracks[i].theta * TMath::RadToDeg();
    //     double phi = gen_tracks.generated_tracks[i].phi * TMath::RadToDeg();
    //     printf("x0 = %f, y0 = %f, theta_real = %f, phi_real = %f\n", gen_tracks.generated_tracks[i].x0, gen_tracks.generated_tracks[i].y0, theta, phi);
    //  }
    //  cout << "-----------------------" << endl;

    // algoritmo di ricostruzione tracce
    TStopwatch t;
    t.Start();
    tracker.computeTracklets();
    // tracker.computeTrackCandidates();
    tracker.new_algo(0.8);
    t.Stop();

    tracker.printRecoTracks_new_alg(real_tracks);

    // cout << "-+-+-+-+-+-+-+-+-+-+-+-+-+-+" << endl;
    // cout << "Real time: " << t.RealTime() << " s\n";
    // cout << "CPU time:  " << t.CpuTime()  << " s\n";
    // cout << "-+-+-+-+-+-+-+-+-+-+-+-+-+-+" << endl;

    // cout << "stats \n" << s << endl;
}

int main(int argc, char **argv)
{
    TApplication app("ROOT Application", &argc, argv);
    TH1::AddDirectory(false);

    // track simulation
    // int *events;
    // int ev = 30;
    // events = &ev;
    // run(events);

    // simulations sim;
    // sim.sim_only_trk_3L(1000);
    // sim.sim_old_algo(100);
    // sim.sim_trk_32L(10000);

    // reco from MUONS
    eventdata e;
    e.analize_data();

    /*
    //compute theta max
    int xtr1, ytr1, ztr1;
    int xtr2, ytr2, ztr2;
    ztr1 = display::TR1CenterZ + display::TR1Size[2];
    ztr2 = display::TR2CenterZ;
    xtr1 = -0.5 * display::TR1Size[0];
    xtr2 = 2 * display::TR2Size[0] + 1.5 * display::TR2GapX;
    ytr1 = -2.5 * display::TR1Size[1] - 2 * display::TR1GapY;
    ytr2 = 0.5 * display::TR2Size[1];
    double delta_y = ytr2 - ytr1;
    double delta_x = xtr2 - xtr1;
    double delta_z = ztr2 - ztr1;
    double delta_r = TMath::Sqrt(delta_y * delta_y + delta_x * delta_x);
    double r = TMath::Sqrt(delta_x * delta_x + delta_y * delta_y + delta_z * delta_z);
    double phi = TMath::ATan2(delta_y, delta_x);
    double theta = TMath::ACos(delta_z / r);
    double radtodeg = TMath::RadToDeg();
    double thera_deg = theta * radtodeg;
    printf("theta max: %f\n", thera_deg);
    // theta max: 75.274958
    */

    // Thread secondario per leggere l'input senza bloccare l'interfaccia grafica
    std::thread inputThread([]()
                            {
                                std::cin.get();   // Attende INVIO
                                gSystem->Exit(0); // Termina il loop di app.Run()
                            });

    // Avvia il loop degli eventi di ROOT (canvas interattivo)
    app.Run(); // Bloccante ma necessario per l'interazione

    inputThread.join(); // Attende la fine del thread di input

    return 0;
}