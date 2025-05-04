#include "./build/display.h"
#include "./build/stats.h"
#include "./build/LTrackerTrack.h"
#include <string>
#include <array> 
#include <map>
#include <iostream>
#include <cmath>
#include "TSystem.h" // Include ROOT's TSystem header
using namespace std;


void printKeysDirectly(const std::unordered_map<int, LCluster>& map) {
    std::cout << "Keys in the unordered_map: ";
    for (const auto& pair : map) {
        std::cout << pair.first << " "; // Print the key directly
    }
    std::cout << std::endl;
}

void run(){

    gSystem->Load("libHist");
    gSystem->Load("libGraf");
    gSystem->Load("libGraf3d");
    gSystem->Load("libMathCore");
    gSystem->Load("libRIO");

    gSystem->AddIncludePath("-Ibuild");
    gSystem->CompileMacro("./src/stats.cpp", "kg", "", "build");
    gSystem->CompileMacro("./src/LTrackerTrack.cpp", "kg", "", "build");
    gSystem->CompileMacro("./src/display.cpp", "kg", "", "build");

    gSystem->Load("build/LTrackerTrack_cpp.so");


    int events = 10;
    stats s(events);
    display d;
    LTrackerTrack tracker;
    d.draw_TR12();
    d.layers();
    d.tracks(events, false, tracker);
 
    //for each generated track build a cluster for each layer
    //use current algorithm to create the track

    cout << "stats \n" << s << endl;

    for (int i = 0; i < events; i++) {
        if (tracker.tidy_clusters_lay2.find(i) != tracker.tidy_clusters_lay2.end()) {
            std::cout << "Cluster: " << tracker.tidy_clusters_lay2.at(i) << std::endl;
        } else {
            std::cout << "No cluster found for key: " << i << std::endl;
        }
    }


}