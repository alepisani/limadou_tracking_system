#include "./build/display.h"
#include "./build/stats.h"
#include <string>
#include <array> 
#include <map>
#include <iostream>
#include <cmath>
#include "TSystem.h" // Include ROOT's TSystem header
using namespace std;

void run(){

    gSystem->Load("libHist");
    gSystem->Load("libGraf");
    gSystem->Load("libGraf3d");
    gSystem->Load("libMathCore");
    gSystem->Load("libRIO");

    gSystem->AddIncludePath("-Ibuild");
    gSystem->CompileMacro("./src/stats.cpp", "kg", "", "build");
    gSystem->CompileMacro("./src/display.cpp", "kg", "", "build");


    stats s;
    display d;
    d.draw_TR12();
    d.layers();

    int events = 10;
    d.tracks(events, true);
 
    //for each generated track build a cluster for each layer
    //use current algorithm to create the track

    cout << "stats \n" << s << endl;



}