#include "./build/display.h"
#include <string>
#include <array> 
#include <map>
#include <iostream>
#include <cmath>
#include "TSystem.h" // Include ROOT's TSystem header
using namespace std;

void run(){

    gSystem->AddIncludePath("-Ibuild");
    gSystem->CompileMacro("./src/display.cxx", "kg", "", "build");

    display d;
    d.draw_TR12();
    d.layers();

    int events = 10;
    d.tracks(events, true);



}