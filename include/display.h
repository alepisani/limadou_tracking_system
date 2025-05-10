#ifndef DISPLAY_H
#define DISPLAY_H

#include <string>
#include <array> 
#include <map>
#include <vector>
#include <iostream>
#include <cmath>
#include <TCanvas.h> // Include ROOT's TCanvas header
#include <TView.h>   // Include ROOT's TView header
#include "TMath.h"
#include "./LTrackerTrack.h"
using namespace std;


class display{
public:

    //all measure in mm
    
    //trigger 1 (bottom trigger)
    static constexpr double TR1CenterZ = 0.;
    static constexpr double TR1Thickness = 2;
    static constexpr std::array<double, 3> TR1Size = {154.6,32.5,TR1Thickness};
    static constexpr double TR1GapY = 1.9;
    //trigger 2 (top trigger)
    static constexpr double TR2CenterZ = TR1CenterZ+60.5;
    static constexpr double TR2Thickness = 8;
    static constexpr std::array<double, 3> TR2Size = {36,150,TR2Thickness};
    static constexpr double TR2GapX = 2; 
    //layers
    static constexpr int PixelNCols = 1024;
    static constexpr int PixelNRows = 512;
    static constexpr double PixelSizeCols = 0.02924; 
    static constexpr double PixelSizeRows = 0.02688;
    static constexpr double ChipSizeX = PixelSizeCols*PixelNCols;
    static constexpr double ChipSizeY = PixelSizeRows*PixelNRows;
    static constexpr double ChipSizeZ = 0.050;     //50microm
    static constexpr double ChipDistanceX = 0.150;
    static constexpr double ChipDistanceY = 0.150;
    static constexpr double ChipStaveDistanceY = 7.22312;
    static constexpr std::array<double, 3> StaveZ = {17.825+TR1CenterZ,17.825+8.5+TR1CenterZ,17.825+17+TR1CenterZ};
    static constexpr float pitch_x = 29.24 / 1000; // [mm] pixel pitch
    static constexpr float pitch_y = 26.88 / 1000; // [mm] pixel pitch
    double err_cl = 10;
    double pi = TMath::Pi();
    static constexpr float dist_z = 8.5;                      // [mm] distance between planes
    static constexpr float shift_z = 17.825;                  // [mm] distance between trigger layer and firts tracking layer
    static constexpr float z_origin_shift = 26.325; // [mm] distance between origin of the reference point for tracks and firts tracking layer

    void draw_TR12(TCanvas* geom);
    void layers(TCanvas* geom);
    void tracks(int, LTrackerTrack&, TCanvas* geom);   //true p-q  //false p-theta-phi
    

    display();

};


#endif




