#include <string>
#include <array> 
#include <map>
#include <iostream>
#include <cmath>
#include <TApplication.h>
#include <TCanvas.h>
#include <TView.h>
#include <TList.h>
#include <TPolyLine3D.h>
#include "TH1F.h"
#include "TLine.h"
#include "TRandom3.h"
#include "TMarker3DBox.h"
#include "TFile.h"
#include "../include/display.h"
#include "../include/stats.h"
#include "../include/chip.h"

using namespace std;

chip::chip(){}

bool chip::check_chip_status(chip& c){
    return c.status;
}

void chip::print_chip(chip& c, TCanvas* canvas){
    //TMarker3DBox *chip = new TMarker3DBox(a,b,c,d,e,f,g,h);
    //a,b,c coordinate del centro della box
    //metà delle dimensioni della scatola --> a -+ d / b -+ e / c -+ f
    //g,h angoli di orientazione della box, sempre 0,0
    canvas->cd();
    if(check_chip_status(c)){
        TMarker3DBox *chip = new TMarker3DBox(c.x, c.y, c.z, chip::x_dim/2, chip::y_dim/2, chip::z_dim/2, 0,0);
        chip->SetLineWidth(2);
        chip->SetLineColor(kCyan);
        chip->Draw();   
    }
    if(!check_chip_status(c)){
        cout << "WARNING! Il seguente chip risulta spento, id: " << c.id << endl;
    }
}

void chip::set_chip_coordinates(chip& c, double x, double y, double z, unsigned short id){
    c.x = x;
    c.y = y;
    c.z = z;
    c.id = id;
}


chips::chips(){}

chips::chips(const std::array<unsigned short, 150> ChipIds, const std::array<TVector3, 150> chip_coordinates){
    //fill che id_coordinates vector with appropriate values
    for (size_t i = 0; i < ChipIds.size(); ++i) {
        unsigned short chipId = ChipIds[i];
        const TVector3& coordinates = chip_coordinates[i];

        // Aggiungi il TVector3 alla lista (vector) associata al ChipId
        id_coordinates[chipId].push_back(coordinates);
    }
}

void chips::print_all_chips(chips& c, TCanvas* canvas){
    for(int i = 0; i< ChipIds.size(); i++){
        const TVector3& coordinates = chip_coordinates[i];
        double x = coordinates[0];
        double y = coordinates[1];
        double z = coordinates[2];
        chip chip_instance;
        chip_instance.id = ChipIds[i];
        chip_instance.set_chip_coordinates(chip_instance,x,y,z, ChipIds[i]);
        if(!chip_instance.is_chip_dead(chip_instance, chips::dead_chip)){
            chip_instance.print_chip(chip_instance,canvas);            //chip is alive --> print it
        }
        if(chip_instance.is_chip_dead(chip_instance, chips::dead_chip)){
            cout << "WARNING! THE FOLLOWING CHIP IS TURNED OFF, id: " << ChipIds[i] << endl;
        }
    }
}

bool chip::is_chip_dead(chip& c, const std::vector<unsigned short> dead_chip){
    for(int i = 0; i < dead_chip.size(); i++){
        if(c.id == dead_chip[i]){
            return true;                //chip dead --> function return true
        }    
    }
    return false;               //chip alive --> function return false, means not dead
}

const std::vector<unsigned short> chips::dead_chip = {0x07c,0x07b,0x07a,0x079,0x078,
    0x070,0x071,0x072,0x073,0x074};

const std::array<unsigned short, 150> chips::ChipIds = {
    //50 chip per piano
    // top layer
    0x07c,0x07b,0x07a,0x079,0x078,
    0x070,0x071,0x072,0x073,0x074,
    0x37c,0x37b,0x37a,0x379,0x378,
    0x370,0x371,0x372,0x373,0x374,
    0x67c,0x67b,0x67a,0x679,0x678,
    0x670,0x671,0x672,0x673,0x674,
    0x97c,0x97b,0x97a,0x979,0x978,
    0x970,0x971,0x972,0x973,0x974,
    0xc7c,0xc7b,0xc7a,0xc79,0xc78,
    0xc70,0xc71,0xc72,0xc73,0xc74,
    // middle layer
    0x17c,0x17b,0x17a,0x179,0x178,
    0x170,0x171,0x172,0x173,0x174,
    0x47c,0x47b,0x47a,0x479,0x478,
    0x470,0x471,0x472,0x473,0x474,
    0x77c,0x77b,0x77a,0x779,0x778,
    0x770,0x771,0x772,0x773,0x774,
    0xa7c,0xa7b,0xa7a,0xa79,0xa78,
    0xa70,0xa71,0xa72,0xa73,0xa74,
    0xd7c,0xd7b,0xd7a,0xd79,0xd78,
    0xd70,0xd71,0xd72,0xd73,0xd74,
    // bottom layer
    0x27c,0x27b,0x27a,0x279,0x278,
    0x270,0x271,0x272,0x273,0x274,
    0x57c,0x57b,0x57a,0x579,0x578,
    0x570,0x571,0x572,0x573,0x574,
    0x87c,0x87b,0x87a,0x879,0x878,
    0x870,0x871,0x872,0x873,0x874,
    0xb7c,0xb7b,0xb7a,0xb79,0xb78,
    0xb70,0xb71,0xb72,0xb73,0xb74,
    0xe7c,0xe7b,0xe7a,0xe79,0xe78,
    0xe70,0xe71,0xe72,0xe73,0xe74
    };

const std::array<TVector3, 150> chips::chip_coordinates = {
    // top layer (layer2)
    // row 0
    TVector3(-2*(display::ChipSizeX + display::ChipDistanceX), -(2*display::ChipStaveDistanceY + 2.5*display::ChipDistanceY + 4.5*display::ChipSizeY), display::StaveZ[2]),
    TVector3(-1*(display::ChipSizeX + display::ChipDistanceX), -(2*display::ChipStaveDistanceY + 2.5*display::ChipDistanceY + 4.5*display::ChipSizeY), display::StaveZ[2]),
    TVector3( 0, -(2*display::ChipStaveDistanceY + 2.5*display::ChipDistanceY + 4.5*display::ChipSizeY), display::StaveZ[2]),
    TVector3( 1*(display::ChipSizeX + display::ChipDistanceX), -(2*display::ChipStaveDistanceY + 2.5*display::ChipDistanceY + 4.5*display::ChipSizeY), display::StaveZ[2]),
    TVector3( 2*(display::ChipSizeX + display::ChipDistanceX), -(2*display::ChipStaveDistanceY + 2.5*display::ChipDistanceY + 4.5*display::ChipSizeY), display::StaveZ[2]),
    // row 1
    TVector3(-2*(display::ChipSizeX + display::ChipDistanceX), -(2*display::ChipStaveDistanceY + 1.5*display::ChipDistanceY + 3.5*display::ChipSizeY), display::StaveZ[2]),
    TVector3(-1*(display::ChipSizeX + display::ChipDistanceX), -(2*display::ChipStaveDistanceY + 1.5*display::ChipDistanceY + 3.5*display::ChipSizeY), display::StaveZ[2]),
    TVector3( 0, -(2*display::ChipStaveDistanceY + 1.5*display::ChipDistanceY + 3.5*display::ChipSizeY), display::StaveZ[2]),
    TVector3( 1*(display::ChipSizeX + display::ChipDistanceX), -(2*display::ChipStaveDistanceY + 1.5*display::ChipDistanceY + 3.5*display::ChipSizeY), display::StaveZ[2]),
    TVector3( 2*(display::ChipSizeX + display::ChipDistanceX), -(2*display::ChipStaveDistanceY + 1.5*display::ChipDistanceY + 3.5*display::ChipSizeY), display::StaveZ[2]),
    // row 2
    TVector3(-2*(display::ChipSizeX + display::ChipDistanceX), -(display::ChipStaveDistanceY + 1.5*display::ChipDistanceY + 2.5*display::ChipSizeY), display::StaveZ[2]),
    TVector3(-1*(display::ChipSizeX + display::ChipDistanceX), -(display::ChipStaveDistanceY + 1.5*display::ChipDistanceY + 2.5*display::ChipSizeY), display::StaveZ[2]),
    TVector3( 0, -(display::ChipStaveDistanceY + 1.5*display::ChipDistanceY + 2.5*display::ChipSizeY), display::StaveZ[2]),
    TVector3( 1*(display::ChipSizeX + display::ChipDistanceX), -(display::ChipStaveDistanceY + 1.5*display::ChipDistanceY + 2.5*display::ChipSizeY), display::StaveZ[2]),
    TVector3( 2*(display::ChipSizeX + display::ChipDistanceX), -(display::ChipStaveDistanceY + 1.5*display::ChipDistanceY + 2.5*display::ChipSizeY), display::StaveZ[2]),
    // row 3
    TVector3(-2*(display::ChipSizeX + display::ChipDistanceX), -(display::ChipStaveDistanceY + 1.5*display::ChipDistanceY + 1.5*display::ChipSizeY), display::StaveZ[2]),
    TVector3(-1*(display::ChipSizeX + display::ChipDistanceX), -(display::ChipStaveDistanceY + 1.5*display::ChipDistanceY + 1.5*display::ChipSizeY), display::StaveZ[2]),
    TVector3( 0, -(display::ChipStaveDistanceY + 1.5*display::ChipDistanceY + 1.5*display::ChipSizeY), display::StaveZ[2]),
    TVector3( 1*(display::ChipSizeX + display::ChipDistanceX), -(display::ChipStaveDistanceY + 1.5*display::ChipDistanceY + 1.5*display::ChipSizeY), display::StaveZ[2]),
    TVector3( 2*(display::ChipSizeX + display::ChipDistanceX), -(display::ChipStaveDistanceY + 1.5*display::ChipDistanceY + 1.5*display::ChipSizeY), display::StaveZ[2]),
    // row 4
    TVector3(-2*(display::ChipSizeX + display::ChipDistanceX), -0.5*(display::ChipDistanceY + display::ChipSizeY), display::StaveZ[2]),
    TVector3(-1*(display::ChipSizeX + display::ChipDistanceX), -0.5*(display::ChipDistanceY + display::ChipSizeY), display::StaveZ[2]),
    TVector3( 0, -0.5*(display::ChipDistanceY + display::ChipSizeY), display::StaveZ[2]),
    TVector3( 1*(display::ChipSizeX + display::ChipDistanceX), -0.5*(display::ChipDistanceY + display::ChipSizeY), display::StaveZ[2]),
    TVector3( 2*(display::ChipSizeX + display::ChipDistanceX), -0.5*(display::ChipDistanceY + display::ChipSizeY), display::StaveZ[2]),
    // row 5
    TVector3(-2*(display::ChipSizeX + display::ChipDistanceX), +0.5*(display::ChipDistanceY + display::ChipSizeY), display::StaveZ[2]),
    TVector3(-1*(display::ChipSizeX + display::ChipDistanceX), +0.5*(display::ChipDistanceY + display::ChipSizeY), display::StaveZ[2]),
    TVector3( 0, +0.5*(display::ChipDistanceY + display::ChipSizeY), display::StaveZ[2]),
    TVector3( 1*(display::ChipSizeX + display::ChipDistanceX), +0.5*(display::ChipDistanceY + display::ChipSizeY), display::StaveZ[2]),
    TVector3( 2*(display::ChipSizeX + display::ChipDistanceX), +0.5*(display::ChipDistanceY + display::ChipSizeY), display::StaveZ[2]),
    // row 6
    TVector3(-2*(display::ChipSizeX + display::ChipDistanceX), +(display::ChipStaveDistanceY + 1.5*display::ChipDistanceY + 1.5*display::ChipSizeY), display::StaveZ[2]),
    TVector3(-1*(display::ChipSizeX + display::ChipDistanceX), +(display::ChipStaveDistanceY + 1.5*display::ChipDistanceY + 1.5*display::ChipSizeY), display::StaveZ[2]),
    TVector3( 0, +(display::ChipStaveDistanceY + 1.5*display::ChipDistanceY + 1.5*display::ChipSizeY), display::StaveZ[2]),
    TVector3( 1*(display::ChipSizeX + display::ChipDistanceX), +(display::ChipStaveDistanceY + 1.5*display::ChipDistanceY + 1.5*display::ChipSizeY), display::StaveZ[2]),
    TVector3( 2*(display::ChipSizeX + display::ChipDistanceX), +(display::ChipStaveDistanceY + 1.5*display::ChipDistanceY + 1.5*display::ChipSizeY), display::StaveZ[2]),
    // row 7
    TVector3(-2*(display::ChipSizeX + display::ChipDistanceX), +(display::ChipStaveDistanceY + 1.5*display::ChipDistanceY + 2.5*display::ChipSizeY), display::StaveZ[2]),
    TVector3(-1*(display::ChipSizeX + display::ChipDistanceX), +(display::ChipStaveDistanceY + 1.5*display::ChipDistanceY + 2.5*display::ChipSizeY), display::StaveZ[2]),
    TVector3( 0, +(display::ChipStaveDistanceY + 1.5*display::ChipDistanceY + 2.5*display::ChipSizeY), display::StaveZ[2]),
    TVector3( 1*(display::ChipSizeX + display::ChipDistanceX), +(display::ChipStaveDistanceY + 1.5*display::ChipDistanceY + 2.5*display::ChipSizeY), display::StaveZ[2]),
    TVector3( 2*(display::ChipSizeX + display::ChipDistanceX), +(display::ChipStaveDistanceY + 1.5*display::ChipDistanceY + 2.5*display::ChipSizeY), display::StaveZ[2]),
    // row 8
    TVector3(-2*(display::ChipSizeX + display::ChipDistanceX), +(2*display::ChipStaveDistanceY + 1.5*display::ChipDistanceY + 3.5*display::ChipSizeY), display::StaveZ[2]),
    TVector3(-1*(display::ChipSizeX + display::ChipDistanceX), +(2*display::ChipStaveDistanceY + 1.5*display::ChipDistanceY + 3.5*display::ChipSizeY), display::StaveZ[2]),
    TVector3( 0, +(2*display::ChipStaveDistanceY + 1.5*display::ChipDistanceY + 3.5*display::ChipSizeY), display::StaveZ[2]),
    TVector3( 1*(display::ChipSizeX + display::ChipDistanceX), +(2*display::ChipStaveDistanceY + 1.5*display::ChipDistanceY + 3.5*display::ChipSizeY), display::StaveZ[2]),
    TVector3( 2*(display::ChipSizeX + display::ChipDistanceX), +(2*display::ChipStaveDistanceY + 1.5*display::ChipDistanceY + 3.5*display::ChipSizeY), display::StaveZ[2]),
    // row 9
    TVector3(-2*(display::ChipSizeX + display::ChipDistanceX), +(2*display::ChipStaveDistanceY + 2.5*display::ChipDistanceY + 4.5*display::ChipSizeY), display::StaveZ[2]),
    TVector3(-1*(display::ChipSizeX + display::ChipDistanceX), +(2*display::ChipStaveDistanceY + 2.5*display::ChipDistanceY + 4.5*display::ChipSizeY), display::StaveZ[2]),
    TVector3( 0, +(2*display::ChipStaveDistanceY + 2.5*display::ChipDistanceY + 4.5*display::ChipSizeY), display::StaveZ[2]),
    TVector3( 1*(display::ChipSizeX + display::ChipDistanceX), +(2*display::ChipStaveDistanceY + 2.5*display::ChipDistanceY + 4.5*display::ChipSizeY), display::StaveZ[2]),
    TVector3( 2*(display::ChipSizeX + display::ChipDistanceX), +(2*display::ChipStaveDistanceY + 2.5*display::ChipDistanceY + 4.5*display::ChipSizeY), display::StaveZ[2]),


    //middle layer
    //row0
    TVector3(-2*(display::ChipSizeX + display::ChipDistanceX), -(2*display::ChipStaveDistanceY + 2.5*display::ChipDistanceY + 4.5*display::ChipSizeY), display::StaveZ[1]),
    TVector3(-1*(display::ChipSizeX + display::ChipDistanceX), -(2*display::ChipStaveDistanceY + 2.5*display::ChipDistanceY + 4.5*display::ChipSizeY), display::StaveZ[1]),
    TVector3(0, -(2*display::ChipStaveDistanceY + 2.5*display::ChipDistanceY + 4.5*display::ChipSizeY), display::StaveZ[1]),
    TVector3(1*(display::ChipSizeX + display::ChipDistanceX), -(2*display::ChipStaveDistanceY + 2.5*display::ChipDistanceY + 4.5*display::ChipSizeY), display::StaveZ[1]),
    TVector3(2*(display::ChipSizeX + display::ChipDistanceX), -(2*display::ChipStaveDistanceY + 2.5*display::ChipDistanceY + 4.5*display::ChipSizeY), display::StaveZ[1]),
    //row1
    TVector3(-2*(display::ChipSizeX + display::ChipDistanceX), -(2*display::ChipStaveDistanceY + 1.5*display::ChipDistanceY + 3.5*display::ChipSizeY), display::StaveZ[1]),
    TVector3(-(display::ChipSizeX + display::ChipDistanceX), -(2*display::ChipStaveDistanceY + 1.5*display::ChipDistanceY + 3.5*display::ChipSizeY), display::StaveZ[1]),
    TVector3(0, -(2*display::ChipStaveDistanceY + 1.5*display::ChipDistanceY + 3.5*display::ChipSizeY), display::StaveZ[1]),
    TVector3((display::ChipSizeX + display::ChipDistanceX), -(2*display::ChipStaveDistanceY + 1.5*display::ChipDistanceY + 3.5*display::ChipSizeY), display::StaveZ[1]),
    TVector3(2*(display::ChipSizeX + display::ChipDistanceX), -(2*display::ChipStaveDistanceY + 1.5*display::ChipDistanceY + 3.5*display::ChipSizeY), display::StaveZ[1]),
    //row2
    TVector3(-2*(display::ChipSizeX + display::ChipDistanceX), -(display::ChipStaveDistanceY + 1.5*display::ChipDistanceY + 2.5*display::ChipSizeY), display::StaveZ[1]),
    TVector3(-(display::ChipSizeX + display::ChipDistanceX), -(display::ChipStaveDistanceY + 1.5*display::ChipDistanceY + 2.5*display::ChipSizeY), display::StaveZ[1]),
    TVector3(0, -(display::ChipStaveDistanceY + 1.5*display::ChipDistanceY + 2.5*display::ChipSizeY), display::StaveZ[1]),
    TVector3((display::ChipSizeX + display::ChipDistanceX), -(display::ChipStaveDistanceY + 1.5*display::ChipDistanceY + 2.5*display::ChipSizeY), display::StaveZ[1]),
    TVector3(2*(display::ChipSizeX + display::ChipDistanceX), -(display::ChipStaveDistanceY + 1.5*display::ChipDistanceY + 2.5*display::ChipSizeY), display::StaveZ[1]),
    //row3
    TVector3(-2*(display::ChipSizeX + display::ChipDistanceX), -(display::ChipStaveDistanceY + 1.5*display::ChipDistanceY + 1.5*display::ChipSizeY), display::StaveZ[1]),
    TVector3(-(display::ChipSizeX + display::ChipDistanceX), -(display::ChipStaveDistanceY + 1.5*display::ChipDistanceY + 1.5*display::ChipSizeY), display::StaveZ[1]),
    TVector3(0, -(display::ChipStaveDistanceY + 1.5*display::ChipDistanceY + 1.5*display::ChipSizeY), display::StaveZ[1]),
    TVector3((display::ChipSizeX + display::ChipDistanceX), -(display::ChipStaveDistanceY + 1.5*display::ChipDistanceY + 1.5*display::ChipSizeY), display::StaveZ[1]),
    TVector3(2*(display::ChipSizeX + display::ChipDistanceX), -(display::ChipStaveDistanceY + 1.5*display::ChipDistanceY + 1.5*display::ChipSizeY), display::StaveZ[1]),
    //row4
    TVector3(-2*(display::ChipSizeX + display::ChipDistanceX), -0.5*(display::ChipDistanceY + display::ChipSizeY), display::StaveZ[1]),
    TVector3(-(display::ChipSizeX + display::ChipDistanceX), -0.5*(display::ChipDistanceY + display::ChipSizeY), display::StaveZ[1]),
    TVector3(0, -0.5*(display::ChipDistanceY + display::ChipSizeY), display::StaveZ[1]),
    TVector3((display::ChipSizeX + display::ChipDistanceX), -0.5*(display::ChipDistanceY + display::ChipSizeY), display::StaveZ[1]),
    TVector3(2*(display::ChipSizeX + display::ChipDistanceX), -0.5*(display::ChipDistanceY + display::ChipSizeY), display::StaveZ[1]),
    //row5
    TVector3(-2*(display::ChipSizeX + display::ChipDistanceX), +0.5*(display::ChipDistanceY + display::ChipSizeY), display::StaveZ[1]),
    TVector3(-(display::ChipSizeX + display::ChipDistanceX), +0.5*(display::ChipDistanceY + display::ChipSizeY), display::StaveZ[1]),
    TVector3(0, +0.5*(display::ChipDistanceY + display::ChipSizeY), display::StaveZ[1]),
    TVector3((display::ChipSizeX + display::ChipDistanceX), +0.5*(display::ChipDistanceY + display::ChipSizeY), display::StaveZ[1]),
    TVector3(2*(display::ChipSizeX + display::ChipDistanceX), +0.5*(display::ChipDistanceY + display::ChipSizeY), display::StaveZ[1]),
    //row6
    TVector3(-2*(display::ChipSizeX + display::ChipDistanceX), +(display::ChipStaveDistanceY + 1.5*display::ChipDistanceY + 1.5*display::ChipSizeY), display::StaveZ[1]),
    TVector3(-(display::ChipSizeX + display::ChipDistanceX), +(display::ChipStaveDistanceY + 1.5*display::ChipDistanceY + 1.5*display::ChipSizeY), display::StaveZ[1]),
    TVector3(0, +(display::ChipStaveDistanceY + 1.5*display::ChipDistanceY + 1.5*display::ChipSizeY), display::StaveZ[1]),
    TVector3((display::ChipSizeX + display::ChipDistanceX), +(display::ChipStaveDistanceY + 1.5*display::ChipDistanceY + 1.5*display::ChipSizeY), display::StaveZ[1]),
    TVector3(2*(display::ChipSizeX + display::ChipDistanceX), +(display::ChipStaveDistanceY + 1.5*display::ChipDistanceY + 1.5*display::ChipSizeY), display::StaveZ[1]),
    //row7
    TVector3(-2*(display::ChipSizeX + display::ChipDistanceX), +(display::ChipStaveDistanceY + 1.5*display::ChipDistanceY + 2.5*display::ChipSizeY), display::StaveZ[1]),
    TVector3(-(display::ChipSizeX + display::ChipDistanceX), +(display::ChipStaveDistanceY + 1.5*display::ChipDistanceY + 2.5*display::ChipSizeY), display::StaveZ[1]),
    TVector3(0, +(display::ChipStaveDistanceY + 1.5*display::ChipDistanceY + 2.5*display::ChipSizeY), display::StaveZ[1]),
    TVector3((display::ChipSizeX + display::ChipDistanceX), +(display::ChipStaveDistanceY + 1.5*display::ChipDistanceY + 2.5*display::ChipSizeY), display::StaveZ[1]),
    TVector3(2*(display::ChipSizeX + display::ChipDistanceX), +(display::ChipStaveDistanceY + 1.5*display::ChipDistanceY + 2.5*display::ChipSizeY), display::StaveZ[1]),
    //row8
    TVector3(-2*(display::ChipSizeX + display::ChipDistanceX), +(2*display::ChipStaveDistanceY + 1.5*display::ChipDistanceY + 3.5*display::ChipSizeY), display::StaveZ[1]),
    TVector3(-(display::ChipSizeX + display::ChipDistanceX), +(2*display::ChipStaveDistanceY + 1.5*display::ChipDistanceY + 3.5*display::ChipSizeY), display::StaveZ[1]),
    TVector3(0, +(2*display::ChipStaveDistanceY + 1.5*display::ChipDistanceY + 3.5*display::ChipSizeY), display::StaveZ[1]),
    TVector3((display::ChipSizeX + display::ChipDistanceX), +(2*display::ChipStaveDistanceY + 1.5*display::ChipDistanceY + 3.5*display::ChipSizeY), display::StaveZ[1]),
    TVector3(2*(display::ChipSizeX + display::ChipDistanceX), +(2*display::ChipStaveDistanceY + 1.5*display::ChipDistanceY + 3.5*display::ChipSizeY), display::StaveZ[1]),
    //row9
    TVector3(-2*(display::ChipSizeX + display::ChipDistanceX), +(2*display::ChipStaveDistanceY + 2.5*display::ChipDistanceY + 4.5*display::ChipSizeY), display::StaveZ[1]),
    TVector3(-1*(display::ChipSizeX + display::ChipDistanceX), +(2*display::ChipStaveDistanceY + 2.5*display::ChipDistanceY + 4.5*display::ChipSizeY), display::StaveZ[1]),
    TVector3(0, +(2*display::ChipStaveDistanceY + 2.5*display::ChipDistanceY + 4.5*display::ChipSizeY), display::StaveZ[1]),
    TVector3(1*(display::ChipSizeX + display::ChipDistanceX), +(2*display::ChipStaveDistanceY + 2.5*display::ChipDistanceY + 4.5*display::ChipSizeY), display::StaveZ[1]),
    TVector3(2*(display::ChipSizeX + display::ChipDistanceX), +(2*display::ChipStaveDistanceY + 2.5*display::ChipDistanceY + 4.5*display::ChipSizeY), display::StaveZ[1]),


    //bottom layer
    //row0
    TVector3(-2*(display::ChipSizeX + display::ChipDistanceX), -(2*display::ChipStaveDistanceY + 2.5*display::ChipDistanceY + 4.5*display::ChipSizeY), display::StaveZ[0]),
    TVector3(-1*(display::ChipSizeX + display::ChipDistanceX), -(2*display::ChipStaveDistanceY + 2.5*display::ChipDistanceY + 4.5*display::ChipSizeY), display::StaveZ[0]),
    TVector3(0, -(2*display::ChipStaveDistanceY + 2.5*display::ChipDistanceY + 4.5*display::ChipSizeY), display::StaveZ[0]),
    TVector3(1*(display::ChipSizeX + display::ChipDistanceX), -(2*display::ChipStaveDistanceY + 2.5*display::ChipDistanceY + 4.5*display::ChipSizeY), display::StaveZ[0]),
    TVector3(2*(display::ChipSizeX + display::ChipDistanceX), -(2*display::ChipStaveDistanceY + 2.5*display::ChipDistanceY + 4.5*display::ChipSizeY), display::StaveZ[0]),
    //row1
    TVector3(-2*(display::ChipSizeX + display::ChipDistanceX), -(2*display::ChipStaveDistanceY + 1.5*display::ChipDistanceY + 3.5*display::ChipSizeY), display::StaveZ[0]),
    TVector3(-(display::ChipSizeX + display::ChipDistanceX), -(2*display::ChipStaveDistanceY + 1.5*display::ChipDistanceY + 3.5*display::ChipSizeY), display::StaveZ[0]),
    TVector3(0, -(2*display::ChipStaveDistanceY + 1.5*display::ChipDistanceY + 3.5*display::ChipSizeY), display::StaveZ[0]),
    TVector3((display::ChipSizeX + display::ChipDistanceX), -(2*display::ChipStaveDistanceY + 1.5*display::ChipDistanceY + 3.5*display::ChipSizeY), display::StaveZ[0]),
    TVector3(2*(display::ChipSizeX + display::ChipDistanceX), -(2*display::ChipStaveDistanceY + 1.5*display::ChipDistanceY + 3.5*display::ChipSizeY), display::StaveZ[0]),
    //row2
    TVector3(-2*(display::ChipSizeX + display::ChipDistanceX), -(display::ChipStaveDistanceY + 1.5*display::ChipDistanceY + 2.5*display::ChipSizeY), display::StaveZ[0]),
    TVector3(-(display::ChipSizeX + display::ChipDistanceX), -(display::ChipStaveDistanceY + 1.5*display::ChipDistanceY + 2.5*display::ChipSizeY), display::StaveZ[0]),
    TVector3(0, -(display::ChipStaveDistanceY + 1.5*display::ChipDistanceY + 2.5*display::ChipSizeY), display::StaveZ[0]),
    TVector3((display::ChipSizeX + display::ChipDistanceX), -(display::ChipStaveDistanceY + 1.5*display::ChipDistanceY + 2.5*display::ChipSizeY), display::StaveZ[0]),
    TVector3(2*(display::ChipSizeX + display::ChipDistanceX), -(display::ChipStaveDistanceY + 1.5*display::ChipDistanceY + 2.5*display::ChipSizeY), display::StaveZ[0]),
    //row3
    TVector3(-2*(display::ChipSizeX + display::ChipDistanceX), -(display::ChipStaveDistanceY + 1.5*display::ChipDistanceY + 1.5*display::ChipSizeY), display::StaveZ[0]),
    TVector3(-(display::ChipSizeX + display::ChipDistanceX), -(display::ChipStaveDistanceY + 1.5*display::ChipDistanceY + 1.5*display::ChipSizeY), display::StaveZ[0]),
    TVector3(0, -(display::ChipStaveDistanceY + 1.5*display::ChipDistanceY + 1.5*display::ChipSizeY), display::StaveZ[0]),
    TVector3((display::ChipSizeX + display::ChipDistanceX), -(display::ChipStaveDistanceY + 1.5*display::ChipDistanceY + 1.5*display::ChipSizeY), display::StaveZ[0]),
    TVector3(2*(display::ChipSizeX + display::ChipDistanceX), -(display::ChipStaveDistanceY + 1.5*display::ChipDistanceY + 1.5*display::ChipSizeY), display::StaveZ[0]),
    //row4
    TVector3(-2*(display::ChipSizeX + display::ChipDistanceX), -0.5*(display::ChipDistanceY + display::ChipSizeY), display::StaveZ[0]),
    TVector3(-(display::ChipSizeX + display::ChipDistanceX), -0.5*(display::ChipDistanceY + display::ChipSizeY), display::StaveZ[0]),
    TVector3(0, -0.5*(display::ChipDistanceY + display::ChipSizeY), display::StaveZ[0]),
    TVector3((display::ChipSizeX + display::ChipDistanceX), -0.5*(display::ChipDistanceY + display::ChipSizeY), display::StaveZ[0]),
    TVector3(2*(display::ChipSizeX + display::ChipDistanceX), -0.5*(display::ChipDistanceY + display::ChipSizeY), display::StaveZ[0]),
    //row5
    TVector3(-2*(display::ChipSizeX + display::ChipDistanceX), +0.5*(display::ChipDistanceY + display::ChipSizeY), display::StaveZ[0]),
    TVector3(-(display::ChipSizeX + display::ChipDistanceX), +0.5*(display::ChipDistanceY + display::ChipSizeY), display::StaveZ[0]),
    TVector3(0, +0.5*(display::ChipDistanceY + display::ChipSizeY), display::StaveZ[0]),
    TVector3((display::ChipSizeX + display::ChipDistanceX), +0.5*(display::ChipDistanceY + display::ChipSizeY), display::StaveZ[0]),
    TVector3(2*(display::ChipSizeX + display::ChipDistanceX), +0.5*(display::ChipDistanceY + display::ChipSizeY), display::StaveZ[0]),
    //row6
    TVector3(-2*(display::ChipSizeX + display::ChipDistanceX), +(display::ChipStaveDistanceY + 1.5*display::ChipDistanceY + 1.5*display::ChipSizeY), display::StaveZ[0]),
    TVector3(-(display::ChipSizeX + display::ChipDistanceX), +(display::ChipStaveDistanceY + 1.5*display::ChipDistanceY + 1.5*display::ChipSizeY), display::StaveZ[0]),
    TVector3(0, +(display::ChipStaveDistanceY + 1.5*display::ChipDistanceY + 1.5*display::ChipSizeY), display::StaveZ[0]),
    TVector3((display::ChipSizeX + display::ChipDistanceX), +(display::ChipStaveDistanceY + 1.5*display::ChipDistanceY + 1.5*display::ChipSizeY), display::StaveZ[0]),
    TVector3(2*(display::ChipSizeX + display::ChipDistanceX), +(display::ChipStaveDistanceY + 1.5*display::ChipDistanceY + 1.5*display::ChipSizeY), display::StaveZ[0]),
    //row7
    TVector3(-2*(display::ChipSizeX + display::ChipDistanceX), +(display::ChipStaveDistanceY + 1.5*display::ChipDistanceY + 2.5*display::ChipSizeY), display::StaveZ[0]),
    TVector3(-(display::ChipSizeX + display::ChipDistanceX), +(display::ChipStaveDistanceY + 1.5*display::ChipDistanceY + 2.5*display::ChipSizeY), display::StaveZ[0]),
    TVector3(0, +(display::ChipStaveDistanceY + 1.5*display::ChipDistanceY + 2.5*display::ChipSizeY), display::StaveZ[0]),
    TVector3((display::ChipSizeX + display::ChipDistanceX), +(display::ChipStaveDistanceY + 1.5*display::ChipDistanceY + 2.5*display::ChipSizeY), display::StaveZ[0]),
    TVector3(2*(display::ChipSizeX + display::ChipDistanceX), +(display::ChipStaveDistanceY + 1.5*display::ChipDistanceY + 2.5*display::ChipSizeY), display::StaveZ[0]),
    //row8
    TVector3(-2*(display::ChipSizeX + display::ChipDistanceX), +(2*display::ChipStaveDistanceY + 1.5*display::ChipDistanceY + 3.5*display::ChipSizeY), display::StaveZ[0]),
    TVector3(-(display::ChipSizeX + display::ChipDistanceX), +(2*display::ChipStaveDistanceY + 1.5*display::ChipDistanceY + 3.5*display::ChipSizeY), display::StaveZ[0]),
    TVector3(0, +(2*display::ChipStaveDistanceY + 1.5*display::ChipDistanceY + 3.5*display::ChipSizeY), display::StaveZ[0]),
    TVector3((display::ChipSizeX + display::ChipDistanceX), +(2*display::ChipStaveDistanceY + 1.5*display::ChipDistanceY + 3.5*display::ChipSizeY), display::StaveZ[0]),
    TVector3(2*(display::ChipSizeX + display::ChipDistanceX), +(2*display::ChipStaveDistanceY + 1.5*display::ChipDistanceY + 3.5*display::ChipSizeY), display::StaveZ[0]),
    //row9
    TVector3(-2*(display::ChipSizeX + display::ChipDistanceX), +(2*display::ChipStaveDistanceY + 2.5*display::ChipDistanceY + 4.5*display::ChipSizeY), display::StaveZ[0]),
    TVector3(-1*(display::ChipSizeX + display::ChipDistanceX), +(2*display::ChipStaveDistanceY + 2.5*display::ChipDistanceY + 4.5*display::ChipSizeY), display::StaveZ[0]),
    TVector3(0, +(2*display::ChipStaveDistanceY + 2.5*display::ChipDistanceY + 4.5*display::ChipSizeY), display::StaveZ[0]),
    TVector3(1*(display::ChipSizeX + display::ChipDistanceX), +(2*display::ChipStaveDistanceY + 2.5*display::ChipDistanceY + 4.5*display::ChipSizeY), display::StaveZ[0]),
    TVector3(2*(display::ChipSizeX + display::ChipDistanceX), +(2*display::ChipStaveDistanceY + 2.5*display::ChipDistanceY + 4.5*display::ChipSizeY), display::StaveZ[0]),
};


