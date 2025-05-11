#ifndef CHIP_H
#define CHIP_H

#include <string>
#include <array> 
#include <map>
#include <vector>
#include "TVector3.h"
#include <iostream>
#include <cmath>
#include <TCanvas.h> // Include ROOT's TCanvas header
#include <TView.h>   // Include ROOT's TView header
#include "TMath.h"
#include "display.h"
using namespace std;

struct chip{
    chip();
    //chip id
    //coordinate del chip (centro del chip)
    double x = -100;
    double y = -100;
    double z = 0;
    double x_dim = display::ChipSizeX;
    double y_dim = display::ChipSizeY;
    double z_dim = display::ChipSizeZ;
    unsigned short id = 0xFFFF;
    bool status = true;         //status (acceso = true/spento = false)

    bool check_chip_status(chip& c);
    void print_chip(chip& c, TCanvas* canvas, int);
    void set_chip_coordinates(chip& c, double x, double y, double z, unsigned short id);
    bool is_chip_dead(chip& c, const std::vector<unsigned short> dead_chip);
    void is_dead_chip_tracking_smt(chip& c, LCluster pL2, LCluster mL1, LCluster qL0);
};

struct chips {
    chips();
    //chips(const std::array<unsigned short, 150>& ChipIds, const std::array<TVector3, 150>& chip_coordinates);
    void print_all_chips(chips& c, TCanvas* canvas);
    void print_id_coordinates();

    static const std::vector<unsigned short> dead_chip;
    static const std::array<unsigned short, 150> ChipIds;
    static const std::array<TVector3, 150> chip_coordinates;
    static std::unordered_map<unsigned short, TVector3> id_coordinates; 
};

#endif



  /* Tracker naming scheme
   *
   *        ix=0    ix=1    ix=2    ix=3    ix=4    
   *     +----------------------------------------> X 
   *     | +-----+ +-----+ +-----+ +-----+ +-----+ \
   * iy=0| | 0x7c| | 0x7b| | 0x7a| | 0x79| | 0x78| |  \             row0
   *     | +-----+ +-----+ +-----+ +-----+ +-----+ |  | 
   *     |  ic=9    ic=8    ic=7    ic=6    ic=5   |  |  \          row1
   *     | +-----+ +-----+ +-----+ +-----+ +-----+ |  |  |
   * iy=1| | 0x70| | 0x71| | 0x72| | 0x73| | 0x74| |  |  |          row2
   *     | +-----+ +-----+ +-----+ +-----+ +-----+ |  |  |
   *     |  ic=0    ic=1    ic=2    ic=3    ic=4   |  |  |           |   
   *     |                                         |iz=0 |
   *     | +-----+ +-----+ +-----+ +-----+ +-----+ |  |  |
   * iy=2| |0x37c| |0x37b| |0x37a| |0x379| |0x738| |  |iz=1          |
   *     | +-----+ +-----+ +-----+ +-----+ +-----+ |  |  |
   *     |  ic=9    ic=8    ic=7    ic=6    ic=5   |  |  |iz=2
   *     | +-----+ +-----+ +-----+ +-----+ +-----+ |  |  |
   * iy=3| |0x370| |0x371| |0x372| |0x373| |0x374| |  |  |
   *     | +-----+ +-----+ +-----+ +-----+ +-----+ |  |  |
   *     |  ic=0    ic=1    ic=2    ic=3    ic=4   |  |  |
   *     |                   ...                  ...
   *   Y V
   */