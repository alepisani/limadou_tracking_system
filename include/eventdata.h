#ifndef FROMDATA_H
#define FROMDATA_H

#include <string>
#include <array> 
#include <map>
#include <vector>
#include <iostream>
#include <cmath>
#include <TCanvas.h> // Include ROOT's TCanvas header
#include <TView.h>   // Include ROOT's TView header
#include "TMath.h"
using namespace std;


//funz che prende i dati e me li disegna sul plot
//estrapola metriche delle stats


struct eventdata{

    std::vector<unsigned char> DIR_turret_idx;
    std::vector<unsigned char> DIR_stave_idx;
    std::vector<unsigned char> DIR_chip_idx0;
    std::vector<unsigned char> DIR_chip_idx1;
    std::vector<unsigned int>  DIR_chip_id;
    std::vector<float> DIR_xpos;
    std::vector<float> DIR_ypos;
    std::vector<float> DIR_zpos;
    std::vector<int> DIR_cls_idx;

    eventdata();       //default constructor
    void takedata();    
    // data();
    // from_int_to_hex();
    // print_data_on_canvas();

};

extern std::unordered_map<int, eventdata> alldata;



#endif