#ifndef FROMDATA_H
#define FROMDATA_H

#include <string>
#include <array> 
#include <map>
#include <vector>
#include <iostream>
#include <cmath>
#include <TCanvas.h> 
#include <TView.h>   
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
    friend std::ostream &operator<<(std::ostream &output, const eventdata &ev);
    void takedata();    
    void print_data_on_canvas(TCanvas* can);

};

extern std::unordered_map<int, eventdata> alldata;



#endif



//17.825
//26.325
//34.825
//175085