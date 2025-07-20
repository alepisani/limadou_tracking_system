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

class LTrackerCluster;

//funz che prende i dati e me li disegna sul plot
//estrapola metriche delle stats


struct eventdata{

    /* std::vector<unsigned char> DIR_turret_idx;
    std::vector<unsigned char> DIR_stave_idx;
    std::vector<unsigned char> DIR_chip_idx0;
    std::vector<unsigned char> DIR_chip_idx1;
    std::vector<unsigned int>  DIR_chip_id;
    std::vector<float> DIR_xpos;
    std::vector<float> DIR_ypos;
    std::vector<float> DIR_zpos;
    std::vector<int> DIR_cls_idx; */

    std::vector<unsigned int> cls_size;
    std::vector<float> cls_mean_x;
    std::vector<float> cls_mean_y;
    std::vector<float> cls_mean_z;
    std::vector<float> cls_mean_x_err;
    std::vector<float> cls_mean_y_err;

    eventdata();       //default constructor
    friend std::ostream &operator<<(std::ostream &output, const eventdata &ev);
    void takedata();    
    void print_data_on_canvas(TCanvas* can);
    void from_ev_to_cluster(LTrackerCluster&, eventdata&);

};

extern std::unordered_map<int, eventdata> alldata;



#endif



//17.825
//26.325
//34.825
//175085