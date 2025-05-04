#ifndef LTRACKERCLUSTER_H
#define LTRACKERCLUSTER_H 

#include <vector>
#include <cstddef>
#include <iostream>
#include "TMath.h"

class LTrackerCluster {
public:
    LTrackerCluster();
    void Reset();
    std::vector<int> GetClusterIdx() {return cls_idx;}

    // residuals for the second tracking methos
    std::vector<float> cls_res_x_m2;
    std::vector<float> cls_res_y_m2;

private:
    std::vector<float> cls_mean_x;
    std::vector<float> cls_mean_y;
    std::vector<float> cls_mean_z;
    std::vector<int> cls_chip_id;
    std::vector<float> cls_mean_err_x;
    std::vector<float> cls_mean_err_y;
    std::vector<float> cls_res_x;
    std::vector<float> cls_res_y;
    std::vector<unsigned int> cls_row_span;
    std::vector<unsigned int> cls_col_span;
    std::vector<unsigned int> cls_size;
    std::vector<unsigned char> cls_pattern;
    std::vector<int> cls_pattern_position;
    std::vector<int> cls_idx;
    std::vector<int> cls_track_idx;


};







#endif