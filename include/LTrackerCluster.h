#ifndef LTRACKERCLUSTER_H
#define LTRACKERCLUSTER_H 

#include <vector>
#include <cstddef>
#include <iostream>
#include "TMath.h"
#include "./eventdata.h"
class LTrackerTrack;

class LTrackerCluster {
public:
    LTrackerCluster();
    void Reset();
    friend std::ostream &operator<<(std::ostream &output, const LTrackerCluster &cluster);
    void AddTrackIdx(std::vector<int> track_idx) {cls_track_idx = track_idx;}
    void AddResiduals(std::vector<float> res_x, std::vector<float> res_y) {cls_res_x = res_x; cls_res_y = res_y;}
    void from_ltt_to_cluster(const LTrackerTrack &ltt);
    //void CalculateClusterPosition(LTrackerSignal signal);
    std::vector<float> GetClusterMeanX() {return cls_mean_x;}
    std::vector<float> GetClusterMeanY() {return cls_mean_y;}
    std::vector<float> GetClusterMeanZ() {return cls_mean_z;}
    std::vector<int> GetClusterChipId() {return cls_chip_id;}
    std::vector<float> GetClusterMeanErrX() {return cls_mean_err_x;}
    std::vector<float> GetClusterMeanErrY() {return cls_mean_err_y;}
    std::vector<float> GetClusterResX() {return cls_res_x;}
    std::vector<float> GetClusterResY() {return cls_res_y;}
    std::vector<float> GetClusterResXm2() {return cls_res_x_m2;}
    std::vector<float> GetClusterResYm2() {return cls_res_y_m2;}
    std::vector<unsigned int> GetClusterRowSpan() {return cls_row_span;}
    std::vector<unsigned int> GetClusterColSpan() {return cls_col_span;}
    std::vector<unsigned int> GetClusterSize() {return cls_size;}
    std::vector<int> GetIndex() {return cls_idx;}
    std::vector<unsigned char> GetClusterPattern() {return cls_pattern;}
    std::vector<int> GetClusterPatternPosition() {return cls_pattern_position;}
    std::vector<int> GetClusterIdx() {return cls_idx;}
    std::vector<int> GetClsTrackIdx() {return cls_track_idx;}

    // residuals for the second tracking methos
    std::vector<float> cls_res_x_m2;
    std::vector<float> cls_res_y_m2;
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