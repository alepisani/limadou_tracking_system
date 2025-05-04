#include "../build/LTrackerCluster.h"
#include <random>
#include <iostream>
#include <vector>
#include <cmath>

LTrackerCluster::LTrackerCluster() {
    Reset();
}

void LTrackerCluster::Reset() {
    cls_mean_x.clear();
    cls_mean_y.clear();
    cls_mean_z.clear();
    cls_chip_id.clear();
    cls_mean_err_x.clear();
    cls_mean_err_y.clear();
    cls_res_x.clear();
    cls_res_y.clear();
    cls_res_x_m2.clear();
    cls_res_y_m2.clear();
    cls_row_span.clear();
    cls_col_span.clear();
    cls_size.clear();
    cls_pattern.clear();
    cls_idx.clear();
    cls_track_idx.clear();
    return;
}