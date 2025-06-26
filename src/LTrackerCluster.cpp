#include "../include/LTrackerCluster.h"
#include "../include/eventdata.h"
#include "../include/stats.h"
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



//sostituisci signal con i dati presi da eventdata
void LTrackerCluster::CalculateClusterPosition(eventdata ev) {
    //std::vector<int> cls_id = GetClusterIndex(signal);
    std::vector<float> x = ev.DIR_xpos;
    std::vector<float> y = ev.DIR_ypos;
    std::vector<float> z = ev.DIR_zpos;
    //std::vector<int> row = signal.GetPixelsRow();
    //std::vector<int> col = signal.GetPixelsCol();
    //std::vector<int> chip_id = ev.DIR_chip_id;
    std::vector<int> chip_id(ev.DIR_chip_id.begin(), ev.DIR_chip_id.end());
    // calculate unique cluster id
    std::vector<int> unique_cls_id;
    std::vector<int> cls_id = ev.DIR_cls_idx;
    for (int i = 0; i < cls_id.size(); ++i) {
        if (std::find(unique_cls_id.begin(), unique_cls_id.end(), cls_id[i]) == unique_cls_id.end()) {
            unique_cls_id.push_back(cls_id[i]);
        }
    }
    // calculate cluster mean position
    int counter = 0;
    for (int i = 0; i < unique_cls_id.size(); ++i)
    {
        float mean_x = 0;
        float mean_y = 0;
        float mean_z = 0;
        int cid = 0;
        int cls_pix_chip_id = 0;
        int size = 0;
        //unsigned int min_row = 0xffffffff;
        //unsigned int min_col = 0xffffffff;
        //unsigned int max_row = 0x0;
        //unsigned int max_col = 0x0;
        //unsigned int row_span = 0;
        //unsigned int col_span = 0;
        // buffer to build cluster pattern
        //std::vector<unsigned int> cls_rows_tmp, cls_cols_tmp;
        for (int j = 0; j < cls_id.size(); ++j)
        {
            if (cls_id[j] == unique_cls_id[i])
            {
                mean_x += x[j];
                mean_y += y[j];
                mean_z = z[j];
                cid = cls_id[j];
                cls_pix_chip_id = chip_id[j];
                size++;
                // evaluate min/max row and col
                /*
                if (row[j] < min_row)
                    min_row = row[j];
                if (col[j] < min_col)
                    min_col = col[j];
                if (row[j] > max_row)
                    max_row = row[j];
                if (col[j] > max_col)
                    max_col = col[j];
                // store (row, col) in buffer
                cls_rows_tmp.push_back(row[j]);
                cls_cols_tmp.push_back(col[j]);
                */
            }
        }
        mean_x /= size;
        mean_y /= size;
        //mean_z /= size;
        // get cluster span
        /*
        row_span = max_row - min_row + 1;
        col_span = max_col - min_col + 1;
         build cluster pattern
        int nbytes = ((int)row_span * (int)col_span) / 8;
        if(((int)row_span * (int)col_span) % 8 != 0){
            nbytes++;
        }
        std::vector<unsigned char> pattern(nbytes, 0);
        for (int i = 0; i < size; i++)
        {
            auto row = cls_rows_tmp[i] - min_row;
            auto col = cls_cols_tmp[i] - min_col;
            int bit = row * col_span + col;
            int element = bit / 8;
            bit = bit % 8;
            pattern[element] |= 128 >> bit;  // 128 = 0b10000000
        }
        */
        cls_mean_x.push_back(mean_x);
        cls_mean_y.push_back(mean_y);
        cls_mean_z.push_back(mean_z);
        cls_chip_id.push_back(cls_pix_chip_id);
        //cls_row_span.push_back(row_span);
        //cls_col_span.push_back(col_span);
        //cls_mean_err_x.push_back(GetClusterErrX(int(col_span)));
        //cls_mean_err_y.push_back(GetClusterErrY(int(row_span)));
        cls_size.push_back(size);
        //int position = std::distance(cls_pattern.begin(), cls_pattern.end());
        //cls_pattern_position.push_back(position);
        //cls_pattern.insert(cls_pattern.end(), pattern.begin(), pattern.end());
        cls_idx.push_back(unique_cls_id[i]);
        if(i==0){stats::hmbh1L++;}
        if(i==1){stats::hmbh2L++; stats::hmbh1L--;}
        if(i==2){stats::hmbh3L++; stats::hmbh2L--;}
    }
}
