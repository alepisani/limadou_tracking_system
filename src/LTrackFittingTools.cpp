#include "../include/LTrackFittingTools.h"
#include "../include/LTrackerCluster.h"
#include "../include/LTrackerTrack.h"
//#include "Line.hh"
//#include "analysis_const.hh"
//#include "Geometry.hh"
#include "TMath.h"
#include "TGraph2D.h"
#include "TRandom2.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TF2.h"
#include "TH1.h"
#include "Math/Functor.h"
#include "Math/Vector3D.h"
#include "TPolyLine3D.h"
#include "Fit/Fitter.h"
#include "TH2F.h"
#include <Eigen/Core>
#include <Eigen/Dense>

using namespace Eigen;
using namespace ROOT::Math;

#include <random>
#include <iostream>
#include <vector>
#include <cmath>

template<class Vector3d>
std::pair < Vector3d, Vector3d > best_line_from_points(const std::vector<Vector3d> & c)
{
	// copy coordinates to  matrix in Eigen format
	size_t num_atoms = c.size();
	Eigen::Matrix< typename Vector3d::Scalar, Eigen::Dynamic, Eigen::Dynamic > centers(num_atoms, 3);
	for (size_t i = 0; i < num_atoms; ++i) centers.row(i) = c[i];

	Vector3d origin = centers.colwise().mean();
	Eigen::MatrixXd centered = centers.rowwise() - origin.transpose();
	Eigen::MatrixXd cov = centered.adjoint() * centered;
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(cov);
	Vector3d axis = eig.eigenvectors().col(2).normalized();
    if (axis.z() < 0) axis = -axis;
	return std::make_pair(origin, axis);
}


double distanceEvPointLine(double x_point, double y_point, double q, double m) {
    return std::abs(y_point - (m * x_point + q)) / std::sqrt(1 + m * m);
}

double point_Rtheta(double x, double y, double theta) {
    return std::cos(theta) * x + std::sin(theta) * y;
}

void fitline(double rho, double theta, double& q, double& m) {
    double x_0 = rho * std::cos(theta);
    double y_0 = rho * std::sin(theta);
    q = y_0 + x_0 / std::tan(theta);
    m = -1 / std::tan(theta);
}

void IterativeHoughTransform(std::vector<std::vector<float>>& points, int i_ev, std::string dim1, std::string dim2) {
    int theta_binning_HT = 1000; // Set the appropriate value
    int rho_binning_HT = 1000; // Set the appropriate value
    float d_max = 5; // Set the appropriate value

    int i_track = 0;
    std::vector<std::vector<float>> df_ev;
    for (const std::vector<float>& point : points) {
        if (point.at(4) == -1) {
            df_ev.push_back(point);
        }
    }
    float theta_step = M_PI/10000;
    int iterations = 0;
    int df_size = df_ev.size();
    std::vector<double> theta_tot, rho_tot;
    while ((df_ev.size() >= 2) && (iterations < 5 )) {
        iterations ++;
        theta_tot.clear();
        rho_tot.clear();
        for (const std::vector<float>& point : df_ev) {
            if (dim1 == "x_pos" && dim2 == "z_pos") {
                for (double theta = 0; theta < M_PI; theta += theta_step) {
                    double rho = point_Rtheta(point.at(0), point.at(2), theta);
                    rho_tot.push_back(rho);
                    theta_tot.push_back(theta);
                }
            }
            else if (dim1 == "y_pos" && dim2 == "z_pos") {
                for (double theta = 0; theta < M_PI; theta += theta_step) {
                    double rho = point_Rtheta(point.at(1), point.at(2), theta);
                    rho_tot.push_back(rho);
                    theta_tot.push_back(theta);
                }
            }
            else if (dim1 == "x_pos" && dim2 == "y_pos") {
                for (double theta = 0; theta < M_PI; theta += theta_step) {
                    double rho = point_Rtheta(point.at(0), point.at(1), theta);
                    rho_tot.push_back(rho);
                    theta_tot.push_back(theta);
                }
            }
        }
        std::vector<std::vector<int>> hist(theta_binning_HT, std::vector<int>(rho_binning_HT, 0));
        for (size_t i = 0; i < theta_tot.size(); ++i) {
            int x_ind = static_cast<int>((theta_tot[i] * theta_binning_HT) / M_PI);
            int y_ind = static_cast<int>(((rho_tot[i] + 300) * rho_binning_HT) / 600);
            ++hist[x_ind][y_ind];
        }
        int max_count = 0;
        int x_ind_max = 0;
        int y_ind_max = 0;   // Find the bin with the maximum count

        for (int i = 0; i < theta_binning_HT; ++i) {
            for (int j = 0; j < rho_binning_HT; ++j) {
                if (hist[i][j] > max_count) {
                    max_count = hist[i][j];
                    x_ind_max = i;
                    y_ind_max = j;
                }
            }
        }

        // Convert the bin indices to theta and rho values
        double theta_max = (static_cast<double>(x_ind_max) * M_PI) / theta_binning_HT + M_PI / (theta_binning_HT*2);
        double rho_max = ((static_cast<double>(y_ind_max) * 600) / rho_binning_HT) - 300 + 600 / (rho_binning_HT*2);

        double q, m;
        fitline(rho_max, theta_max, q, m);
        for (std::vector<float>& point : df_ev) {
            if (dim1 == "x_pos" && dim2 == "z_pos") {
                double d = distanceEvPointLine(point.at(0), point.at(2), q, m);
                if (d < d_max) { // Set the appropriate threshold value
                    point.at(4) = i_track;
                    points[point.at(3)].at(4) = i_track;
                    
                }
            }
            else if (dim1 == "y_pos" && dim2 == "z_pos") {
                double d = distanceEvPointLine(point.at(1), point.at(2), q, m);
                if (d < d_max) { // Set the appropriate threshold value
                    point.at(4) = i_track;
                    points[point.at(3)].at(4) = i_track;
                }
            }
            else if (dim1 == "x_pos" && dim2 == "y_pos") {
                double d = distanceEvPointLine(point.at(0), point.at(1), q, m);
                if (d < d_max) { // Set the appropriate threshold value
                    point.at(4) = i_track;
                    points[point.at(3)].at(4) = i_track;
                }
            }
        }
        ++i_track;
        // Remove the fitted points from the input points
        df_ev.erase(std::remove_if(df_ev.begin(), df_ev.end(), [&](const std::vector<float>& point) {
            return point.at(4) != -1;
        }), df_ev.end());
        df_size = df_ev.size();
        
    }
}

float distance3dpoints(Vector3d line_point, Vector3d line_direction, const std::vector<Vector3d> & points3d){
    // loop over 3d points
    float distance = 0; 
    for (int i = 0 ; i < points3d.size() ; i++){ 
        // calculate distance
        Vector3d point_line = points3d[i] - line_point;
        Vector3d point_line_cross = point_line.cross(line_direction);
        distance += point_line_cross.norm() / line_direction.norm();
    }
    return distance;
}

std::vector<std::vector<float>> ChooseBestHough(std::vector<std::vector<float>> pointsxy, std::vector<std::vector<float>> pointsxz, std::vector<std::vector<float>> pointsyz){
    //find unique track number
    std::vector<int> trk_nr_xy, trk_nr_xz, trk_nr_yz;
    for (std::vector<float>& point : pointsxz) {
        trk_nr_xz.push_back(point.at(4));
    }
    for (std::vector<float>& point : pointsyz) {
        trk_nr_yz.push_back(point.at(4));
    }

    std::sort(trk_nr_xz.begin(), trk_nr_xz.end());
    std::sort(trk_nr_yz.begin(), trk_nr_yz.end());

    auto last_xz = std::unique(trk_nr_xz.begin(), trk_nr_xz.end());
    auto last_yz = std::unique(trk_nr_yz.begin(), trk_nr_yz.end());

    trk_nr_xz.erase(last_xz, trk_nr_xz.end());
    trk_nr_yz.erase(last_yz, trk_nr_yz.end());

    double dist_xz = 0.;
    double dist_yz = 0.;

    int npoints_xz = 0;
    int npoints_yz = 0;

    for (int i = 0; i < trk_nr_xz.size(); ++i) {
        if (trk_nr_xz[i] == -1) {
            continue;
        }
        std::vector<std::vector<float>> pointsxz_trk;
        for (const std::vector<float>& point : pointsxz) {
            if (point.at(4) == trk_nr_xz[i]) {
                npoints_xz ++;
                pointsxz_trk.push_back(point);
            }
        }
        std::vector<Vector3d> points3d;
        for (const std::vector<float>& point : pointsxz_trk) {
            points3d.push_back(Vector3d(point.at(0), point.at(1), point.at(2)));
        }
        std::pair < Vector3d, Vector3d > best_fiteigen = best_line_from_points(points3d);
        Vector3d line_point = best_fiteigen.first;
        Vector3d line_direction = best_fiteigen.second;

        dist_xz += distance3dpoints(line_point,line_direction,points3d);
    }
    for (int i = 0; i < trk_nr_yz.size(); ++i) {
        if (trk_nr_yz[i] == -1) {
            continue;
        }
        std::vector<std::vector<float>> pointsyz_trk;
        for (const std::vector<float>& point : pointsyz) {
            if (point.at(4) == trk_nr_yz[i]) {
                npoints_yz ++;
                pointsyz_trk.push_back(point);
            }
        }
        std::vector<Vector3d> points3d;
        for (const std::vector<float>& point : pointsyz_trk) {
            points3d.push_back(Vector3d(point.at(0), point.at(1), point.at(2)));
        }
        std::pair < Vector3d, Vector3d > best_fiteigen = best_line_from_points(points3d);
        Vector3d line_point = best_fiteigen.first;
        Vector3d line_direction = best_fiteigen.second;

        dist_yz += distance3dpoints(line_point,line_direction,points3d);
    }
    std::vector<std::vector<float>> best_points;
    if (npoints_xz == 0) dist_xz = 1000000;
    if (npoints_yz == 0) dist_yz = 1000000;
    if (dist_xz <= dist_yz) {
        best_points = pointsxz;
    }
    else {
        best_points = pointsyz;
    }
    return best_points;
}

std::vector<std::vector<float>> RemoveHorizontal(std::vector<std::vector<float>> points){
    // select unique track numbers
    std::vector<int> trk_nr;
    for (const std::vector<float>& point : points) {
        trk_nr.push_back(point.at(4));
    }
    std::sort(trk_nr.begin(), trk_nr.end());
    auto last = std::unique(trk_nr.begin(), trk_nr.end());
    trk_nr.erase(last, trk_nr.end());

    //check if track have at least two z coordinates
    std::vector<int> trk_nr_z;
    for (int i = 0; i < trk_nr.size(); ++i) {
        if (trk_nr[i] == -1) {
            continue;
        }
        int trk_nr_i = trk_nr[i];
        std::vector<float> z_pos;
        for (int j = 0; j < points.size(); ++j) {
            if (points.at(j).at(4) == trk_nr_i) {
                z_pos.push_back(points.at(j).at(2));
            }
        }
        std::sort(z_pos.begin(), z_pos.end());
        auto last_z = std::unique(z_pos.begin(), z_pos.end());
        z_pos.erase(last_z, z_pos.end());
        if (z_pos.size() < 2) {
            for (int j = 0; j < points.size(); ++j) {
                if (points.at(j).at(4) == trk_nr_i) {
                    points.at(j).at(4) = -1;
                }
            }
        }
    }
    return points;
}

void CalculateResiduals(LTrackerCluster &cluster, LTrack &track){
    std::vector<float> res_x_vec, res_y_vec;
    std::vector<float> chi2;

    std::vector<int> cls_trk_idx = cluster.GetClsTrackIdx();
    std::vector<float> cls_x_pos = cluster.GetClusterMeanX();
    std::vector<float> cls_y_pos = cluster.GetClusterMeanY();
    std::vector<float> cls_z_pos = cluster.GetClusterMeanZ();
    std::vector<float> cls_err_x = cluster.GetClusterMeanErrX();
    std::vector<float> cls_err_y = cluster.GetClusterMeanErrY();

    std::vector<float> x0 = track.GetX0();
    std::vector<float> y0 = track.GetY0();
    std::vector<float> z0 = track.GetZ0();
    std::vector<float> theta = track.GetTheta();     //degree
    std::vector<float> phi = track.GetPhi();         //degree

    for (int i = 0; i < cls_trk_idx.size(); ++i) {
        if (cls_trk_idx[i] < 0){
            res_x_vec.push_back(-999.);
            res_y_vec.push_back(-999.);
        }
        else{
            //track coordinates
            float x = cls_x_pos[i];
            //cluster coordinates
            float y = cls_y_pos[i];
            float z = cls_z_pos[i];
            float x_trk = x0[cls_trk_idx[i]]  + std::tan(theta[cls_trk_idx[i]] *TMath::Pi()/180) * std::cos(phi[cls_trk_idx[i]]  *TMath::Pi()/180) * (z - z0[cls_trk_idx[i]] );
            float y_trk = y0[cls_trk_idx[i]]  + std::tan(theta[cls_trk_idx[i]]  *TMath::Pi()/180) * std::sin(phi[cls_trk_idx[i]]  *TMath::Pi()/180) * (z - z0[cls_trk_idx[i]] );
            float res_x = x - x_trk;
            float res_y = y - y_trk;
            res_x_vec.push_back(res_x);
            res_y_vec.push_back(res_y);
        }
    }

    for (int i = 0; i < phi.size(); ++i) {
        float chi2_tot = 0;
        float chi2_x = 0;
        float chi2_y = 0;
        for (int j = 0; j < cls_trk_idx.size(); ++j) {
            if (cls_trk_idx[j] == i) {
                chi2_x += std::pow(res_x_vec[j]/cls_err_x[j],2);
                chi2_y += std::pow(res_y_vec[j]/cls_err_y[j],2);
                chi2_tot += std::sqrt(chi2_x + chi2_y);
            }
        }
        chi2.push_back(chi2_tot);
    }

    cluster.AddResiduals(res_x_vec, res_y_vec);
    track.AddChi2(chi2);
}

void HoughTransform3D(LTrackerCluster &cluster, LTrack &track) {
    std::vector<std::vector<float>> pointsxy, pointsxz, pointsyz;
    std::vector<int> noise_cls_idx;
    std::vector<int> cls_idx;
    std::vector<float> cls_x_pos, cls_y_pos, cls_z_pos;
    cls_x_pos = cluster.GetClusterMeanX();
    cls_y_pos = cluster.GetClusterMeanY();
    cls_z_pos = cluster.GetClusterMeanZ();
    cls_idx = cluster.GetClusterIdx();
    std::vector<int> cls_trk_idx(cls_idx.size(), -1);

    int counter = 0;
    for (int i = 0; i < cls_x_pos.size(); ++i) {
        std::vector<float> point; // Create separate vectors in each iteration
        if (cls_idx[i] == -999) {
            cls_trk_idx[i] = -999;
            continue;
        }
        point.push_back(cls_x_pos[i]); // cls x
        point.push_back(cls_y_pos[i]); // cls y
        point.push_back(cls_z_pos[i]); // cls z
        point.push_back(counter); // cls idx
        point.push_back(-1); // track idx
        point.push_back(i); // old order

        pointsxy.push_back(point);
        pointsxz.push_back(point);
        pointsyz.push_back(point);

        counter ++;
    }

    //cout << pointsxy.size() << endl;
    //cout << pointsxz.size() << endl;
    //cout << pointsyz.size() << endl;



    if (pointsxy.size() < 2) {
        for (int i = 0; i < cls_trk_idx.size(); ++i) {
            cls_trk_idx[i] = -1;
        }
        cluster.AddTrackIdx(cls_trk_idx);
        return;
    }

    IterativeHoughTransform(pointsxy, 0, "x_pos", "y_pos");
    pointsxy = RemoveHorizontal(pointsxy);
    IterativeHoughTransform(pointsxz, 0, "x_pos", "z_pos");
    pointsxz = RemoveHorizontal(pointsxz);
    IterativeHoughTransform(pointsyz, 0, "y_pos", "z_pos");
    pointsyz = RemoveHorizontal(pointsyz);
    std::vector<std::vector<float>> best_ht = ChooseBestHough(pointsxy, pointsxz, pointsyz);
    std::vector<int> point_track_nr;
    // loop over ordered 
    cout << "quante best_ht? " << best_ht.size() << endl;
    for (const std::vector<float>& point : best_ht) {
        point_track_nr.push_back(point.at(4));
    }
    

    std::vector<int> unique_trk_nr;
    for (const std::vector<float>& point : best_ht) {
        if (point.at(4) == -1) {
            cout << "A" << endl;
            continue;
        }
        if (std::find(unique_trk_nr.begin(), unique_trk_nr.end(), point.at(4)) == unique_trk_nr.end()) {
            unique_trk_nr.push_back(point.at(4));
        }
    }
    int counter_trk = 0;
    for (int i = 0; i < unique_trk_nr.size(); ++i) {
        std::vector<std::vector<float>> track_points;
        std::vector<int> point_track_order;
        for (const std::vector<float>& point : best_ht) {
            if (point.at(4) == unique_trk_nr[i]) {
                track_points.push_back(point);
                point_track_order.push_back(point.at(5));
            }
        }
        if (track_points.size() < 2) {
            for (int j = 0; j < point_track_order.size(); ++j) {
                cls_trk_idx[point_track_order[j]] = -1;
            }
        }
        else{
            for (int j = 0; j < point_track_order.size(); ++j) {
                cls_trk_idx[point_track_order[j]] = i;
            }
        }
        //fill points in Vector3D
        std::vector<Vector3d> points3d;
        int npoints_trk = 0;
        for (const std::vector<float>& point : track_points) {
            points3d.push_back(Vector3d(point.at(0), point.at(1), point.at(2)));
            npoints_trk ++;
        }
        std::pair < Vector3d, Vector3d > best_fiteigen = best_line_from_points(points3d);
        track.AddTrack(counter_trk,npoints_trk,best_fiteigen.first.x(), best_fiteigen.first.y(), best_fiteigen.first.z(), std::acos(best_fiteigen.second.z()) * 180/TMath::Pi(), std::atan2(best_fiteigen.second.y(),best_fiteigen.second.x())*180/TMath::Pi());
        cout << "z_hough? " << best_fiteigen.first.z() << endl;
        LTrackerTrack ltt;
        counter_trk ++;
    }
    //cluster.AddTrackIdx(cls_trk_idx);
}