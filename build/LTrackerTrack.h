#ifndef LTRACKERTRACK_H
#define LTRACKERTRACK_H

#include <vector>
#include <unordered_map>
#include <ostream>
#include "TObject.h"
#include "./LTrackerCluster.h"
#include <TCanvas.h>
#include <TView.h>

struct LTrackCandidate
{
  int id = -1;
  int n_clus = -1;
  std::vector<int> clus_id;
  std::vector<int> tracklet_id;
  float theta = -999.;
  float phi = -999.;
  float x0 = -999.;
  float y0 = -999.;
  float z0 = -999.;
  float err_theta = -1.;
  float err_phi = -1.;
  float err_x0 = -1.;
  float err_y0 = -1.;
  float chi2 = -1.;

  void print_trackcandidate(LTrackCandidate& tr, TCanvas* reco);

  friend std::ostream &operator<<(std::ostream &output, const LTrackCandidate &tr)
  {
    output << "track id : " << tr.id << std::endl;
    output << "n_clus : " << tr.n_clus << std::endl;
    output << "clus_id: ";
    for (auto p : tr.clus_id)
    {
      output << p << ", ";
    }
    output << std::endl;
    output << "tracklet_id: ";
    for (auto p : tr.tracklet_id)
    {
      output << p << ", ";
    }
    output << std::endl;
    output << "theta : " << tr.theta << " +- " << tr.err_theta << std::endl;
    output << "phi : " << tr.phi << " +- " << tr.err_phi << std::endl;
    output << "x0 : " << tr.x0 << " +- " << tr.err_x0 << std::endl;
    output << "y0 : " << tr.y0 << " +- " << tr.err_y0 << std::endl;
    output << "z0 : " << tr.z0 << std::endl;
    output << "chi2: " << tr.chi2 << std::endl;
    return output;
  }
};

struct LCluster
{
  LCluster();
  void fill_cluster(LCluster&, double, double, double, double, double, double, int);

  int id = -1;
  float x = -999.;
  float errx = -999.;
  float y = -999.;
  float erry = -999.;
  float z = -999.;
  float errz = -999.;

  friend std::ostream &operator<<(std::ostream &output, const LCluster &cl)
  {
      output << "cluster id : " << cl.id << std::endl;
      output << "x : " << cl.x << " +- " << cl.errx << std::endl;
      output << "y : " << cl.y << " +- " << cl.erry << std::endl;
      output << "z : " << cl.z << " +- " << cl.errz << std::endl;
      return output;
  }
};

struct LTracklet
{
  
  LTracklet();
    int id = -1;
  int firstClusterId = -1;
  int secondClusterId = -1;

  friend std::ostream &operator<<(std::ostream &output, const LTracklet &trkl)
  {
    output << "tracklet id : " << trkl.id << std::endl;
    output << "firstClusterId: " << trkl.firstClusterId << std::endl;
    output << "secondClusterId: " << trkl.secondClusterId << std::endl;
    return output;
  }
};

class LTrackerTrack 
{

public:
  LTrackerTrack();
  void Reset();
  void createTracklet(std::pair<int, LCluster> cl_l0, std::pair<int, LCluster> cl_l1, std::vector<LTracklet> &tracklet_vector, int &tracklet_counter);
  void computeTracklets();
  static double fct(const std::vector<LCluster> &clusters, const double *par);
  void fitStraightLine(const std::vector<LCluster> &clusters, LTrackCandidate &trkCand);
  void addSpuriousTracks(std::vector<int> &used_tracklets, std::vector<int> &used_clusters, std::vector<LTracklet> &tracklets, std::unordered_map<int, LCluster> &cluster_map_first_layer, std::unordered_map<int, LCluster> &cluster_map_second_layer);
  void computeTrackCandidates();


  
  // clusters separated by layer
  std::unordered_map<int, LCluster> tidy_clusters_lay0;
  std::unordered_map<int, LCluster> tidy_clusters_lay1;
  std::unordered_map<int, LCluster> tidy_clusters_lay2;
  // tracklets separated by layer
  std::vector<LTracklet> tracklet_lay01;
  std::vector<LTracklet> tracklet_lay12;
  std::vector<LTracklet> tracklet_lay02;
  // track candidates
  std::vector<LTrackCandidate> track_candidates;
  // used clusters in the process
  std::vector<int> used_clusters_lay0;
  std::vector<int> used_clusters_lay1;
  std::vector<int> used_clusters_lay2;
  // final tracks
  std::vector<LTrackCandidate> tracks;
  
  float chi2_cut = 100;

};






#endif  