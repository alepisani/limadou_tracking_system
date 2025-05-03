#ifndef LTRACKERTRACK_H
#define LTRACKERTRACK_H

#include <vector>
#include <unordered_map>
#include <ostream>
#include "TObject.h"

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

class LTrackerTrack{

public:
  LTrackerTrack();
  void Reset();



  // clusters separated by layer
  std::unordered_map<int, LCluster> tidy_clusters_lay0;
  std::unordered_map<int, LCluster> tidy_clusters_lay1;
  std::unordered_map<int, LCluster> tidy_clusters_lay2;
  // tracklets separated by layer
  std::vector<LTracklet> tracklet_lay01;
  std::vector<LTracklet> tracklet_lay12;
  std::vector<LTracklet> tracklet_lay02;
  // track candidates
  //std::vector<LTrackCandidate> track_candidates;
  // used clusters in the process
  std::vector<int> used_clusters_lay0;
  std::vector<int> used_clusters_lay1;
  std::vector<int> used_clusters_lay2;
  // final tracks
  //std::vector<LTrackCandidate> tracks;

};




#endif