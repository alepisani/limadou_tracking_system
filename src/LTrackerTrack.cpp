#include "TMath.h"
#include "TSystem.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "../build/stats.h"
#include "../build/display.h"
#include "../build/LTrackerTrack.h"

#include <cmath>
#include <algorithm>
#include <functional>

LCluster::LCluster(){}
LTracklet::LTracklet(){}
LTrackerTrack::LTrackerTrack(){}

void LTrackerTrack::Reset()
{
  tidy_clusters_lay0.clear();
  tidy_clusters_lay1.clear();
  tidy_clusters_lay2.clear();
  tracklet_lay01.clear();
  tracklet_lay12.clear();
  tracklet_lay02.clear();
  //track_candidates.clear();
  used_clusters_lay0.clear();
  used_clusters_lay1.clear();
  used_clusters_lay2.clear();
  //tracks.clear();
}

void LCluster::fill_cluster(LCluster& cl, double x, double y, double z, double errx, double erry, double errz, int id){
    cl.x = x;
    cl.errx = errx;
    cl.y = y;
    cl.erry = erry;
    cl.z = z;
    cl.errz = errz;
    cl.id = id;
}

void LTrackerTrack::createTracklet(std::pair<int, LCluster> cl_l0, std::pair<int, LCluster> cl_l1, std::vector<LTracklet> &tracklet_vector, int &tracklet_counter)
{
  LTracklet tracklet;
  tracklet.firstClusterId = cl_l0.first;
  tracklet.secondClusterId = cl_l1.first;
  tracklet.id = tracklet_counter++;
  tracklet_vector.push_back(tracklet);
}

void LTrackerTrack::computeTracklets()
{
  int tracklet_counter = 0;
  // layer 0 - layer 1
  for (auto &cl_l0 : tidy_clusters_lay0)
  {
    for (auto &cl_l1 : tidy_clusters_lay1)
    {
      createTracklet(cl_l0, cl_l1, tracklet_lay01, tracklet_counter);
    }
  }
  // layer 1 - layer 2
  for (auto &cl_l1 : tidy_clusters_lay1)
  {
    for (auto &cl_l2 : tidy_clusters_lay2)
    {
      createTracklet(cl_l1, cl_l2, tracklet_lay12, tracklet_counter);
    }
  }
  // layer 0 - layer 2
  for (auto &cl_l0 : tidy_clusters_lay0)
  {
    for (auto &cl_l2 : tidy_clusters_lay2)
    {
      createTracklet(cl_l0, cl_l2, tracklet_lay02, tracklet_counter);
    }
  }
}