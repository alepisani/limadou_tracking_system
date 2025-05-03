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


