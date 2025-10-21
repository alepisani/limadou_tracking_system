#ifndef __LTRACKFITTINGTOOLS__
#define __LTRACKFITTINGTOOLS__ 1

//#include "LTrackerSignal.hh"
#include "LTrackerCluster.h"
#include "LTrack.h"
//#include "LTrackerMask.hh"

#include <vector>
#include <cstddef>

struct EvPoint {
    double x;
    double y;
    double z;
    int trk_nr;
    int idx;
    EvPoint(double x, double y, double z, int idx) : x(x), y(y), z(z), idx(idx), trk_nr(-1) {}
};

void HoughTransform3D(LTrackerCluster &cluster,LTrack &track);
void CalculateResiduals(LTrackerCluster &cluster,LTrack &track);
void ChooseBestProjection();

#endif