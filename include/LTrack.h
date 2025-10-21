#ifndef __LTRACK__
#define __LTRACK__ 1

#include <vector>
#include <cstddef>
#include "./LTrackerTrack.h"

class LTrack {
public:
    LTrack();
    void Reset();
    void Dump();
    void AddTrack(int trk_idx, int trk_npoints, float x0, float y0, float z0, float theta, float phi);
    void AddChi2(std::vector<float> chi2_track) {chi2 = chi2_track;}
    void from_LTrack_to_TrackCand(LTrack &track, LTrackerTrack &ltt);

    std::vector<int> GetTrkIdx() {return trk_idx;}
    std::vector<int> GetNPoints() {return npoints;}
    std::vector<float> GetX0() {return x0;}
    std::vector<float> GetY0() {return y0;}
    std::vector<float> GetZ0() {return z0;}
    std::vector<float> GetTheta() {return theta;}
    std::vector<float> GetPhi() {return phi;}
    std::vector<float> GetChi2() {return chi2;}
private:
    std::vector<int> trk_idx;
    std::vector<int> npoints;
    std::vector<float> x0;
    std::vector<float> y0;
    std::vector<float> z0;
    std::vector<float> theta;
    std::vector<float> phi;
    std::vector<float> chi2;
};

#endif
