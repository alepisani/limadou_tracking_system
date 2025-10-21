 #include "LTrack.h"
//#include "detector_const.hh"
#include <random>
#include <chrono>
#include <iostream>

void LTrack::Reset() {
    trk_idx.clear();
    npoints.clear();
    x0.clear();
    y0.clear();
    z0.clear();
    theta.clear();
    phi.clear();
    chi2.clear();
  return;
}

LTrack::LTrack() {
  Reset();
}

void LTrack::Dump() {
    std::cout << "LTrack::Dump()" << std::endl;
    for(int i=0; i<x0.size(); i++){
        std::cout << "idx: " << trk_idx[i] << std::endl;
        std::cout << "npoints: " << npoints[i] << std::endl;
        std::cout << "x0: " << x0[i] << std::endl;
        std::cout << "y0: " << y0[i] << std::endl;
        std::cout << "z0: " << z0[i] << std::endl;
        std::cout << "theta: " << theta[i] << std::endl;
        std::cout << "phi: " << phi[i] << std::endl;
        std::cout << "chi2: " << chi2[i] << std::endl;
    std::cout << "-------------------------" << std::endl;
    }
}

void LTrack::AddTrack(int idx, int trk_npoints, float trk_x0, float trk_y0, float trk_z0, float trk_theta, float trk_phi) {
    trk_idx.push_back(idx);
    npoints.push_back(trk_npoints);
    x0.push_back(trk_x0);
    y0.push_back(trk_y0);
    z0.push_back(trk_z0);
    theta.push_back(trk_theta);
    phi.push_back(trk_phi);
    cout << trk_z0 << endl;

    return;
}

void LTrack::from_LTrack_to_TrackCand(LTrack &track, LTrackerTrack &ltt){

  for(int i = 0; i < track.x0.size(); ++i){
    LTrackCandidate trkCand;
    trkCand.x0 = track.x0[i];
    trkCand.y0 = track.y0[i];
    trkCand.z0 = track.z0[i];
    cout << "z_trkcand? " << trkCand.z0 << endl;
    trkCand.theta = track.theta[i] * TMath::Pi() / 180;
    trkCand.phi = track.phi[i] * TMath::Pi() / 180;
    ltt.tracks.push_back(trkCand);
    //cout << "AAAAAAAAAAAAA" << endl;
    //cout << trkCand << endl;
  }
  


}
