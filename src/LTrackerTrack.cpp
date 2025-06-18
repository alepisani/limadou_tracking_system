#include "TMath.h"
#include "TSystem.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include <TPolyLine3D.h>
#include <TCanvas.h>
#include <TView.h>
#include "../include/stats.h"
#include "../include/display.h"
#include "../include/LTrackerTrack.h"
#include "../include/LTrackerCluster.h"

#include <cmath>
#include <algorithm>
#include <functional>
#include "LTrackerTrack.h"

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
  cout << "gg" << endl;
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

// Distance function to be minimised
double LTrackerTrack::fct(const std::vector<LCluster> &clusters, const double *par)
{
  double x0 = par[0];
  double y0 = par[1];
  double theta = par[2];
  double phi = par[3];

  double chi2 = 0.0;
  for (int i = 0; i < 3; i++)
  {
    double x_fit = x0 + (clusters[i].z - display::z_origin_shift) * TMath::Tan(theta) * TMath::Cos(phi);
    double y_fit = y0 + (clusters[i].z - display::z_origin_shift) * TMath::Tan(theta) * TMath::Sin(phi);
    double dx = (clusters[i].x - x_fit) / clusters[i].errx;
    double dy = (clusters[i].y - y_fit) / clusters[i].erry;
    chi2 += dx * dx + dy * dy;
  }

  return chi2;
};

void LTrackerTrack::fitStraightLine(const std::vector<LCluster> &clusters, LTrackCandidate &trkCand)
{
  if (clusters.size() < 3)
  {
    return;
  }

  ROOT::Math::Minimizer *min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "");

  min->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
  min->SetMaxIterations(10000);      // for GSL
  min->SetTolerance(0.01);
  min->SetPrintLevel(0);

  // Create a functor putting together the static function fct and the ThreeClusters object
  std::function<double(const std::vector<LCluster> &, const double *)> unboundFct = fct;
  std::function<double(const double *)> boundFct = std::bind(unboundFct, clusters, std::placeholders::_1);
  ROOT::Math::Functor fct(boundFct, 4);

  // Assign fct to the minimiser
  min->SetFunction(fct);

  // Set the initial parameter values and step sizes
  double variable[4] = {0.0, 0.0, 0.0, 0.0};
  double step[4] = {0.01, 0.01, 0.01, 0.01};

  min->SetVariable(0, "x0", variable[0], step[0]);
  min->SetVariable(1, "y0", variable[1], step[1]);
  min->SetLimitedVariable(2, "theta", variable[2], step[2], 0., TMath::Pi() / 2);
  min->SetLimitedVariable(3, "phi", variable[3], step[3], -1 * TMath::Pi(), TMath::Pi());

  // Perform minimisation
  min->Minimize();

  const double *params = min->X();
  const double *errors = min->Errors();

  // Set the fitted values, errors, and chi-square in the LTrackCandidate object
  trkCand.x0 = params[0];
  trkCand.y0 = params[1];
  trkCand.z0 = display::z_origin_shift;
  trkCand.theta = params[2] * 180 / TMath::Pi();
  trkCand.phi = params[3] * 180 / TMath::Pi();
  trkCand.err_x0 = errors[0];
  trkCand.err_y0 = errors[1];
  trkCand.err_theta = errors[2] * 180 / TMath::Pi();
  trkCand.err_phi = errors[3] * 180 / TMath::Pi();
  trkCand.chi2 = min->MinValue();
}

// Add unused tracklets to without used clusters to final tracks (chi2 = 1 by definition)
void LTrackerTrack::addSpuriousTracks(std::vector<int> &used_tracklets, std::vector<int> &used_clusters, std::vector<LTracklet> &tracklets, std::unordered_map<int, LCluster> &cluster_map_first_layer, std::unordered_map<int, LCluster> &cluster_map_second_layer)
{
  for (auto &trkl : tracklets)
  {
    // Reject tracklets already used in full tracks
    if (std::find_if(used_tracklets.begin(), used_tracklets.end(), [&](int id)
                     { return id == trkl.id; }) != used_tracklets.end())
    {
      continue;
    }
    // Reject tracklets with clusters already used in full tracks
    if (std::find_if(used_clusters.begin(), used_clusters.end(), [&](int id)
                     { return (id == trkl.firstClusterId || id == trkl.secondClusterId); }) != used_clusters.end())
    {
      continue;
    }

    // Consider remaining tracklets
    auto first_cluster = cluster_map_first_layer[trkl.firstClusterId];
    auto second_cluster = cluster_map_second_layer[trkl.secondClusterId];
    double tan_theta = std::hypot((double)first_cluster.x - (double)second_cluster.x, (double)first_cluster.y - (double)second_cluster.y) / ((double)second_cluster.z - (double)first_cluster.z);
    double delta_y = (double)second_cluster.y - (double)first_cluster.y;
    double delta_x = (double)second_cluster.x - (double)first_cluster.x;
    double delta_r = TMath::Sqrt(delta_y * delta_y + delta_x * delta_x);
    double cos_phi = delta_x / delta_r;
    double sin_phi = delta_y / delta_r;
    double x0 = (double)first_cluster.x - ((double)first_cluster.z - display::z_origin_shift) * tan_theta * cos_phi;
    double y0 = (double)first_cluster.y - ((double)first_cluster.z - display::z_origin_shift) * tan_theta * sin_phi;
    double z0 = (double)first_cluster.z;

    LTrackCandidate spurious;
    spurious.n_clus = 2;
    spurious.clus_id.push_back(first_cluster.id);
    spurious.clus_id.push_back(second_cluster.id);
    spurious.tracklet_id.push_back(trkl.id);
    spurious.theta = TMath::ATan(tan_theta) * 180 / TMath::Pi();
    spurious.phi = TMath::ACos(cos_phi) * 180 / TMath::Pi();
    if (delta_y < 0)
    {
      spurious.phi *= -1.;
    }
    spurious.x0 = x0;
    spurious.y0 = y0;
    spurious.z0 = z0;
    spurious.err_x0 = -1.;
    spurious.err_y0 = -1.;
    spurious.err_theta = -1.;
    spurious.err_phi = -1.;
    spurious.chi2 = -1.;
    tracks.push_back(spurious);
  }
}

void LTrackerTrack::computeTrackCandidates(TCanvas* reco)
{
  //inside () should have LTrackerCluster &clusterer
  int candidateCounter = 0;

  //std::vector<int> clusterer_indices = clusterer.GetClusterIdx();           
  //int n_cluster_clusterer = (int)clusterer_indices.size();
  //clusterer.cls_res_x_m2.resize(n_cluster_clusterer);
  //clusterer.cls_res_y_m2.resize(n_cluster_clusterer);

  cout << "trackelt_lay01: " << tracklet_lay01.size() << endl;
  cout << "trackelt_lay12: " << tracklet_lay12.size() << endl;
  cout << "trackelt_lay02: " << tracklet_lay02.size() << endl;

  for (auto &trkl01 : tracklet_lay01)
  {
    if (candidateCounter > 1000)
    {
      std::cout << "Too many track candidates" << std::endl;
      break;
    }
    for (auto &trkl12 : tracklet_lay12)
    {
      if (trkl01.secondClusterId == trkl12.firstClusterId)
      {
        std::vector<LCluster> clus_vec = {tidy_clusters_lay0[trkl01.firstClusterId],
                                          tidy_clusters_lay1[trkl01.secondClusterId],
                                          tidy_clusters_lay2[trkl12.secondClusterId]};
        LTrackCandidate trkCand;
        fitStraightLine(clus_vec, trkCand);
        trkCand.id = candidateCounter++;
        trkCand.n_clus = 3;
        trkCand.clus_id = {trkl01.firstClusterId, trkl01.secondClusterId, trkl12.secondClusterId};
        // compute residuals
        for (auto &clus : clus_vec)
        {
          float l_res_x = clus.x - (trkCand.x0 + (clus.z - display::z_origin_shift) * TMath::Tan(trkCand.theta) * TMath::Cos(trkCand.phi));
          float l_res_y = clus.y - (trkCand.y0 + (clus.z - display::z_origin_shift) * TMath::Tan(trkCand.theta) * TMath::Sin(trkCand.phi));

          //int position = std::distance(clusterer_indices.begin(), std::find(clusterer_indices.begin(), clusterer_indices.end(), clus.id));
          //clusterer.cls_res_x_m2[position] = l_res_x;
          //clusterer.cls_res_y_m2[position] = l_res_y;
        }
        trkCand.tracklet_id = {trkl01.id, trkl12.id};
        track_candidates.push_back(trkCand);
      }
    }
  }
  // Sort track candidates by descending chi2
  std::sort(track_candidates.begin(), track_candidates.end(), [](LTrackCandidate &a, LTrackCandidate &b)
            { return a.chi2 < b.chi2; });
  // Remove candidates with large chi2
  auto new_end = std::remove_if(track_candidates.begin(), track_candidates.end(), [&](LTrackCandidate &trk)
                                { return trk.chi2 > chi2_cut; });
  track_candidates.erase(new_end, track_candidates.end());

  // Remove track candidates sharing clusters
  auto share_clusters = [](LTrackCandidate &trk1, LTrackCandidate &trk2)
  { return (trk1.clus_id[0] == trk2.clus_id[0] ||
            trk1.clus_id[1] == trk2.clus_id[1] ||
            trk1.clus_id[2] == trk2.clus_id[2]); };

  track_candidates.erase(std::unique(track_candidates.begin(), track_candidates.end(), share_clusters), track_candidates.end());

  // Record used tracklets and clusters
  std::vector<int> used_tracklets;
  std::vector<int> used_clusters;

  for (auto &trk : track_candidates)
  {
    tracks.push_back(trk);
    used_tracklets.insert(used_tracklets.end(), trk.tracklet_id.begin(), trk.tracklet_id.end());
    used_clusters.insert(used_clusters.end(), trk.clus_id.begin(), trk.clus_id.end());
  }

  addSpuriousTracks(used_tracklets, used_clusters, tracklet_lay01, tidy_clusters_lay0, tidy_clusters_lay1);
  addSpuriousTracks(used_tracklets, used_clusters, tracklet_lay12, tidy_clusters_lay1, tidy_clusters_lay2);
  addSpuriousTracks(used_tracklets, used_clusters, tracklet_lay02, tidy_clusters_lay0, tidy_clusters_lay2);

  // Reassigning track id
  for (int i = 0; i < tracks.size(); i++)
  {
    tracks[i].id = i;
  }


  cout << "tracce ricostruite: " << tracks.size() << endl;


}


void LTrackerTrack::printRecoTracks(TCanvas* reco) {

  for (auto &trk : tracks){
    float x1, y1, dz, x2, y2, t, p;
    dz = 100;         
    t = (trk.theta/90)*2*TMath::Pi();
    p = (trk.phi/90)*2*TMath::Pi();
    x1 = trk.x0 + (trk.z0 - dz)*(TMath::Tan(t))*(TMath::Cos(p));
    y1 = trk.y0 + (trk.z0 - dz)*(TMath::Tan(t))*(TMath::Sin(p));
    x2 = trk.x0 + (trk.z0 + dz)*(TMath::Tan(t))*(TMath::Cos(p));
    y2 = trk.y0 + (trk.z0 + dz)*(TMath::Tan(t))*(TMath::Sin(p));
    Double_t x_line[3] = {x2, trk.x0, x1};
    Double_t y_line[3] = {y2, trk.y0, y1};
    Double_t z_line[3] = {-dz, trk.z0, dz};
    TPolyLine3D* line_track = new TPolyLine3D(3, x_line, y_line, z_line);
    line_track->SetLineWidth(2);
    line_track->SetLineColor(kRed);
    line_track->Draw();
    reco->Update();
    
    cout << "theta: " << trk.theta << endl;
    cout << "phi: " << trk.phi << endl;

  }

}

