#include "TMath.h"
#include "TSystem.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include <TPolyLine3D.h>
#include <TCanvas.h>
#include <TView.h>
#include "TMarker3DBox.h"
#include "../include/stats.h"
#include "../include/display.h"
#include "../include/LTrackerTrack.h"
#include "../include/LTrackerCluster.h"

#include <cmath>
#include <algorithm>
#include <functional>
#include "LTrackerTrack.h"

LCluster::LCluster() {}
LTracklet::LTracklet() {}
LTrackerTrack::LTrackerTrack() {}

void LTrackerTrack::Reset()
{
  tidy_clusters_lay0.clear();
  tidy_clusters_lay1.clear();
  tidy_clusters_lay2.clear();
  tracklet_lay01.clear();
  tracklet_lay12.clear();
  tracklet_lay02.clear();
  track_candidates.clear();
  used_clusters_lay0.clear();
  used_clusters_lay1.clear();
  used_clusters_lay2.clear();
  tracks.clear();
}

void LCluster::fill_cluster(LCluster &cl, double x, double y, double z, double errx, double erry, double errz, int id)
{
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
  /* cout << "~~~~~~~~~~~~~~~~~~~~~" << endl;
  cout << "tracklet_counter: " << tracklet_counter << endl;
  cout << "trackelt_lay01: " << tracklet_lay01.size() << endl;
  cout << "trackelt_lay12: " << tracklet_lay12.size() << endl;
  cout << "trackelt_lay02: " << tracklet_lay02.size() << endl;
  cout << "~~~~~~~~~~~~~~~~~~~~~" << endl; */
}

void LTrackerTrack::print_tracklet(const LCluster cl_0, const LCluster cl_2)
{

  double x0, y0, z0, x2, y2, z2;
  x0 = cl_0.x;
  y0 = cl_0.y;
  z0 = cl_0.z;
  x2 = cl_2.x;
  y2 = cl_2.y;
  z2 = cl_2.z;
  Double_t x_line[2] = {x0, x2};
  Double_t y_line[2] = {y0, y2};
  Double_t z_line[2] = {z0, z2};
  TPolyLine3D *trk = new TPolyLine3D(2, x_line, y_line, z_line);
  trk->SetLineWidth(2);
  trk->SetLineColor(kGreen);
  trk->Draw();
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

/*
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
  //min->SetLimitedVariable(2, "theta", variable[2], step[2], 0., TMath::Pi() / 2);
  min->SetLimitedVariable(2, "theta", variable[2], step[2], 0, 4*TMath::Pi());
  //min->SetLimitedVariable(3, "phi", variable[3], step[3], -1 * TMath::Pi(), TMath::Pi());
  min->SetLimitedVariable(3, "phi", variable[3], step[3], -2 * TMath::Pi(), 2 * TMath::Pi());

  // Perform minimisation
  min->Minimize();

  const double *params = min->X();
  const double *errors = min->Errors();

  // Set the fitted values, errors, and chi-square in the LTrackCandidate object
  trkCand.x0 = params[0];
  trkCand.y0 = params[1];
  trkCand.z0 = display::z_origin_shift;
  trkCand.theta = params[2] * 180 / TMath::Pi();        //from rad -> gradi
  trkCand.phi = params[3] * 180 / TMath::Pi();          //from rad -> gradi
  trkCand.err_x0 = errors[0];
  trkCand.err_y0 = errors[1];
  trkCand.err_theta = errors[2] * 180 / TMath::Pi();
  trkCand.err_phi = errors[3] * 180 / TMath::Pi();
  trkCand.chi2 = min->MinValue();
}
*/

void LTrackerTrack::fitStraightLine(const std::vector<LCluster> &clusters, LTrackCandidate &trkCand)
{
  if (clusters.size() < 3)
    return;

  auto minimize = [&](double *initialVars, double *steps, const char *tag) -> std::tuple<const double *, const double *, double>
  {
    ROOT::Math::Minimizer *min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "");
    min->SetMaxFunctionCalls(1000000);
    min->SetMaxIterations(10000);
    min->SetTolerance(0.01);
    min->SetPrintLevel(0);

    std::function<double(const std::vector<LCluster> &, const double *)> unboundFct = fct;
    std::function<double(const double *)> boundFct = std::bind(unboundFct, clusters, std::placeholders::_1);
    ROOT::Math::Functor f(boundFct, 4);
    min->SetFunction(f);

    min->SetVariable(0, "x0", initialVars[0], steps[0]);
    min->SetVariable(1, "y0", initialVars[1], steps[1]);
    min->SetLimitedVariable(2, "theta", initialVars[2], steps[2], -TMath::Pi(), 4 * TMath::Pi());
    min->SetLimitedVariable(3, "phi", initialVars[3], steps[3], -2 * TMath::Pi(), 2 * TMath::Pi());

    min->Minimize();

    return {min->X(), min->Errors(), min->MinValue()};
  };

  double vars[4] = {0.0, 0.0, 0.0, 0.0};
  double step_coarse[4] = {0.01, 0.01, 0.01, 0.01};
  double step_fine[4] = {0.001, 0.001, 0.001, 0.001};

  const double *params, *errors;
  double chi2;
  std::tie(params, errors, chi2) = minimize(vars, step_coarse, "coarse");

  if (chi2 <= 100.0)
  {
    double refined_vars[4] = {params[0], params[1], params[2], params[3]};
    std::tie(params, errors, chi2) = minimize(refined_vars, step_fine, "fine");
  }

  trkCand.x0 = params[0];
  trkCand.y0 = params[1];
  trkCand.z0 = display::z_origin_shift;
  trkCand.theta = params[2] * 180 / TMath::Pi();
  trkCand.phi = params[3] * 180 / TMath::Pi();
  trkCand.err_x0 = errors[0];
  trkCand.err_y0 = errors[1];
  trkCand.err_theta = errors[2] * 180 / TMath::Pi();
  trkCand.err_phi = errors[3] * 180 / TMath::Pi();
  trkCand.chi2 = chi2;
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
    double x0 = (double)first_cluster.x - ((double)first_cluster.z - display::z_origin_shift) * tan_theta * cos_phi; // theta e phi in rad
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

void LTrackerTrack::New_addSpuriousTracks(std::vector<int> &used_tracklets, std::vector<int> &used_clusters)
{
  computeTracklets();
  LTrackerTrack t;

  for (auto &trkl01 : tracklet_lay01)
  {
    // Reject tracklets already used in full tracks
    if (std::find_if(used_tracklets.begin(), used_tracklets.end(), [&](int id)
                     { return id == trkl01.id; }) != used_tracklets.end())
    {
      continue;
    }
    // Reject tracklets with clusters already used in full tracks
    if (std::find_if(used_clusters.begin(), used_clusters.end(), [&](int id)
                     { return (id == trkl01.firstClusterId || id == trkl01.secondClusterId); }) != used_clusters.end())
    {
      continue;
    }
    auto cls_lay0 = tidy_clusters_lay0[trkl01.firstClusterId];
    auto cls_lay1 = tidy_clusters_lay1[trkl01.secondClusterId];
    // t.print_tracklet(cls_lay0, cls_lay1);
    double delta_y = (double)cls_lay1.y - (double)cls_lay0.y;
    double delta_x = (double)cls_lay1.x - (double)cls_lay0.x;
    double delta_z = (double)cls_lay1.z - (double)cls_lay0.z;
    double r = TMath::Sqrt(pow(delta_x, 2) + pow(delta_y, 2) + pow(delta_z, 2));
    double phi = TMath::ATan2(delta_y, delta_x);
    double theta = TMath::ACos(delta_z / r);
    double x2 = (double)cls_lay0.x + 2 * display::dist_z * TMath::Tan(theta) * TMath::Cos(phi);
    double y2 = (double)cls_lay0.y + 2 * display::dist_z * TMath::Tan(theta) * TMath::Sin(phi);
    double z2 = display::StaveZ[2];

    LTrackCandidate spurious;
    spurious.n_clus = 2;
    spurious.clus_id.push_back(cls_lay0.id);
    spurious.clus_id.push_back(cls_lay1.id);
    spurious.tracklet_id.push_back(trkl01.id);
    spurious.theta = theta * 180 / TMath::Pi();
    spurious.phi = phi * 180 / TMath::Pi();
    spurious.err_x0 = -1.;
    spurious.err_y0 = -1.;
    spurious.err_theta = -1.;
    spurious.err_phi = -1.;
    spurious.chi2 = -1.;

    double x1 = (double)cls_lay0.x + display::dist_z * TMath::Tan(theta) * TMath::Cos(phi);
    double y1 = (double)cls_lay0.y + display::dist_z * TMath::Tan(theta) * TMath::Sin(phi);

    if (!display::is_inside_the_layers(x2, y2) && t.track_hit_TR(x1, y1, theta, phi))
    {
      spurious.x0 = cls_lay1.x;
      spurious.y0 = cls_lay1.y;
      spurious.z0 = cls_lay1.z;
      if (cls_lay0.id == cls_lay1.id)
      {
        stats::hmrtar++;
      }
      stats::hmrt++;
      //cout << "010101010100101" << endl;
      tracks.push_back(spurious);
      used_tracklets.push_back(trkl01.id);
      used_clusters.push_back(cls_lay0.id);
      used_clusters.push_back(cls_lay1.id);
    }
  }

  for (auto &trkl12 : tracklet_lay12)
  {
    // Reject tracklets already used in full tracks
    if (std::find_if(used_tracklets.begin(), used_tracklets.end(), [&](int id)
                     { return id == trkl12.id; }) != used_tracklets.end())
    {
      continue;
    }
    // Reject tracklets with clusters already used in full tracks
    if (std::find_if(used_clusters.begin(), used_clusters.end(), [&](int id)
                     { return (id == trkl12.firstClusterId || id == trkl12.secondClusterId); }) != used_clusters.end())
    {
      continue;
    }
    auto cls_lay1 = tidy_clusters_lay1[trkl12.firstClusterId];
    auto cls_lay2 = tidy_clusters_lay2[trkl12.secondClusterId];
    // t.print_tracklet(cls_lay0, cls_lay1);
    double delta_y = (double)cls_lay2.y - (double)cls_lay1.y;
    double delta_x = (double)cls_lay2.x - (double)cls_lay1.x;
    double delta_z = (double)cls_lay2.z - (double)cls_lay1.z;
    double delta_r = TMath::Sqrt(delta_y * delta_y + delta_x * delta_x);
    double r = TMath::Sqrt(pow(delta_x, 2) + pow(delta_y, 2) + pow(delta_z, 2));
    double phi = TMath::ATan2(delta_y, delta_x);
    double theta = TMath::ACos(delta_z / r);
    double x0 = (double)cls_lay1.x - display::dist_z * TMath::Tan(theta) * TMath::Cos(phi);
    double y0 = (double)cls_lay1.y - display::dist_z * TMath::Tan(theta) * TMath::Sin(phi);
    double z0 = display::StaveZ[0];

    LTrackCandidate spurious;
    spurious.n_clus = 2;
    spurious.clus_id.push_back(cls_lay1.id);
    spurious.clus_id.push_back(cls_lay2.id);
    spurious.tracklet_id.push_back(trkl12.id);
    spurious.theta = theta * 180 / TMath::Pi();
    spurious.phi = phi * 180 / TMath::Pi();
    spurious.err_x0 = -1.;
    spurious.err_y0 = -1.;
    spurious.err_theta = -1.;
    spurious.err_phi = -1.;
    spurious.chi2 = -1.;

    if (!display::is_inside_the_layers(x0, y0) && t.track_hit_TR(cls_lay1.x, cls_lay1.y, theta, phi))
    {
      spurious.x0 = cls_lay1.x;
      spurious.y0 = cls_lay1.y;
      spurious.z0 = cls_lay1.z;
      if (cls_lay1.id == cls_lay2.id)
      {
        stats::hmrtar++;
      }
      stats::hmrt++;
      //cout << "212121212121221122121" << endl;
      tracks.push_back(spurious);
      used_tracklets.push_back(trkl12.id);
      used_clusters.push_back(cls_lay1.id);
      used_clusters.push_back(cls_lay2.id);
    }
  }

  for (auto &trkl02 : tracklet_lay02)
  {
    // Reject tracklets already used in full tracks
    if (std::find_if(used_tracklets.begin(), used_tracklets.end(), [&](int id)
                     { return id == trkl02.id; }) != used_tracklets.end())
    {
      continue;
    }
    // Reject tracklets with clusters already used in full tracks
    if (std::find_if(used_clusters.begin(), used_clusters.end(), [&](int id)
                     { return (id == trkl02.firstClusterId || id == trkl02.secondClusterId); }) != used_clusters.end())
    {
      continue;
    }
    auto cls_lay0 = tidy_clusters_lay0[trkl02.firstClusterId];
    auto cls_lay2 = tidy_clusters_lay2[trkl02.secondClusterId];
    // t.print_tracklet(cls_lay0, cls_lay1);
    double delta_y = (double)cls_lay0.y - (double)cls_lay2.y;
    double delta_x = (double)cls_lay0.x - (double)cls_lay2.x;
    double delta_z = (double)cls_lay0.z - (double)cls_lay2.z;
    double delta_r = TMath::Sqrt(delta_y * delta_y + delta_x * delta_x);
    double r = TMath::Sqrt(pow(delta_x, 2) + pow(delta_y, 2) + pow(delta_z, 2));
    double phi = TMath::ATan2(delta_y, delta_x);
    double theta = TMath::ACos(delta_z / r);
    double x1 = (double)cls_lay0.x + display::dist_z * TMath::Tan(theta) * TMath::Cos(phi);
    double y1 = (double)cls_lay0.y + display::dist_z * TMath::Tan(theta) * TMath::Sin(phi);
    double z1 = display::StaveZ[1];

    LTrackCandidate spurious;
    spurious.n_clus = 2;
    spurious.clus_id.push_back(cls_lay0.id);
    spurious.clus_id.push_back(cls_lay2.id);
    spurious.tracklet_id.push_back(trkl02.id);
    spurious.theta = theta * 180 / TMath::Pi();
    spurious.phi = phi * 180 / TMath::Pi();
    spurious.err_x0 = -1.;
    spurious.err_y0 = -1.;
    spurious.err_theta = -1.;
    spurious.err_phi = -1.;
    spurious.chi2 = -1.;

    if (!display::is_inside_the_layers(x1, y1) && t.track_hit_TR(x1, y1, theta, phi))
    {
      spurious.x0 = x1;
      spurious.y0 = y1;
      spurious.z0 = z1;
      if (cls_lay0.id == cls_lay2.id)
      {
        stats::hmrtar++;
      }
      stats::hmrt++;
      //cout << "020202020200202" << endl;
      tracks.push_back(spurious);
      used_tracklets.push_back(trkl02.id);
      used_clusters.push_back(cls_lay0.id);
      used_clusters.push_back(cls_lay2.id);
    }
  }
}

void LTrackerTrack::computeTrackCandidates()
{
  // inside () should have LTrackerCluster &clusterer
  int candidateCounter = 0;

  /*   std::vector<int> clusterer_indices = clusterer.GetClusterIdx();
    int n_cluster_clusterer = (int)clusterer_indices.size();
    clusterer.cls_res_x_m2.resize(n_cluster_clusterer);
    clusterer.cls_res_y_m2.resize(n_cluster_clusterer); */

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
        if (trkCand.clus_id[0] == trkCand.clus_id[1] && trkCand.clus_id[1] == trkCand.clus_id[2] && trkCand.clus_id[0] == trkCand.clus_id[2])
        {
          stats::hmrtar++;
        }

        // compute residuals
        for (auto &clus : clus_vec)
        {
          float l_res_x = clus.x - (trkCand.x0 + (clus.z - display::z_origin_shift) * TMath::Tan(trkCand.theta) * TMath::Cos(trkCand.phi)); // trkCand esce in gradi e qua viene usato in rad
          float l_res_y = clus.y - (trkCand.y0 + (clus.z - display::z_origin_shift) * TMath::Tan(trkCand.theta) * TMath::Sin(trkCand.phi));

          /*           int position = std::distance(clusterer_indices.begin(), std::find(clusterer_indices.begin(), clusterer_indices.end(), clus.id));
                    clusterer.cls_res_x_m2[position] = l_res_x;
                    clusterer.cls_res_y_m2[position] = l_res_y; */
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
    // cout << trk.chi2 << endl;
  }

  // addSpuriousTracks(used_tracklets, used_clusters, tracklet_lay01, tidy_clusters_lay0, tidy_clusters_lay1);
  // addSpuriousTracks(used_tracklets, used_clusters, tracklet_lay12, tidy_clusters_lay1, tidy_clusters_lay2);
  // addSpuriousTracks(used_tracklets, used_clusters, tracklet_lay02, tidy_clusters_lay0, tidy_clusters_lay2);

  // Reassigning track id
  for (int i = 0; i < tracks.size(); i++)
  {
    tracks[i].id = i;
  }

  stats::hmrt = tracks.size();
}

void LTrackerTrack::new_algo(double radius)
{

  double r;
  int candidateCounter = 0;
  cout << "trkl size " << tracklet_lay02.size() << endl;

  for (auto &trkl02 : tracklet_lay02)
  {
    // calcola punto sul layer 1
    LCluster clus_0 = tidy_clusters_lay0[trkl02.firstClusterId];
    LCluster clus_2 = tidy_clusters_lay2[trkl02.secondClusterId];

    // intorno del punto nel quale si cerca un cluster
    double x1, y1, z1, r;
    x1 = (clus_2.x + clus_0.x) / 2;
    y1 = (clus_2.y + clus_0.y) / 2;
    z1 = display::StaveZ[1];
    r = radius;

    for (const auto &tcl1 : tidy_clusters_lay1)
    {
      const LCluster &clus_1 = tcl1.second;
      int key_cls1 = tcl1.first;
      if (clus_1.x < x1 + r + clus_1.errx && clus_1.x > x1 - r - clus_1.errx &&
          clus_1.y < y1 + r + clus_1.erry && clus_1.y > y1 - r - clus_1.erry)
      {
        LTrackerTrack t;
        std::vector<LCluster> clus_vec = {clus_0, clus_1, clus_2};
        LTrackCandidate trkCand;

        fitStraightLine(clus_vec, trkCand);

        trkCand.id = candidateCounter++;
        trkCand.n_clus = clus_vec.size();
        trkCand.tracklet_id = {trkl02.id};
        trkCand.clus_id = {trkl02.firstClusterId, clus_1.id, trkl02.secondClusterId};
        double theta, phi;
        theta = trkCand.theta * TMath::DegToRad();
        phi = trkCand.phi * TMath::DegToRad();

        // check if recotrk passa dai trigger
        if (t.track_hit_TR((double)trkCand.x0, (double)trkCand.y0, theta, phi) && trkCand.chi2 < chi2_cut)
        {
          track_candidates.push_back(trkCand);
          if (clus_0.id == clus_1.id && clus_1.id == clus_2.id && clus_0.id == clus_2.id && trkCand.chi2 < chi2_cut)
          {
            stats::hmrtar++;
          }
        }
      }
    }
  }

  // Sort track candidates by descending chi2
  std::sort(track_candidates.begin(), track_candidates.end(), [](LTrackCandidate &a, LTrackCandidate &b)
            { return a.chi2 < b.chi2; });
  // Remove candidates with large chi2
  double init_size_trkCand = track_candidates.size();
  chi2_cut = 50;
  auto new_end = std::remove_if(track_candidates.begin(), track_candidates.end(), [&](LTrackCandidate &trk)
                                { return trk.chi2 > chi2_cut; });
  track_candidates.erase(new_end, track_candidates.end());

  // Record used tracklets and clusters
  std::vector<int> used_tracklets;
  std::vector<int> used_clusters;

  for (auto &trk : track_candidates)
  {
    tracks.push_back(trk);
    used_tracklets.insert(used_tracklets.end(), trk.tracklet_id.begin(), trk.tracklet_id.end());
    used_clusters.insert(used_clusters.end(), trk.clus_id.begin(), trk.clus_id.end());
  }

  // Reassigning track id
  for (int i = 0; i < tracks.size(); i++)
  {
    tracks[i].id = i;
  }

  stats::hmrt = tracks.size();
  
  New_addSpuriousTracks(used_tracklets, used_clusters);


}

bool LTrackerTrack::track_hit_TR(double x1, double y1, double theta, double phi)
{
  //  +----------+   xTR2t, yTR2t
  //  |   TR2    |
  //  +----------+   xTR2b, yTR2b
  //
  //  +----------+   xTR1t, yTR1t
  //  |   TR1    |
  //  +----------+   xTR1b, yTR1b

  double xTR1b = x1 - (display::StaveZ[1] + display::TR1Thickness / 2) * TMath::Tan(theta) * TMath::Cos(phi);
  double yTR1b = y1 - (display::StaveZ[1] + display::TR1Thickness / 2) * TMath::Tan(theta) * TMath::Sin(phi);
  double xTR1t = x1 - (display::StaveZ[1] - display::TR1Thickness / 2) * TMath::Tan(theta) * TMath::Cos(phi);
  double yTR1t = y1 - (display::StaveZ[1] - display::TR1Thickness / 2) * TMath::Tan(theta) * TMath::Sin(phi);
  double xTR2b = x1 + (-display::StaveZ[1] + display::TR2CenterZ - display::TR2Thickness / 2) * TMath::Tan(theta) * TMath::Cos(phi);
  double yTR2b = y1 + (-display::StaveZ[1] + display::TR2CenterZ - display::TR2Thickness / 2) * TMath::Tan(theta) * TMath::Sin(phi);
  double xTR2t = x1 + (-display::StaveZ[1] + display::TR2CenterZ + display::TR2Thickness / 2) * TMath::Tan(theta) * TMath::Cos(phi);
  double yTR2t = y1 + (-display::StaveZ[1] + display::TR2CenterZ + display::TR2Thickness / 2) * TMath::Tan(theta) * TMath::Sin(phi);

  if ((xTR1b < display::TR1Size[0] / 2 && xTR1b > -display::TR1Size[0] / 2 &&
       ((yTR1b < (2.5 * display::TR1Size[1] + 2 * display::TR1GapY) && yTR1b > (1.5 * display::TR1Size[1] + 2 * display::TR1GapY)) ||
        (yTR1b < (1.5 * display::TR1Size[1] + 1 * display::TR1GapY) && yTR1b > (0.5 * display::TR1Size[1] + 1 * display::TR1GapY)) ||
        (yTR1b < (0.5 * display::TR1Size[1] + 0 * display::TR1GapY) && yTR1b > -(0.5 * display::TR1Size[1] + 0 * display::TR1GapY)) ||
        (yTR1b < -(0.5 * display::TR1Size[1] + 1 * display::TR1GapY) && yTR1b > -(1.5 * display::TR1Size[1] + 1 * display::TR1GapY)) ||
        (yTR1b < -(1.5 * display::TR1Size[1] + 2 * display::TR1GapY) && yTR1b > -(2.5 * display::TR1Size[1] + 2 * display::display::TR1GapY))) &&
       (yTR2b < display::TR2Size[1] / 2 && yTR2b > -display::TR2Size[1] / 2 &&
            (xTR2b < (2 * display::TR2Size[0] + 1.5 * display::TR2GapX) && xTR2b > (1 * display::TR2Size[0] + 1.5 * display::TR2GapX)) ||
        (xTR2b < (1 * display::TR2Size[0] + 0.5 * display::TR2GapX) && xTR2b > (0 * display::TR2Size[0] + 0.5 * display::TR2GapX)) ||
        (xTR2b < -(0 * display::TR2Size[0] + 0.5 * display::TR2GapX) && xTR2b > -(1 * display::TR2Size[0] + 0.5 * display::TR2GapX)) ||
        (xTR2b < -(1 * display::TR2Size[0] + 1.5 * display::TR2GapX) && xTR2b > -(2 * display::TR2Size[0] + 1.5 * display::display::TR2GapX)))) ||
      ((xTR1t < display::TR1Size[0] / 2 && xTR1t > -display::TR1Size[0] / 2 &&
        ((yTR1t < (2.5 * display::TR1Size[1] + 2 * display::TR1GapY) && yTR1t > (1.5 * display::TR1Size[1] + 2 * display::TR1GapY)) ||
         (yTR1t < (1.5 * display::TR1Size[1] + 1 * display::TR1GapY) && yTR1t > (0.5 * display::TR1Size[1] + 1 * display::TR1GapY)) ||
         (yTR1t < (0.5 * display::TR1Size[1] + 0 * display::TR1GapY) && yTR1t > -(0.5 * display::TR1Size[1] + 0 * display::TR1GapY)) ||
         (yTR1t < -(0.5 * display::TR1Size[1] + 1 * display::TR1GapY) && yTR1t > -(1.5 * display::TR1Size[1] + 1 * display::TR1GapY)) ||
         (yTR1t < -(1.5 * display::TR1Size[1] + 2 * display::TR1GapY) && yTR1t > -(2.5 * display::TR1Size[1] + 2 * display::display::TR1GapY))) &&
        (yTR2t < display::TR2Size[1] / 2 && yTR2t > -display::TR2Size[1] / 2 &&
             (xTR2t < (2 * display::TR2Size[0] + 1.5 * display::TR2GapX) && xTR2t > (1 * display::TR2Size[0] + 1.5 * display::TR2GapX)) ||
         (xTR2t < (1 * display::TR2Size[0] + 0.5 * display::TR2GapX) && xTR2t > (0 * display::TR2Size[0] + 0.5 * display::TR2GapX)) ||
         (xTR2t < -(0 * display::TR2Size[0] + 0.5 * display::TR2GapX) && xTR2t > -(1 * display::TR2Size[0] + 0.5 * display::TR2GapX)) ||
         (xTR2t < -(1 * display::TR2Size[0] + 1.5 * display::TR2GapX) && xTR2t > -(2 * display::TR2Size[0] + 1.5 * display::display::TR2GapX))))))
  {
    return true;
  }
  return false;
}

void LTrackerTrack::printRecoTracks_old_alg(TCanvas *reco, int events)
{

  int i = 0;
  for (auto &trk : tracks)
  {
    if (i >= events)
    {
      break;
    }
    ++i;
    float x1, y1, z1, dz, x2, y2, z2, t, p;
    // dz = display::dist_z;
    dz = 100;
    t = trk.theta * TMath::DegToRad();
    p = trk.phi * TMath::DegToRad();
    x2 = trk.x0 + dz * (TMath::Tan(t)) * (TMath::Cos(p));
    y2 = trk.y0 + dz * (TMath::Tan(t)) * (TMath::Sin(p));
    z2 = trk.z0 + dz;
    x1 = trk.x0 - dz * (TMath::Tan(t)) * (TMath::Cos(p));
    y1 = trk.y0 - dz * (TMath::Tan(t)) * (TMath::Sin(p));
    z1 = trk.z0 - dz;
    Double_t x_line[3] = {x1, trk.x0, x2};
    Double_t y_line[3] = {y1, trk.y0, y2};
    Double_t z_line[3] = {z1, trk.z0, z2};
    TPolyLine3D *line_track = new TPolyLine3D(3, x_line, y_line, z_line);
    line_track->SetLineWidth(2);
    line_track->SetLineColor(kRed);
    line_track->Draw();

    TMarker3DBox *g = new TMarker3DBox(x2, y2, z2, 0, 0, 0, 0, 0);
    g->Draw();
    TMarker3DBox *m = new TMarker3DBox(trk.x0, trk.y0, trk.z0, 0, 0, 0, 0, 0);
    m->Draw();
    TMarker3DBox *f = new TMarker3DBox(x1, y1, z1, 0, 0, 0, 0, 0);
    f->Draw();

    /*     cout << "L2 ------ x:" << x2 << ",   y: " << y2 << ",    z: " << z2 << endl;
        cout << "L1 ------ x:" << trk.x0 << "+- " << trk.err_x0 << ",   y: " << trk.y0 << "+- " << trk.err_y0 << ",    z: " << trk.z0 << endl;
        cout << "L0 ------ x:" << x1 << ",   y: " << y1 << ",    z: " << z1 << endl;
        cout << "(gradi) theta: " << trk.theta << "     phi: " << trk.phi << endl;
        cout << "(rad)   theta: " << t         << "     phi: " <<     p   << endl;
        cout << "~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;  */
  }

  reco->Update();
}

void LTrackerTrack::printRecoTracks_new_alg(TCanvas *reco)
{

  int i = 0;
  cout << "trakcs size" << tracks.size() << endl;
  for (auto &trk : tracks)
  {
    float x1, y1, z1, dz, x2, y2, z2, t, p;

    dz = 100;

    t = trk.theta * TMath::DegToRad();
    p = trk.phi * TMath::DegToRad();
    x2 = trk.x0 + dz * (TMath::Tan(t)) * (TMath::Cos(p));
    y2 = trk.y0 + dz * (TMath::Tan(t)) * (TMath::Sin(p));
    z2 = trk.z0 + dz;
    x1 = trk.x0 - dz * (TMath::Tan(t)) * (TMath::Cos(p));
    y1 = trk.y0 - dz * (TMath::Tan(t)) * (TMath::Sin(p));
    z1 = trk.z0 - dz;
    Double_t x_line[3] = {x1, trk.x0, x2};
    Double_t y_line[3] = {y1, trk.y0, y2};
    Double_t z_line[3] = {z1, trk.z0, z2};
    TPolyLine3D *line_track = new TPolyLine3D(3, x_line, y_line, z_line);
    line_track->SetLineWidth(2);
    line_track->SetLineColor(kRed);
    line_track->Draw();

    TMarker3DBox *g = new TMarker3DBox(x2, y2, z2, 0, 0, 0, 0, 0);
    g->Draw();
    TMarker3DBox *m = new TMarker3DBox(trk.x0, trk.y0, trk.z0, 0, 0, 0, 0, 0);
    m->Draw();
    TMarker3DBox *f = new TMarker3DBox(x1, y1, z1, 0, 0, 0, 0, 0);
    f->Draw();

    /*     cout << "x: " << trk.x0 << "+- " << trk.err_x0 << "| y: " << trk.y0 << "+- " << trk.err_y0 << "| z: " << trk.z0 << endl;
        cout << "theta  " << t << "phi " << p << endl;
        cout << "theta  " << trk.theta << "phi " << trk.phi << endl;
        cout << "bottom " << x1 << y1 << z1 << endl;
        cout << "top " << x2 << y2 << z2 << endl; */
  }

  reco->Update();
}

std::ostream &operator<<(std::ostream &output, const LTrackerTrack &tracker)
{
  output << "CLUSTERS" << std::endl
         << std::endl;
  output << "-------------------" << std::endl;
  output << "layer 0 " << std::endl
         << std::endl;
  for (auto &cl : tracker.tidy_clusters_lay0)
  {
    output << cl.second << std::endl;
  }

  output << "-------------------" << std::endl;
  output << "layer 1 " << std::endl
         << std::endl;
  for (auto &cl : tracker.tidy_clusters_lay1)
  {
    output << cl.second << std::endl;
  }
  output << "-------------------" << std::endl;
  output << "layer 2 " << std::endl
         << std::endl;
  for (auto &cl : tracker.tidy_clusters_lay2)
  {
    output << cl.second << std::endl;
  }
  /*
  output << std::endl
         << "TRACKLETS" << std::endl
         << std::endl;
  output << "-------------------" << std::endl;
  output << "layer 0 - layer 1" << std::endl
         << std::endl;
  for (auto &trkl : tracker.tracklet_lay01)
  {
    output << trkl << std::endl;
  }
  output << "-------------------" << std::endl;
  output << "layer 1 - layer 2" << std::endl
         << std::endl;
  for (auto &trkl : tracker.tracklet_lay12)
  {
    output << trkl << std::endl;
  }
  output << "-------------------" << std::endl;
  output << "layer 0 - layer 2 " << std::endl
         << std::endl;
  for (auto &trkl : tracker.tracklet_lay02)
  {
    output << trkl << std::endl;
  }
    */
  return output;
}
