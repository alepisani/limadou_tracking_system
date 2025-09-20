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
  std::unordered_map<int, LCluster>().swap(tidy_clusters_lay0);
  std::unordered_map<int, LCluster>().swap(tidy_clusters_lay1);
  std::unordered_map<int, LCluster>().swap(tidy_clusters_lay2);
  std::vector<LTracklet>().swap(tracklet_lay01);
  std::vector<LTracklet>().swap(tracklet_lay12);
  std::vector<LTracklet>().swap(tracklet_lay02);
  std::vector<LTrackCandidate>().swap(track_candidates);
  std::vector<LTrackCandidate>().swap(tracks);
  std::vector<int>().swap(used_clusters_lay0);
  std::vector<int>().swap(used_clusters_lay1);
  std::vector<int>().swap(used_clusters_lay2);
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
}

void LTrackerTrack::print_tracklet(const LCluster cl_0, const LCluster cl_2)
{

  float x0, y0, z0, x2, y2, z2;
  float fx0, fy0, fz0, fx2, fy2, fz2;
  x0 = cl_0.x;
  y0 = cl_0.y;
  z0 = cl_0.z;
  x2 = cl_2.x;
  y2 = cl_2.y;
  z2 = cl_2.z;
  // float dx = x2 - x0;
  // float dy = y2 - y0;
  // fx0 = x0 - dx;
  // fy0 = y0 - dy;
  // fz0 = z0 - display::dist_z;
  // fx2 = x2 + dx;
  // fy2 = y2 + dy;
  // fz2 = z2 + display::dist_z;

  // Double_t x_line[4] = {fx0, x0, x2, fx2};
  // Double_t y_line[4] = {fy0, y0, y2, fy2};
  // Double_t z_line[4] = {fz0, z0, z2, fz2};
  Double_t x_line[2] = {x0, x2};
  Double_t y_line[2] = {y0, y2};
  Double_t z_line[2] = {z0, z2};
  TPolyLine3D *trk = new TPolyLine3D(2, x_line, y_line, z_line);
  trk->SetLineWidth(1);
  trk->SetLineColor(kBlue);
  trk->Draw();
}

void LTrackerTrack::print_all_tracklet(LTrackerTrack ltt)
{

  for (auto &trkl01 : ltt.tracklet_lay01)
  {
    auto cls_lay0 = tidy_clusters_lay0[trkl01.firstClusterId];
    auto cls_lay1 = tidy_clusters_lay1[trkl01.secondClusterId];
    print_tracklet(cls_lay0, cls_lay1);
  }

  for (auto &trkl12 : ltt.tracklet_lay12)
  {
    auto cls_lay1 = tidy_clusters_lay1[trkl12.firstClusterId];
    auto cls_lay2 = tidy_clusters_lay2[trkl12.secondClusterId];
    print_tracklet(cls_lay1, cls_lay2);
  }

  for (auto &trkl02 : ltt.tracklet_lay02)
  {
    auto cls_lay0 = tidy_clusters_lay0[trkl02.firstClusterId];
    auto cls_lay2 = tidy_clusters_lay2[trkl02.secondClusterId];
    print_tracklet(cls_lay0, cls_lay2);
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
  float pi = TMath::Pi();

  if (clusters.size() < 3)
    return;

  auto minimize = [&](double *initialVars, double *steps, const char *tag) -> std::tuple<const double *, const double *, double>
  {
    // could change Minuit2 with other minimizer such as minuit or pleanty other
    ROOT::Math::Minimizer *min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "");
    min->SetMaxFunctionCalls(1000000);
    min->SetMaxIterations(10000);
    min->SetTolerance(0.0001);
    min->SetPrintLevel(0);

    std::function<double(const std::vector<LCluster> &, const double *)> unboundFct = fct;
    std::function<double(const double *)> boundFct = std::bind(unboundFct, clusters, std::placeholders::_1);
    ROOT::Math::Functor f(boundFct, 4);
    min->SetFunction(f);

    min->SetVariable(0, "x0", initialVars[0], steps[0]);
    min->SetVariable(1, "y0", initialVars[1], steps[1]);

    // we're using a wider range because we are trying to optimise a circular function.
    // maybe the function goes to one end but the optimal fit value could be just on the edge of the other end.
    // the real theta distribution is between (0,pi/2), now is (-pi/2, pi);
    // the real  phi distribution is between (-pi,+pi), now is (-2pi, 2pi);

    // what real intervals should look like
    min->SetLimitedVariable(2, "theta", initialVars[2], steps[2], 0., pi / 2);
    min->SetLimitedVariable(3, "phi", initialVars[3], steps[3], -1 * pi, 1 * pi);
    min->Minimize();
    return {min->X(), min->Errors(), min->MinValue()};
  };

  // inizialise the variable for the fit function gives better accuracy in the algorithm

  double dx = clusters[2].x - clusters[0].x;
  double dy = clusters[2].y - clusters[0].y;
  double dz = clusters[2].z - clusters[0].z;
  double r = sqrt(dx * dx + dy * dy + dz * dz);
  double guess_theta = acos(dz / r);
  double guess_phi = atan2(dy, dx);

  double vars[4] = {clusters[0].x, clusters[0].y, guess_theta, guess_phi};
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
  trkCand.err_x0 = errors[0];
  trkCand.err_y0 = errors[1];
  trkCand.err_theta = errors[2];
  trkCand.err_phi = errors[3];
  trkCand.chi2 = chi2;
  trkCand.theta = params[2];
  trkCand.phi = params[3];

  // if(trkCand.phi < -pi) trkCand.phi += 2 * pi;
  // if(trkCand.phi > +pi) trkCand.phi -= 2 * pi;
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
    // spurious.theta = TMath::ATan(tan_theta) * 180 / TMath::Pi();
    // spurious.phi = TMath::ACos(cos_phi) * 180 / TMath::Pi();
    spurious.theta = TMath::ATan(tan_theta);
    spurious.phi = TMath::ACos(cos_phi);
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
  tracks.reserve(tracklet_lay01.size() + tracklet_lay02.size() + tracklet_lay12.size());

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
    double r = TMath::Sqrt(delta_x * delta_x + delta_y * delta_y + delta_z * delta_z);
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
    // spurious.theta = theta * 180 / TMath::Pi();
    // spurious.phi = phi * 180 / TMath::Pi();
    spurious.theta = theta;
    spurious.phi = phi;
    spurious.err_x0 = -1.;
    spurious.err_y0 = -1.;
    spurious.err_theta = -1.;
    spurious.err_phi = -1.;
    spurious.chi2 = -1.;

    double x1 = (double)cls_lay0.x + display::dist_z * TMath::Tan(theta) * TMath::Cos(phi);
    double y1 = (double)cls_lay0.y + display::dist_z * TMath::Tan(theta) * TMath::Sin(phi);

    if (!display::is_inside_the_layers(&x2, &y2) && t.track_hit_TR(x1, y1, theta, phi))
    {
      spurious.x0 = cls_lay1.x;
      spurious.y0 = cls_lay1.y;
      spurious.z0 = cls_lay1.z;
      if (cls_lay0.id == cls_lay1.id)
      {
        stats::hmrtar++;
        stats::hmrtar2++;
      }
      else
        stats::hmrtaf2++;
      stats::hmrt++;
      // cout << "010101010100101" << endl;
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
    double r = TMath::Sqrt(delta_x * delta_x + delta_y * delta_y + delta_z * delta_z);
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
    // spurious.theta = theta * 180 / TMath::Pi();
    // spurious.phi = phi * 180 / TMath::Pi();
    spurious.theta = theta;
    spurious.phi = phi;
    spurious.err_x0 = -1.;
    spurious.err_y0 = -1.;
    spurious.err_theta = -1.;
    spurious.err_phi = -1.;
    spurious.chi2 = -1.;

    if (!display::is_inside_the_layers(&x0, &y0) && t.track_hit_TR(cls_lay1.x, cls_lay1.y, theta, phi))
    {
      spurious.x0 = cls_lay1.x;
      spurious.y0 = cls_lay1.y;
      spurious.z0 = cls_lay1.z;
      if (cls_lay1.id == cls_lay2.id)
      {
        stats::hmrtar++;
        stats::hmrtar2++;
      }
      else
        stats::hmrtaf2++;
      stats::hmrt++;
      // cout << "212121212121221122121" << endl;
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
    double delta_y = (double)cls_lay2.y - (double)cls_lay0.y;
    double delta_x = (double)cls_lay2.x - (double)cls_lay0.x;
    double delta_z = (double)cls_lay2.z - (double)cls_lay0.z;
    double delta_r = TMath::Sqrt(delta_y * delta_y + delta_x * delta_x);
    double r = TMath::Sqrt(delta_x * delta_x + delta_y * delta_y + delta_z * delta_z);
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
    // spurious.theta = theta * 180 / TMath::Pi();
    // spurious.phi = phi * 180 / TMath::Pi();
    spurious.theta = theta;
    spurious.phi = phi;
    spurious.err_x0 = -1.;
    spurious.err_y0 = -1.;
    spurious.err_theta = -1.;
    spurious.err_phi = -1.;
    spurious.chi2 = -1.;

    if (!display::is_inside_the_layers(&x1, &y1) && t.track_hit_TR(x1, y1, theta, phi))
    {
      spurious.x0 = x1;
      spurious.y0 = y1;
      spurious.z0 = z1;
      if (cls_lay0.id == cls_lay2.id)
      {
        stats::hmrtar++;
        stats::hmrtar2++;
      }
      else
        stats::hmrtaf2++;
      stats::hmrt++;
      // cout << "020202020200202" << endl;
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
        else
          stats::hmrtaf3++;

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
  addSpuriousTracks(used_tracklets, used_clusters, tracklet_lay02, tidy_clusters_lay0, tidy_clusters_lay2);

  // Reassigning track id
  for (int i = 0; i < tracks.size(); i++)
  {
    tracks[i].id = i;
  }

  stats::hmrt = tracks.size();
}

void LTrackerTrack::new_algo(double radius)
{
  chi2_cut = 5000;
  float degtorad = TMath::DegToRad();
  float pi = TMath::Pi();
  int candidateCounter = 0;
  std::vector<LCluster> clus_vec;
  clus_vec.reserve(3);
  track_candidates.clear();
  tracks.clear();

  for (const auto &trkl02 : tracklet_lay02)
  {
    // calcola punto sul layer 1
    // define the object inside the for cicle to improve memory allocation, at the end of each iteration will be naturally deleated
    const LCluster &clus_0 = tidy_clusters_lay0[trkl02.firstClusterId];
    const LCluster &clus_2 = tidy_clusters_lay2[trkl02.secondClusterId];

    // intorno del punto nel quale si cerca un cluster
    double x1 = (clus_2.x + clus_0.x) * 0.5; // moltiplications are faster than divisions
    double y1 = (clus_2.y + clus_0.y) * 0.5;
    double r = radius;

    for (const auto &tcl1 : tidy_clusters_lay1)
    {
      const LCluster &clus_1 = tcl1.second;
      const int &key_cls1 = tcl1.first;
      double dx = std::abs(clus_1.x - x1);
      double dy = std::abs(clus_1.y - y1);
      if (dx < radius && dy < radius)
      {
        LTrackerTrack t;
        LTrackCandidate trkCand;
        clus_vec.clear();
        clus_vec.assign({clus_0, clus_1, clus_2});

        fitStraightLine(clus_vec, trkCand);

        trkCand.id = candidateCounter++;
        trkCand.n_clus = clus_vec.size();
        trkCand.tracklet_id = {trkl02.id};
        trkCand.clus_id = {trkl02.firstClusterId, clus_1.id, trkl02.secondClusterId};

        // check if recotrk passa dai trigger
        if (t.track_hit_TR((double)trkCand.x0, (double)trkCand.y0, trkCand.theta, trkCand.phi) && trkCand.chi2 < chi2_cut)
        {
          track_candidates.push_back(trkCand);
          if (clus_0.id == clus_1.id && clus_1.id == clus_2.id && clus_0.id == clus_2.id)
          {
            stats::hmrtar++;
            stats::hmrtar3++;
          }
          else
            stats::hmrtaf3++;
        }
      }
    }
  }

  // Record used tracklets and clusters
  std::vector<int> used_tracklets;
  std::vector<int> used_clusters;
  used_tracklets.reserve(tracklet_lay01.size() + tracklet_lay02.size() + tracklet_lay12.size());
  used_clusters.reserve(used_clusters_lay0.size() + used_clusters_lay1.size() + used_clusters_lay2.size());

  for (auto &trk : track_candidates)
  {
    tracks.push_back(trk);
    used_tracklets.insert(used_tracklets.end(), trk.tracklet_id.begin(), trk.tracklet_id.end());
    used_clusters.insert(used_clusters.end(), trk.clus_id.begin(), trk.clus_id.end());
  }

  New_addSpuriousTracks(used_tracklets, used_clusters);

  // Reassigning track id
  for (int i = 0; i < tracks.size(); i++)
  {
    tracks[i].id = i;
    // printf("chi2: %f\n", tracks[i].chi2);
  }

  stats::hmrt = tracks.size();
  // printf("tracks %ld\n", tracks.size());
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

  if (
      // --- Bottom face ---
      ((xTR1b < display::TR1Size[0] / 2 && xTR1b > -display::TR1Size[0] / 2 &&
        ((yTR1b < (2.5 * display::TR1Size[1] + 2 * display::TR1GapY) && yTR1b > (1.5 * display::TR1Size[1] + 2 * display::TR1GapY)) ||
         (yTR1b < (1.5 * display::TR1Size[1] + 1 * display::TR1GapY) && yTR1b > (0.5 * display::TR1Size[1] + 1 * display::TR1GapY)) ||
         (yTR1b < (0.5 * display::TR1Size[1] + 0 * display::TR1GapY) && yTR1b > -(0.5 * display::TR1Size[1] + 0 * display::TR1GapY)) ||
         (yTR1b < -(0.5 * display::TR1Size[1] + 1 * display::TR1GapY) && yTR1b > -(1.5 * display::TR1Size[1] + 1 * display::TR1GapY)) ||
         (yTR1b < -(1.5 * display::TR1Size[1] + 2 * display::TR1GapY) && yTR1b > -(2.5 * display::TR1Size[1] + 2 * display::TR1GapY)))) &&
       (yTR2b < display::TR2Size[1] / 2 && yTR2b > -display::TR2Size[1] / 2 &&
        ((xTR2b < (2 * display::TR2Size[0] + 1.5 * display::TR2GapX) && xTR2b > (1 * display::TR2Size[0] + 1.5 * display::TR2GapX)) ||
         (xTR2b < (1 * display::TR2Size[0] + 0.5 * display::TR2GapX) && xTR2b > (0 * display::TR2Size[0] + 0.5 * display::TR2GapX)) ||
         (xTR2b < -(0 * display::TR2Size[0] + 0.5 * display::TR2GapX) && xTR2b > -(1 * display::TR2Size[0] + 0.5 * display::TR2GapX)) ||
         (xTR2b < -(1 * display::TR2Size[0] + 1.5 * display::TR2GapX) && xTR2b > -(2 * display::TR2Size[0] + 1.5 * display::TR2GapX))))) ||
      // --- Top face ---
      ((xTR1t < display::TR1Size[0] / 2 && xTR1t > -display::TR1Size[0] / 2 &&
        ((yTR1t < (2.5 * display::TR1Size[1] + 2 * display::TR1GapY) && yTR1t > (1.5 * display::TR1Size[1] + 2 * display::TR1GapY)) ||
         (yTR1t < (1.5 * display::TR1Size[1] + 1 * display::TR1GapY) && yTR1t > (0.5 * display::TR1Size[1] + 1 * display::TR1GapY)) ||
         (yTR1t < (0.5 * display::TR1Size[1] + 0 * display::TR1GapY) && yTR1t > -(0.5 * display::TR1Size[1] + 0 * display::TR1GapY)) ||
         (yTR1t < -(0.5 * display::TR1Size[1] + 1 * display::TR1GapY) && yTR1t > -(1.5 * display::TR1Size[1] + 1 * display::TR1GapY)) ||
         (yTR1t < -(1.5 * display::TR1Size[1] + 2 * display::TR1GapY) && yTR1t > -(2.5 * display::TR1Size[1] + 2 * display::TR1GapY)))) &&
       (yTR2t < display::TR2Size[1] / 2 && yTR2t > -display::TR2Size[1] / 2 &&
        ((xTR2t < (2 * display::TR2Size[0] + 1.5 * display::TR2GapX) && xTR2t > (1 * display::TR2Size[0] + 1.5 * display::TR2GapX)) ||
         (xTR2t < (1 * display::TR2Size[0] + 0.5 * display::TR2GapX) && xTR2t > (0 * display::TR2Size[0] + 0.5 * display::TR2GapX)) ||
         (xTR2t < -(0 * display::TR2Size[0] + 0.5 * display::TR2GapX) && xTR2t > -(1 * display::TR2Size[0] + 0.5 * display::TR2GapX)) ||
         (xTR2t < -(1 * display::TR2Size[0] + 1.5 * display::TR2GapX) && xTR2t > -(2 * display::TR2Size[0] + 1.5 * display::TR2GapX))))))
  {
    return true;
  }
  return false;
}

void LTrackerTrack::printRecoTracks_new_alg(TCanvas *reco)
{

  int i = 0;
  for (auto &trk : tracks)
  {
    float x1, y1, z1, dz, x2, y2, z2;

    dz = 100;

    printf("x0 = %f, y0 = %f, theta_reco = %f, phi_reco = %f\n", trk.x0, trk.y0, trk.theta * TMath::RadToDeg(), trk.phi * TMath::RadToDeg());
    x2 = trk.x0 + dz * (TMath::Tan(trk.theta)) * (TMath::Cos(trk.phi));
    y2 = trk.y0 + dz * (TMath::Tan(trk.theta)) * (TMath::Sin(trk.phi));
    z2 = trk.z0 + dz;
    x1 = trk.x0 - dz * (TMath::Tan(trk.theta)) * (TMath::Cos(trk.phi));
    y1 = trk.y0 - dz * (TMath::Tan(trk.theta)) * (TMath::Sin(trk.phi));
    z1 = trk.z0 - dz;
    Double_t x_line[3] = {x1, trk.x0, x2};
    Double_t y_line[3] = {y1, trk.y0, y2};
    Double_t z_line[3] = {z1, trk.z0, z2};
    TPolyLine3D *line_track = new TPolyLine3D(3, x_line, y_line, z_line);
    line_track->SetLineWidth(2);
    if (trk.chi2 > 10)
      line_track->SetLineColor(kGreen);
    else
      line_track->SetLineColor(kRed);
    line_track->Draw();

    TMarker3DBox *g = new TMarker3DBox(x2, y2, z2, 0, 0, 0, 0, 0);
    g->Draw();
    TMarker3DBox *m = new TMarker3DBox(trk.x0, trk.y0, trk.z0, 0, 0, 0, 0, 0);
    m->Draw();
    TMarker3DBox *f = new TMarker3DBox(x1, y1, z1, 0, 0, 0, 0, 0);
    f->Draw();

    if (trk.chi2 > 0.)
    {
      double R = 5;      // radius of the circle
      const int N = 100; // number of points to make circle smooth
      Double_t x_circ[N], y_circ[N], z_circ[N];

      for (int j = 0; j < N; j++)
      {
        double phi = 2 * TMath::Pi() * j / (N - 1); // angle
        x_circ[j] = trk.x0 + R * TMath::Cos(phi);
        y_circ[j] = trk.y0 + R * TMath::Sin(phi);
        z_circ[j] = trk.z0; // circle in XY plane
      }

      // make polyline
      TPolyLine3D *circle = new TPolyLine3D(N, x_circ, y_circ, z_circ);
      circle->SetLineColor(kBlue);
      circle->SetLineWidth(2);
      circle->Draw();
    }
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
