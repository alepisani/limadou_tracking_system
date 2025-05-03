#ifndef LTRACKERTRACK_H
#define LTRACKERTRACK_H

#include <vector>
#include <unordered_map>
#include <ostream>
#include "TObject.h"

struct LCluster
{
    LCluster();

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

#endif