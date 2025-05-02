#include <string>
#include <array> 
#include <map>
#include <iostream>
#include <cmath>
#include <TApplication.h>
#include <TCanvas.h>
#include <TView.h>
#include <TList.h>
#include <TPolyLine3D.h>
#include "TH1F.h"
#include "TLine.h"
#include "TRandom3.h"
#include "TMarker3DBox.h"
#include "../build/stats.h"
using namespace std;

static double stats::hmgt = 0;
static double stats::hmgthTR1 = 0;
static double stats::hmgthL2 = 0;
static double stats::hmgthL1 = 0;
static double stats::hmgthL0 = 0;
static double stats::hmrt = 0;
static double stats::hmgthL012 = 0;