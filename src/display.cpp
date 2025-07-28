#include <string>
#include <array> 
#include <map>
#include <iostream>
#include <cmath>
#include <TApplication.h>
#include <TCanvas.h>
#include <TView.h>
#include <TList.h>
#include <TTree.h>
#include <TPolyLine3D.h>
#include "TH1F.h"
#include "TRandom3.h"
#include "TMarker3DBox.h"
#include "TFile.h"
#include "../include/display.h"
#include "../include/stats.h"
#include "../include/chip.h"
#include "../include/LTrackerTrack.h"
#include "display.h"
using namespace std;

display::display() {}


void display::draw_TR12(TCanvas* geom){

    geom->cd();
    //draw trigger 1
    TMarker3DBox *TR1_0 = new TMarker3DBox(0,0,0,TR1Size[0]/2,TR1Size[1]/2,TR1Size[2]/2,0,0);
    TR1_0->SetLineColor(kRed);
    TR1_0->SetLineWidth(3);
    TR1_0->Draw();
    TMarker3DBox *TR1_1 = new TMarker3DBox(0,TR1Size[1]+TR1GapY,0,TR1Size[0]/2,TR1Size[1]/2,TR1Size[2]/2,0,0);
    TR1_1->SetLineColor(kRed);
    TR1_1->SetLineWidth(3);
    TR1_1->Draw();
    TMarker3DBox *TR1_2 = new TMarker3DBox(0,2*(TR1Size[1]+TR1GapY),0,TR1Size[0]/2,TR1Size[1]/2,TR1Size[2]/2,0,0);
    TR1_2->SetLineColor(kRed);
    TR1_2->SetLineWidth(3);
    TR1_2->Draw();
    TMarker3DBox *TR1_3 = new TMarker3DBox(0,-(TR1Size[1]+TR1GapY),0,TR1Size[0]/2,TR1Size[1]/2,TR1Size[2]/2,0,0);
    TR1_3->SetLineColor(kRed);
    TR1_3->SetLineWidth(3);
    TR1_3->Draw();
    TMarker3DBox *TR1_4 = new TMarker3DBox(0,-2*(TR1Size[1]+TR1GapY),0,TR1Size[0]/2,TR1Size[1]/2,TR1Size[2]/2,0,0);
    TR1_4->SetLineColor(kRed);
    TR1_4->SetLineWidth(3);
    TR1_4->Draw();

    //draw trigger 2
    TMarker3DBox *TR2_0 = new TMarker3DBox(-1.5*(TR2Size[0]+TR2GapX),0,TR2CenterZ, TR2Size[0]/2, TR2Size[1]/2, TR2Size[2]/2, 0, 0);
    TR2_0->SetLineColor(kRed);
    TR2_0->SetLineWidth(3);
    TR2_0->Draw();
    TMarker3DBox *TR2_1 = new TMarker3DBox(-0.5*(TR2Size[0]+TR2GapX),0,TR2CenterZ, TR2Size[0]/2, TR2Size[1]/2, TR2Size[2]/2, 0, 0);
    TR2_1->SetLineColor(kRed);
    TR2_1->SetLineWidth(3);
    TR2_1->Draw();
    TMarker3DBox *TR2_2 = new TMarker3DBox(0.5*(TR2Size[0]+TR2GapX),0,TR2CenterZ, TR2Size[0]/2, TR2Size[1]/2, TR2Size[2]/2, 0, 0);
    TR2_2->SetLineColor(kRed);
    TR2_2->SetLineWidth(3);
    TR2_2->Draw();
    TMarker3DBox *TR2_3 = new TMarker3DBox(1.5*(TR2Size[0]+TR2GapX),0,TR2CenterZ, TR2Size[0]/2, TR2Size[1]/2, TR2Size[2]/2, 0, 0);
    TR2_3->SetLineColor(kRed);
    TR2_3->SetLineWidth(3);
    TR2_3->Draw();
    geom->Update();
    
}

void display::take_distributions(){
    
    TFile *file = TFile::Open("../data_beam_test/TEST_MUONS_m_MAIN_1000.0MeV_-999.0deg_-0.05V_boot207_run510_L2.root");
    if (!file || file->IsZombie()) {
        std::cerr << "Errore nell'aprire il file ROOT\n";
    }

    // Carica il TTree (sostituisci con il nome corretto se diverso)
    TTree *tree = (TTree*)file->Get("L2;1");
    if (!tree) {
        std::cerr << "TTree non trovato nel file\n";
    }

    tree->SetBranchAddress("theta", &theta);
    tree->SetBranchAddress("phi", &phi);
    tree->SetBranchAddress("cls_size", &cls_size);
    Long64_t n = tree->GetEntries();

    for (Long64_t i = 0; i < n; ++i) {
        tree->GetEntry(i);
        if (theta) {
            allTheta.insert(allTheta.end(), theta->begin(), theta->end());
        }
        if (phi) {
            allPhi.insert(allPhi.end(), phi->begin(), phi->end());
        }
        if (cls_size) {
            all_cls_size.insert(all_cls_size.end(), cls_size->begin(), cls_size->end());
        }
    }

    /*
    cout << "size theta: " << allTheta.size() << "size phi: " << allPhi.size();
    for (int i = 0; i < allTheta.size(); i++){
        cout << allTheta[i] << ", " << allPhi[i] << "           ";

    }
    */
}

void display::tracks(int events, LTrackerTrack& tracker, TCanvas* geom){

stats::hmgt = events;
int nbins = events;

//TH1F(name, title, nbins, xlow, xup)
TH1F* hxTR2 = new TH1F("hxTR2", "x_trigger2;x_trigger1;counts", nbins, -TR2Size[0]*2.5, TR2Size[0]*2.5);
TH1F* hyTR2 = new TH1F("hyTR2", "y_trigger2;y_trigger1;counts", nbins, -TR2Size[1]/2, TR2Size[1]/2);
TH1F* hzTR2 = new TH1F("hzTR2", "z_trigger2;z_trigger2;counts", nbins, -TR2Size[2]/2+TR2CenterZ, TR2Size[2]/2+TR2CenterZ);
TH1F* hxTR1 = new TH1F("hxTR1", "x_trigger1;x_trigger1;counts", nbins, -TR1Size[0]/2, TR1Size[0]/2);
TH1F* hyTR1 = new TH1F("hyTR1", "y_trigger1;y_trigger1;counts", nbins, -TR1Size[1]*3, TR1Size[1]*3);
TH1F* hzTR1 = new TH1F("hzTR1", "z_trigger1;z_trigger1;counts", nbins, -TR1Size[2]/2+TR1CenterZ, TR1Size[2]/2+TR1CenterZ);
TH1F* hphi = new TH1F("hphi", "phi;#phi;counts", nbins, -pi, pi);
TH1F* htheta = new TH1F("htheta", "theta;#theta;counts", nbins, 0, pi/2);
TH1F* htheta3L = new TH1F("htheta3L", "theta;#theta;counts", nbins, 0, pi/2);
TH1F* hphi3L = new TH1F("hphi3L", "phi;#phi;counts", nbins, -pi, pi);
TH1F* htheta2L = new TH1F("htheta2L", "theta;#theta;counts", nbins, 0, pi/2);
TH1F* hphi2L = new TH1F("hphi2L", "phi;#phi;counts", nbins, -pi, pi);
TH1F* htheta1L = new TH1F("htheta1L", "theta;#theta;counts", nbins, 0, pi/2);
TH1F* hphi1L = new TH1F("hphi1L", "phi;#phi;counts", nbins, -pi, pi);
TH1F* htheta0L = new TH1F("htheta0L", "theta;#theta;counts", nbins, 0, pi/2);
TH1F* hphi0L = new TH1F("hphi0L", "phi;#phi;counts", nbins, -pi, pi);
TH1F* hx = new TH1F("x", "x;x;counts", nbins, -TR2Size[0]*2.5, TR2Size[0]*2.5);
TH1F* hy = new TH1F("y", "y;y;counts", nbins, -100, 100);

//puoi settare il seed for reproducibility
//TRandom3 *rnd = new TRandom3(45133); 
TRandom3 *rnd = new TRandom3(0); 



//MC
for (int i=0; i < events; i++){
    double xTR2, yTR2, zTR2, xTR1, yTR1, zTR1;
    double xL2, xL1, xL0, yL2, yL1, yL0;
    float phi, theta, cls_size, cls_size_x, cls_size_y; 
    double smiring_xL2, smiring_xL1, smiring_xL0, smiring_yL2, smiring_yL1, smiring_yL0;
    stats::hitL0=false;
    stats::hitL1=false;
    stats::hitL2=false;

    //i layer acquisiscono il segnale quando arriva un segnale AND dagli scintillatori
    //genero traccie misurabili come traccie che passano nei due scintillatori TR2, TR1
    
    do {
        double xTR2_fake = rnd->Uniform(-TR2Size[0]*2,TR2Size[0]*2);     
        if(xTR2_fake>0 && xTR2_fake<TR2Size[0]){xTR2 = xTR2_fake + 0.5*TR2GapX;}
        if(xTR2_fake<0 && xTR2_fake>-TR2Size[0]){xTR2 = xTR2_fake - 0.5*TR2GapX;}
        if(xTR2_fake<2*TR2Size[0] && xTR2_fake>TR2Size[0]){xTR2 = xTR2_fake + 1.5*TR2GapX;}
        if(xTR2_fake>-2*TR2Size[0] && xTR2_fake<-TR2Size[0]){xTR2 = xTR2_fake - 1.5*TR2GapX;}
        yTR2 = rnd->Uniform(-TR2Size[1]/2,TR2Size[1]/2);
        zTR2 = rnd->Uniform(TR2CenterZ-TR2Thickness/2,TR2CenterZ+TR2Thickness/2);
        zTR1 = rnd->Uniform(TR1CenterZ-TR1Thickness/2,TR1CenterZ+TR1Thickness/2);

        //phi = rnd->Uniform(-pi,pi);
        
        /*
        double THETA; double y;
        do {
            THETA = gRandom->Uniform(0, TMath::Pi()/2);
            y = gRandom->Uniform(0, 1);
        } while (y > TMath::Sin(THETA) * TMath::Cos(THETA) * TMath::Cos(THETA));
        theta = THETA;
        */
        //theta = rnd->Uniform(0,pi/2);

        //prendo gli angoli dalla distribuzione di muoni
        int ind_theta = rnd->Uniform(0, allTheta.size()-1);
        int ind_phi = rnd->Uniform(0, allPhi.size()-1);
        int ind_cls_size = rnd->Uniform(0, all_cls_size.size()-1);
        theta = (allTheta[ind_theta]/180)*TMath::Pi();
        phi = (allPhi[ind_phi]/180)*TMath::Pi();
        cls_size = (all_cls_size[ind_cls_size]);
        cls_size_x = cls_size*PixelSizeRows;
        cls_size_y = cls_size*PixelSizeCols;

        xL2 = xTR2 + (zTR2-StaveZ[2])*(TMath::Tan(theta))*(TMath::Cos(phi));
        yL2 = yTR2 + (zTR2-StaveZ[2])*(TMath::Tan(theta))*(TMath::Sin(phi));
        xL1 = xTR2 + (zTR2-StaveZ[1])*(TMath::Tan(theta))*(TMath::Cos(phi));
        yL1 = yTR2 + (zTR2-StaveZ[1])*(TMath::Tan(theta))*(TMath::Sin(phi));
        xL0 = xTR2 + (zTR2-StaveZ[0])*(TMath::Tan(theta))*(TMath::Cos(phi));
        yL0 = yTR2 + (zTR2-StaveZ[0])*(TMath::Tan(theta))*(TMath::Sin(phi));

        //applica smiring
        smiring_xL2 = rnd->Uniform(-cls_size_x, +cls_size_x);
        smiring_xL1 = rnd->Uniform(-cls_size_x, +cls_size_x);
        smiring_xL0 = rnd->Uniform(-cls_size_x, +cls_size_x);
        smiring_yL2 = rnd->Uniform(-cls_size_y, +cls_size_y);
        smiring_yL1 = rnd->Uniform(-cls_size_y, +cls_size_y);
        smiring_yL0 = rnd->Uniform(-cls_size_y, +cls_size_y);
        xL2 = xL2 + smiring_xL2;
        xL1 = xL1 + smiring_xL1;
        xL0 = xL0 + smiring_xL0;
        yL2 = yL2 + smiring_yL2;
        yL1 = yL1 + smiring_yL1;
        yL0 = yL0 + smiring_yL0;

        xTR1 = xTR2 + (zTR2-TR1CenterZ)*(TMath::Tan(theta))*(TMath::Cos(phi));
        yTR1 = yTR2 + (zTR2-TR1CenterZ)*(TMath::Tan(theta))*(TMath::Sin(phi));

        

    } while (!(xTR1 < TR1Size[0]/2 && xTR1 > -TR1Size[0]/2 &&
        (
            (yTR1 < (2.5*TR1Size[1]+2*TR1GapY) && yTR1 > (1.5*TR1Size[1]+2*TR1GapY)) ||
            (yTR1 < (1.5*TR1Size[1]+1*TR1GapY) && yTR1 > (0.5*TR1Size[1]+1*TR1GapY)) ||
            (yTR1 < (0.5*TR1Size[1]+0*TR1GapY) && yTR1 > -(0.5*TR1Size[1]+0*TR1GapY)) ||
            (yTR1 < -(0.5*TR1Size[1]+1*TR1GapY) && yTR1 > -(1.5*TR1Size[1]+1*TR1GapY)) ||
            (yTR1 < -(1.5*TR1Size[1]+2*TR1GapY) && yTR1 > -(2.5*TR1Size[1]+2*TR1GapY)) 
        ) 
    ));
    
    


    if(xTR1 < TR1Size[0]/2 && xTR1 > -TR1Size[0]/2 &&
        (
            (yTR1 < (2.5*TR1Size[1]+2*TR1GapY) && yTR1 > (1.5*TR1Size[1]+2*TR1GapY)) ||
            (yTR1 < (1.5*TR1Size[1]+1*TR1GapY) && yTR1 > (0.5*TR1Size[1]+1*TR1GapY)) ||
            (yTR1 < (0.5*TR1Size[1]+0*TR1GapY) && yTR1 > -(0.5*TR1Size[1]+0*TR1GapY)) ||
            (yTR1 < -(0.5*TR1Size[1]+1*TR1GapY) && yTR1 > -(1.5*TR1Size[1]+1*TR1GapY)) ||
            (yTR1 < -(1.5*TR1Size[1]+2*TR1GapY) && yTR1 > -(2.5*TR1Size[1]+2*TR1GapY)) 
        )
      ){
        stats::hmgthTR1++;
    }
    
    


    //fake hit rate (rate = 10^-6 per event)
    double rate = rnd->Uniform(0,1);
    if (rate < 1./1000.){
        //x
        double fakehit_x = rnd->Uniform(-(2.5*display::ChipSizeX+2*display::ChipDistanceX),+(2.5*display::ChipSizeX+2*display::ChipDistanceX));
        //y
        double FAKEhit_y = rnd->Uniform(-(5*display::ChipSizeY+2.5*display::ChipDistanceY),+(5*display::ChipSizeY+2.5*display::ChipDistanceY));
        double fakehit_y;
        if(FAKEhit_y < (0.5*display::ChipDistanceY+1*display::ChipSizeY) && FAKEhit_y > -(0.5*display::ChipDistanceY+1*display::ChipSizeY)){
            fakehit_y = FAKEhit_y;
        }
        if((FAKEhit_y < (1.5*display::ChipDistanceY+3*display::ChipSizeY) && FAKEhit_y > (0.5*display::ChipDistanceY+1*display::ChipSizeY)) ||
          (FAKEhit_y > -(1.5*display::ChipDistanceY+3*display::ChipSizeY) && FAKEhit_y < -(0.5*display::ChipDistanceY+1*display::ChipSizeY))){
            fakehit_y = FAKEhit_y + display::ChipStaveDistanceY;
          }
        if((FAKEhit_y < (2.5*display::ChipDistanceY+5*display::ChipSizeY) && FAKEhit_y > (1.5*display::ChipDistanceY+3*display::ChipSizeY)) ||
          (FAKEhit_y > -(2.5*display::ChipDistanceY+5*display::ChipSizeY) && FAKEhit_y < -(1.5*display::ChipDistanceY+3*display::ChipSizeY))){
            fakehit_y = FAKEhit_y + 2*display::ChipStaveDistanceY;
          }
        //z
        int r = 1 + (int)(rnd->Uniform(0,3));
        double fakehit_z;
        if(r == 1){fakehit_z = display::StaveZ[0];}
        if(r == 2){fakehit_z = display::StaveZ[1];}
        if(r == 3){fakehit_z = display::StaveZ[2];}

        TMarker3DBox *fake_hit = new TMarker3DBox(fakehit_x, fakehit_y, fakehit_z, err_cl, err_cl, 0, 0, 0);
        fake_hit->Draw();

        stats::fakehit++;
    }


    //plotting tracks
    
    Double_t x_line[2] = {xTR2, xTR1};
    Double_t y_line[2] = {yTR2, yTR1};
    Double_t z_line[2] = {zTR2, zTR1};
    TPolyLine3D* line_track = new TPolyLine3D(2, x_line, y_line, z_line);
    line_track->SetLineColor(kBlue);
    line_track->SetLineWidth(4);
    geom->cd();
    line_track->Draw();
    
    
    /*
    cout << "real data" << endl;
    cout << "L2 ------ x:" << xL2 << ",   y: " << yL2 << ",    z: " << StaveZ[2] << endl; 
    cout << "L1 ------ x:" << xL1 << ",   y: " << yL1 << ",    z: " << StaveZ[1] << endl; 
    cout << "L0 ------ x:" << xL0 << ",   y: " << yL0 << ",    z: " << StaveZ[0] << endl; 
    cout << "(gra)theta" << (theta*180)/TMath::Pi() << ",     phi" << (phi*180)/TMath::Pi() << endl;
    cout << "(rad) theta: " << theta << "     phi: " << phi << endl;
    cout << "-----------------------------------------------------------------------------"; 
    */

    LCluster pL2; LCluster mL1; LCluster qL0;
    //check if the track hitted the staves in layer 2
    if((xL2 < ChipSizeX*2.5 + ChipDistanceX && xL2 > -(ChipSizeX*2.5 + ChipDistanceX)) &&
    ((yL2 < ChipSizeY*5 + ChipStaveDistanceY*2 + ChipDistanceY*2.5 && yL2 > ChipSizeY*4 + ChipStaveDistanceY*2 + ChipDistanceY*2.5) ||
    (yL2 < ChipSizeY*4 + ChipStaveDistanceY*2 + ChipDistanceY*1.5 && yL2 > ChipSizeY*3 + ChipStaveDistanceY*2 + ChipDistanceY*1.5) ||
    (yL2 < ChipSizeY*3 + ChipStaveDistanceY*1 + ChipDistanceY*1.5 && yL2 > ChipSizeY*2 + ChipStaveDistanceY*1 + ChipDistanceY*1.5) ||
    (yL2 < ChipSizeY*2 + ChipStaveDistanceY*1 + ChipDistanceY*0.5 && yL2 > ChipSizeY*1 + ChipStaveDistanceY*1 + ChipDistanceY*0.5) ||
    (yL2 < ChipSizeY*1 + ChipStaveDistanceY*0 + ChipDistanceY*0.5 && yL2 > ChipSizeY*0 + ChipStaveDistanceY*0 + ChipDistanceY*0.5) ||
    (yL2 > -(ChipSizeY*1 + ChipStaveDistanceY*0 + ChipDistanceY*0.5) && yL2 < -(ChipSizeY*0 + ChipStaveDistanceY*0 + ChipDistanceY*0.5)) ||
    (yL2 > -(ChipSizeY*2 + ChipStaveDistanceY*1 + ChipDistanceY*0.5) && yL2 < -(ChipSizeY*1 + ChipStaveDistanceY*1 + ChipDistanceY*0.5)) ||
    (yL2 > -(ChipSizeY*3 + ChipStaveDistanceY*1 + ChipDistanceY*1.5) && yL2 < -(ChipSizeY*2 + ChipStaveDistanceY*1 + ChipDistanceY*1.5)) ||
    (yL2 > -(ChipSizeY*4 + ChipStaveDistanceY*2 + ChipDistanceY*1.5) && yL2 < -(ChipSizeY*3 + ChipStaveDistanceY*2 + ChipDistanceY*1.5)) ||
    (yL2 > -(ChipSizeY*5 + ChipStaveDistanceY*2 + ChipDistanceY*2.5) && yL2 < -(ChipSizeY*4 + ChipStaveDistanceY*2 + ChipDistanceY*2.5)))
    ){
    stats::hmgthL2++;
    stats::hitL2 = true;
    //only if hitted the chips goes to tidy_clusters
    pL2.fill_cluster(pL2, xL2, yL2, StaveZ[2], cls_size_x, cls_size_y, 0, i);
    tracker.tidy_clusters_lay2.try_emplace(i,pL2);
    TMarker3DBox *p = new TMarker3DBox(xL2, yL2, StaveZ[2], err_cl, err_cl,0,0,0);
    p->Draw();
    }
    //check if the track hitted the staves in layer 1
    if((xL1 < ChipSizeX*2.5 + ChipDistanceX && xL1 > -(ChipSizeX*2.5 + ChipDistanceX)) &&
    ((yL1 < ChipSizeY*5 + ChipStaveDistanceY*2 + ChipDistanceY*2.5 && yL1 > ChipSizeY*4 + ChipStaveDistanceY*2 + ChipDistanceY*2.5) ||
    (yL1 < ChipSizeY*4 + ChipStaveDistanceY*2 + ChipDistanceY*1.5 && yL1 > ChipSizeY*3 + ChipStaveDistanceY*2 + ChipDistanceY*1.5) ||
    (yL1 < ChipSizeY*3 + ChipStaveDistanceY*1 + ChipDistanceY*1.5 && yL1 > ChipSizeY*2 + ChipStaveDistanceY*1 + ChipDistanceY*1.5) ||
    (yL1 < ChipSizeY*2 + ChipStaveDistanceY*1 + ChipDistanceY*0.5 && yL1 > ChipSizeY*1 + ChipStaveDistanceY*1 + ChipDistanceY*0.5) ||
    (yL1 < ChipSizeY*1 + ChipStaveDistanceY*0 + ChipDistanceY*0.5 && yL1 > ChipSizeY*0 + ChipStaveDistanceY*0 + ChipDistanceY*0.5) ||
    (yL1 > -(ChipSizeY*1 + ChipStaveDistanceY*0 + ChipDistanceY*0.5) && yL1 < -(ChipSizeY*0 + ChipStaveDistanceY*0 + ChipDistanceY*0.5)) ||
    (yL1 > -(ChipSizeY*2 + ChipStaveDistanceY*1 + ChipDistanceY*0.5) && yL1 < -(ChipSizeY*1 + ChipStaveDistanceY*1 + ChipDistanceY*0.5)) ||
    (yL1 > -(ChipSizeY*3 + ChipStaveDistanceY*1 + ChipDistanceY*1.5) && yL1 < -(ChipSizeY*2 + ChipStaveDistanceY*1 + ChipDistanceY*1.5)) ||
    (yL1 > -(ChipSizeY*4 + ChipStaveDistanceY*2 + ChipDistanceY*1.5) && yL1 < -(ChipSizeY*3 + ChipStaveDistanceY*2 + ChipDistanceY*1.5)) ||
    (yL1 > -(ChipSizeY*5 + ChipStaveDistanceY*2 + ChipDistanceY*2.5) && yL1 < -(ChipSizeY*4 + ChipStaveDistanceY*2 + ChipDistanceY*2.5)))
    ){
    stats::hmgthL1++;
    stats::hitL1 = true;
    mL1.fill_cluster(mL1, xL1, yL1, StaveZ[1], cls_size_x, cls_size_y, 0, i);
    tracker.tidy_clusters_lay1.try_emplace(i,mL1);
    TMarker3DBox *m = new TMarker3DBox(xL1, yL1, StaveZ[1], err_cl, err_cl,0,0,0);
    m->Draw();
    }
    //check if the track hitted the staves in layer 0
    if((xL0 < ChipSizeX*2.5 + ChipDistanceX && xL0 > -(ChipSizeX*2.5 + ChipDistanceX)) &&
    ((yL0 < ChipSizeY*5 + ChipStaveDistanceY*2 + ChipDistanceY*2.5 && yL0 > ChipSizeY*4 + ChipStaveDistanceY*2 + ChipDistanceY*2.5) ||
    (yL0 < ChipSizeY*4 + ChipStaveDistanceY*2 + ChipDistanceY*1.5 && yL0 > ChipSizeY*3 + ChipStaveDistanceY*2 + ChipDistanceY*1.5) ||
    (yL0 < ChipSizeY*3 + ChipStaveDistanceY*1 + ChipDistanceY*1.5 && yL0 > ChipSizeY*2 + ChipStaveDistanceY*1 + ChipDistanceY*1.5) ||
    (yL0 < ChipSizeY*2 + ChipStaveDistanceY*1 + ChipDistanceY*0.5 && yL0 > ChipSizeY*1 + ChipStaveDistanceY*1 + ChipDistanceY*0.5) ||
    (yL0 < ChipSizeY*1 + ChipStaveDistanceY*0 + ChipDistanceY*0.5 && yL0 > ChipSizeY*0 + ChipStaveDistanceY*0 + ChipDistanceY*0.5) ||
    (yL0 > -(ChipSizeY*1 + ChipStaveDistanceY*0 + ChipDistanceY*0.5) && yL0 < -(ChipSizeY*0 + ChipStaveDistanceY*0 + ChipDistanceY*0.5)) ||
    (yL0 > -(ChipSizeY*2 + ChipStaveDistanceY*1 + ChipDistanceY*0.5) && yL0 < -(ChipSizeY*1 + ChipStaveDistanceY*1 + ChipDistanceY*0.5)) ||
    (yL0 > -(ChipSizeY*3 + ChipStaveDistanceY*1 + ChipDistanceY*1.5) && yL0 < -(ChipSizeY*2 + ChipStaveDistanceY*1 + ChipDistanceY*1.5)) ||
    (yL0 > -(ChipSizeY*4 + ChipStaveDistanceY*2 + ChipDistanceY*1.5) && yL0 < -(ChipSizeY*3 + ChipStaveDistanceY*2 + ChipDistanceY*1.5)) ||
    (yL0 > -(ChipSizeY*5 + ChipStaveDistanceY*2 + ChipDistanceY*2.5) && yL0 < -(ChipSizeY*4 + ChipStaveDistanceY*2 + ChipDistanceY*2.5)))
    ){
    stats::hmgthL0++;
    stats::hitL0 = true;
    qL0.fill_cluster(qL0, xL0, yL0, StaveZ[0], cls_size_x, cls_size_y, 0, i);
    tracker.tidy_clusters_lay0.try_emplace(i,qL0);
    TMarker3DBox *q = new TMarker3DBox(xL0, yL0, StaveZ[0], err_cl, err_cl, 0,0,0);
    q->Draw();
    }

    chip c;
    c.is_dead_chip_tracking_smt(c, pL2, mL1, qL0);

    
    //contiamo quante traccie hanno colpito quanti layer ciascuna (3layer, 2layer, 1layer, 0layer)
    if(stats::hitL2 && stats::hitL1 && stats::hitL0 && stats::hmgthTR1){
        stats::hmgthL012++;
        htheta3L->Fill(theta);
        hphi3L->Fill(phi);
    }
    if((stats::hitL0 && stats::hitL1 && !stats::hitL2)||(stats::hitL1 && stats::hitL2 && !stats::hitL0)||(stats::hitL0 && stats::hitL2 && !stats::hitL1)){
        stats::hmgth2L++;
        htheta2L->Fill(theta);
        hphi2L->Fill(phi);
    }
    if((stats::hitL0 && !stats::hitL1 && !stats::hitL2)||(stats::hitL1 && !stats::hitL2 && !stats::hitL0)||(!stats::hitL0 && stats::hitL2 && !stats::hitL1)){
        stats::hmgth1L++;
        htheta1L->Fill(theta);
        hphi1L->Fill(phi);
    }
    if(!stats::hitL0 && !stats::hitL1 && !stats::hitL2){
        stats::hmgth0L++;
        htheta0L->Fill(theta);
        hphi0L->Fill(phi);
    }
    //SOLO SE VUOI FARE TRACKING SULLE TRACCIE CHE HANNO COLPITO I 3 LAYER
    //creo solo le tracklet di traccie con hit su L012
/*     if(stats::hitL0 && stats::hitL1 && stats::hitL2){
        tracker.tidy_clusters_lay2.try_emplace(i,pL2);
        tracker.tidy_clusters_lay1.try_emplace(i,mL1);
        tracker.tidy_clusters_lay0.try_emplace(i,qL0);
    }  */
    
    
    
    
    
    LTrackCandidate real_track;
    real_track.x0 = xL1;
    real_track.y0 = yL1;
    real_track.z0 = StaveZ[1];
    real_track.err_x0 = cls_size_x;
    real_track.err_y0 = cls_size_y;
    real_track.theta = theta;
    real_track.phi = phi;
    generated_tracks.push_back(real_track);

    hxTR2->Fill(xTR2);
    hyTR2->Fill(yTR2);
    hzTR2->Fill(zTR2);
    hxTR1->Fill(xTR1);
    hyTR1->Fill(yTR1);
    hzTR1->Fill(zTR1);
    hphi->Fill(phi);
    htheta->Fill(theta);

    if(stats::hitL2 || stats::hitL1 || stats::hitL0){
        hx->Fill(xL2);
        hx->Fill(xL1);
        hx->Fill(xL0);
        hy->Fill(yL2);
        hy->Fill(yL1);
        hy->Fill(yL0);
    }


}

char file[200];
//output mc distributions
//geom->SaveAs("../data/plot.png");
sprintf(file,"../data/statsHisttracks.root");
TFile f(file,"RECREATE");
hxTR2->Write();
hyTR2->Write();
hzTR2->Write();
hxTR1->Write();
hyTR1->Write();
hzTR1->Write();
htheta->Write();
hphi->Write();
hphi3L->Write();
hphi2L->Write();
hphi1L->Write();
hphi0L->Write();
htheta3L->Write();
htheta2L->Write();
htheta1L->Write();
htheta0L->Write();
hx->Write();
hy->Write();
f.Close();

}



void display::tracks_no_print_hist(int events, LTrackerTrack& tracker){

stats::hmgt = events;

//puoi settare il seed for reproducibility
TRandom3 *rnd = new TRandom3(0); 

//MC
for (int i=0; i < events; i++){
    double xTR2, yTR2, zTR2, xTR1, yTR1, zTR1;
    double xL2, xL1, xL0, yL2, yL1, yL0;
    float phi, theta, cls_size, cls_size_x, cls_size_y; 
    double smiring_xL2, smiring_xL1, smiring_xL0, smiring_yL2, smiring_yL1, smiring_yL0;
    stats::hitL0=false;
    stats::hitL1=false;
    stats::hitL2=false;

    //i layer acquisiscono il segnale quando arriva un segnale AND dagli scintillatori
    //genero traccie misurabili come traccie che passano nei due scintillatori TR2, TR1
    
    do {
        double xTR2_fake = rnd->Uniform(-TR2Size[0]*2,TR2Size[0]*2);     
        if(xTR2_fake>0 && xTR2_fake<TR2Size[0]){xTR2 = xTR2_fake + 0.5*TR2GapX;}
        if(xTR2_fake<0 && xTR2_fake>-TR2Size[0]){xTR2 = xTR2_fake - 0.5*TR2GapX;}
        if(xTR2_fake<2*TR2Size[0] && xTR2_fake>TR2Size[0]){xTR2 = xTR2_fake + 1.5*TR2GapX;}
        if(xTR2_fake>-2*TR2Size[0] && xTR2_fake<-TR2Size[0]){xTR2 = xTR2_fake - 1.5*TR2GapX;}
        yTR2 = rnd->Uniform(-TR2Size[1]/2,TR2Size[1]/2);
        zTR2 = rnd->Uniform(TR2CenterZ-TR2Thickness/2,TR2CenterZ+TR2Thickness/2);
        zTR1 = rnd->Uniform(TR1CenterZ-TR1Thickness/2,TR1CenterZ+TR1Thickness/2);

        //prendo gli angoli dalla distribuzione di muoni
        int ind_theta = rnd->Uniform(0, allTheta.size()-1);
        int ind_phi = rnd->Uniform(0, allPhi.size()-1);
        int ind_cls_size = rnd->Uniform(0, all_cls_size.size()-1);
        theta = (allTheta[ind_theta]/180)*TMath::Pi();
        phi = (allPhi[ind_phi]/180)*TMath::Pi();
        cls_size = (all_cls_size[ind_cls_size]);
        cls_size_x = cls_size*PixelSizeRows;
        cls_size_y = cls_size*PixelSizeCols;

        //cout << "cls size xy " << cls_size_x << "   " << cls_size_y << endl; 


        xL2 = xTR2 + (zTR2-StaveZ[2])*(TMath::Tan(theta))*(TMath::Cos(phi));
        yL2 = yTR2 + (zTR2-StaveZ[2])*(TMath::Tan(theta))*(TMath::Sin(phi));
        xL1 = xTR2 + (zTR2-StaveZ[1])*(TMath::Tan(theta))*(TMath::Cos(phi));
        yL1 = yTR2 + (zTR2-StaveZ[1])*(TMath::Tan(theta))*(TMath::Sin(phi));
        xL0 = xTR2 + (zTR2-StaveZ[0])*(TMath::Tan(theta))*(TMath::Cos(phi));
        yL0 = yTR2 + (zTR2-StaveZ[0])*(TMath::Tan(theta))*(TMath::Sin(phi));

        //applica smiring
        smiring_xL2 = rnd->Uniform(-cls_size_x, +cls_size_x);
        smiring_xL1 = rnd->Uniform(-cls_size_x, +cls_size_x);
        smiring_xL0 = rnd->Uniform(-cls_size_x, +cls_size_x);
        smiring_yL2 = rnd->Uniform(-cls_size_y, +cls_size_y);
        smiring_yL1 = rnd->Uniform(-cls_size_y, +cls_size_y);
        smiring_yL0 = rnd->Uniform(-cls_size_y, +cls_size_y);
        xL2 = xL2 + smiring_xL2;
        xL1 = xL1 + smiring_xL1;
        xL0 = xL0 + smiring_xL0;
        yL2 = yL2 + smiring_yL2;
        yL1 = yL1 + smiring_yL1;
        yL0 = yL0 + smiring_yL0;


        xTR1 = xTR2 + (zTR2-TR1CenterZ)*(TMath::Tan(theta))*(TMath::Cos(phi));
        yTR1 = yTR2 + (zTR2-TR1CenterZ)*(TMath::Tan(theta))*(TMath::Sin(phi));
    } while (!(xTR1 < TR1Size[0]/2 && xTR1 > -TR1Size[0]/2 &&
        (
            (yTR1 < (2.5*TR1Size[1]+2*TR1GapY) && yTR1 > (1.5*TR1Size[1]+2*TR1GapY)) ||
            (yTR1 < (1.5*TR1Size[1]+1*TR1GapY) && yTR1 > (0.5*TR1Size[1]+1*TR1GapY)) ||
            (yTR1 < (0.5*TR1Size[1]+0*TR1GapY) && yTR1 > -(0.5*TR1Size[1]+0*TR1GapY)) ||
            (yTR1 < -(0.5*TR1Size[1]+1*TR1GapY) && yTR1 > -(1.5*TR1Size[1]+1*TR1GapY)) ||
            (yTR1 < -(1.5*TR1Size[1]+2*TR1GapY) && yTR1 > -(2.5*TR1Size[1]+2*TR1GapY)) 
        )
    ));
    

    if(xTR1 < TR1Size[0]/2 && xTR1 > -TR1Size[0]/2 &&
        (
            (yTR1 < (2.5*TR1Size[1]+2*TR1GapY) && yTR1 > (1.5*TR1Size[1]+2*TR1GapY)) ||
            (yTR1 < (1.5*TR1Size[1]+1*TR1GapY) && yTR1 > (0.5*TR1Size[1]+1*TR1GapY)) ||
            (yTR1 < (0.5*TR1Size[1]+0*TR1GapY) && yTR1 > -(0.5*TR1Size[1]+0*TR1GapY)) ||
            (yTR1 < -(0.5*TR1Size[1]+1*TR1GapY) && yTR1 > -(1.5*TR1Size[1]+1*TR1GapY)) ||
            (yTR1 < -(1.5*TR1Size[1]+2*TR1GapY) && yTR1 > -(2.5*TR1Size[1]+2*TR1GapY)) 
        )
      ){
        stats::hmgthTR1++;
    }
    
    
    //fake hit rate (rate = 10^-6 per event)
    /* double rate = rnd->Uniform(0,1);
    LCluster fakehit;
    if (rate < 1./1000.){
        //x
        double fakehit_x = rnd->Uniform(-(2.5*display::ChipSizeX+2*display::ChipDistanceX),+(2.5*display::ChipSizeX+2*display::ChipDistanceX));
        //y
        double FAKEhit_y = rnd->Uniform(-(5*display::ChipSizeY+2.5*display::ChipDistanceY),+(5*display::ChipSizeY+2.5*display::ChipDistanceY));
        double fakehit_y;
        if(FAKEhit_y < (0.5*display::ChipDistanceY+1*display::ChipSizeY) && FAKEhit_y > -(0.5*display::ChipDistanceY+1*display::ChipSizeY)){
            fakehit_y = FAKEhit_y;
        }
        if((FAKEhit_y < (1.5*display::ChipDistanceY+3*display::ChipSizeY) && FAKEhit_y > (0.5*display::ChipDistanceY+1*display::ChipSizeY)) ||
          (FAKEhit_y > -(1.5*display::ChipDistanceY+3*display::ChipSizeY) && FAKEhit_y < -(0.5*display::ChipDistanceY+1*display::ChipSizeY))){
            fakehit_y = FAKEhit_y + display::ChipStaveDistanceY;
          }
        if((FAKEhit_y < (2.5*display::ChipDistanceY+5*display::ChipSizeY) && FAKEhit_y > (1.5*display::ChipDistanceY+3*display::ChipSizeY)) ||
          (FAKEhit_y > -(2.5*display::ChipDistanceY+5*display::ChipSizeY) && FAKEhit_y < -(1.5*display::ChipDistanceY+3*display::ChipSizeY))){
            fakehit_y = FAKEhit_y + 2*display::ChipStaveDistanceY;
          }
        //z
        int r = 1 + (int)(rnd->Uniform(0,3));
        double fakehit_z;
        if(r == 1){
            fakehit_z = display::StaveZ[0];
            fakehit.fill_cluster(fakehit, fakehit_x, fakehit_y, fakehit_z, cls_size_x, cls_size_y, 0, -1);
            tracker.tidy_clusters_lay0.try_emplace(-1, fakehit);
        }
        if(r == 2){
            fakehit_z = display::StaveZ[1];
            fakehit.fill_cluster(fakehit, fakehit_x, fakehit_y, fakehit_z, cls_size_x, cls_size_y, 0, -1);
            tracker.tidy_clusters_lay1.try_emplace(-1, fakehit);
        }
        if(r == 3){
            fakehit_z = display::StaveZ[2];
            fakehit.fill_cluster(fakehit, fakehit_x, fakehit_y, fakehit_z, cls_size_x, cls_size_y, 0, -1);
            tracker.tidy_clusters_lay2.try_emplace(-1, fakehit);
        }

        stats::fakehit++;
    }
 */

    LCluster pL2; LCluster mL1; LCluster qL0;
    //check if the track hitted the staves in layer 2
    if((xL2 < ChipSizeX*2.5 + ChipDistanceX && xL2 > -(ChipSizeX*2.5 + ChipDistanceX)) &&
    ((yL2 < ChipSizeY*5 + ChipStaveDistanceY*2 + ChipDistanceY*2.5 && yL2 > ChipSizeY*4 + ChipStaveDistanceY*2 + ChipDistanceY*2.5) ||
    (yL2 < ChipSizeY*4 + ChipStaveDistanceY*2 + ChipDistanceY*1.5 && yL2 > ChipSizeY*3 + ChipStaveDistanceY*2 + ChipDistanceY*1.5) ||
    (yL2 < ChipSizeY*3 + ChipStaveDistanceY*1 + ChipDistanceY*1.5 && yL2 > ChipSizeY*2 + ChipStaveDistanceY*1 + ChipDistanceY*1.5) ||
    (yL2 < ChipSizeY*2 + ChipStaveDistanceY*1 + ChipDistanceY*0.5 && yL2 > ChipSizeY*1 + ChipStaveDistanceY*1 + ChipDistanceY*0.5) ||
    (yL2 < ChipSizeY*1 + ChipStaveDistanceY*0 + ChipDistanceY*0.5 && yL2 > ChipSizeY*0 + ChipStaveDistanceY*0 + ChipDistanceY*0.5) ||
    (yL2 > -(ChipSizeY*1 + ChipStaveDistanceY*0 + ChipDistanceY*0.5) && yL2 < -(ChipSizeY*0 + ChipStaveDistanceY*0 + ChipDistanceY*0.5)) ||
    (yL2 > -(ChipSizeY*2 + ChipStaveDistanceY*1 + ChipDistanceY*0.5) && yL2 < -(ChipSizeY*1 + ChipStaveDistanceY*1 + ChipDistanceY*0.5)) ||
    (yL2 > -(ChipSizeY*3 + ChipStaveDistanceY*1 + ChipDistanceY*1.5) && yL2 < -(ChipSizeY*2 + ChipStaveDistanceY*1 + ChipDistanceY*1.5)) ||
    (yL2 > -(ChipSizeY*4 + ChipStaveDistanceY*2 + ChipDistanceY*1.5) && yL2 < -(ChipSizeY*3 + ChipStaveDistanceY*2 + ChipDistanceY*1.5)) ||
    (yL2 > -(ChipSizeY*5 + ChipStaveDistanceY*2 + ChipDistanceY*2.5) && yL2 < -(ChipSizeY*4 + ChipStaveDistanceY*2 + ChipDistanceY*2.5)))
    ){
    stats::hmgthL2++;
    stats::hitL2 = true;
    //only if hitted the chips goes to tidy_clusters
    pL2.fill_cluster(pL2, xL2, yL2, StaveZ[2], cls_size_x, cls_size_y, 0, i);
    tracker.tidy_clusters_lay2.try_emplace(i,pL2);
    }
    //check if the track hitted the staves in layer 1
    if((xL1 < ChipSizeX*2.5 + ChipDistanceX && xL1 > -(ChipSizeX*2.5 + ChipDistanceX)) &&
    ((yL1 < ChipSizeY*5 + ChipStaveDistanceY*2 + ChipDistanceY*2.5 && yL1 > ChipSizeY*4 + ChipStaveDistanceY*2 + ChipDistanceY*2.5) ||
    (yL1 < ChipSizeY*4 + ChipStaveDistanceY*2 + ChipDistanceY*1.5 && yL1 > ChipSizeY*3 + ChipStaveDistanceY*2 + ChipDistanceY*1.5) ||
    (yL1 < ChipSizeY*3 + ChipStaveDistanceY*1 + ChipDistanceY*1.5 && yL1 > ChipSizeY*2 + ChipStaveDistanceY*1 + ChipDistanceY*1.5) ||
    (yL1 < ChipSizeY*2 + ChipStaveDistanceY*1 + ChipDistanceY*0.5 && yL1 > ChipSizeY*1 + ChipStaveDistanceY*1 + ChipDistanceY*0.5) ||
    (yL1 < ChipSizeY*1 + ChipStaveDistanceY*0 + ChipDistanceY*0.5 && yL1 > ChipSizeY*0 + ChipStaveDistanceY*0 + ChipDistanceY*0.5) ||
    (yL1 > -(ChipSizeY*1 + ChipStaveDistanceY*0 + ChipDistanceY*0.5) && yL1 < -(ChipSizeY*0 + ChipStaveDistanceY*0 + ChipDistanceY*0.5)) ||
    (yL1 > -(ChipSizeY*2 + ChipStaveDistanceY*1 + ChipDistanceY*0.5) && yL1 < -(ChipSizeY*1 + ChipStaveDistanceY*1 + ChipDistanceY*0.5)) ||
    (yL1 > -(ChipSizeY*3 + ChipStaveDistanceY*1 + ChipDistanceY*1.5) && yL1 < -(ChipSizeY*2 + ChipStaveDistanceY*1 + ChipDistanceY*1.5)) ||
    (yL1 > -(ChipSizeY*4 + ChipStaveDistanceY*2 + ChipDistanceY*1.5) && yL1 < -(ChipSizeY*3 + ChipStaveDistanceY*2 + ChipDistanceY*1.5)) ||
    (yL1 > -(ChipSizeY*5 + ChipStaveDistanceY*2 + ChipDistanceY*2.5) && yL1 < -(ChipSizeY*4 + ChipStaveDistanceY*2 + ChipDistanceY*2.5)))
    ){
    stats::hmgthL1++;
    stats::hitL1 = true;
    mL1.fill_cluster(mL1, xL1, yL1, StaveZ[1], cls_size_x, cls_size_y, 0, i);
    tracker.tidy_clusters_lay1.try_emplace(i,mL1);
    }
    //check if the track hitted the staves in layer 0
    if((xL0 < ChipSizeX*2.5 + ChipDistanceX && xL0 > -(ChipSizeX*2.5 + ChipDistanceX)) &&
    ((yL0 < ChipSizeY*5 + ChipStaveDistanceY*2 + ChipDistanceY*2.5 && yL0 > ChipSizeY*4 + ChipStaveDistanceY*2 + ChipDistanceY*2.5) ||
    (yL0 < ChipSizeY*4 + ChipStaveDistanceY*2 + ChipDistanceY*1.5 && yL0 > ChipSizeY*3 + ChipStaveDistanceY*2 + ChipDistanceY*1.5) ||
    (yL0 < ChipSizeY*3 + ChipStaveDistanceY*1 + ChipDistanceY*1.5 && yL0 > ChipSizeY*2 + ChipStaveDistanceY*1 + ChipDistanceY*1.5) ||
    (yL0 < ChipSizeY*2 + ChipStaveDistanceY*1 + ChipDistanceY*0.5 && yL0 > ChipSizeY*1 + ChipStaveDistanceY*1 + ChipDistanceY*0.5) ||
    (yL0 < ChipSizeY*1 + ChipStaveDistanceY*0 + ChipDistanceY*0.5 && yL0 > ChipSizeY*0 + ChipStaveDistanceY*0 + ChipDistanceY*0.5) ||
    (yL0 > -(ChipSizeY*1 + ChipStaveDistanceY*0 + ChipDistanceY*0.5) && yL0 < -(ChipSizeY*0 + ChipStaveDistanceY*0 + ChipDistanceY*0.5)) ||
    (yL0 > -(ChipSizeY*2 + ChipStaveDistanceY*1 + ChipDistanceY*0.5) && yL0 < -(ChipSizeY*1 + ChipStaveDistanceY*1 + ChipDistanceY*0.5)) ||
    (yL0 > -(ChipSizeY*3 + ChipStaveDistanceY*1 + ChipDistanceY*1.5) && yL0 < -(ChipSizeY*2 + ChipStaveDistanceY*1 + ChipDistanceY*1.5)) ||
    (yL0 > -(ChipSizeY*4 + ChipStaveDistanceY*2 + ChipDistanceY*1.5) && yL0 < -(ChipSizeY*3 + ChipStaveDistanceY*2 + ChipDistanceY*1.5)) ||
    (yL0 > -(ChipSizeY*5 + ChipStaveDistanceY*2 + ChipDistanceY*2.5) && yL0 < -(ChipSizeY*4 + ChipStaveDistanceY*2 + ChipDistanceY*2.5)))
    ){
    stats::hmgthL0++;
    stats::hitL0 = true;
    qL0.fill_cluster(qL0, xL0, yL0, StaveZ[0], cls_size_x, cls_size_y, 0, i);
    tracker.tidy_clusters_lay0.try_emplace(i,qL0);
    }

    //chip c;
    //c.is_dead_chip_tracking_smt(c, pL2, mL1, qL0);

    
    //contiamo quante traccie hanno colpito quanti layer ciascuna (3layer, 2layer, 1layer, 0layer)
    if(stats::hitL2 && stats::hitL1 && stats::hitL0 && stats::hmgthTR1){
        stats::hmgthL012++;
    }
    if((stats::hitL0 && stats::hitL1 && !stats::hitL2)||(stats::hitL1 && stats::hitL2 && !stats::hitL0)||(stats::hitL0 && stats::hitL2 && !stats::hitL1)){
        stats::hmgth2L++;
    }
    if((stats::hitL0 && !stats::hitL1 && !stats::hitL2)||(stats::hitL1 && !stats::hitL2 && !stats::hitL0)||(!stats::hitL0 && stats::hitL2 && !stats::hitL1)){
        stats::hmgth1L++;
    }
    if(!stats::hitL0 && !stats::hitL1 && !stats::hitL2){
        stats::hmgth0L++;
    }
    //SOLO SE VUOI FARE TRACKING SULLE TRACCIE CHE HANNO COLPITO I 3 LAYER
    //creo solo le tracklet di traccie con hit su L012
    /* if(stats::hitL0 && stats::hitL1 && stats::hitL2){
        tracker.tidy_clusters_lay2.try_emplace(i,pL2);
        tracker.tidy_clusters_lay1.try_emplace(i,mL1);
        tracker.tidy_clusters_lay0.try_emplace(i,qL0);
        //cout << "STAI FACENDO RECO SOLO CON TRACCIE CHE HANNO COLPITO TUTTI E TRE I LAYER" << endl;
    } */
    
    //CONSIDERA TUTTI I CLUSTER
    
    
    
    


    LTrackCandidate real_track;
    real_track.x0 = xL1;
    real_track.y0 = yL1;
    real_track.z0 = StaveZ[1];
    real_track.err_x0 = cls_size_x;
    real_track.err_y0 = cls_size_y;
    real_track.theta = theta;
    real_track.phi = phi;
    generated_tracks.push_back(real_track);
    
}




}