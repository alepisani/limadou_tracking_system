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
#include "../build/display.h"
#include "../build/stats.h"
#include "../build/LTrackerTrack.h"
using namespace std;

display::display() {
    //geom = new TCanvas("geom", "3D geom", 800, 600);
    //geometry = TView::CreateView(1);
    //geometry->SetRange(-100, -100, 0, 100, 100, 70);
    //geometry->ShowAxis();
}


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

void display::layers(TCanvas* geom){

    geom->cd();
    //draw layer 1
    //row, cols (00) 
    //row0
    TMarker3DBox *stave1_00 = new TMarker3DBox(-2*(ChipSizeX+ChipDistanceX),-(2*ChipStaveDistanceY+2.5*ChipDistanceY+4.5*ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_00->Draw();
    TMarker3DBox *stave1_01 = new TMarker3DBox(-1*(ChipSizeX+ChipDistanceX),-(2*ChipStaveDistanceY+2.5*ChipDistanceY+4.5*ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_01->Draw();
    TMarker3DBox *stave1_02 = new TMarker3DBox(0,-(2*ChipStaveDistanceY+2.5*ChipDistanceY+4.5*ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_02->Draw();
    TMarker3DBox *stave1_03 = new TMarker3DBox(1*(ChipSizeX+ChipDistanceX),-(2*ChipStaveDistanceY+2.5*ChipDistanceY+4.5*ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_03->Draw();
    TMarker3DBox *stave1_04 = new TMarker3DBox(2*(ChipSizeX+ChipDistanceX),-(2*ChipStaveDistanceY+2.5*ChipDistanceY+4.5*ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_04->Draw();
    //row1
    TMarker3DBox *stave1_10 = new TMarker3DBox(-2*(ChipSizeX+ChipDistanceX),-(2*ChipStaveDistanceY+1.5*ChipDistanceY+3.5*ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_10->Draw();
    TMarker3DBox *stave1_11 = new TMarker3DBox(-(ChipSizeX+ChipDistanceX),-(2*ChipStaveDistanceY+1.5*ChipDistanceY+3.5*ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_11->Draw();
    TMarker3DBox *stave1_12 = new TMarker3DBox(0,-(2*ChipStaveDistanceY+1.5*ChipDistanceY+3.5*ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_12->Draw();
    TMarker3DBox *stave1_13 = new TMarker3DBox((ChipSizeX+ChipDistanceX),-(2*ChipStaveDistanceY+1.5*ChipDistanceY+3.5*ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_13->Draw();
    TMarker3DBox *stave1_14 = new TMarker3DBox(2*(ChipSizeX+ChipDistanceX),-(2*ChipStaveDistanceY+1.5*ChipDistanceY+3.5*ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_14->Draw();
    //row2
    TMarker3DBox *stave1_20 = new TMarker3DBox(-2*(ChipSizeX+ChipDistanceX),-(ChipStaveDistanceY+1.5*ChipDistanceY+2.5*ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_20->Draw();
    TMarker3DBox *stave1_21 = new TMarker3DBox(-(ChipSizeX+ChipDistanceX),-(ChipStaveDistanceY+1.5*ChipDistanceY+2.5*ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_21->Draw();
    TMarker3DBox *stave1_22 = new TMarker3DBox(0,-(ChipStaveDistanceY+1.5*ChipDistanceY+2.5*ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_22->Draw();
    TMarker3DBox *stave1_23 = new TMarker3DBox((ChipSizeX+ChipDistanceX),-(ChipStaveDistanceY+1.5*ChipDistanceY+2.5*ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_23->Draw();
    TMarker3DBox *stave1_24 = new TMarker3DBox(2*(ChipSizeX+ChipDistanceX),-(ChipStaveDistanceY+1.5*ChipDistanceY+2.5*ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_24->Draw();
    //row3
    TMarker3DBox *stave1_30 = new TMarker3DBox(-2*(ChipSizeX+ChipDistanceX),-(ChipStaveDistanceY+1.5*ChipDistanceY+1.5*ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_30->Draw();
    TMarker3DBox *stave1_31 = new TMarker3DBox(-(ChipSizeX+ChipDistanceX),-(ChipStaveDistanceY+1.5*ChipDistanceY+1.5*ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_31->Draw();
    TMarker3DBox *stave1_32 = new TMarker3DBox(0,-(ChipStaveDistanceY+1.5*ChipDistanceY+1.5*ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_32->Draw();
    TMarker3DBox *stave1_33 = new TMarker3DBox((ChipSizeX+ChipDistanceX),-(ChipStaveDistanceY+1.5*ChipDistanceY+1.5*ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_33->Draw();
    TMarker3DBox *stave1_34 = new TMarker3DBox(2*(ChipSizeX+ChipDistanceX),-(ChipStaveDistanceY+1.5*ChipDistanceY+1.5*ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_34->Draw();
    //row4
    TMarker3DBox *stave1_40 = new TMarker3DBox(-2*(ChipSizeX+ChipDistanceX),-0.5*(ChipDistanceY+ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_40->Draw();
    TMarker3DBox *stave1_41 = new TMarker3DBox(-(ChipSizeX+ChipDistanceX),-0.5*(ChipDistanceY+ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_41->Draw();
    TMarker3DBox *stave1_42 = new TMarker3DBox(0,-0.5*(ChipDistanceY+ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_42->Draw();
    TMarker3DBox *stave1_43 = new TMarker3DBox((ChipSizeX+ChipDistanceX),-0.5*(ChipDistanceY+ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_43->Draw();
    TMarker3DBox *stave1_44 = new TMarker3DBox(2*(ChipSizeX+ChipDistanceX),-0.5*(ChipDistanceY+ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_44->Draw();
    //row5
    TMarker3DBox *stave1_50 = new TMarker3DBox(-2*(ChipSizeX+ChipDistanceX),+0.5*(ChipDistanceY+ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_50->Draw();
    TMarker3DBox *stave1_51 = new TMarker3DBox(-(ChipSizeX+ChipDistanceX),+0.5*(ChipDistanceY+ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_51->Draw();
    TMarker3DBox *stave1_52 = new TMarker3DBox(0,+0.5*(ChipDistanceY+ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_52->Draw();
    TMarker3DBox *stave1_53 = new TMarker3DBox((ChipSizeX+ChipDistanceX),+0.5*(ChipDistanceY+ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_53->Draw();
    TMarker3DBox *stave1_54 = new TMarker3DBox(2*(ChipSizeX+ChipDistanceX),+0.5*(ChipDistanceY+ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_54->Draw();
    //row6
    TMarker3DBox *stave1_60 = new TMarker3DBox(-2*(ChipSizeX+ChipDistanceX),+(ChipStaveDistanceY+1.5*ChipDistanceY+1.5*ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_60->Draw();
    TMarker3DBox *stave1_61 = new TMarker3DBox(-(ChipSizeX+ChipDistanceX),+(ChipStaveDistanceY+1.5*ChipDistanceY+1.5*ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_61->Draw();
    TMarker3DBox *stave1_62 = new TMarker3DBox(0,+(ChipStaveDistanceY+1.5*ChipDistanceY+1.5*ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_62->Draw();
    TMarker3DBox *stave1_63 = new TMarker3DBox((ChipSizeX+ChipDistanceX),+(ChipStaveDistanceY+1.5*ChipDistanceY+1.5*ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_63->Draw();
    TMarker3DBox *stave1_64 = new TMarker3DBox(2*(ChipSizeX+ChipDistanceX),+(ChipStaveDistanceY+1.5*ChipDistanceY+1.5*ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_64->Draw();
    //row7
    TMarker3DBox *stave1_70 = new TMarker3DBox(-2*(ChipSizeX+ChipDistanceX),+(ChipStaveDistanceY+1.5*ChipDistanceY+2.5*ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_70->Draw();
    TMarker3DBox *stave1_71 = new TMarker3DBox(-(ChipSizeX+ChipDistanceX),+(ChipStaveDistanceY+1.5*ChipDistanceY+2.5*ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_71->Draw();
    TMarker3DBox *stave1_72 = new TMarker3DBox(0,+(ChipStaveDistanceY+1.5*ChipDistanceY+2.5*ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_72->Draw();
    TMarker3DBox *stave1_73 = new TMarker3DBox((ChipSizeX+ChipDistanceX),+(ChipStaveDistanceY+1.5*ChipDistanceY+2.5*ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_73->Draw();
    TMarker3DBox *stave1_74 = new TMarker3DBox(2*(ChipSizeX+ChipDistanceX),+(ChipStaveDistanceY+1.5*ChipDistanceY+2.5*ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_74->Draw();
    //row8
    TMarker3DBox *stave1_80 = new TMarker3DBox(-2*(ChipSizeX+ChipDistanceX),+(2*ChipStaveDistanceY+1.5*ChipDistanceY+3.5*ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_80->Draw();
    TMarker3DBox *stave1_81 = new TMarker3DBox(-(ChipSizeX+ChipDistanceX),+(2*ChipStaveDistanceY+1.5*ChipDistanceY+3.5*ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_81->Draw();
    TMarker3DBox *stave1_82 = new TMarker3DBox(0,+(2*ChipStaveDistanceY+1.5*ChipDistanceY+3.5*ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_82->Draw();
    TMarker3DBox *stave1_83 = new TMarker3DBox((ChipSizeX+ChipDistanceX),+(2*ChipStaveDistanceY+1.5*ChipDistanceY+3.5*ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_83->Draw();
    TMarker3DBox *stave1_84 = new TMarker3DBox(2*(ChipSizeX+ChipDistanceX),+(2*ChipStaveDistanceY+1.5*ChipDistanceY+3.5*ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_84->Draw();
    //row9
    TMarker3DBox *stave1_90 = new TMarker3DBox(-2*(ChipSizeX+ChipDistanceX),+(2*ChipStaveDistanceY+2.5*ChipDistanceY+4.5*ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_90->Draw();
    TMarker3DBox *stave1_91 = new TMarker3DBox(-1*(ChipSizeX+ChipDistanceX),+(2*ChipStaveDistanceY+2.5*ChipDistanceY+4.5*ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_91->Draw();
    TMarker3DBox *stave1_92 = new TMarker3DBox(0,+(2*ChipStaveDistanceY+2.5*ChipDistanceY+4.5*ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_92->Draw();
    TMarker3DBox *stave1_93 = new TMarker3DBox(1*(ChipSizeX+ChipDistanceX),+(2*ChipStaveDistanceY+2.5*ChipDistanceY+4.5*ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_93->Draw();
    TMarker3DBox *stave1_94 = new TMarker3DBox(2*(ChipSizeX+ChipDistanceX),+(2*ChipStaveDistanceY+2.5*ChipDistanceY+4.5*ChipSizeY),StaveZ[0],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave1_94->Draw();

    //layer2
    //row0
    TMarker3DBox *stave2_00 = new TMarker3DBox(-2*(ChipSizeX+ChipDistanceX),-(2*ChipStaveDistanceY+2.5*ChipDistanceY+4.5*ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_00->Draw();
    TMarker3DBox *stave2_01 = new TMarker3DBox(-1*(ChipSizeX+ChipDistanceX),-(2*ChipStaveDistanceY+2.5*ChipDistanceY+4.5*ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_01->Draw();
    TMarker3DBox *stave2_02 = new TMarker3DBox(0,-(2*ChipStaveDistanceY+2.5*ChipDistanceY+4.5*ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_02->Draw();
    TMarker3DBox *stave2_03 = new TMarker3DBox(1*(ChipSizeX+ChipDistanceX),-(2*ChipStaveDistanceY+2.5*ChipDistanceY+4.5*ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_03->Draw();
    TMarker3DBox *stave2_04 = new TMarker3DBox(2*(ChipSizeX+ChipDistanceX),-(2*ChipStaveDistanceY+2.5*ChipDistanceY+4.5*ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_04->Draw();
    //row1
    TMarker3DBox *stave2_10 = new TMarker3DBox(-2*(ChipSizeX+ChipDistanceX),-(2*ChipStaveDistanceY+1.5*ChipDistanceY+3.5*ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_10->Draw();
    TMarker3DBox *stave2_11 = new TMarker3DBox(-(ChipSizeX+ChipDistanceX),-(2*ChipStaveDistanceY+1.5*ChipDistanceY+3.5*ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_11->Draw();
    TMarker3DBox *stave2_12 = new TMarker3DBox(0,-(2*ChipStaveDistanceY+1.5*ChipDistanceY+3.5*ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_12->Draw();
    TMarker3DBox *stave2_13 = new TMarker3DBox((ChipSizeX+ChipDistanceX),-(2*ChipStaveDistanceY+1.5*ChipDistanceY+3.5*ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_13->Draw();
    TMarker3DBox *stave2_14 = new TMarker3DBox(2*(ChipSizeX+ChipDistanceX),-(2*ChipStaveDistanceY+1.5*ChipDistanceY+3.5*ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_14->Draw();
    //row2
    TMarker3DBox *stave2_20 = new TMarker3DBox(-2*(ChipSizeX+ChipDistanceX),-(ChipStaveDistanceY+1.5*ChipDistanceY+2.5*ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_20->Draw();
    TMarker3DBox *stave2_21 = new TMarker3DBox(-(ChipSizeX+ChipDistanceX),-(ChipStaveDistanceY+1.5*ChipDistanceY+2.5*ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_21->Draw();
    TMarker3DBox *stave2_22 = new TMarker3DBox(0,-(ChipStaveDistanceY+1.5*ChipDistanceY+2.5*ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_22->Draw();
    TMarker3DBox *stave2_23 = new TMarker3DBox((ChipSizeX+ChipDistanceX),-(ChipStaveDistanceY+1.5*ChipDistanceY+2.5*ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_23->Draw();
    TMarker3DBox *stave2_24 = new TMarker3DBox(2*(ChipSizeX+ChipDistanceX),-(ChipStaveDistanceY+1.5*ChipDistanceY+2.5*ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_24->Draw();
    //row3
    TMarker3DBox *stave2_30 = new TMarker3DBox(-2*(ChipSizeX+ChipDistanceX),-(ChipStaveDistanceY+1.5*ChipDistanceY+1.5*ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_30->Draw();
    TMarker3DBox *stave2_31 = new TMarker3DBox(-(ChipSizeX+ChipDistanceX),-(ChipStaveDistanceY+1.5*ChipDistanceY+1.5*ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_31->Draw();
    TMarker3DBox *stave2_32 = new TMarker3DBox(0,-(ChipStaveDistanceY+1.5*ChipDistanceY+1.5*ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_32->Draw();
    TMarker3DBox *stave2_33 = new TMarker3DBox((ChipSizeX+ChipDistanceX),-(ChipStaveDistanceY+1.5*ChipDistanceY+1.5*ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_33->Draw();
    TMarker3DBox *stave2_34 = new TMarker3DBox(2*(ChipSizeX+ChipDistanceX),-(ChipStaveDistanceY+1.5*ChipDistanceY+1.5*ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_34->Draw();
    //row4
    TMarker3DBox *stave2_40 = new TMarker3DBox(-2*(ChipSizeX+ChipDistanceX),-0.5*(ChipDistanceY+ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_40->Draw();
    TMarker3DBox *stave2_41 = new TMarker3DBox(-(ChipSizeX+ChipDistanceX),-0.5*(ChipDistanceY+ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_41->Draw();
    TMarker3DBox *stave2_42 = new TMarker3DBox(0,-0.5*(ChipDistanceY+ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_42->Draw();
    TMarker3DBox *stave2_43 = new TMarker3DBox((ChipSizeX+ChipDistanceX),-0.5*(ChipDistanceY+ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_43->Draw();
    TMarker3DBox *stave2_44 = new TMarker3DBox(2*(ChipSizeX+ChipDistanceX),-0.5*(ChipDistanceY+ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_44->Draw();
    //row5
    TMarker3DBox *stave2_50 = new TMarker3DBox(-2*(ChipSizeX+ChipDistanceX),+0.5*(ChipDistanceY+ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_50->Draw();
    TMarker3DBox *stave2_51 = new TMarker3DBox(-(ChipSizeX+ChipDistanceX),+0.5*(ChipDistanceY+ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_51->Draw();
    TMarker3DBox *stave2_52 = new TMarker3DBox(0,+0.5*(ChipDistanceY+ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_52->Draw();
    TMarker3DBox *stave2_53 = new TMarker3DBox((ChipSizeX+ChipDistanceX),+0.5*(ChipDistanceY+ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_53->Draw();
    TMarker3DBox *stave2_54 = new TMarker3DBox(2*(ChipSizeX+ChipDistanceX),+0.5*(ChipDistanceY+ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_54->Draw();
    //row6
    TMarker3DBox *stave2_60 = new TMarker3DBox(-2*(ChipSizeX+ChipDistanceX),+(ChipStaveDistanceY+1.5*ChipDistanceY+1.5*ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_60->Draw();
    TMarker3DBox *stave2_61 = new TMarker3DBox(-(ChipSizeX+ChipDistanceX),+(ChipStaveDistanceY+1.5*ChipDistanceY+1.5*ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_61->Draw();
    TMarker3DBox *stave2_62 = new TMarker3DBox(0,+(ChipStaveDistanceY+1.5*ChipDistanceY+1.5*ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_62->Draw();
    TMarker3DBox *stave2_63 = new TMarker3DBox((ChipSizeX+ChipDistanceX),+(ChipStaveDistanceY+1.5*ChipDistanceY+1.5*ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_63->Draw();
    TMarker3DBox *stave2_64 = new TMarker3DBox(2*(ChipSizeX+ChipDistanceX),+(ChipStaveDistanceY+1.5*ChipDistanceY+1.5*ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_64->Draw();
    //row7
    TMarker3DBox *stave2_70 = new TMarker3DBox(-2*(ChipSizeX+ChipDistanceX),+(ChipStaveDistanceY+1.5*ChipDistanceY+2.5*ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_70->Draw();
    TMarker3DBox *stave2_71 = new TMarker3DBox(-(ChipSizeX+ChipDistanceX),+(ChipStaveDistanceY+1.5*ChipDistanceY+2.5*ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_71->Draw();
    TMarker3DBox *stave2_72 = new TMarker3DBox(0,+(ChipStaveDistanceY+1.5*ChipDistanceY+2.5*ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_72->Draw();
    TMarker3DBox *stave2_73 = new TMarker3DBox((ChipSizeX+ChipDistanceX),+(ChipStaveDistanceY+1.5*ChipDistanceY+2.5*ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_73->Draw();
    TMarker3DBox *stave2_74 = new TMarker3DBox(2*(ChipSizeX+ChipDistanceX),+(ChipStaveDistanceY+1.5*ChipDistanceY+2.5*ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_74->Draw();
    //row8
    TMarker3DBox *stave2_80 = new TMarker3DBox(-2*(ChipSizeX+ChipDistanceX),+(2*ChipStaveDistanceY+1.5*ChipDistanceY+3.5*ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_80->Draw();
    TMarker3DBox *stave2_81 = new TMarker3DBox(-(ChipSizeX+ChipDistanceX),+(2*ChipStaveDistanceY+1.5*ChipDistanceY+3.5*ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_81->Draw();
    TMarker3DBox *stave2_82 = new TMarker3DBox(0,+(2*ChipStaveDistanceY+1.5*ChipDistanceY+3.5*ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_82->Draw();
    TMarker3DBox *stave2_83 = new TMarker3DBox((ChipSizeX+ChipDistanceX),+(2*ChipStaveDistanceY+1.5*ChipDistanceY+3.5*ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_83->Draw();
    TMarker3DBox *stave2_84 = new TMarker3DBox(2*(ChipSizeX+ChipDistanceX),+(2*ChipStaveDistanceY+1.5*ChipDistanceY+3.5*ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_84->Draw();
    //row9
    TMarker3DBox *stave2_90 = new TMarker3DBox(-2*(ChipSizeX+ChipDistanceX),+(2*ChipStaveDistanceY+2.5*ChipDistanceY+4.5*ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_90->Draw();
    TMarker3DBox *stave2_91 = new TMarker3DBox(-1*(ChipSizeX+ChipDistanceX),+(2*ChipStaveDistanceY+2.5*ChipDistanceY+4.5*ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_91->Draw();
    TMarker3DBox *stave2_92 = new TMarker3DBox(0,+(2*ChipStaveDistanceY+2.5*ChipDistanceY+4.5*ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_92->Draw();
    TMarker3DBox *stave2_93 = new TMarker3DBox(1*(ChipSizeX+ChipDistanceX),+(2*ChipStaveDistanceY+2.5*ChipDistanceY+4.5*ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_93->Draw();
    TMarker3DBox *stave2_94 = new TMarker3DBox(2*(ChipSizeX+ChipDistanceX),+(2*ChipStaveDistanceY+2.5*ChipDistanceY+4.5*ChipSizeY),StaveZ[1],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave2_94->Draw();

    //layer3
    //row0
    TMarker3DBox *stave3_00 = new TMarker3DBox(-2*(ChipSizeX+ChipDistanceX),-(2*ChipStaveDistanceY+2.5*ChipDistanceY+4.5*ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_00->Draw();
    TMarker3DBox *stave3_01 = new TMarker3DBox(-1*(ChipSizeX+ChipDistanceX),-(2*ChipStaveDistanceY+2.5*ChipDistanceY+4.5*ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_01->Draw();
    TMarker3DBox *stave3_02 = new TMarker3DBox(0,-(2*ChipStaveDistanceY+2.5*ChipDistanceY+4.5*ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_02->Draw();
    TMarker3DBox *stave3_03 = new TMarker3DBox(1*(ChipSizeX+ChipDistanceX),-(2*ChipStaveDistanceY+2.5*ChipDistanceY+4.5*ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_03->Draw();
    TMarker3DBox *stave3_04 = new TMarker3DBox(2*(ChipSizeX+ChipDistanceX),-(2*ChipStaveDistanceY+2.5*ChipDistanceY+4.5*ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_04->Draw();
    //row1
    TMarker3DBox *stave3_10 = new TMarker3DBox(-2*(ChipSizeX+ChipDistanceX),-(2*ChipStaveDistanceY+1.5*ChipDistanceY+3.5*ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_10->Draw();
    TMarker3DBox *stave3_11 = new TMarker3DBox(-(ChipSizeX+ChipDistanceX),-(2*ChipStaveDistanceY+1.5*ChipDistanceY+3.5*ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_11->Draw();
    TMarker3DBox *stave3_12 = new TMarker3DBox(0,-(2*ChipStaveDistanceY+1.5*ChipDistanceY+3.5*ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_12->Draw();
    TMarker3DBox *stave3_13 = new TMarker3DBox((ChipSizeX+ChipDistanceX),-(2*ChipStaveDistanceY+1.5*ChipDistanceY+3.5*ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_13->Draw();
    TMarker3DBox *stave3_14 = new TMarker3DBox(2*(ChipSizeX+ChipDistanceX),-(2*ChipStaveDistanceY+1.5*ChipDistanceY+3.5*ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_14->Draw();
    //row2
    TMarker3DBox *stave3_20 = new TMarker3DBox(-2*(ChipSizeX+ChipDistanceX),-(ChipStaveDistanceY+1.5*ChipDistanceY+2.5*ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_20->Draw();
    TMarker3DBox *stave3_21 = new TMarker3DBox(-(ChipSizeX+ChipDistanceX),-(ChipStaveDistanceY+1.5*ChipDistanceY+2.5*ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_21->Draw();
    TMarker3DBox *stave3_22 = new TMarker3DBox(0,-(ChipStaveDistanceY+1.5*ChipDistanceY+2.5*ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_22->Draw();
    TMarker3DBox *stave3_23 = new TMarker3DBox((ChipSizeX+ChipDistanceX),-(ChipStaveDistanceY+1.5*ChipDistanceY+2.5*ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_23->Draw();
    TMarker3DBox *stave3_24 = new TMarker3DBox(2*(ChipSizeX+ChipDistanceX),-(ChipStaveDistanceY+1.5*ChipDistanceY+2.5*ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_24->Draw();
    //row3
    TMarker3DBox *stave3_30 = new TMarker3DBox(-2*(ChipSizeX+ChipDistanceX),-(ChipStaveDistanceY+1.5*ChipDistanceY+1.5*ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_30->Draw();
    TMarker3DBox *stave3_31 = new TMarker3DBox(-(ChipSizeX+ChipDistanceX),-(ChipStaveDistanceY+1.5*ChipDistanceY+1.5*ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_31->Draw();
    TMarker3DBox *stave3_32 = new TMarker3DBox(0,-(ChipStaveDistanceY+1.5*ChipDistanceY+1.5*ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_32->Draw();
    TMarker3DBox *stave3_33 = new TMarker3DBox((ChipSizeX+ChipDistanceX),-(ChipStaveDistanceY+1.5*ChipDistanceY+1.5*ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_33->Draw();
    TMarker3DBox *stave3_34 = new TMarker3DBox(2*(ChipSizeX+ChipDistanceX),-(ChipStaveDistanceY+1.5*ChipDistanceY+1.5*ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_34->Draw();
    //row4
    TMarker3DBox *stave3_40 = new TMarker3DBox(-2*(ChipSizeX+ChipDistanceX),-0.5*(ChipDistanceY+ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_40->Draw();
    TMarker3DBox *stave3_41 = new TMarker3DBox(-(ChipSizeX+ChipDistanceX),-0.5*(ChipDistanceY+ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_41->Draw();
    TMarker3DBox *stave3_42 = new TMarker3DBox(0,-0.5*(ChipDistanceY+ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_42->Draw();
    TMarker3DBox *stave3_43 = new TMarker3DBox((ChipSizeX+ChipDistanceX),-0.5*(ChipDistanceY+ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_43->Draw();
    TMarker3DBox *stave3_44 = new TMarker3DBox(2*(ChipSizeX+ChipDistanceX),-0.5*(ChipDistanceY+ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_44->Draw();
    //row5
    TMarker3DBox *stave3_50 = new TMarker3DBox(-2*(ChipSizeX+ChipDistanceX),+0.5*(ChipDistanceY+ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_50->Draw();
    TMarker3DBox *stave3_51 = new TMarker3DBox(-(ChipSizeX+ChipDistanceX),+0.5*(ChipDistanceY+ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_51->Draw();
    TMarker3DBox *stave3_52 = new TMarker3DBox(0,+0.5*(ChipDistanceY+ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_52->Draw();
    TMarker3DBox *stave3_53 = new TMarker3DBox((ChipSizeX+ChipDistanceX),+0.5*(ChipDistanceY+ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_53->Draw();
    TMarker3DBox *stave3_54 = new TMarker3DBox(2*(ChipSizeX+ChipDistanceX),+0.5*(ChipDistanceY+ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_54->Draw();
    //row6
    TMarker3DBox *stave3_60 = new TMarker3DBox(-2*(ChipSizeX+ChipDistanceX),+(ChipStaveDistanceY+1.5*ChipDistanceY+1.5*ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_60->Draw();
    TMarker3DBox *stave3_61 = new TMarker3DBox(-(ChipSizeX+ChipDistanceX),+(ChipStaveDistanceY+1.5*ChipDistanceY+1.5*ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_61->Draw();
    TMarker3DBox *stave3_62 = new TMarker3DBox(0,+(ChipStaveDistanceY+1.5*ChipDistanceY+1.5*ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_62->Draw();
    TMarker3DBox *stave3_63 = new TMarker3DBox((ChipSizeX+ChipDistanceX),+(ChipStaveDistanceY+1.5*ChipDistanceY+1.5*ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_63->Draw();
    TMarker3DBox *stave3_64 = new TMarker3DBox(2*(ChipSizeX+ChipDistanceX),+(ChipStaveDistanceY+1.5*ChipDistanceY+1.5*ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_64->Draw();
    //row7
    TMarker3DBox *stave3_70 = new TMarker3DBox(-2*(ChipSizeX+ChipDistanceX),+(ChipStaveDistanceY+1.5*ChipDistanceY+2.5*ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_70->Draw();
    TMarker3DBox *stave3_71 = new TMarker3DBox(-(ChipSizeX+ChipDistanceX),+(ChipStaveDistanceY+1.5*ChipDistanceY+2.5*ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_71->Draw();
    TMarker3DBox *stave3_72 = new TMarker3DBox(0,+(ChipStaveDistanceY+1.5*ChipDistanceY+2.5*ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_72->Draw();
    TMarker3DBox *stave3_73 = new TMarker3DBox((ChipSizeX+ChipDistanceX),+(ChipStaveDistanceY+1.5*ChipDistanceY+2.5*ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_73->Draw();
    TMarker3DBox *stave3_74 = new TMarker3DBox(2*(ChipSizeX+ChipDistanceX),+(ChipStaveDistanceY+1.5*ChipDistanceY+2.5*ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_74->Draw();
    //row8
    TMarker3DBox *stave3_80 = new TMarker3DBox(-2*(ChipSizeX+ChipDistanceX),+(2*ChipStaveDistanceY+1.5*ChipDistanceY+3.5*ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_80->Draw();
    TMarker3DBox *stave3_81 = new TMarker3DBox(-(ChipSizeX+ChipDistanceX),+(2*ChipStaveDistanceY+1.5*ChipDistanceY+3.5*ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_81->Draw();
    TMarker3DBox *stave3_82 = new TMarker3DBox(0,+(2*ChipStaveDistanceY+1.5*ChipDistanceY+3.5*ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_82->Draw();
    TMarker3DBox *stave3_83 = new TMarker3DBox((ChipSizeX+ChipDistanceX),+(2*ChipStaveDistanceY+1.5*ChipDistanceY+3.5*ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_83->Draw();
    TMarker3DBox *stave3_84 = new TMarker3DBox(2*(ChipSizeX+ChipDistanceX),+(2*ChipStaveDistanceY+1.5*ChipDistanceY+3.5*ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_84->Draw();
    //row9
    TMarker3DBox *stave3_90 = new TMarker3DBox(-2*(ChipSizeX+ChipDistanceX),+(2*ChipStaveDistanceY+2.5*ChipDistanceY+4.5*ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_90->Draw();
    TMarker3DBox *stave3_91 = new TMarker3DBox(-1*(ChipSizeX+ChipDistanceX),+(2*ChipStaveDistanceY+2.5*ChipDistanceY+4.5*ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_91->Draw();
    TMarker3DBox *stave3_92 = new TMarker3DBox(0,+(2*ChipStaveDistanceY+2.5*ChipDistanceY+4.5*ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_92->Draw();
    TMarker3DBox *stave3_93 = new TMarker3DBox(1*(ChipSizeX+ChipDistanceX),+(2*ChipStaveDistanceY+2.5*ChipDistanceY+4.5*ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_93->Draw();
    TMarker3DBox *stave3_94 = new TMarker3DBox(2*(ChipSizeX+ChipDistanceX),+(2*ChipStaveDistanceY+2.5*ChipDistanceY+4.5*ChipSizeY),StaveZ[2],ChipSizeX/2,ChipSizeY/2,0,0,0);
    stave3_94->Draw();
    
    geom->Update();

}

void display::tracks(int events, bool track_generation, LTrackerTrack& tracker, TCanvas* geom){

//TH1F(name, title, nbins, xlow, xup)
TH1F* hxTR2 = new TH1F("hxTR2", "x_trigger2;x_trigger1;counts", events, -TR2Size[0]*2.5, TR2Size[0]*2.5);
TH1F* hyTR2 = new TH1F("hyTR2", "y_trigger2;y_trigger1;counts", events, -TR2Size[1]/2, TR2Size[1]/2);
TH1F* hzTR2 = new TH1F("hzTR2", "z_trigger2;z_trigger2;counts", events, -TR2Size[2]/2+TR2CenterZ, TR2Size[2]/2+TR2CenterZ);
TH1F* hxTR1 = new TH1F("hxTR1", "x_trigger1;x_trigger1;counts", events, -TR1Size[0]/2, TR1Size[0]/2);
TH1F* hyTR1 = new TH1F("hyTR1", "y_trigger1;y_trigger1;counts", events, -TR1Size[1]*3, TR1Size[1]*3);
TH1F* hzTR1 = new TH1F("hzTR1", "z_trigger1;z_trigger1;counts", events, -TR1Size[2]/2+TR1CenterZ, TR1Size[2]/2+TR1CenterZ);
TH1F* hphi = new TH1F("hphi", "phi;#phi;counts", events, -pi, pi);
TH1F* htheta = new TH1F("htheta", "theta;#theta;counts", events, -pi/2, pi/2);

//puoi settare il seed for reproducibility
TRandom3 *rnd = new TRandom3(0); 

//MC
for (int i=0; i < events; i++){
    double xTR2, yTR2, zTR2, xTR1, yTR1, zTR1, phi, theta; 

    //i layer acquisiscono il segnale quando arriva un segnale AND dagli scintillatori
    //genero traccie misurabili come traccie che passano nei due scintillatori TR2, TR1
    if(track_generation){
        double xTR2_fake = rnd->Uniform(-TR2Size[0]*2,TR2Size[0]*2);     
        if(xTR2_fake>0 && xTR2_fake<TR2Size[0]){xTR2 = xTR2_fake + 0.5*TR2GapX;}
        if(xTR2_fake<0 && xTR2_fake>-TR2Size[0]){xTR2 = xTR2_fake - 0.5*TR2GapX;}
        if(xTR2_fake<2*TR2Size[0] && xTR2_fake>TR2Size[0]){xTR2 = xTR2_fake + 1.5*TR2GapX;}
        if(xTR2_fake>-2*TR2Size[0] && xTR2_fake<-TR2Size[0]){xTR2 = xTR2_fake - 1.5*TR2GapX;}
        yTR2 = rnd->Uniform(-TR2Size[1]/2,TR2Size[1]/2);
        zTR2 = rnd->Uniform(TR2CenterZ-TR2Thickness/2,TR2CenterZ+TR2Thickness/2);

        xTR1 = rnd->Uniform(-TR1Size[0]/2,TR1Size[0]/2);
        double yTR1_fake = rnd->Uniform(-TR1Size[1]*2.5,TR1Size[1]*2.5);     
        if(yTR1_fake<TR1Size[1]/2 && yTR1_fake>-TR1Size[1]/2){yTR1 = yTR1_fake;}
        if(yTR1_fake<1.5*TR1Size[1] && yTR1_fake>0.5*TR1Size[1]){yTR1 = yTR1_fake + TR1GapY;}
        if(yTR1_fake<2.5*TR1Size[1] && yTR1_fake>1.5*TR1Size[1]){yTR1 = yTR1_fake + 2*TR1GapY;}
        if(yTR1_fake>-1.5*TR1Size[1] && yTR1_fake<-0.5*TR1Size[1]){yTR1 = yTR1_fake - TR1GapY;}
        if(yTR1_fake>-2.5*TR1Size[1] && yTR1_fake<-1.5*TR1Size[1]){yTR1 = yTR1_fake - 2*TR1GapY;}
        zTR1 = rnd->Uniform(TR1CenterZ-TR1Thickness/2,TR1CenterZ+TR1Thickness/2);

        //uso ATan2 instead of ATan to avoid losing information
        phi = TMath::ATan2((yTR2-yTR1),(xTR2-xTR1));   
        theta = TMath::ATan2(TMath::Hypot(xTR1-xTR2,yTR1-yTR2),zTR2-zTR1);

        stats::hmgthTR1++;

    }  
    if(!track_generation){
        double xTR2_fake = rnd->Uniform(-TR2Size[0]*2,TR2Size[0]*2);     
        if(xTR2_fake>0 && xTR2_fake<TR2Size[0]){xTR2 = xTR2_fake + 0.5*TR2GapX;}
        if(xTR2_fake<0 && xTR2_fake>-TR2Size[0]){xTR2 = xTR2_fake - 0.5*TR2GapX;}
        if(xTR2_fake<2*TR2Size[0] && xTR2_fake>TR2Size[0]){xTR2 = xTR2_fake + 1.5*TR2GapX;}
        if(xTR2_fake>-2*TR2Size[0] && xTR2_fake<-TR2Size[0]){xTR2 = xTR2_fake - 1.5*TR2GapX;}
        yTR2 = rnd->Uniform(-TR2Size[1]/2,TR2Size[1]/2);
        zTR2 = rnd->Uniform(TR2CenterZ-TR2Thickness/2,TR2CenterZ+TR2Thickness/2);

        phi = rnd->Uniform(2.*pi)-pi;
        double THETA = rnd->Uniform(pi)-pi/2;
        theta = pow(TMath::Cos(THETA),2); 

        xTR1 = xTR2 + (zTR2-TR1CenterZ)*(TMath::Sin(theta))*(TMath::Cos(phi))*(1/(TMath::Cos(theta)));
        yTR1 = yTR2 + (zTR2-TR1CenterZ)*(TMath::Sin(theta))*(TMath::Sin(phi)*(1/(TMath::Cos(theta))));
    }
    hxTR2->Fill(xTR2);
    hyTR2->Fill(yTR2);
    hzTR2->Fill(zTR2);
    hxTR1->Fill(xTR1);
    hyTR1->Fill(yTR1);
    hzTR1->Fill(zTR1);
    hphi->Fill(phi);
    htheta->Fill(theta);

    double xL2 = xTR2 + (zTR2-StaveZ[2])*(TMath::Tan(theta))*(TMath::Cos(phi));
    double yL2 = yTR2 + (zTR2-StaveZ[2])*(TMath::Tan(theta))*(TMath::Sin(phi));
    double xL1 = xTR2 + (zTR2-StaveZ[1])*(TMath::Tan(theta))*(TMath::Cos(phi));
    double yL1 = yTR2 + (zTR2-StaveZ[1])*(TMath::Tan(theta))*(TMath::Sin(phi));
    double xL0 = xTR2 + (zTR2-StaveZ[0])*(TMath::Tan(theta))*(TMath::Cos(phi));
    double yL0 = yTR2 + (zTR2-StaveZ[0])*(TMath::Tan(theta))*(TMath::Sin(phi));

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
    
    LCluster pL2;
    pL2.fill_cluster(pL2, xL2, yL2, StaveZ[2], err_cl, err_cl, err_cl, i);
    tracker.tidy_clusters_lay2.try_emplace(i,pL2);

    LCluster mL1;
    mL1.fill_cluster(mL1, xL1, yL1, StaveZ[1], err_cl, err_cl, err_cl, i);
    tracker.tidy_clusters_lay1.try_emplace(i,mL1);

    LCluster qL0;
    qL0.fill_cluster(qL0, xL0, yL0, StaveZ[0], err_cl, err_cl, err_cl, i);
    tracker.tidy_clusters_lay0.try_emplace(i,qL0);

    //plotting tracks
    Double_t x_line[2] = {xTR2, xTR1};
    Double_t y_line[2] = {yTR2, yTR1};
    Double_t z_line[2] = {zTR2, zTR1};
    TPolyLine3D* line_track = new TPolyLine3D(2, x_line, y_line, z_line);
    line_track->SetLineColor(kBlue);
    line_track->SetLineWidth(2);
    geom->cd();
    line_track->Draw();
    geom->Update();

}

}


