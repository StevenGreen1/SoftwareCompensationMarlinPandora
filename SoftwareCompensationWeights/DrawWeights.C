#include <TGraph.h>
#include <TLegend.h>
#include <TAxis.h>

#include <iostream>
#include <fstream>

using namespace std;

float mip2gev = 0.0225;

void DrawWeights(int Cell = 30, bool drawOldParameters = false){

  float m_dim = Cell/100.;

  const int nE = 14;
  double Ebeam[nE] = {1, 3, 5, 7, 10, 20, 30, 40, 50, 60, 70, 80, 90, 95};

  const int nbinWei = 10;
  float lowMIP[nbinWei]  = {0.3,  2, 5.5,  8, 10, 14, 17, 21, 25, 30};//First bin has to start at 0.3 now as the cut in MarlinPandora changed
  float highMIP[nbinWei] = {  2, 5.5,   8, 10, 14, 17, 21, 25, 30, 1e6};
  
  float lowDensity[nbinWei] , highDensity[nbinWei];
  float volume = m_dim*m_dim*0.265;
  
  for (int i = 0; i < nbinWei; i++){
    lowDensity[i]  = lowMIP[i]*mip2gev/volume;
    highDensity[i] = highMIP[i]*mip2gev/volume;
  }
  
  float rho[nbinWei];
  for (int i = 0; i < nbinWei; i++){
    rho[i] = (lowDensity[i]+highDensity[i])/2;
    if (i==(nbinWei-1))
      rho[i] = 40*mip2gev/volume;
  }

  //Old definition for 3x3 cm*cm - CALICE meeting
  //double rho[12] = {3.25, 7, 11.5, 17, 23.5, 31, 39.5, 49, 59.5, 71, 83.5, 97};
    
  TGraph *gWei[nE];

  //float par1[9] = {2.4515, -0.0253871, 0.000409333, -0.0238832, -0.000260914, -3.16977e-06, 0.472525, 0.648817, -0.0589877};
  //float par1[9] = {2.47822, -0.0273163, 0.000431585, -0.128046, -0.00122036, -1.88328e-05, 0.460827, 0.633081, -0.0597235};
  float par1[9] = {2.41985,-0.0902137,0.00133613,-0.118089,0.00373349,-0.00014658,0.126489,0.16879,-0.123544};
  //float par2[9] = {1.78842, -0.0390278, 0.000733125, -0.0710077, -0.0017393, -5.59103e-05, 0.365751, 0.485933, -0.0896917};
  float par2[9] = {2.43157, -0.0689246, 0.000814196, -0.0953147, 0.00371114, -8.79431e-05, 0.0542545, 0.071707,3 -0.102446};
/*
Prameter 0 is 2.41985
Parameter 1 is -0.0902137
Parameter 2 is 0.00133613
Parameter 3 is -0.118089
Parameter 4 is 0.00373349
Parameter 5 is -0.00014658
Parameter 6 is 0.126489
Parameter 7 is 0.16879
Parameter 8 is -0.123544
*/
/*
Parameter 0 is 2.64238
Parameter 1 is -0.105172
Parameter 2 is 0.00149956
Parameter 3 is -0.138385
Parameter 4 is 0.00689485
Parameter 5 is -0.000185867
Parameter 6 is 0.0729758
Parameter 7 is 0.0973938
Parameter 8 is -0.128571
*/

  for (int ie = 0; ie < nE; ie++){
    gWei[ie] = new TGraph();
    gWei[ie]->SetMarkerStyle(20);
    gWei[ie]->SetMarkerSize(0.8);
    
    cout << "Ebeam " << Ebeam[ie] << endl;
    
    //Now make the weight array

    for (int ibin = 0; ibin < nbinWei; ibin++){
      float p1 = 0, p2 = 0, p3 = 0;

      if (drawOldParameters){
	p1 = par1[0] + par1[1]*Ebeam[ie] + par1[2]*Ebeam[ie]*Ebeam[ie];
	p2 = par1[3] + par1[4]*Ebeam[ie] + par1[5]*Ebeam[ie]*Ebeam[ie];
	p3 = par1[6]/(par1[7]+exp(par1[8]*Ebeam[ie]));
      } else {
	p1 = par2[0] + par2[1]*Ebeam[ie] + par2[2]*Ebeam[ie]*Ebeam[ie];
	p2 = par2[3] + par2[4]*Ebeam[ie] + par2[5]*Ebeam[ie]*Ebeam[ie];
	p3 = par2[6]/(par2[7]+exp(par2[8]*Ebeam[ie]));
      }

      double weight = p1*exp(p2*rho[ibin]) + p3;
      cout << "ibin " << ibin << " weight " << weight << endl;           
      
      gWei[ie]->SetPoint(ibin,rho[ibin],weight);
      
    }
  }

  gWei[4]->SetMarkerColor(4);
  gWei[6]->SetMarkerColor(kOrange+1);
  gWei[8]->SetMarkerColor(2);
  gWei[10]->SetMarkerColor(3);

  gWei[0]->SetMinimum(0);
  gWei[0]->SetMaximum(5);
  gWei[4]->SetMinimum(0);
  gWei[4]->SetMaximum(6);
  gWei[6]->SetMinimum(0);
  gWei[6]->SetMaximum(6);

  TLegend *leg = new TLegend(0.6,0.6,0.85,0.85,"Weights");
  leg->SetFillColor(0);
  leg->SetLineColor(0);
  TString legname = "Beam Energy: "; legname += Ebeam[4]; legname += " GeV";
  leg->AddEntry(gWei[4],legname,"p");
  legname = "Beam Energy: "; legname += Ebeam[6]; legname += " GeV";
  leg->AddEntry(gWei[6],legname,"p");
  legname = "Beam Energy: "; legname += Ebeam[8]; legname += " GeV";
  leg->AddEntry(gWei[8],legname,"p");
  legname = "Beam Energy: "; legname += Ebeam[10]; legname += " GeV";
  leg->AddEntry(gWei[10],legname,"p");

  gWei[4]->GetXaxis()->SetTitle("Energy Density [GeV/dm^{3}]");
  gWei[4]->GetXaxis()->SetTitleSize(0.055);
  gWei[4]->GetXaxis()->SetLabelSize(0.055);
  gWei[4]->GetYaxis()->SetTitle("Weight values");
  gWei[4]->GetYaxis()->SetTitleSize(0.055);
  gWei[4]->GetYaxis()->SetLabelSize(0.055);  

  gWei[4]->Draw("AP");
  gWei[6]->Draw("PS");
  gWei[8]->Draw("PS");
  gWei[10]->Draw("PS");
//  gWei[9]->Draw("PS");
  leg->Draw();

  double Ebeam[nE] = {1, 3, 5, 7, 10, 20, 30, 40, 50, 60, 70, 80, 90, 95};
}
