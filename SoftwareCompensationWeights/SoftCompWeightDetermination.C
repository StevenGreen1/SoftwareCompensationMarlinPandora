#include <iostream>
#include <fstream>

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TAxis.h>
#include <TH2F.h>
#include <TString.h>
#include <TLeaf.h>
#include <TRandom.h>
#include <TMinuit.h>
#include <TF1.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <vector>
#include <TROOT.h>
#include <TStyle.h>
#include <TPaveStats.h>
#include <TLegend.h>
#include <TGraph.h>
#include <TStyle.h>

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

#include <iostream>

using namespace std;

bool useEini = false;
bool useSumEhit = false;

const int nE = 10;
double Ebeam[nE] = {10,20,30,40,50,60,70,80,90,100};

const int NBIN = 10;
const int NPAR = 9;

float lowDensity[NBIN]  = {0, 2,   5, 7.5, 9.5, 13, 16,   20, 23.5,  28};
float highDensity[NBIN] = {2, 5, 7.5, 9.5,  13, 16, 20, 23.5,   28,  1e6};
float binDensity[NBIN];

bool m_Debug = false;

typedef std::vector<float> FloatVector;
typedef std::vector<int> IntVector;

double m_p[NBIN*3];


void MakeBinDensity(){
  for (int ibin = 0; ibin < NBIN; ibin++){
    binDensity[ibin] = (lowDensity[ibin]+highDensity[ibin])/2;
    if (ibin==(NBIN-1))
      binDensity[ibin] = 30;
  }
}

/*-------------------------------------------------------------------*/

int FindHitBin(float hitEnergy, float cellvolume){  
  bool m_debugDensity = false;
  
  float hitEnergyDensity = hitEnergy/cellvolume;
  
  int hitBin = -1;
  for (int ibin = 0; ibin < NBIN; ibin++){
    if (hitEnergyDensity>=lowDensity[ibin] && hitEnergyDensity<highDensity[ibin])
      hitBin = ibin;
  }
  
  if (hitBin<0)
    cout << "Can't find hit bin" << endl;
  else
    return hitBin;   
}

/*-------------------------------------------------------------------*/

float FindHitBinDensity(float hitEnergy, float cellvolume){  
  bool m_debugDensity = false;
  
  float hitEnergyDensity = hitEnergy/cellvolume;
  
  float hitBinDensity = 0;
  for (int ibin = 0; ibin < NBIN; ibin++){
    if (hitEnergyDensity>=lowDensity[ibin] && hitEnergyDensity<highDensity[ibin])
      hitBinDensity = binDensity[ibin];
  }
  
  if (hitBinDensity<=0)
    cout << "Can't fin hit bin density" << endl;
  else
    return hitBinDensity;   
}


//========================================================================

class EventC 
{
    public:
        /**
          *  @brief  Constructor
          */
        EventC() 
        { 
            m_HCalBinEnergy->clear();
            m_HCalHitEnergy->clear();
            m_ECalEnergy = 0;
            m_Ereco = 0;
            m_Ebeam = 0;
        }
 
        /**
          *  @brief  Destructor
          */
        ~EventC(){}

        /**
          *  @brief  Constructor
          */
        EventC(double Ebeam, double Ereco, FloatVector *pHitEnergies, FloatVector *pCellSize0, FloatVector *pCellSize1, FloatVector *pCellThickness, IntVector *pHitType)
        { 


Modify this to use the hit type to determine whether it's a hcal or ecal hit.  ecal hits have hit type 1, hcal have hit type 2 and others have hit type 3.


            m_HCalBinEnergy = new FloatVector;
            m_HCalBinEnergy->reserve(NBIN);
            m_HCalHitEnergy = new FloatVector;
            m_HCalHitEnergy->reserve(1000);
            m_Ebeam = Ebeam;
            m_Ereco = (double)Ereco;    
            m_ECalEnergy = 0;
    
            double EnergyPerBin[NBIN];

            for (int i = 0; i < NBIN; i++)
            {
                EnergyPerBin[i] = 0;
            }
    
            for (unsigned int ihit = 0; ihit < HitEnergies->size(); ihit++)
            {
                if (m_Debug) cout << "ihit " << ihit << " energy " << HitEnergies->at(ihit) << endl;
                if (HitEnergies->at(ihit)<=0) continue;

                //cout << "hit " << ihit << " cellsize1 " << CellSize1->at(ihit) << endl;
      
                if (CellSize1->at(ihit)==30)//HCAL hits      
                {
                    m_HCalHitEnergy->push_back(HitEnergies->at(ihit));
                    float cellvolume = (CellSize0->at(ihit))*(CellSize1->at(ihit))*(CellThickness->at(ihit))/1000000.;
                    int hitBin = FindHitBin(HitEnergies->at(ihit),cellvolume);
                    EnergyPerBin[hitBin] += HitEnergies->at(ihit);
                }
                else //ECAL hits 
                {
                    m_ECalEnergy += HitEnergies->at(ihit);
                }
            }

            for (int i = 0; i < NBIN; i++)
            {
                if (m_Debug) cout << "i " << i << " EnergyPerBin " << EnergyPerBin[i] << " inMIP " << EnergyPerBin[i] << endl;
                m_HCalBinEnergy->push_back(EnergyPerBin[i]);
            }
            if (m_Debug) cout << "HCalBinEnergy size " << m_HCalBinEnergy->size() << endl;
        }

        double Ebeam() { return m_Ebeam; }
        double Ereco() { return m_Ereco; }
        double ECalEnergy() { return m_ECalEnergy; }
        FloatVector *HCalBinEnergy() { return m_HCalBinEnergy; }
        FloatVector *HCalHitEnergy() { return m_HCalHitEnergy; }

    private:
        double m_Ebeam;
        double m_Ereco;
        double m_ECalEnergy;
        FloatVector *m_HCalBinEnergy;
        FloatVector *m_HCalHitEnergy;
};

/*----------------------------------------------------------------------------------------------------------- Get list of events -----------------------------------------------------------------------------------------------------------*/
std::vector<EventC> m_EventList;

void GetListOfEvents()
{
    m_EventList.clear();

    TString path = "/r06/lc/sg568/DESYCollaboration/Training/RootFiles";

    for (int ifi = 0; ifi < nE; ifi++)
    {    
        for (int serialNumber = 0; serialNumber <= 10; serialNumber++)
        {
        TString fname = path + "/RegisterHitsForSCAlgorithm_PandoraSettingsDefault_NoECorr_Training_ILD_o1_v06_Kaon0L_pdg_";
        fname += Ebeam[ifi];
        fname += "GeV_";
        fname += serialNumber;
        fname += ".root";

        std::ifstream infile(fname);
        if (!infile.good()) continue;

        cout << "input file " << fname << endl;

        TFile *f = TFile::Open(fname);

        TTree *m_tree = (TTree*)f->Get("HitEnergyTree");

        float energyOfPfo(-1.f);
        float rawEnergyOfCluster(-1.f);
        FloatVector *pHitEnergies(NULL);
        FloatVector *pCellSize0(NULL);
        FloatVector *pCellSize1(NULL);
        FloatVector *pCellThickness(NULL);
        IntVector *pHitType(NULL);

        TBranch *b_EnergyOfPfo;
        TBranch *b_RawEnergyOfCluster;
        TBranch *b_HitEnergies;
        TBranch *b_CellSize0;
        TBranch *b_CellSize1;
        TBranch *b_CellThickness;
        TBranch *b_HitType;

        m_tree->SetBranchAddress("EnergyOfPfo", energyOfPfo, &b_EnergyOfPfo);
        m_tree->SetBranchAddress("RawEnergyOfCluster", rawEnergyOfCluster, &b_RawEnergyOfCluster);
        m_tree->SetBranchAddress("HitEnergies", &pHitEnergies, &b_HitEnergies);
        m_tree->SetBranchAddress("CellSize0", &pCellSize0, &b_CellSize0);
        m_tree->SetBranchAddress("CellSize1", &pCellSize1, &b_CellSize1);
        m_tree->SetBranchAddress("CellThickness", &pCellThickness, &b_CellThickness);
        m_tree->SetBranchAddress("HitType", &pHitType, &b_HitType);

        for (long int iev = 0; iev < m_tree->GetEntries(); iev++)
        {
            m_tree->GetEntry(iev);
      
            if (pHitEnergies->size()==0) continue;
      
            EventC eventC(Ebeam[ifi],pfoE->at(0), pHitEnergies, pCellSize0, pCellSize1, pCellThickness, pHitType);
            m_EventList.push_back(eventC);
        }
        }
    }

    for (int ifi = 0; ifi < nE; ifi++)
    {
        for (int serialNumber = 0; serialNumber <= 10; serialNumber++)
        {
        TString fname = path + "/RegisterHitsForSCAlgorithm_PandoraSettingsDefault_NoECorr_Training_ILD_o1_v06_Neutron_pdg_";
        fname += Ebeam[ifi];
        fname += "GeV_";
        fname += serialNumber;
        fname += ".root";

        std::ifstream infile(fname);
        if (!infile.good()) continue;

        cout << "input file " << fname << endl;

        TFile *f = TFile::Open(fname);

        TTree *m_tree = (TTree*)f->Get("HitEnergyTree");

        float energyOfPfo(-1.f);
        float rawEnergyOfCluster(-1.f);
        FloatVector *pHitEnergies(NULL);
        FloatVector *pCellSize0(NULL);
        FloatVector *pCellSize1(NULL);
        FloatVector *pCellThickness(NULL);
        IntVector *pHitType(NULL);

        TBranch *b_EnergyOfPfo;
        TBranch *b_RawEnergyOfCluster;
        TBranch *b_HitEnergies;
        TBranch *b_CellSize0;
        TBranch *b_CellSize1;
        TBranch *b_CellThickness;
        TBranch *b_HitType;

        m_tree->SetBranchAddress("EnergyOfPfo", energyOfPfo, &b_EnergyOfPfo);
        m_tree->SetBranchAddress("RawEnergyOfCluster", rawEnergyOfCluster, &b_RawEnergyOfCluster);
        m_tree->SetBranchAddress("HitEnergies", &pHitEnergies, &b_HitEnergies);
        m_tree->SetBranchAddress("CellSize0", &pCellSize0, &b_CellSize0);
        m_tree->SetBranchAddress("CellSize1", &pCellSize1, &b_CellSize1);
        m_tree->SetBranchAddress("CellThickness", &pCellThickness, &b_CellThickness);
        m_tree->SetBranchAddress("HitType", &pHitType, &b_HitType);

        for (long int iev = 0; iev < m_tree->GetEntries(); iev++)
        {
            m_tree->GetEntry(iev);

            if (pHitEnergies->size()==0) continue;

            EventC eventC(Ebeam[ifi],pfoE->at(0), pHitEnergies, pCellSize0, pCellSize1, pCellThickness, pHitType);
            m_EventList.push_back(eventC);
        }
        }
    }

    cout << "m_EventList size " << m_EventList.size() << endl;
}

/*----------------------------------------------------------------------------------------------------------- MakeHitDensityHisto -----------------------------------------------------------------------------------------------------------*/
/*
void MakeBinDensityHisto(double cellsize = 30){

  GetListOfEvents(cellsize);

  //gStyle->SetOptStat(0);

  TH1F *hECal[10];
  TProfile *hHCal[10];
  TH1F *hHCalHitEnergy[10];
  TH1F *hsumEnergyHCal[10];
  TString hname = ""; 
  
  for (int index = 0; index < nE; index++){
    hname = "EnergyEcal"; hname += Ebeam[index]; hname += "GeV";
    hECal[index] = new TH1F(hname, hname, 100, 0, 100);
    
    hname = "EnergyHcal"; hname += Ebeam[index]; hname += "GeV";
    hHCal[index] = new TProfile(hname, hname, NBIN, 0, NBIN);  

    hname = "HcalHitEnergy"; hname += Ebeam[index]; hname += "GeV";
    hHCalHitEnergy[index] = new TH1F(hname, hname, 120, 0, 120);  

    hname = "sumEnergyHCal"; hname += Ebeam[index]; hname += "GeV";
    hsumEnergyHCal[index] = new TH1F(hname, hname, 100, 0, 100);
  }

  TLegend *leg = new TLegend(0.6, 0.65, 0.9, 0.85);

  for (unsigned int iev = 0; iev < m_EventList.size(); iev++){
    EventC event = m_EventList[iev];
    double ebeam = event.Ebeam();
    int index = -1;
    for (int i = 0; i < nE; i++){
      if (ebeam == Ebeam[i])
	index = i;
    }

    double ECalEnergy = event.ECalEnergy();
    FloatVector *HCalBinEnergy = event.HCalBinEnergy();
    FloatVector *HCalHitEnergy = event.HCalHitEnergy();

    if (ECalEnergy>0)
      hECal[index]->Fill(ECalEnergy);

    for (unsigned int i = 0; i < HCalHitEnergy->size(); i++){
      if (HCalHitEnergy->at(i)<=0) continue;
      hHCalHitEnergy[index]->Fill(HCalHitEnergy->at(i)/0.0225);
    }
    
    double sumAllEHCal = 0;
    for (int i = 0; i < NBIN; i++){
      hHCal[index]->Fill(i,HCalBinEnergy->at(i));      
      sumAllEHCal += HCalBinEnergy->at(i);
    }
    sumAllEHCal += ECalEnergy;

    hsumEnergyHCal[index]->Fill(sumAllEHCal);
    
    //if (m_Debug){
    cout << "ECalEnergy " << ECalEnergy << endl;
    cout << "ereco " << event.Ereco() << " sumAllEHCal " << sumAllEHCal << endl;
      //}
  }  
  
  TCanvas *c1 = new TCanvas("c1", "", 1100, 600);
  c1->Divide(2);
  c1->cd(1);
  hECal[4]->SetLineColor(kBlue+2);
  hECal[4]->Draw(); 
  for (int i = 0; i < nE; i++){
    hECal[i]->SetLineWidth(2);
    //hECal[i]->SetLineColor(i+1);
    if (i==9)
      hECal[i]->SetLineColor(kOrange+1);
    TString legname = ""; legname += Ebeam[i]; legname += "GeV";
    leg->AddEntry(hECal[i],legname,"l");
    hECal[i]->Draw("same");      
  }
  c1->cd(2);
  hsumEnergyHCal[4]->SetLineColor(kBlue+2);
  hsumEnergyHCal[4]->Draw(); 
  for (int i = 0; i < nE; i++){
    hsumEnergyHCal[i]->SetLineWidth(2);
    //hECal[i]->SetLineColor(i+1);
    if (i==9)
      hsumEnergyHCal[i]->SetLineColor(kOrange+1);
    TString legname = ""; legname += Ebeam[i]; legname += "GeV";
    leg->AddEntry(hsumEnergyHCal[i],legname,"l");
    hsumEnergyHCal[i]->Draw("same");      
  }
  //leg->Draw();

  TCanvas *c2 = new TCanvas("c2", "", 1100, 600);
  c2->Divide(2);
  c2->cd(1);
  hHCalHitEnergy[0]->SetMinimum(0);
  //hHCalHitEnergy[0]->SetMaximum(500000);
  hHCalHitEnergy[4]->SetLineColor(kBlue+2);
  hHCalHitEnergy[4]->SetFillColor(kAzure-5);
  hHCalHitEnergy[4]->GetXaxis()->SetTitle("Hit energy [MIP]");
  hHCalHitEnergy[4]->GetXaxis()->SetTitleSize(0.045);
  hHCalHitEnergy[4]->Draw();
  for (int i = 0; i < nE; i++){
    hHCalHitEnergy[i]->SetLineWidth(2);
    //hHCalHitEnergy[i]->SetLineColor(i+1);
    if (i==9)
      hHCalHitEnergy[i]->SetLineColor(kOrange+1);
    hHCalHitEnergy[i]->Draw("same");      
  }
  //leg->Draw();
  c2->cd(2);
  hHCal[0]->SetMinimum(0);
  hHCal[0]->SetMaximum(2);
  //hHCal[0]->DrawNormalized();
  hHCal[4]->SetLineColor(kBlue+2);
  hHCal[4]->SetFillColor(kAzure-5);
  hHCal[4]->GetXaxis()->SetTitle("Bin energy");
  hHCal[4]->GetXaxis()->SetTitleSize(0.045);
  hHCal[4]->Draw();
  for (int i = 0; i < nE; i++){
    hHCal[i]->SetLineWidth(2);
    //hHCal[i]->SetLineColor(i+1);
    if (i==9)
      hHCal[i]->SetLineColor(kOrange+1);
    hHCal[i]->Draw("same");      
  }
  //leg->Draw();

}
*/
/*----------------------------------------------------------------------------------------------------------- Fitting  -----------------------------------------------------------------------------------------------------------*/

double chi2(const double *par){
  bool m_DebugChi2 = false;

  double Chi2 = 0;

  for (unsigned int iev = 0; iev < m_EventList.size(); iev++){
    EventC event = m_EventList[iev];

    double ebeam = event.Ebeam();
    double ereco = event.Ereco();

    double ECalEnergy = event.ECalEnergy();
    FloatVector *HCalBinEnergy = event.HCalBinEnergy();

    if (m_DebugChi2) cout << "ebeam " << ebeam << " ereco " << ereco << endl;
    if (m_DebugChi2) cout << "HCalBinEnergy size " << HCalBinEnergy->size() << endl;

    double E_SC = 0;

    //Add energy from ECal
    E_SC += ECalEnergy;

    if (m_DebugChi2) cout << "ECalEnergy " << ECalEnergy << endl;

    //Sum bins in HCal
    for (unsigned int ibin = 0; ibin < HCalBinEnergy->size(); ibin++){
      double BinEnergy = HCalBinEnergy->at(ibin);
      
      //if (BinEnergy<=0) continue;

      if (m_DebugChi2) cout << "ibin " << ibin << endl;
      if (m_DebugChi2) cout << "HCalBinEnergy " << HCalBinEnergy->at(ibin) << endl;
      
      double Weight = 0;
      if (NPAR==3)
	Weight = par[0]*exp(par[1]*binDensity[ibin])+par[2];
      else if (NPAR==9){
	float p1 = par[0] + par[1]*ebeam + par[2]*ebeam*ebeam;
	float p2 = par[3] + par[4]*ebeam + par[5]*ebeam*ebeam;
	float p3 = par[6]/(par[7]+exp(par[8]*ebeam));

	Weight = p1*exp(p2*binDensity[ibin])+p3;
      }

      if (m_DebugChi2) cout << "binDensity " << binDensity[ibin] << " Weight = " << Weight << endl;      
      
      double BinEnergyWeighted = BinEnergy*Weight;

      //cout << "bin " << ibin << " BinEnergy " << BinEnergy << " Weight " << Weight << " BinEnergyWeighted " << BinEnergyWeighted << endl;

      E_SC += BinEnergyWeighted;
    }

    if (m_DebugChi2) cout << "E_SC " << E_SC << endl;
    //cout << "E_SC " << E_SC << endl;

    double diff = (E_SC - ebeam)*(E_SC - ebeam)/(0.5*ebeam);
    Chi2 += diff;
  }

  return Chi2;
}

//========================================================================================

void Fit()
{
    //Make the list of density values
    MakeBinDensity();

    //Get list of events
    GetListOfEvents();

    int npar = NPAR;
    double step = 0.000001;

    double fitPar0 = 2.5;
    double fitPar1 = -0.02;
    double fitPar2 = 0.5;
    double fitPar3, fitPar4, fitPar5, fitPar6, fitPar7, fitPar8;

    const char *minName = "Minuit2";
    const char *algoName = "Minos";
  
    //ROOT::Minuit2::Minuit2Minimizer
    ROOT::Math::Minimizer* myMinuit = ROOT::Math::Factory::CreateMinimizer(minName,algoName);

    // set tolerance , etc...
    myMinuit->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
    //myMinuit->SetMaxIterations(10000);  // for GSL
    myMinuit->SetTolerance(1.);
    myMinuit->SetPrintLevel(1);
  
    // create function wrapper for minmizer
    ROOT::Math::Functor f(&chi2,npar);

    myMinuit->SetFunction(f);

    if (NPAR==3)
    {
        myMinuit->SetVariable(0, "p0", fitPar0, step);
        myMinuit->SetVariable(1, "p1", fitPar1, step);
        myMinuit->SetVariable(2, "p2", fitPar2, step);
    }

    else if (NPAR==9)
    {
        fitPar0 = 1.0000000;
        myMinuit->SetVariable(0, "p0", fitPar0, step);
        fitPar1 = -0.02999999;
        myMinuit->SetVariable(1, "p1", fitPar1, step);
        fitPar2 = 0.0002999999;
        myMinuit->SetVariable(2, "p2", fitPar2, step);    
        fitPar3 = -0.01999999;
        myMinuit->SetVariable(3, "p3", fitPar3, step);    
        fitPar4 = -0.0003999999;
        myMinuit->SetVariable(4, "p4", fitPar4, step);    
        fitPar5 = -0.000001999999;
        myMinuit->SetVariable(5, "p5", fitPar5, step);    
        fitPar6 = 0.1999999;
        myMinuit->SetVariable(6, "p6", fitPar6, step);    
        fitPar7 = 0.2999999;
        myMinuit->SetVariable(7, "p7", fitPar7, step);    
        fitPar8 = -0.0999999;
        myMinuit->SetVariable(8, "p8", fitPar8, step);    
    }

    // do the minimization
    myMinuit->Minimize();

    // Printouts

    const double *xs = myMinuit->X();
    std::cout << "Minimum: chi2 = " << myMinuit->MinValue()  << std::endl;

    // expected minimum is 0
    if ( myMinuit->MinValue()  < 1.E-4  && f(xs) < 1.E-4)
    {
        std::cout << "Minimizer " << minName << " - " << algoName << "   converged to the right minimum" << std::endl;
    }

    else 
    {
        std::cout << "Minimizer " << minName << " - " << algoName << "   failed to converge !!!" << std::endl;
        std::cout << "MinValue " << myMinuit->MinValue() << std::endl;
    }
  
    for (int i = 0; i < npar; i++)
    {
        cout << "Parameter " << i << " is " << xs[i] << endl;
        m_p[i] = xs[i];
    }

    for (int i = 0; i < npar; i++)
    {
        cout << xs[i] << " ";
    }
}

