#include <iostream>
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

double m_dim = 0.3;
double mip2gev = 0.0225;

const int NBIN = 10;
const int NPAR = 9;

float lowMIP[NBIN]  = {0.3,  2, 5.5,  8, 10, 14, 17, 21, 25, 30};//First bin has to start at 0.3 now as the cut in MarlinPandora changed
float highMIP[NBIN] = {  2, 5.5,   8, 10, 14, 17, 21, 25, 30, 1e6};
float rho[NBIN];

void MakeRho(){
  double cellthickness = 0.265;//26.5mm according to simulation
  double volume = m_dim*m_dim*cellthickness;
  
  for (int ibin = 0; ibin < NBIN; ibin++){
    rho[ibin] = (lowMIP[ibin]+highMIP[ibin])/2;
    if (ibin==(NBIN-1))
      rho[ibin] = 40;

    rho[ibin] *= mip2gev;
    rho[ibin] /= volume;
  }
}

const int nE = 10;
//double Ebeam[nE] = {1, 3, 5, 7};
double Ebeam[nE] = {10, 20, 30, 40, 50, 60, 70, 80, 90, 95};

bool m_Debug = false;

typedef std::vector<float> FloatVector;

double m_p[NBIN*3];


/*----------------------------------------------------------------------------------------------------------- Energy density bin definition -----------------------------------------------------------------------------------------------------------*/

int FindHitEnergyBin(double hitEnergy){  
  bool m_debugDensity = false;

  //double volume = m_dim*m_dim*0.05;
  
  double rho = 0;

  double hitEnergyMIP = hitEnergy/mip2gev;

  if (m_Debug) cout << "hitEnergy " << hitEnergy << " hitEnergyMIP " << hitEnergyMIP << endl;

  int index = -1;

  for (int ibin = 0; ibin < NBIN; ibin++){
    if (m_Debug) cout << "ibin " << ibin << " lowMIP " << lowMIP[ibin] << " highMIP " << highMIP[ibin] << endl;

    if (hitEnergyMIP>=lowMIP[ibin] && hitEnergyMIP<highMIP[ibin]){
      if (m_Debug) cout << "hit in bin " << ibin << endl;
      index = ibin;
    }
  }

  if (m_Debug) cout << "index " << index << endl;

  if (index<0){
    cout << "In FindBinDensity: something is wrong!" << endl;
    cout << "hitEnergy " << hitEnergy << " hitEnergyMIP " << hitEnergyMIP << endl;
  }
  return index;   
}


/*----------------------------------------------------------------------------------------------------------- Define class EventC -----------------------------------------------------------------------------------------------------------*/

class EventC {

public:
  EventC() { 
    m_HCalBinEnergy->clear();
    m_HCalHitEnergy->clear();
    m_ECalEnergy = 0;
    m_Ereco = 0;
    m_Ebeam = 0;
  }
  ~EventC(){}

  EventC(double Ebeam, double Ereco, FloatVector *hitECALenergyvec, FloatVector *hitHCALenergyvec){ 
    m_HCalBinEnergy = new FloatVector;
    m_HCalBinEnergy->reserve(NBIN);
    m_HCalHitEnergy = new FloatVector;
    m_HCalHitEnergy->reserve(1000);
    m_Ebeam = Ebeam;
    m_Ereco = (double)Ereco;    
    m_ECalEnergy = 0;

    for (unsigned int ihit = 0; ihit < hitECALenergyvec->size(); ihit++){
      if (m_Debug) cout << "ihit " << ihit << " energy " << hitECALenergyvec->at(ihit) << endl;
      if (hitECALenergyvec->at(ihit)<=0) continue;
      m_ECalEnergy += hitECALenergyvec->at(ihit);
    }
    
    if (m_Debug) cout << "m_ECalEnergy " << m_ECalEnergy << endl;

    double EnergyPerBin[NBIN];
    for (int i = 0; i < NBIN; i++){
      EnergyPerBin[i] = 0;
    }

    for (unsigned int ihit = 0; ihit < hitHCALenergyvec->size(); ihit++){
      if (m_Debug) cout << "ihit " << ihit << " energy " << hitHCALenergyvec->at(ihit) << endl;
      if (hitHCALenergyvec->at(ihit)<=0) continue;
      m_HCalHitEnergy->push_back(hitHCALenergyvec->at(ihit));
      int indexBin = FindHitEnergyBin(hitHCALenergyvec->at(ihit));
      if (indexBin<0)
	cout << "Something is wrong!" << endl;
      else
	EnergyPerBin[indexBin] += hitHCALenergyvec->at(ihit);
    }

    if (m_Debug) cout << "NBIN " << NBIN << endl;

    for (int i = 0; i < NBIN; i++){
      if (m_Debug) cout << "i " << i << " EnergyPerBin " << EnergyPerBin[i] 
			<< " inMIP " << EnergyPerBin[i] << endl;
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

void GetListOfEvents(int cellsize = 30){
  //m_dim = (double)cellsize/100.;
  m_EventList.clear();

  //Get tree
  TString path = "/nfs/dust/ilc/user/huonglan/HCAL_Optimisation/SoftwareCompensation/RECO/Output/HcalTimingCut100/";

  for (int ifi = 0; ifi < nE; ifi++){    
    //if (ifi!=11) continue;

    TString fname = path; fname += Ebeam[ifi]; 
    fname += "GeV_Cell"; fname += cellsize; fname += ".root";

    cout << "input file " << fname << endl;

    TFile *f = TFile::Open(fname);

    TTree *m_tree = (TTree*)f->Get("PfoAnalysisTree");

    float pfoE = 0;
    FloatVector *hitECALEnergyVec = 0;
    FloatVector *hitHCALEnergyVec = 0;
    
    TBranch *b_hitECALEnergyVec;
    TBranch *b_hitHCALEnergyVec;

    TBranch *b_pfoEnergyTotal;

    m_tree->SetBranchAddress("hitECALEnergyVec", &hitECALEnergyVec, &b_hitECALEnergyVec);
    m_tree->SetBranchAddress("hitHCALEnergyVec", &hitHCALEnergyVec, &b_hitHCALEnergyVec);
    m_tree->SetBranchAddress("pfoEnergyTotal", &pfoE, &b_pfoEnergyTotal);

    for (long int iev = 0; iev < m_tree->GetEntries(); iev++){
      m_tree->GetEntry(iev);

      if (hitHCALEnergyVec->size()==0 && hitECALEnergyVec->size()==0) continue;
      if (pfoE<=0) continue;

      EventC eventC(Ebeam[ifi],pfoE,hitECALEnergyVec,hitHCALEnergyVec);
    
      m_EventList.push_back(eventC);
    }
  }

  cout << "m_EventList size " << m_EventList.size() << endl;
  
}

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
      
      if (m_DebugChi2) cout << "ibin " << ibin << endl;
      if (m_DebugChi2) cout << "HCalBinEnergy " << HCalBinEnergy->at(ibin) << endl;
      
      double Weight = 0;
      if (NPAR==3)
	Weight = par[0]*exp(par[1]*rho[ibin])+par[2];
      else if (NPAR==9){
	float p1 = par[0] + par[1]*ebeam + par[2]*ebeam*ebeam;
	float p2 = par[3] + par[4]*ebeam + par[5]*ebeam*ebeam;
	float p3 = par[6]/(par[7]+exp(par[8]*ebeam));

	Weight = p1*exp(p2*rho[ibin])+p3;
      }

      if (m_DebugChi2) cout << "rho " << rho[ibin] << " Weight = " << Weight << endl;      
      
      double BinEnergyWeighted = BinEnergy*Weight;

      E_SC += BinEnergyWeighted;
    }

    if (m_DebugChi2) cout << "E_SC " << E_SC << endl;
    //cout << "E_SC " << E_SC << endl;

    double diff = (E_SC - ebeam)*(E_SC - ebeam)/(0.5*ebeam);
    Chi2 += diff;
  }

  return Chi2;
}


void Fit(int cellsize = 30){
  //Make the list of density values
  m_dim = cellsize/100.;
  MakeRho();

  //Get list of events
  GetListOfEvents(cellsize);

  int npar = NPAR;
  double step = 0.000001;

  double fitPar0 = 2.5;
  double fitPar1 = -0.02;
  double fitPar2 = 0.5;
  double fitPar3, fitPar4, fitPar5, fitPar6, fitPar7, fitPar8;

  const char *minName = "Minuit2";
  const char *algoName = "Migrad";
  
  //ROOT::Minuit2::Minuit2Minimizer
  ROOT::Math::Minimizer* myMinuit = ROOT::Math::Factory::CreateMinimizer(minName,algoName);

  // set tolerance , etc...
  myMinuit->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
  myMinuit->SetMaxIterations(10000);  // for GSL
  //myMinuit->SetTolerance(1.);
  myMinuit->SetPrintLevel(1);
  
  // create function wrapper for minmizer
  ROOT::Math::Functor f(&chi2,npar);

  myMinuit->SetFunction(f);

  if (NPAR==3){
    myMinuit->SetVariable(0, "p0", fitPar0, step);
    myMinuit->SetVariable(1, "p1", fitPar1, step);
    myMinuit->SetVariable(2, "p2", fitPar2, step);
  } else if (NPAR==9){
    fitPar0 = 1.999999;
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

  const double *xs = myMinuit->X();
  std::cout << "Minimum: chi2 = " << myMinuit->MinValue()  << std::endl;

  // expected minimum is 0
  if ( myMinuit->MinValue()  < 1.E-4  && f(xs) < 1.E-4)
    std::cout << "Minimizer " << minName << " - " << algoName
	      << "   converged to the right minimum" << std::endl;
  else {
    std::cout << "Minimizer " << minName << " - " << algoName
	      << "   failed to converge !!!" << std::endl;
    std::cout << "MinValue " << myMinuit->MinValue() << std::endl;
  }
  
  for (int i = 0; i < npar; i++){
    cout << "Parameter " << i << " is " << xs[i] << endl;
    m_p[i] = xs[i];
  }

  for (int i = 0; i < npar; i++){
    cout << xs[i] << " ";
    //if (i%3==2) cout << endl;
  }
  cout << endl;
}

