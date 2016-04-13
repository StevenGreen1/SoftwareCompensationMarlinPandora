#include "SoftCompWeightDetermination.h"

//==============================================================================

int main(int argc, char **argv)
{
    gStyle->SetOptStat(kFALSE); 

    SoftCompWeightDetermination(argv[1],argv[2],std::atoi(argv[3]));
}

//==============================================================================

SoftCompWeightDetermination::SoftCompWeightDetermination(TString rootFiles) : 
    m_RootFiles(rootFiles)
{
    this->LoadResults(rootFiles);
}

//==============================================================================

SoftCompWeightDetermination::~SoftCompWeightDetermination()
{
}

//==============================================================================

void SoftCompWeightDetermination::LoadResults(TString rootFiles)
{
    for (std::vector<int>::const_iterator it = m_Energies.begin(); it != m_Energies.end(); ++it)
    {
        TString energy;
        energy.Form("%i",*it);

        TChain *pTChain = new TChain("B4");
        pTChain->Add(rootFiles);

        std::vector<double> energyGap(m_NumLayers,0.f);
        std::vector<double> energyAbs(m_NumLayers,0.f);
 
        m_AveGapE.push_back(aveGapE);
        m_AveAbsE.push_back(aveAbsE);
    }
}

//==============================================================================

