#include "SoftCompWeightDetermination.h"

//==============================================================================

int main(int argc, char **argv)
{
    gStyle->SetOptStat(kFALSE); 
    std::string argument1 = argv[1];
    SoftCompWeightDetermination *pSoftCompWeightDetermination = new SoftCompWeightDetermination(argument1);
}

//==============================================================================

SoftCompWeightDetermination::SoftCompWeightDetermination(std::string rootFiles) : 
    m_RootFiles(rootFiles)
{
    int energyArray[10] = { 10, 20, 30, 40, 50, 60, 70, 80, 90, 100 };
    IntVector energies(energyArray, energyArray + sizeof(energyArray)/sizeof(energyArray[0]));
    m_Energies.insert(m_Energies.end(), energies.begin(), energies.end());

    this->MakeBinDensities();
    this->LoadEvents();
    this->Fit();
}

//==============================================================================

SoftCompWeightDetermination::~SoftCompWeightDetermination()
{
}

//==============================================================================

void SoftCompWeightDetermination::Fit()
{
    const char *minName = "Minuit2";
    const char *algoName = "Minos";

    ROOT::Math::Minimizer* pMinimizer = ROOT::Math::Factory::CreateMinimizer(minName, algoName);
    pMinimizer->SetMaxFunctionCalls(1000000);
    pMinimizer->SetTolerance(1.f);
    pMinimizer->SetPrintLevel(1);

    ROOT::Math::Functor functorChi2(this, &SoftCompWeightDetermination::Chi2,9);
    pMinimizer->SetFunction(functorChi2);

    const float step = 0.000001;
    double fitPar0 = 1.0000000;
    pMinimizer->SetVariable(0, "p0", fitPar0, step);
    double fitPar1(-0.02999999);
    pMinimizer->SetVariable(1, "p1", fitPar1, step);
    double fitPar2(0.0002999999);
    pMinimizer->SetVariable(2, "p2", fitPar2, step);
    double fitPar3(-0.01999999);
    pMinimizer->SetVariable(3, "p3", fitPar3, step);
    double fitPar4(-0.0003999999);
    pMinimizer->SetVariable(4, "p4", fitPar4, step);
    double fitPar5(-0.000001999999);
    pMinimizer->SetVariable(5, "p5", fitPar5, step);
    double fitPar6(0.1999999);
    pMinimizer->SetVariable(6, "p6", fitPar6, step);
    double fitPar7(0.2999999);
    pMinimizer->SetVariable(7, "p7", fitPar7, step);
    double fitPar8(-0.0999999);
    pMinimizer->SetVariable(8, "p8", fitPar8, step);

    // do the minimization
    pMinimizer->Minimize();

    // Printouts
    const double *xs = pMinimizer->X();
    std::cout << "Minimum: chi2 = " << pMinimizer->MinValue()  << std::endl;

    // expected minimum is 0

    if ( pMinimizer->MinValue()  < 1.E-4  && this->Chi2(xs) < 1.E-4)
    {
        std::cout << "Minimizer " << minName << " - " << algoName << "   converged to the right minimum" << std::endl;
    }

    else
    {
        std::cout << "Minimizer " << minName << " - " << algoName << "   failed to converge !!!" << std::endl;
        std::cout << "MinValue " << pMinimizer->MinValue() << std::endl;
    }

    for (int i = 0; i < 9; i++)
    {
        std::cout << "Parameter " << i << " is " << xs[i] << std::endl;
    }

    for (int i = 0; i < 9; i++)
    {
        std::cout << xs[i] << " ";
    }
}

//==============================================================================

void SoftCompWeightDetermination::MakeBinDensities()
{
    const float lowDensityBins[10]  = {0, 2,   5, 7.5, 9.5, 13, 16,   20, 23.5,  28};
    FloatVector lowEdgeDensityBins(lowDensityBins, lowDensityBins + sizeof(lowDensityBins)/sizeof(lowDensityBins[0]));
    m_LowEdgeDensityBins.insert(m_LowEdgeDensityBins.end(), lowEdgeDensityBins.begin(), lowEdgeDensityBins.end());

    const float highDensityBins[10] = {2, 5, 7.5, 9.5,  13, 16, 20, 23.5,   28,  1e6};
    FloatVector highEdgeDensityBins(highDensityBins, highDensityBins + sizeof(highDensityBins)/sizeof(highDensityBins[0]));
    m_HighEdgeDensityBins.insert(m_HighEdgeDensityBins.end(), highEdgeDensityBins.begin(), highEdgeDensityBins.end());

    for (FloatVector::iterator it = m_LowEdgeDensityBins.begin(); it != m_LowEdgeDensityBins.end(); it++)
    {
        const int position(it - m_LowEdgeDensityBins.begin());
        const float lowBin(*it);
        const float highBin(m_HighEdgeDensityBins.at(position));
        const float centralBin((lowBin + highBin) / 2);

        if (position == m_LowEdgeDensityBins.size() - 1)
        {
            m_CentreDensityBin.push_back(30);
        }
        else
        {
            m_CentreDensityBin.push_back(centralBin);
        }
    }
}

//==============================================================================

void SoftCompWeightDetermination::LoadEvents()
{
    for (std::vector<int>::const_iterator it = m_Energies.begin(); it != m_Energies.end(); ++it)
    {
        std::string energy = this->NumberToString(*it);
        std::string rootFile = m_RootFiles;
        rootFile.replace(rootFile.find("ENERGY"), std::string("ENERGY").length(), energy);

        TChain *pTChain = new TChain("HitEnergyTree");
        pTChain->Add(rootFile.c_str());

        float mcEnergy(*it);
        float pfoEnergy(-1.f);
        float rawClusterEnergy(-1.f);
        FloatVector *pHitEnergies(NULL);
        FloatVector *pCellSize0(NULL);
        FloatVector *pCellSize1(NULL);
        FloatVector *pCellThickness(NULL);
        IntVector *pHitType(NULL);

        pTChain->SetBranchAddress("EnergyOfPfo",&pfoEnergy);
        pTChain->SetBranchAddress("RawEnergyOfCluster",&rawClusterEnergy);
        pTChain->SetBranchAddress("HitEnergies",&pHitEnergies);
        pTChain->SetBranchAddress("CellSize0",&pCellSize0);
        pTChain->SetBranchAddress("CellSize1",&pCellSize1);
        pTChain->SetBranchAddress("CellThickness",&pCellThickness);
        pTChain->SetBranchAddress("HitType",&pHitType);

        if (0 == pTChain->GetEntries()) continue;

        for (unsigned int i = 0; i < pTChain->GetEntries(); i++)
        {
            pTChain->GetEntry(i);
            EventClass *pEventClass = new EventClass(mcEnergy, pfoEnergy, rawClusterEnergy, pHitEnergies, pCellSize0, pCellSize1, pCellThickness, pHitType);
            m_EventVector.push_back(pEventClass);
        }
    }
}

//==============================================================================

double SoftCompWeightDetermination::Chi2(const double *par)
{
    float chi2(0.f);
   
    for (EventVector::iterator it = m_EventVector.begin(); it != m_EventVector.end(); it++)
    {
        EventClass *pEventClass = *it;
        const float softCompEnergy = this->SoftwareCompensatedEnergy(pEventClass, par);
        const float mcEnergy = pEventClass->GetMCEnergy();
        float diff = (softCompEnergy - mcEnergy)*(softCompEnergy - mcEnergy)/(0.5*mcEnergy);
        chi2 += diff;
    }
    return chi2;
}
//==============================================================================

float SoftCompWeightDetermination::SoftwareCompensatedEnergy(EventClass *pEventClass, const double *par)
{
    float softCompEnergy(0.f);
    const float mcEnergy = pEventClass->GetMCEnergy();

    softCompEnergy += pEventClass->GetECalEnergy();

    const float p1 = par[0] + par[1]*mcEnergy + par[2]*mcEnergy*mcEnergy;
    const float p2 = par[3] + par[4]*mcEnergy + par[5]*mcEnergy*mcEnergy;
    const float p3 = par[6] / (par[7] + exp(par[8]*mcEnergy));

    FloatVector hcalHitEnergies = pEventClass->GetHCalHitEnergies();
    FloatVector hcalHitEnergyDensities = pEventClass->GetHCalHitEnergyDensities();

    for (FloatVector::iterator it = hcalHitEnergies.begin(); it != hcalHitEnergies.end(); it++) 
    {
        unsigned int position = it - hcalHitEnergies.begin();
        const float binEnergy = *it;
        const float binEnergyDensity = hcalHitEnergyDensities.at(position);
        const float weight = p1*exp(p2*binEnergyDensity) + p3;
        softCompEnergy += weight * binEnergy;
    }
    return softCompEnergy; 
}

//==============================================================================

template <class T>
std::string SoftCompWeightDetermination::NumberToString ( T Number )
{
    std::ostringstream ss;
    ss << Number;
    return ss.str();
}

//==============================================================================
