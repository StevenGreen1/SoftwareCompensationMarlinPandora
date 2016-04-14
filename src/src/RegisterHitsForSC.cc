/**
 *  @file   MarlinPandora/src/RegisterHitsForSC.cc
 * 
 *  @brief  Implementation of the register hits for SC algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"
#include "PandoraMonitoringApi.h"

#include "RegisterHitsForSC.h"

using namespace pandora;

//------------------------------------------------------------------------------------------------------------------------------------------

RegisterHitsForSC::RegisterHitsForSC() :
    m_myRootFileName("HitsForSC.root")
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

RegisterHitsForSC::~RegisterHitsForSC()
{   
    std::cout << "Exiting RegisterHitsForSC algorithm and saving root file : " << m_myRootFileName << std::endl;
    PANDORA_MONITORING_API(SaveTree(this->GetPandora(), "HitEnergyTree", m_myRootFileName, "UPDATE"));
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode RegisterHitsForSC::Run()
{
    const pandora::PfoList *pPfoList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pPfoList));
    const int numberOfPfos(pPfoList->size());

    if (numberOfPfos != 1)
    {
        return STATUS_CODE_SUCCESS;
    }

    const pandora::Pfo *const pPfo = *pPfoList->begin();
    const float pfoEnergy(pPfo->GetEnergy());
    const pandora::ClusterList *pClusterList = &pPfo->GetClusterList();
    const int numberOfClusters(pClusterList->size());

    if (numberOfClusters != 1)
    {
        return STATUS_CODE_SUCCESS;
    }

    const Cluster *const pCluster = *pClusterList->begin();

    const pandora::OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());
    pandora::CaloHitList nonIsolatedCaloHitList;
    orderedCaloHitList.GetCaloHitList(nonIsolatedCaloHitList);
    const pandora::CaloHitList &isolatedCaloHitList(pCluster->GetIsolatedCaloHitList());
    pandora::CaloHitList clusterCaloHitList;
    clusterCaloHitList.insert(nonIsolatedCaloHitList.begin(), nonIsolatedCaloHitList.end());
    clusterCaloHitList.insert(isolatedCaloHitList.begin(), isolatedCaloHitList.end());

    const float rawEnergyOfCluster(pCluster->GetHadronicEnergy());

    std::vector<float> cellSize0;// = new std::vector<float>;
    std::vector<float> cellSize1;// = new std::vector<float>;
    std::vector<float> cellThickness;// = new std::vector<float>;
    std::vector<float> hitEnergies;// = new std::vector<float>;
    std::vector<int> hitType;// = new std::vector<int>;

    for(pandora::CaloHitList::const_iterator hitIter = clusterCaloHitList.begin() , endhitIter = clusterCaloHitList.end() ; endhitIter != hitIter ; ++hitIter)
    {
        const pandora::CaloHit *pCaloHit = *hitIter;
        const float cellSize0ToAdd(pCaloHit->GetCellSize0());
        const float cellSize1ToAdd(pCaloHit->GetCellSize1());
        const float cellThicknessToAdd(pCaloHit->GetCellThickness());
        const float cellHadronicEnergy(pCaloHit->GetHadronicEnergy());

        cellSize0.push_back(cellSize0ToAdd);
        cellSize1.push_back(cellSize1ToAdd);
        cellThickness.push_back(cellThicknessToAdd);
        hitEnergies.push_back(cellHadronicEnergy);

        if (HCAL == pCaloHit->GetHitType())
        {
            hitType.push_back(2);
        }

        else if (ECAL == pCaloHit->GetHitType())
        {
            hitType.push_back(1);
        }

        else
        {
            hitType.push_back(3);
        }
    }

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "HitEnergyTree", "EnergyOfPfo", pfoEnergy));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "HitEnergyTree", "RawEnergyOfCluster", rawEnergyOfCluster));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "HitEnergyTree", "HitEnergies", &hitEnergies));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "HitEnergyTree", "CellSize0", &cellSize0));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "HitEnergyTree", "CellSize1", &cellSize1));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "HitEnergyTree", "CellThickness", &cellThickness));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "HitEnergyTree", "HitType", &hitType));
    PANDORA_MONITORING_API(FillTree(this->GetPandora(), "HitEnergyTree"));

   return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode RegisterHitsForSC::ReadSettings(const TiXmlHandle xmlHandle)
{
    // Read settings from xml file here

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MyRootFileName", m_myRootFileName));

    return STATUS_CODE_SUCCESS;
}
