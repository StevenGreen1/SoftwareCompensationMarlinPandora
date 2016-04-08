/**
 *  @file   MarlinPandora/src/CellEnergyAlgorithm.cc
 * 
 *  @brief  Implementation of the cell energy algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "CellEnergyAlgorithm.h"

using namespace pandora;

//------------------------------------------------------------------------------------------------------------------------------------------

CellEnergyAlgorithm::CellEnergyAlgorithm() :
    m_RootFileName("Test.root")
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

CellEnergyAlgorithm::~CellEnergyAlgorithm()
{   
    std::cout << "Exiting Cell Energy Algorithm and saving root file : " << m_RootFileName << std::endl;
    PANDORA_MONITORING_API(SaveTree(this->GetPandora(), "ClusterEnergyTree", m_RootFileName, "UPDATE"));
    //PandoraMonitoringApi::SaveTree(this->GetPandora(), )
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CellEnergyAlgorithm::Run()
{
    std::cout << "===== Running Cell Energy Algorithm =====" << std::endl;

    // Algorithm code here
    const pandora::PfoList *pPfoList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pPfoList));

    std::vector<float> pfoMomentum;
    std::vector<float> pfoCosTheta;

    for (pandora::PfoList::const_iterator pfoIter = pPfoList->begin(), pfoIterEnd = pPfoList->end(); pfoIter != pfoIterEnd; ++pfoIter)
    {
        const pandora::Pfo *const pPfo = *pfoIter;
        const CartesianVector cartesianVector = pPfo->GetMomentum();
        const float px = cartesianVector.GetX();
        const float py = cartesianVector.GetY();
        const float pz = cartesianVector.GetZ();
        std::cout << px << " " << py << " " << pz << std::endl; 
        const float momentum(std::sqrt(px * px + py * py + pz * pz));
        const float cosTheta((momentum > std::numeric_limits<float>::epsilon()) ? pz / momentum : -999.f);
        pfoMomentum.push_back(momentum);
        pfoCosTheta.push_back(cosTheta);
    }

    int numberOfPfos(pPfoList->size());

    const pandora::ClusterList *pCurrentClusterList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pCurrentClusterList));

    std::vector<int> numberOfHitsInCluster;
    std::vector<float> rawEnergyOfCluster;
    std::vector<float> correctedEnergyOfCluster;
    std::vector<float> clusterDistributionMetrics;
    std::vector<int> nECalHits;
    std::vector<int> nHCalHits;
    std::vector<int> isEMShower;
    std::vector<int> nHitsAboveMean;

    for (pandora::ClusterList::const_iterator clusterIter = pCurrentClusterList->begin(), clusterIterEnd = pCurrentClusterList->end(); clusterIter != clusterIterEnd; ++clusterIter)
    {
        const Cluster *const pCluster = *clusterIter;
        const bool emShower(PandoraContentApi::GetPlugins(*this)->GetParticleId()->IsEmShower(pCluster));
        int necalHits(0);
        int nhcalHits(0);

        int emshower(0);
        if (emShower) emshower = 1;

        pandora::CaloHitList clusterCaloHitList;
        pCluster->GetOrderedCaloHitList().GetCaloHitList(clusterCaloHitList);

        this->ClusterType(clusterCaloHitList,necalHits,nhcalHits);

        int nHitsAbvMean(0);
        float clusterDistributionMetric(0.f);
        float clusterHadEnergy(pCluster->GetHadronicEnergy());
        this->ClusterDistributionMetric(clusterCaloHitList,clusterDistributionMetric,clusterHadEnergy,nHitsAbvMean);

//        if (!isEMShower)
        std::cout << "Number of hits in cluster : " << pCluster->GetNCaloHits() << std::endl;
        numberOfHitsInCluster.push_back(pCluster->GetNCaloHits());
        rawEnergyOfCluster.push_back(pCluster->GetHadronicEnergy());
        correctedEnergyOfCluster.push_back(pCluster->GetCorrectedHadronicEnergy(this->GetPandora()));
        clusterDistributionMetrics.push_back(clusterDistributionMetric);
        nECalHits.push_back(necalHits);
        nHCalHits.push_back(nhcalHits);
        isEMShower.push_back(emshower);
        nHitsAboveMean.push_back(nHitsAbvMean);

        this->Display(clusterCaloHitList);
    }

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ClusterEnergyTree", "RawEnergyOfCluster", &rawEnergyOfCluster));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ClusterEnergyTree", "CorrectedEnergyOfCluster", &correctedEnergyOfCluster));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ClusterEnergyTree", "ClusterDistributionMetrics", &clusterDistributionMetrics));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ClusterEnergyTree", "NumberOfHitsInCluster", &numberOfHitsInCluster));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ClusterEnergyTree", "NumberOfPFOsInEvent", numberOfPfos));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ClusterEnergyTree", "PFOMomentum", &pfoMomentum));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ClusterEnergyTree", "PFOCosTheta", &pfoCosTheta));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ClusterEnergyTree", "nECalHits", &nECalHits));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ClusterEnergyTree", "nHCalHits", &nHCalHits));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ClusterEnergyTree", "IsEMShower", &isEMShower));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ClusterEnergyTree", "nHitsAboveMean", &nHitsAboveMean));
    PANDORA_MONITORING_API(FillTree(this->GetPandora(), "ClusterEnergyTree"));

   return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CellEnergyAlgorithm::ClusterType(const pandora::CaloHitList &caloHitList, int &nECalHits, int &nHCalHits) const
{
    for(pandora::CaloHitList::const_iterator iter = caloHitList.begin() , endIter = caloHitList.end() ; endIter != iter ; ++iter)
    {   
        const pandora::CaloHit *pCaloHit = *iter;
        
        if (HCAL == pCaloHit->GetHitType())
            nHCalHits++;
        
        else if (ECAL == pCaloHit->GetHitType())
            nECalHits++;
    }
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CellEnergyAlgorithm::ClusterDistributionMetric(const pandora::CaloHitList &caloHitList, float &clusterDistributionMetric, float &clusterHadronicEnergy, int &nHitsAbvMean) const
{
    std::vector<float> cellEnergyDensities;
    float averageEnergyDensity(0.f);
    for(pandora::CaloHitList::const_iterator iter = caloHitList.begin() , endIter = caloHitList.end() ; endIter != iter ; ++iter)
    {
        const pandora::CaloHit *pCaloHit = *iter;
        const float cellVolume(pCaloHit->GetCellSize0() * pCaloHit->GetCellSize1() * pCaloHit->GetCellThickness());
        const float cellEnergy(pCaloHit->GetHadronicEnergy()/clusterHadronicEnergy);
        const float cellEnergyDensity(cellEnergy/cellVolume);
        cellEnergyDensities.push_back(cellEnergyDensity);
        averageEnergyDensity += cellEnergyDensity;
    }
    averageEnergyDensity /= caloHitList.size();

    float rms(0.f);
    float minCellEnergyDensity(std::numeric_limits<float>::max());
    float maxCellEnergyDensity(-std::numeric_limits<float>::max());

    for(std::vector<float>::iterator it = cellEnergyDensities.begin(); it != cellEnergyDensities.end(); ++it)
    {
        float energyDensity = *it;
        rms += pow(energyDensity-averageEnergyDensity,2);
        if (energyDensity > averageEnergyDensity) nHitsAbvMean++;
        if (energyDensity > maxCellEnergyDensity) maxCellEnergyDensity = energyDensity;
        if (energyDensity < minCellEnergyDensity) minCellEnergyDensity = energyDensity;
    }
    clusterDistributionMetric = pow(rms,0.5);
/*    float lastPre90PercentPoint(0.f);

    sort(cellEnergyDensities.begin(), cellEnergyDensities.end());

    for(std::vector<float>::iterator it = cellEnergyDensities.begin(); it != cellEnergyDensities.end(); ++it)
    {   
        float energyDensity = *it;
        if (energyDensity < minCellEnergyDensity + 0.9 * (maxCellEnergyDensity - minCellEnergyDensity) ) lastPre90PercentPoint = energyDensity;
    }

    //float energyDensityRange(maxCellEnergyDensity-minCellEnergyDensity);

    clusterDistributionMetric = (lastPre90PercentPoint - minCellEnergyDensity) / ( maxCellEnergyDensity - minCellEnergyDensity );
*/
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CellEnergyAlgorithm::Display(const CaloHitList &caloHitList) const
{
    std::string detectorView = "default";
    PandoraMonitoringApi::SetEveDisplayParameters(this->GetPandora(), true, (detectorView.find("xz") != std::string::npos) ? DETECTOR_VIEW_XZ : (detectorView.find("xy") != std::string::npos) ? DETECTOR_VIEW_XY : DETECTOR_VIEW_DEFAULT, -1.f, 0.5);

    PandoraMonitoringApi::VisualizeCaloHits(this->GetPandora(), &caloHitList, "Cluster", AUTOENERGY);
    PandoraMonitoringApi::ViewEvent(this->GetPandora());

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CellEnergyAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    // Read settings from xml file here

//    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "HotHadCellEnergy", m_HotHadCellEnergy));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "RootFileName", m_RootFileName));

    return STATUS_CODE_SUCCESS;
}
