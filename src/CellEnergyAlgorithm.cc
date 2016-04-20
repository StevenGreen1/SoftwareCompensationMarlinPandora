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
    std::vector<int> nECalHits;
    std::vector<int> nHCalHits;
    std::vector<int> isEMShower;
    std::vector<int> numberTrackAssociations;

    for (pandora::ClusterList::const_iterator clusterIter = pCurrentClusterList->begin(), clusterIterEnd = pCurrentClusterList->end(); clusterIter != clusterIterEnd; ++clusterIter)
    {
        const Cluster *const pCluster = *clusterIter;
        const bool emShower(PandoraContentApi::GetPlugins(*this)->GetParticleId()->IsEmShower(pCluster));
        int necalHits(0);
        int nhcalHits(0);
        const float clusterHadEnergy = pCluster->GetHadronicEnergy();
        const float clusterCorHadEnergy = pCluster->GetCorrectedHadronicEnergy(this->GetPandora());
        const pandora::TrackList &trackList = pCluster->GetAssociatedTrackList();
        const unsigned int nTrackAssociations = trackList.size();

        int emshower(0);
        if (emShower) emshower = 1;

        // PFO only has the nonIsolated calo hits associated to it.  Topologically this is right, but the energy of the
        // isolated hits is used in PFO creation.
        const pandora::OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());
        // This converts the OrderedCaloHitList into a standard CaloHitList
        pandora::CaloHitList nonIsolatedCaloHitList;
        orderedCaloHitList.GetCaloHitList(nonIsolatedCaloHitList);
        // Isolated calo hit list associated to the cluster
        const pandora::CaloHitList &isolatedCaloHitList(pCluster->GetIsolatedCaloHitList());
        pandora::CaloHitList clusterCaloHitList;
        clusterCaloHitList.insert(nonIsolatedCaloHitList.begin(), nonIsolatedCaloHitList.end());
        clusterCaloHitList.insert(isolatedCaloHitList.begin(), isolatedCaloHitList.end());

        this->ClusterType(clusterCaloHitList,necalHits,nhcalHits);

        const bool isSoftwareCompOn(std::fabs(clusterCorHadEnergy - clusterHadEnergy) > std::numeric_limits<float>::min());

        if ( isSoftwareCompOn )
        {
            std::cout << "NTracks : " << nTrackAssociations << std::endl;
            std::cout << "Software compensation changing energy : " << isSoftwareCompOn << std::endl;
            std::cout << "Number of hits in cluster   : " << pCluster->GetNCaloHits() << std::endl;
            std::cout << "Hadronic Energy of Cluster  : " << pCluster->GetHadronicEnergy() << std::endl;
            std::cout << "Corrected Energy of Cluster : " << clusterCorHadEnergy << std::endl;
            std::cout << "NECalHits                   : " << necalHits << std::endl;
            std::cout << "NHCalHits                   : " << nhcalHits << std::endl;

            std::string detectorView = "default";
            PandoraMonitoringApi::SetEveDisplayParameters(this->GetPandora(), true, (detectorView.find("xz") != std::string::npos) ? DETECTOR_VIEW_XZ : (detectorView.find("xy") != std::string::npos) ? DETECTOR_VIEW_XY : DETECTOR_VIEW_DEFAULT, -1.f, clusterHadEnergy);
            PandoraMonitoringApi::VisualizeCaloHits(this->GetPandora(), &clusterCaloHitList, "SoftCompCaloHits", SOFTCOMPWEIGHT);
            PandoraMonitoringApi::VisualizeCaloHits(this->GetPandora(), &clusterCaloHitList, "SoftCompCaloHits", CLUSTERHADE);
            PandoraMonitoringApi::VisualizeCaloHits(this->GetPandora(), &clusterCaloHitList, "SoftCompCaloHits", CELLENERGYDENSITY);
            PandoraMonitoringApi::VisualizeParticleFlowObjects(this->GetPandora(), pPfoList, "AllPFOs", AUTO, true, true);
            PandoraMonitoringApi::ViewEvent(this->GetPandora());
        }

        numberTrackAssociations.push_back(nTrackAssociations);
        numberOfHitsInCluster.push_back(pCluster->GetNCaloHits());
        rawEnergyOfCluster.push_back(pCluster->GetHadronicEnergy());
        correctedEnergyOfCluster.push_back(clusterCorHadEnergy);
        nECalHits.push_back(necalHits);
        nHCalHits.push_back(nhcalHits);
        isEMShower.push_back(emshower);
    }

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ClusterEnergyTree", "RawEnergyOfCluster", &rawEnergyOfCluster));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ClusterEnergyTree", "CorrectedEnergyOfCluster", &correctedEnergyOfCluster));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ClusterEnergyTree", "NumberOfHitsInCluster", &numberOfHitsInCluster));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ClusterEnergyTree", "NumberOfPFOsInEvent", numberOfPfos));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ClusterEnergyTree", "PFOMomentum", &pfoMomentum));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ClusterEnergyTree", "PFOCosTheta", &pfoCosTheta));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ClusterEnergyTree", "nECalHits", &nECalHits));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ClusterEnergyTree", "nHCalHits", &nHCalHits));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ClusterEnergyTree", "IsEMShower", &isEMShower));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ClusterEnergyTree", "NumberTrackAssociations", &numberTrackAssociations));
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

StatusCode CellEnergyAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    // Read settings from xml file here

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "RootFileName", m_RootFileName));

    return STATUS_CODE_SUCCESS;
}
