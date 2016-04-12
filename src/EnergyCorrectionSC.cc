/**
 *  @file   MarlinPandoraSC/src/EnergyCorrectionSC.cc
 * 
 *  @brief  Implementation of the SC energy correction function algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"
#include "EnergyCorrectionSC.h"
#include "PandoraPFANewProcessor.h"

using namespace pandora;

EnergyCorrectionSC::EnergyCorrectionSC() :
  m_clusterMinHadEnergy(10.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EnergyCorrectionSC::MakeEnergyCorrections(const pandora::Cluster *const pCluster, float &correctedHadronicEnergy) const
{
    if (NULL == pCluster)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);

    if (0 == pCluster->GetNCaloHits())
    {
        correctedHadronicEnergy = 0.f;
        return STATUS_CODE_SUCCESS;
    }

    const float clusterHadEnergy = pCluster->GetHadronicEnergy();

    if (clusterHadEnergy < m_clusterMinHadEnergy)
    {
        return STATUS_CODE_SUCCESS;
    }

//    const pandora::TrackList &trackList = pCluster->GetAssociatedTrackList();
//    const unsigned int nTrackAssociations = trackList.size();

//    if (0 == nTrackAssociations)
//    {
//        return STATUS_CODE_SUCCESS;
//    }

    //Here to check the cluster type so that the corresponding correction function is called
    bool isHCalCluster(false), isECalCluster(false);
    
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
/*
    int N_lim(0), N_av(0);
    float maxDiff(std::numeric_limits<float>::min());
    float maxHitEnergy(std::numeric_limits<float>::min());
    float minHitEnergy(std::numeric_limits<float>::max());
    const int numberCaloHits(clusterCaloHitList.size());
    const float meanHitEnergy(clusterHadEnergy/numberCaloHits);
    std::vector<float> clusterHitEnergy;

    for(pandora::CaloHitList::const_iterator iter = clusterCaloHitList.begin() , endIter = clusterCaloHitList.end() ; endIter != iter ; ++iter)
    {
        const pandora::CaloHit *pCaloHit = *iter;
        const float hitEnergy = pCaloHit->GetHadronicEnergy();
        const float hitEnergyinMIP = pCaloHit->GetMipEquivalentEnergy();
        clusterHitEnergy.push_back(hitEnergy);
   
        if (fabs(hitEnergy-meanHitEnergy)>maxDiff)
        {
            maxDiff = fabs(hitEnergy-meanHitEnergy);
        }

        if (hitEnergy>maxHitEnergy)
        {
            maxHitEnergy = hitEnergy;
        }

        if (hitEnergy<minHitEnergy)
        {
            minHitEnergy = hitEnergy;
        }

        if (hitEnergy<meanHitEnergy)
        {
            N_av++;
        }

        if (hitEnergyinMIP<5.)
        {
            N_lim++;
        }
    }
   
    float C_global(0.f);

    if (N_av==0)
    {
        return STATUS_CODE_SUCCESS;
    }
    else
    {
        C_global = (float)N_lim/N_av;
    }

    if (C_global>=1.2 || maxDiff/(maxHitEnergy-minHitEnergy)<0.8 )
    {
        return STATUS_CODE_SUCCESS;
    }

    std::cout << "Software Compensation Doing Something At Clustering" << std::endl;
*/
    this->clusterType(clusterCaloHitList, isECalCluster, isHCalCluster);

    if (isECalCluster)
    {
        std::cout << "ECal Cluster" << std::endl;
        return STATUS_CODE_SUCCESS;
    }
    else if (isHCalCluster)
    {
        std::cout << "HCal Cluster" << std::endl;
        EnergyCorrectionSC::SCClusterEnergyCorrectionFunction(clusterHadEnergy,clusterCaloHitList, correctedHadronicEnergy);
    }
    else
    {
        std::cout << "Split Cluster" << std::endl;
        EnergyCorrectionSC::SCHCalECalSplitClusterEnergyCorrectionFunction(clusterHadEnergy,clusterCaloHitList, correctedHadronicEnergy);
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EnergyCorrectionSC::FindDensity(const pandora::CaloHit *const pCaloHit, float &energyDensity) const
{
    const int NBIN = 10;
    float lowMIP[NBIN]  = {0.3,  2, 5.5,  8, 10, 14, 17, 21, 25, 30};
    float highMIP[NBIN] = {  2, 5.5,   8, 10, 14, 17, 21, 25, 30, 1e6};

    const float cellVolume = pCaloHit->GetCellSize0() * pCaloHit->GetCellSize1() * pCaloHit->GetCellThickness() / 1000000;
    const float mipEquivalentEnergy = pCaloHit->GetMipEquivalentEnergy();
    const float hitEnergyHadronic(pCaloHit->GetHadronicEnergy());

    for (int ibin = 0; ibin < NBIN; ibin++)
    {
        if (mipEquivalentEnergy >= lowMIP[ibin] && mipEquivalentEnergy < highMIP[ibin])
        {
            energyDensity = (lowMIP[ibin]+highMIP[ibin])/2;
            if (ibin==(NBIN-1))
            {
                energyDensity = 40;
            }
            const float mip2gev = hitEnergyHadronic / mipEquivalentEnergy;
            energyDensity = energyDensity * mip2gev;
            energyDensity /= cellVolume;
        }
    }
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EnergyCorrectionSC::clusterType(const pandora::CaloHitList &caloHitList, bool &isECalCluster, bool &isHCalCluster) const
{
    int nECalHits(0), nHCalHits(0);

    for(pandora::CaloHitList::const_iterator iter = caloHitList.begin() , endIter = caloHitList.end() ; endIter != iter ; ++iter)
    {
        const pandora::CaloHit *pCaloHit = *iter;

        if (HCAL == pCaloHit->GetHitType())
            nHCalHits++;

        else if (ECAL == pCaloHit->GetHitType())
            nECalHits++;
    }

    if (nECalHits != 0 && nHCalHits == 0)
        isECalCluster = true;

    else if (nHCalHits != 0 && nECalHits == 0)
        isHCalCluster = true;

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EnergyCorrectionSC::ECalClusterEnergyCorrectionFunction(const pandora::CaloHitList &caloHitList, float &energyCorrection) const
{
    //std::cout << "========ECAL CONTAINED CLUSTER========" << std::endl;

    for(pandora::CaloHitList::const_iterator iter = caloHitList.begin() , endIter = caloHitList.end() ; endIter != iter ; ++iter)
    {
        const pandora::CaloHit *pCaloHit = *iter;
        energyCorrection += pCaloHit->GetHadronicEnergy();
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EnergyCorrectionSC::SCClusterEnergyCorrectionFunction(float clusterEnergyEstimation, const pandora::CaloHitList &caloHitList, float &energyCorrection) const
{
    std::cout << "========HCAL CONTAINED CLUSTER========" << std::endl;

    float E_SC(0.f);

    float p10 = m_SCEnergyConstants1.at(0);
    float p11 = m_SCEnergyConstants1.at(1);
    float p12 = m_SCEnergyConstants1.at(2);

    float p20 = m_SCEnergyConstants1.at(3);
    float p21 = m_SCEnergyConstants1.at(4);
    float p22 = m_SCEnergyConstants1.at(5);

    float p30 = m_SCEnergyConstants1.at(6);
    float p31 = m_SCEnergyConstants1.at(7);
    float p32 = m_SCEnergyConstants1.at(8);

    if (m_cheating)
        clusterEnergyEstimation = m_trueEnergy;

    float p1 = p10 + p11*clusterEnergyEstimation + p12*clusterEnergyEstimation*clusterEnergyEstimation;
    float p2 = p20 + p21*clusterEnergyEstimation + p22*clusterEnergyEstimation*clusterEnergyEstimation;
    float p3 = p30/(p31 + exp(p32*clusterEnergyEstimation));

    for(pandora::CaloHitList::const_iterator iter = caloHitList.begin() , endIter = caloHitList.end() ; endIter != iter ; ++iter)
    {
        const pandora::CaloHit *pCaloHit = *iter;
        const float hitEnergy = pCaloHit->GetHadronicEnergy();
        float rho(0.f);
        this->FindDensity(pCaloHit,rho);
        double weight(p1*exp(p2*rho) + p3);
        float hitEcorr(hitEnergy*weight);
        E_SC += hitEcorr;
    }

    energyCorrection = E_SC;
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EnergyCorrectionSC::SCHCalECalSplitClusterEnergyCorrectionFunction(float clusterEnergyEstimation, const pandora::CaloHitList &caloHitList, float &energyCorrection) const
{
    //std::cout << "========SPLIT CLUSTER========" << std::endl;
    float E_SC(0.f);

    float p10 = m_SCEnergyConstants1.at(0);
    float p11 = m_SCEnergyConstants1.at(1);
    float p12 = m_SCEnergyConstants1.at(2);

    float p20 = m_SCEnergyConstants1.at(3);
    float p21 = m_SCEnergyConstants1.at(4);
    float p22 = m_SCEnergyConstants1.at(5);

    float p30 = m_SCEnergyConstants1.at(6);
    float p31 = m_SCEnergyConstants1.at(7);
    float p32 = m_SCEnergyConstants1.at(8);

    if (m_cheating)
        clusterEnergyEstimation = m_trueEnergy;

    float p1 = p10 + p11*clusterEnergyEstimation + p12*clusterEnergyEstimation*clusterEnergyEstimation;
    float p2 = p20 + p21*clusterEnergyEstimation + p22*clusterEnergyEstimation*clusterEnergyEstimation;
    float p3 = p30/(p31 + exp(p32*clusterEnergyEstimation));

    for(pandora::CaloHitList::const_iterator iter = caloHitList.begin() , endIter = caloHitList.end() ; endIter != iter ; ++iter)
    {
        const pandora::CaloHit *pCaloHit = *iter;

        if (HCAL == pCaloHit->GetHitType())
        {
            const float hitEnergy = pCaloHit->GetHadronicEnergy();
            float rho(0.f);
            this->FindDensity(pCaloHit,rho);
            double weight(p1*exp(p2*rho) + p3);
            float hitEcorr(hitEnergy*weight);
            E_SC += hitEcorr;
        }	
        else
        {           
	  E_SC += pCaloHit->GetHadronicEnergy();
        }
    }
    energyCorrection = E_SC;
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EnergyCorrectionSC::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
        "SCEnergyConstantsHighE", m_SCEnergyConstants1));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
        "SCEnergyConstantsSmallE", m_SCEnergyConstants2));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "Cheating", m_cheating));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "TrueEnergy", m_trueEnergy));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ClusterMinHadEnergy", m_clusterMinHadEnergy));

    return STATUS_CODE_SUCCESS;
}

