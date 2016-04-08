/**
 *  @file   MarlinPandora/src/PfoCreator.cc
 * 
 *  @brief  Implementation of the pfo creator class.
 * 
 *  $Log: $
 */

#include "marlin/Global.h"
#include "marlin/Processor.h"

#include "EVENT/LCCollection.h"

#include "IMPL/ClusterImpl.h"
#include "IMPL/LCCollectionVec.h"
#include "IMPL/LCFlagImpl.h"
#include "IMPL/LCGenericObjectImpl.h"
#include "IMPL/LCRelationImpl.h"
#include "IMPL/ReconstructedParticleImpl.h"

#include "CalorimeterHitType.h"
#include "ClusterShapes.h"

#include "Api/PandoraApi.h"

#include "Objects/Cluster.h"
#include "Objects/ParticleFlowObject.h"
#include "Objects/Track.h"

#include "PandoraPFANewProcessor.h"
#include "PfoCreator.h"

#include "CaloHitCreator.h"

#include <cmath>

PfoCreator::PfoCreator(const Settings &settings, const pandora::Pandora *const pPandora) :
    m_settings(settings),
    m_pPandora(pPandora)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

PfoCreator::~PfoCreator()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

float PfoCreator::FindHitDensityBin(float hitEnergy, float cellsize){
  cellsize /= 100;
  //std::cout << "cellsize1 " << cellsize << std::endl;
  
  float volume = cellsize*cellsize*0.05;
  float mip2gev = 0.0225;
  
  float hitEnergyInMIP = hitEnergy/mip2gev;

  const int NBIN = 10;
  //float lowMIP[NBIN]  = {0.49,  2, 3.5, 5.5,  8, 10, 14, 17, 21, 25, 30, 35};
  //float highMIP[NBIN] = {  2, 3.5, 5.5,   8, 10, 14, 17, 21, 25, 30, 35, 1e6};
  float lowMIP[NBIN]  = {0.3,  2, 5.5,  8, 10, 14, 17, 21, 25, 30};
  float highMIP[NBIN] = {  2, 5.5,   8, 10, 14, 17, 21, 25, 30, 1e6};

  float rho = 0;
  for (int ibin = 0; ibin < NBIN; ibin++){

    if (hitEnergyInMIP>=lowMIP[ibin] && hitEnergyInMIP<highMIP[ibin]){
      rho = (lowMIP[ibin]+highMIP[ibin])/2;
      if (ibin==(NBIN-1))
	rho = 40;
    
      rho *= mip2gev;
      rho /= volume;
    }
  }

  return rho;
}

float PfoCreator::SCEnergyCorrection(const pandora::ParticleFlowObject *const pPfo)
{
  std::cout << "PfoCreator::SCEnergyCorrection being called" << std::endl;
  float PFOenergyEstimation = pPfo->GetEnergy();

  float E_SC(0.f);

  float p10 = m_settings.m_SCparameters.at(0);
  float p11 = m_settings.m_SCparameters.at(1);
  float p12 = m_settings.m_SCparameters.at(2);

  float p20 = m_settings.m_SCparameters.at(3);
  float p21 = m_settings.m_SCparameters.at(4);
  float p22 = m_settings.m_SCparameters.at(5);

  float p30 = m_settings.m_SCparameters.at(6);
  float p31 = m_settings.m_SCparameters.at(7);
  float p32 = m_settings.m_SCparameters.at(8);

  float p1 = p10 + p11*PFOenergyEstimation + p12*PFOenergyEstimation*PFOenergyEstimation;
  float p2 = p20 + p21*PFOenergyEstimation + p22*PFOenergyEstimation*PFOenergyEstimation;
  float p3 = p30/(p31 + exp(p32*PFOenergyEstimation));

  const pandora::ClusterList &clusterList(pPfo->GetClusterList());
  
  for (pandora::ClusterList::const_iterator cIter = clusterList.begin(), cIterEnd = clusterList.end(); cIter != cIterEnd; ++cIter)
    {
      const pandora::Cluster *pPandoraCluster = *cIter;

      pandora::CaloHitList pandoraCaloHitList;
      pPandoraCluster->GetOrderedCaloHitList().GetCaloHitList(pandoraCaloHitList);

      for (pandora::CaloHitList::const_iterator hIter = pandoraCaloHitList.begin(), hIterEnd = pandoraCaloHitList.end(); hIter != hIterEnd; ++hIter)
	{

	  const pandora::CaloHit *pPandoraCaloHit= *hIter;
	  EVENT::CalorimeterHit *pCalorimeterHit = (EVENT::CalorimeterHit*)(pPandoraCaloHit->GetParentCaloHitAddress());
	  
	  const float hitEnergy(pCalorimeterHit->getEnergy());
	  const CHT cht(pCalorimeterHit->getType());
	  
	  const float cellsize = pPandoraCaloHit->GetCellSize1();
	  
	  if (cht.is(CHT::hcal)) {
	    float rho = FindHitDensityBin(hitEnergy,cellsize);
	    float weight = p1*exp(p2*rho)+p3;
	    E_SC += hitEnergy*weight;	  
	  } else {
	    E_SC += hitEnergy;
	  }
	}//end of loop on calo hit list

    //const EVENT::CalorimeterHitVec &calorimeterHitVec((*iter)->getCalorimeterHits());	  
    //for (EVENT::CalorimeterHitVec::const_iterator hitIter = calorimeterHitVec.begin(), hitIterEnd = calorimeterHitVec.end(); hitIter != hitIterEnd; ++hitIter) {
    //const EVENT::CalorimeterHit *pCalorimeterHit = *hitIter;
	    
    //const float hitEnergy(pCalorimeterHit->getEnergy());
    //const CHT cht(pCalorimeterHit->getType());

    //if (cht.is(CHT::hcal)) {
    //float rho = FindHitDensityBin(hitEnergy);
    //float weight = p1*exp(p2*rho)+p3;
    //E_SC += hitEnergy*weight;	  
    //} else {
    //E_SC += hitEnergy;
    //}
    //}
  }

  return E_SC;
}

pandora::StatusCode PfoCreator::CreateParticleFlowObjects(EVENT::LCEvent *pLCEvent)
{
    const pandora::PfoList *pPfoList = NULL;
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::GetCurrentPfoList(*m_pPandora, pPfoList));

    IMPL::LCCollectionVec *pClusterCollection = new IMPL::LCCollectionVec(LCIO::CLUSTER);
    IMPL::LCCollectionVec *pReconstructedParticleCollection = new IMPL::LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);

    IMPL::LCFlagImpl lcFlagImpl(pClusterCollection->getFlag());
    lcFlagImpl.setBit(LCIO::CLBIT_HITS);
    pClusterCollection->setFlag(lcFlagImpl.getFlag());

    std::vector<std::string> subDetectorNames ;
    subDetectorNames.push_back("ecal") ; const unsigned int ecal_Index(0) ;
    subDetectorNames.push_back("hcal") ; const unsigned int hcal_Index(1) ;
    subDetectorNames.push_back("yoke") ; const unsigned int yoke_Index(2) ;
    subDetectorNames.push_back("lcal") ; const unsigned int lcal_Index(3) ;
    subDetectorNames.push_back("lhcal"); const unsigned int lhcal_Index(4);
    subDetectorNames.push_back("bcal") ; const unsigned int bcal_Index(5) ;

    pClusterCollection->parameters().setValues("ClusterSubdetectorNames", subDetectorNames);

    int countPfo = 0;

    // Create lcio "reconstructed particles" from the pandora "particle flow objects"
    for (pandora::PfoList::const_iterator itPFO = pPfoList->begin(), itPFOEnd = pPfoList->end(); itPFO != itPFOEnd; ++itPFO)
    {
        IMPL::ReconstructedParticleImpl *pReconstructedParticle = new ReconstructedParticleImpl();

        const pandora::ClusterAddressList clusterAddressList((*itPFO)->GetClusterAddressList());
        const pandora::TrackAddressList trackAddressList((*itPFO)->GetTrackAddressList());
	countPfo++;

        // Create lcio clusters

	int count = 0;

        for (pandora::ClusterAddressList::const_iterator itCluster = clusterAddressList.begin(), itClusterEnd = clusterAddressList.end();
            itCluster != itClusterEnd; ++itCluster)
        {
            IMPL::ClusterImpl *pCluster = new ClusterImpl();

            const unsigned int nHitsInCluster((*itCluster).size());
	    
	    count++;

            float clusterEnergy(0.);
            float *pHitE = new float[nHitsInCluster];
            float *pHitX = new float[nHitsInCluster];
            float *pHitY = new float[nHitsInCluster];
            float *pHitZ = new float[nHitsInCluster];

            for (unsigned int iHit = 0; iHit < nHitsInCluster; ++iHit)
            {
                EVENT::CalorimeterHit *pCalorimeterHit = (CalorimeterHit*)((*itCluster)[iHit]);
                pCluster->addHit(pCalorimeterHit, 1.0);

                const float caloHitEnergy(pCalorimeterHit->getEnergy());
                clusterEnergy += caloHitEnergy;

                pHitE[iHit] = caloHitEnergy;
                pHitX[iHit] = pCalorimeterHit->getPosition()[0];
                pHitY[iHit] = pCalorimeterHit->getPosition()[1];
                pHitZ[iHit] = pCalorimeterHit->getPosition()[2];

                std::vector<float> &subDetectorEnergies = pCluster->subdetectorEnergies();
                subDetectorEnergies.resize(subDetectorNames.size());

                switch (CHT(pCalorimeterHit->getType()).caloID())
                {
                    case CHT::ecal:  subDetectorEnergies[ecal_Index ] += caloHitEnergy; break;
                    case CHT::hcal:  subDetectorEnergies[hcal_Index ] += caloHitEnergy; break;
                    case CHT::yoke:  subDetectorEnergies[yoke_Index ] += caloHitEnergy; break;
                    case CHT::lcal:  subDetectorEnergies[lcal_Index ] += caloHitEnergy; break;
                    case CHT::lhcal: subDetectorEnergies[lhcal_Index] += caloHitEnergy; break;
                    case CHT::bcal:  subDetectorEnergies[bcal_Index ] += caloHitEnergy; break;
                    default: streamlog_out(DEBUG) << " no subdetector found for hit with type: " << pCalorimeterHit->getType() << std::endl;
                }
            }

            pCluster->setEnergy(clusterEnergy);

            ClusterShapes *pClusterShapes = new ClusterShapes(nHitsInCluster, pHitE, pHitX, pHitY, pHitZ);
            pCluster->setPosition(pClusterShapes->getCentreOfGravity());
            pCluster->setIPhi(std::atan2(pClusterShapes->getEigenVecInertia()[1], pClusterShapes->getEigenVecInertia()[0]));
            pCluster->setITheta(std::acos(pClusterShapes->getEigenVecInertia()[2]));

            pClusterCollection->addElement(pCluster);
            pReconstructedParticle->addCluster(pCluster);

            delete pClusterShapes;
            delete[] pHitE; delete[] pHitX; delete[] pHitY; delete[] pHitZ;
        }

        // Add tracks to the lcio reconstructed particles
        for (pandora::TrackAddressList::const_iterator itTrack = trackAddressList.begin(), itTrackEnd = trackAddressList.end();
            itTrack != itTrackEnd; ++itTrack)
        {
            pReconstructedParticle->addTrack((Track*)(*itTrack));
        }

        float momentum[3] = {(*itPFO)->GetMomentum().GetX(), (*itPFO)->GetMomentum().GetY(), (*itPFO)->GetMomentum().GetZ()};
        pReconstructedParticle->setType((*itPFO)->GetParticleId());
        pReconstructedParticle->setMomentum(momentum);
        pReconstructedParticle->setEnergy((*itPFO)->GetEnergy());
	if ( m_settings.m_applySoftwareCompensation && (pReconstructedParticle->getTracks().empty()) && (pReconstructedParticle->getType()!=22) )
	  pReconstructedParticle->setEnergy(SCEnergyCorrection(*itPFO));
        pReconstructedParticle->setMass((*itPFO)->GetMass());
        pReconstructedParticle->setCharge((*itPFO)->GetCharge());

        pReconstructedParticleCollection->addElement(pReconstructedParticle);
    }

    pLCEvent->addCollection(pClusterCollection, m_settings.m_clusterCollectionName.c_str());
    pLCEvent->addCollection(pReconstructedParticleCollection, m_settings.m_pfoCollectionName.c_str());
    //pLCEvent->add

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

PfoCreator::Settings::Settings()
{

}
