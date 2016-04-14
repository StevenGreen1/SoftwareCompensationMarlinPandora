/**
 *  @file   MarlinPandora/include/PfoCreator.h
 * 
 *  @brief  Header file for the pfo creator class.
 * 
 *  $Log: $
 */

#ifndef PFO_CREATOR_H
#define PFO_CREATOR_H 1

#include "EVENT/LCEvent.h"
#include "EVENT/ReconstructedParticle.h"

#include "Api/PandoraApi.h"

/**
 *  @brief  PfoCreator class
 */

class PfoCreator
{
public:
    /**
     *  @brief  Settings class
     */
    class Settings
    {
    public:
        /**
         *  @brief  Default constructor
         */
        Settings();

        std::string     m_clusterCollectionName;                ///< The name of the cluster output collection
        std::string     m_pfoCollectionName;                    ///< The name of the pfo output collection
	bool            m_applySoftwareCompensation;            ///< The flag used to apply software compensation

	typedef std::vector<float> FloatVector;
	FloatVector     m_SCparameters;                         ///< Parameters used to determine SC weights
    };

    /**
     *  @brief  Constructor
     * 
     *  @param  settings the creator settings
     *  @param  pPandora address of the relevant pandora instance
     */
     PfoCreator(const Settings &settings, const pandora::Pandora *const pPandora);

    /**
     *  @brief  Destructor
     */
     ~PfoCreator();

    /**
     *  @brief  Create particle flow objects
     * 
     *  @param  pLCEvent the lcio event
     */    
    pandora::StatusCode CreateParticleFlowObjects(EVENT::LCEvent *pLCEvent);

    /**
     *  @breif  Find energy desnity of CaloHit in mips per cell
     *
     *  @param  pCaloHit calo hit to find energy density of
     *  @param  energyDensity energy density to set
     */ 

    pandora::StatusCode FindDensity(const pandora::CaloHit *const pCaloHit, float &energyDensity) const;

    /**
     *  @breif  Apply software compensation to PFO
     *
     *  @param  pPfo pfo to apply software compensation to
     *  @param  pfoEnergyEstimator software compensation pfo energy estimator
     */ 
    pandora::StatusCode SCEnergyCorrection(const pandora::ParticleFlowObject *const pPfo, float &pfoEnergyEstimator) const;

private:
    const Settings          m_settings;                         ///< The pfo creator settings
    const pandora::Pandora *m_pPandora;                         ///< Address of the pandora object from which to extract the pfos
};

#endif // #ifndef PFO_CREATOR_H
