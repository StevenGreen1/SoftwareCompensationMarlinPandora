#ifndef EVENT_CLASS_H
#define EVENT_CLASS_H

#include <algorithm>
#include <fstream>
#include <iostream>
#include <list>
#include <math.h>
#include <sstream>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TChain.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TStyle.h"

typedef std::vector<float> FloatVector;
typedef std::vector<int> IntVector;

class EventClass
{
    private:
        // Input Data
        const float     m_MCEnergy;
        const float     m_PFOEnergy;
        const float     m_RawClusterEnergy;
        FloatVector    *m_HitEnergies;
        FloatVector    *m_CellSize0;
        FloatVector    *m_CellSize1;
        FloatVector    *m_CellThickness;
        IntVector      *m_HitType;

        // Derived Data
        float           m_ECalEnergy;
        FloatVector     m_ECalHitEnergies;
        FloatVector     m_ECalHitEnergyDensities;
        float           m_HCalEnergy;
        FloatVector     m_HCalHitEnergies;
        FloatVector     m_HCalHitEnergyDensities;

        /*
         * Create hit energy densities for ECal and HCal hits in same cluster
         */
        void LoadData();

        /*
         * Sum ECal and HCal energies
         */
        void CalorimeterEnergies();

    public:
        /*
         * Default Constructor
         */
        EventClass(float mcEnergy, float pfoEnergy, float rawClusterEnergy, FloatVector *pHitEnergies, FloatVector *pCellSize0, FloatVector *pCellSize1, FloatVector *pCellThickness, IntVector *pHitType);

        /*
         * Default Destructor
         */
        ~EventClass();

        inline float GetMCEnergy() {return m_MCEnergy;}
        inline float GetPFOEnergy() {return m_PFOEnergy;}
        inline float GetRawClusterEnergy() {return m_RawClusterEnergy;}
        inline float GetECalEnergy() {return m_ECalEnergy;}
        inline FloatVector GetECalHitEnergies() {return m_ECalHitEnergies;}
        inline FloatVector GetECalHitEnergyDensities() {return m_ECalHitEnergyDensities;}
        inline float GetHCalEnergy() {return m_HCalEnergy;}
        inline FloatVector GetHCalHitEnergies() {return m_HCalHitEnergies;}
        inline FloatVector GetHCalHitEnergyDensities() {return m_HCalHitEnergyDensities;}
};

#endif
