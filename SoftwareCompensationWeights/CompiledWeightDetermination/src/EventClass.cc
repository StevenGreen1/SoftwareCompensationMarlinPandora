#include "EventClass.h"

EventClass::EventClass(float mcEnergy, float pfoEnergy, float rawClusterEnergy, FloatVector *pHitEnergies, FloatVector *pCellSize0, FloatVector *pCellSize1, FloatVector *pCellThickness, IntVector *pHitType) :
    m_MCEnergy(mcEnergy),
    m_PFOEnergy(pfoEnergy),
    m_RawClusterEnergy(rawClusterEnergy),
    m_HitEnergies(pHitEnergies),
    m_CellSize0(pCellSize0),
    m_CellSize1(pCellSize1),
    m_CellThickness(pCellThickness),
    m_HitType(pHitType)
{
    this->LoadData();
    this->CalorimeterEnergies();
}

//==============================================================================

EventClass::~EventClass()
{
}

//==============================================================================

void EventClass::LoadData()
{
    for (FloatVector::iterator it = m_HitEnergies->begin(); it != m_HitEnergies->end(); it++)
    {
        const int position = it - m_HitEnergies->begin();
        const float hitEnergy(*it);
        const float cellSize0(m_CellSize0->at(position));
        const float cellSize1(m_CellSize1->at(position));
        const float cellThickness(m_CellThickness->at(position));
        const float cellVolume(cellSize0 * cellSize1 * cellThickness);
        const float hitEnergyDensity(hitEnergy/cellVolume);
        const int hitType = m_HitType->at(position);

        if (1 == hitType)
        {
            m_ECalHitEnergies.push_back(hitEnergy);
            m_ECalHitEnergyDensities.push_back(hitEnergyDensity);
        }

        else if (2 == hitType)
        {
            m_HCalHitEnergies.push_back(hitEnergy);
            m_HCalHitEnergyDensities.push_back(hitEnergyDensity);
        }
    }
}

//==============================================================================

void EventClass::CalorimeterEnergies()
{
    float ecalEnergy(0.f), hcalEnergy(0.f);

    for (FloatVector::iterator it = m_ECalHitEnergies.begin(); it != m_ECalHitEnergies.end(); it++)
    {
        ecalEnergy += *it;
    }

    for (FloatVector::iterator it = m_HCalHitEnergies.begin(); it != m_HCalHitEnergies.end(); it++)
    {
        hcalEnergy += *it;
    }

    m_ECalEnergy = ecalEnergy;
    m_HCalEnergy = hcalEnergy;
}

