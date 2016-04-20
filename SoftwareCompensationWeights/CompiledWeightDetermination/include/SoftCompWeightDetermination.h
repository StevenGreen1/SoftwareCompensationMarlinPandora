#ifndef SOFT_COMP_WEIGHT_DETERMINATION_H
#define SOFT_COMP_WEIGHT_DETERMINATION_H

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
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TROOT.h"
#include "TStyle.h"

#include "EventClass.h"

typedef std::vector<EventClass*> EventVector;

class SoftCompWeightDetermination
{
    private:
        std::string     m_RootFiles;
        EventVector     m_EventVector;
        IntVector       m_Energies;
        FloatVector     m_LowEdgeDensityBins;
        FloatVector     m_HighEdgeDensityBins;
        FloatVector     m_CentreDensityBin;

    public:
        /*
         * Default Constructor
         */
        SoftCompWeightDetermination(std::string rootFiles);

        /*
         * Default Destructor
         */
        ~SoftCompWeightDetermination();

        /*
         * Read data 
         */
        void LoadEvents();

        /*
         * Calculate mean bin positions for software compensation weight binning 
         */
        void MakeBinDensities();

        /*
         * Perform TMinuit fit to find software compensation weights
         */
        void Fit();

        /*
         * Work out chi squared for fit
         */
        double Chi2(const double *par);

        /*
         * Return the software compensated energy for the event given using the given parameters
         */
        float SoftwareCompensatedEnergy(EventClass *pEventClass, const double *par);

        template <class T>
        std::string NumberToString(T Number);
};

#endif
