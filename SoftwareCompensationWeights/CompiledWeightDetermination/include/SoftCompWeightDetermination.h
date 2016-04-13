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
#include "TROOT.h"
#include "TStyle.h"

typedef std::vector< std::vector<float> > FloatMatrix;
typedef std::vector<float> FloatRow;
typedef std::vector< std::vector<int> > IntMatrix;
typedef std::vector<int> IntRow;

class SoftCompWeightDetermination
{
    private:
        TString     m_RootFiles;

    public:
        /*
         * Default Constructor
         */
        SoftCompWeightDetermination(TString m_RootFiles);

        /*
         * Default Destructor
         */
        ~SoftCompWeightDetermination();

        /*
         * Read data 
         */
        void LoadResults(TString rootFiles);
};

#endif
