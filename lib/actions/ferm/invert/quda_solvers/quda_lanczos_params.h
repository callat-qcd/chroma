#ifndef __QUDA_LANCZOS_PARAMS_H__
#define __QUDA_LANCZOS_PARAMS_H__

#include "chromabase.h"
#include "actions/ferm/invert/quda_solvers/enum_quda_io.h"


namespace Chroma 
{
    struct LanczosParams 
    {
        Real EigAmin;
        int  EigPolyDeg;
        int  EigNEv;
        int  EigNConv;
        int  EigNKr;
        Real EigTol;
        int  EigBatchedRotate;
        int  EigMaxRestarts;
        bool EigRequireConvergence;
        bool EigUsePolyAcc;
        QudaEigSpectrumType EigSpectrum;
        bool EigUseNormOp;
        /*maybe need these
          QudaSolveType
          QudaSolutionType
        */
        std::string EigLoadPath;
        std::string EigSavePath;

        LanczosParams(XMLReader& in, const std::string& path);
        LanczosParams() {
            EigAmin               = 0;     //0.1;
            EigPolyDeg            = 0;     //100;
            EigNEv                = 0;     //16;
            EigNConv              = 0;     //16;
            EigNKr                = 0;     //32;
            EigTol                = 0.;    //1e-6;
            EigBatchedRotate      = 0;     //10;
            EigMaxRestarts        = 0;     //100;
            EigRequireConvergence = false; //true;
            EigUsePolyAcc         = false; //true;
            EigSpectrum           = SR;
            EigUseNormOp          = false; //true;
            EigLoadPath           = "";
            EigSavePath           = "";
        };

    };

    // Reader
    void read(XMLReader& xml, const std::string& path, LanczosParams& param);
    // Writer
    void write(XMLWriter& xml, const std::string& path, const LanczosParams& param);


}
#endif
