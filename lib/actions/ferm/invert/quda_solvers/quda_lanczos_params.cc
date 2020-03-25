#include "chromabase.h"
#include "actions/ferm/invert/quda_solvers/quda_lanczos_params.h"

using namespace QDP;

namespace Chroma {
    LanczosParams::LanczosParams(XMLReader& xml, const std::string& path)
    {
        XMLReader paramtop(xml, path);
        // Lanczos Eigen Deflation Parameters
        if ( paramtop.count("EigAmin") > 0 )               { read(paramtop, "EigAmin", EigAmin); }
        if ( paramtop.count("EigPolyDeg") > 0 )            { read(paramtop, "EigPolyDeg", EigPolyDeg); }
        if ( paramtop.count("EigNEv") > 0 )                { read(paramtop, "EigNEv", EigNEv); }
        if ( paramtop.count("EigNConv") > 0 )              { read(paramtop, "EigNConv", EigNConv); }
        if ( paramtop.count("EigNKr") > 0 )                { read(paramtop, "EigNKr", EigNKr); }
        if ( paramtop.count("EigTol") > 0 )                { read(paramtop, "EigTol", EigTol); }
        if ( paramtop.count("EigBatchedRotate") > 0 )      { read(paramtop, "EigBatchedRotate", EigBatchedRotate); }
        if ( paramtop.count("EigMaxRestarts") > 0 )        { read(paramtop, "EigMaxRestarts", EigMaxRestarts); }
        if ( paramtop.count("EigRequireConvergence") > 0 ) { read(paramtop, "EigRequireConvergence", EigRequireConvergence); }
        if ( paramtop.count("EigUsePolyAcc") > 0 )         { read(paramtop, "EigUsePolyAcc", EigUsePolyAcc); }
        if ( paramtop.count("EigSpectrum") > 0 )           { read(paramtop, "EigSpectrum", EigSpectrum); }
        if ( paramtop.count("EigUseNormOp") > 0 )          { read(paramtop, "EigUseNormOp", EigUseNormOp); }
        if ( paramtop.count("EigLoadPath") > 0 )           { read(paramtop, "EigLoadPath", EigLoadPath); }
        if ( paramtop.count("EigSavePath") > 0 )           { read(paramtop, "EigSavePath", EigSavePath); }
    };

    void read(XMLReader& xml, const std::string& path, LanczosParams& p)
    {
        LanczosParams tmp(xml, path);
        p = tmp;
    }

    void write(XMLWriter& xml, const std::string& path, const LanczosParams& param)
    {
        push(xml, path);

        write(xml, "EigAmin"              , param.EigAmin);
        write(xml, "EigPolyDeg"           , param.EigPolyDeg);
        write(xml, "EigNEv"               , param.EigNEv);
        write(xml, "EigNConv"             , param.EigNConv);
        write(xml, "EigNKr"               , param.EigNKr);
        write(xml, "EigTol"               , param.EigTol);
        write(xml, "EigBatchedRotate"     , param.EigBatchedRotate);
        write(xml, "EigMaxRestarts"       , param.EigMaxRestarts);
        write(xml, "EigRequireConvergence", param.EigRequireConvergence);
        write(xml, "EigUsePolyAcc"        , param.EigUsePolyAcc);
        write(xml, "EigSpectrum"          , param.EigSpectrum);
        write(xml, "EigUseNormOp"         , param.EigUseNormOp);
        write(xml, "EigLoadPath"          , param.EigLoadPath);
        write(xml, "EigSavePath"          , param.EigSavePath);

        pop(xml);
    }
}
