#ifndef __SYSSOLVER_QUDA_CLOVER_MDWF_PARAMS_H__
#define __SYSSOLVER_QUDA_CLOVER_MDWF_PARAMS_H__

#include "chromabase.h"
#include "io/xml_group_reader.h"

#include "actions/ferm/fermacts/nef_fermact_params_w.h"
#include "actions/ferm/invert/quda_solvers/enum_quda_io.h"
#include "actions/ferm/invert/quda_solvers/quda_gcr_params.h"
#include "actions/ferm/invert/quda_solvers/quda_lanczos_params.h"
#include <string>
#include "handle.h"

namespace Chroma 
{
    struct SysSolverQUDANEFParams { 
        SysSolverQUDANEFParams(XMLReader& xml, const std::string& path);
        SysSolverQUDANEFParams() {
            solverType                  = CG;
            cudaPrecision               = DEFAULT;
            cudaReconstruct             = RECONS_12;
            cudaSloppyPrecision         = DEFAULT;
            cudaSloppyReconstruct       = RECONS_12;
            cudaPreconditionPrecision   = DEFAULT;
            cudaPreconditionReconstruct = RECONS_12;
            RsdToleranceFactor          = Real(10); //< Tolerate if the solution achived is better (less) than rsdToleranceFactor*RsdTarget
            asymmetricP     = false; //< Use asymmetric version of the linear operator
            axialGaugeP     = false; //< Fix Axial Gauge?
            SilentFailP     = false; //< If set to true ignore lack of convergence. Default is 'loud' 
            tuneDslashP     = false; //< v0.3 autotune feature
            verboseP        = false;
            innerParamsP    = false;
            backup_invP     = false;
            dump_on_failP   = false;
            cgnrP           = false;
            InvDeflate      = false;
            checkSolution   = false;
            MatPCType       = ODD_ODD_ASYM; // change me to ODD_ODD when supported
            MatSolutionType = MATPC;        // change me to MAT when supported
        };

        SysSolverQUDANEFParams( const SysSolverQUDANEFParams& p) {
            NEFParams = p.NEFParams;
            AntiPeriodicT = p.AntiPeriodicT;
            MaxIter = p.MaxIter;
            RsdTarget = p.RsdTarget;
            Delta = p.Delta;
            solverType = p.solverType;
            verboseP = p.verboseP;
            asymmetricP = p.asymmetricP;
            cudaPrecision = p.cudaPrecision;
            cudaReconstruct = p.cudaReconstruct;
            cudaSloppyPrecision = p.cudaSloppyPrecision;
            cudaSloppyReconstruct = p.cudaSloppyReconstruct;
            cudaPreconditionPrecision=p.cudaPreconditionPrecision;
            cudaPreconditionReconstruct=p.cudaPreconditionReconstruct;
            axialGaugeP = p.axialGaugeP;
            SilentFailP = p.SilentFailP;
            RsdToleranceFactor = p.RsdToleranceFactor;
            tuneDslashP = p.tuneDslashP;
            innerParamsP = p.innerParamsP;
            innerParams = p.innerParams;
            backup_invP = p.backup_invP;
            backup_inv_param = p.backup_inv_param;
            dump_on_failP = p.dump_on_failP;
            cgnrP= p.cgnrP;
            checkSolution   = p.checkSolution;
            MatPCType       = p.MatPCType;
            MatSolutionType = p.MatSolutionType;
            // Do deflation?
            InvDeflate= p.InvDeflate;
            if ( InvDeflate ) {
                // Lanczos Params
                NEFLanczosParams = p.NEFLanczosParams;
            }
        }

        // Original NEF params
        NEFFermActParams NEFParams;
        bool checkSolution;
        bool AntiPeriodicT;
        int  MaxIter;
        Real RsdTarget;
        Real Delta;
        bool verboseP;
        bool asymmetricP;
        bool axialGaugeP;
        bool SilentFailP;
        bool tuneDslashP;
        bool innerParamsP;
        Real RsdToleranceFactor;
        QudaSolverType    solverType;
        QudaPrecisionType cudaPrecision;
        QudaReconsType    cudaReconstruct;
        QudaPrecisionType cudaSloppyPrecision;
        QudaReconsType    cudaSloppyReconstruct;
        QudaPrecisionType cudaPreconditionPrecision;
        QudaReconsType    cudaPreconditionReconstruct;

        //
        ChromaQudaMatPCType       MatPCType;
        ChromaQudaMatSolutionType MatSolutionType;

        // Deflate?
        bool InvDeflate;
        //Lanczos Params
        //Handle<LanczosParams> NEFLanczosParams;
        LanczosParams NEFLanczosParams;
	
        // GCR Specific params
        // Params for the preconditioner
        Handle<GCRInnerSolverParams> innerParams;

        // XML for Backup Solver
        bool backup_invP;
        GroupXML_t backup_inv_param;
        bool dump_on_failP;
        bool cgnrP;    

    };

    void read(XMLReader& xml, const std::string& path, SysSolverQUDANEFParams& p);

    void write(XMLWriter& xml, const std::string& path, 
               const SysSolverQUDANEFParams& param);



}

#endif


