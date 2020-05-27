// -*- C++ -*-
/*! \file
 *  \brief Interface to QUDA for MDWF solves - only supports const b5 and c5 for now

    May 2020: Andre Walker-Loud and Dean Howarth 
              with help/insight from Balint Joo and Kate Clark
              have enhanced the interface to support:
              1a) Having QUDA do reconstruct of solution
               b) Enabling support for all four ODD/EVEN and ASYM/SYMM type solutions
              2a) Enabling use of the Lanczos deflation in QUDA to accelerate solves
 */
#include <iomanip>


#ifndef __syssolver_linop_quda_nef_h__
#define __syssolver_linop_quda_nef_h__

#include "chroma_config.h"

#ifdef BUILD_QUDA
#include <quda.h>

#include "handle.h"
#include "state.h"
#include "syssolver.h"
#include "linearop.h"
#include "actions/ferm/fermbcs/simple_fermbc.h"
#include "actions/ferm/fermstates/periodic_fermstate.h"
#include "actions/ferm/invert/quda_solvers/syssolver_quda_nef_params.h"
#include "actions/ferm/linop/eoprec_nef_linop_array_w.h"
#include "meas/gfix/temporal_gauge.h"
#include <string>

#include "util/gauge/reunit.h"

//#include <util_quda.h>
#ifdef QDP_IS_QDPJIT
#include "actions/ferm/invert/quda_solvers/qdpjit_memory_wrapper.h"
#endif

namespace Chroma
{

    //! Richardson system solver namespace
    namespace LinOpSysSolverQUDANEFEnv
    {
        //! Register the syssolver
        bool registerAll();
    }



    //! Solve a Clover Fermion System using the QUDA inverter
    /*! \ingroup invert
     * WARNING THIS SOLVER WORKS FOR MOEBIUS/DWF ONLY ***
     */

    class LinOpSysSolverQUDANEF : public LinOpSystemSolverArray<LatticeFermion>
    {
    public:
        typedef LatticeFermion T;
        typedef LatticeColorMatrix U;
        typedef multi1d<LatticeColorMatrix> Q;

        typedef LatticeFermionF TF;
        typedef LatticeColorMatrixF UF;
        typedef multi1d<LatticeColorMatrixF> QF;

        typedef LatticeFermionF TD;
        typedef LatticeColorMatrixF UD;
        typedef multi1d<LatticeColorMatrixF> QD;

        size_t start_site;// Use this to handle MAT or MATPC type solves
        bool   quda_returns_mat;

        typedef WordType<T>::Type_t REALT;
        //! Constructor
        /*!
         * \param M_        Linear operator ( Read )
         * \param invParam  inverter parameters ( Read )
         */
        LinOpSysSolverQUDANEF(Handle< LinearOperatorArray<T> > A_,
                              Handle< FermState<T,Q,Q> > state_,
                              const SysSolverQUDANEFParams& invParam_) :
            A(A_), invParam(invParam_)
        {
            START_CODE();
            // Declare sub_domain which is either all for MAT or rb[1] for MATPC
            const Subset& sub_domain = A->subset();

            QDPIO::cout << "LinOpSysSolverQUDANEF:" << std::endl;

            // 0) Determine if user is using a PREC or UNPREC solver
            //    and if asking QUDA to return a full MAT or MATPC solution
            //    and what type of solve we want QUDA to do

            // move these initializations to the top
            // 2) pull 'new; GAUGE and Invert params
            q_gauge_param  = newQudaGaugeParam();
            quda_inv_param = newQudaInvertParam();
            quda_eig_param = newQudaEigParam();

            // SOLUTION and SOLVER TYPE
            switch( invParam.MatSolutionType ) {
            case MATPC:
                //only symmetric DWF supported at the moment:
                if ( invParam.MatPCType == EVEN_EVEN || invParam.MatPCType == EVEN_EVEN_ASYM ){
                    QDPIO::cerr << "EVEN_EVEN and EVEN_EVEN_ASYM solutions are not supported by Chroma "
                                << "if Chroma is doing the reconstruction (MATPC solution type)" << std::endl;
                    QDP_abort(1);
                }
                QDPIO::cout << "we got QUDA_MATPC_SOLUTION" << std::endl;
                start_site = rb[1].start();
                quda_inv_param.solution_type = QUDA_MATPC_SOLUTION;
                quda_returns_mat = false;
                break;
            case MAT:
                QDPIO::cout << "we got QUDA_MAT_SOLUTION" << std::endl;
                start_site = 0;
                quda_inv_param.solution_type = QUDA_MAT_SOLUTION;
                quda_returns_mat = true;
                break;
            default:// we verified the MAT solution agrees - so default will become the faster one
                if ( sub_domain.numSiteTable() == Layout::sitesOnNode() ) {
                    QDPIO::cout << "we got QUDA_MAT_SOLUTION" << std::endl;
                    start_site = 0;
                    quda_inv_param.solution_type = QUDA_MAT_SOLUTION;
                    quda_returns_mat = true;
                }
                // else - relapse to QUDA returning a MATPC solution and chroma reconstructs
                else {
                    QDPIO::cout << "we got QUDA_MATPC_SOLUTION" << std::endl;
                    start_site = rb[1].start();
                    quda_inv_param.solution_type = QUDA_MATPC_SOLUTION;
                    quda_returns_mat = false;
                }
                break;
            }

            // SOLVER
            switch( invParam.MatPCType ) {
            case ODD_ODD_ASYM:
                QDPIO::cout << "Using ODD_ODD_Asymmetric Linop: A_oo - D_oe A^{-1}_ee D_eo" << std::endl;
                quda_inv_param.matpc_type  =  QUDA_MATPC_ODD_ODD_ASYMMETRIC;
                break;
            case ODD_ODD:
                QDPIO::cout << "Using ODD_ODD_Symmetric Linop: 1 - A_oo ^{-1} Doe Aee^{-1} Deo" << std::endl;
                quda_inv_param.matpc_type  =  QUDA_MATPC_ODD_ODD;
                break;
            case EVEN_EVEN_ASYM:
                QDPIO::cout << "Using EVEN_EVEN_Asymmetric Linop: A_ee - D_eo A^{-1}_oo D_oe" << std::endl;
                quda_inv_param.matpc_type  =  QUDA_MATPC_EVEN_EVEN_ASYMMETRIC;
                break;
            case EVEN_EVEN:
                QDPIO::cout << "Using EVEN_EVEN_Symmetric Linop: 1 - A_ee ^{-1} Deo Aoo^{-1} Doe" << std::endl;
                quda_inv_param.matpc_type  =  QUDA_MATPC_EVEN_EVEN;
                break;
            default:
                // The new interface supports MAT solutions for which E_E_S is probably the fastest
                if ( quda_returns_mat ) {
                    QDPIO::cout << "Using EVEN_EVEN_Symmetric Linop: 1 - A_ee ^{-1} Deo Aoo^{-1} Doe" << std::endl;
                    quda_inv_param.matpc_type  =  QUDA_MATPC_EVEN_EVEN;
                }
                // The default Chroma call was with NEF which only supports O_O_A
                else {
                    QDPIO::cout << "Using ODD_ODD_Asymmetric Linop: A_oo - D_oe A^{-1}_ee D_eo" << std::endl;
                    quda_inv_param.matpc_type  =  QUDA_MATPC_ODD_ODD_ASYMMETRIC;
                }
                break;
            }
            if(quda_returns_mat){
                QDPIO::cout << "We are requesting a MAT solution from QUDA" << std::endl;
            }
            else {
                QDPIO::cout << "We are requesting a MATPC solution from QUDA" << std::endl;
            }

            // FOLLOWING INITIALIZATION in test QUDA program

            // 1) work out cpu_prec, cuda_prec, cuda_prec_sloppy
            int s = sizeof( WordType<T>::Type_t );

            if (s == 4) {
                cpu_prec = QUDA_SINGLE_PRECISION;
            }
            else {
                cpu_prec = QUDA_DOUBLE_PRECISION;
            }


            // Work out GPU precision
            switch( invParam.cudaPrecision ) {
            case HALF:
                gpu_prec = QUDA_HALF_PRECISION;
                break;
            case SINGLE:
                gpu_prec = QUDA_SINGLE_PRECISION;
                break;
            case DOUBLE:
                gpu_prec = QUDA_DOUBLE_PRECISION;
                break;
            default:
                gpu_prec = cpu_prec;
                break;
            }

            // Work out GPU Sloppy precision
            // Default: No Sloppy
            switch( invParam.cudaSloppyPrecision ) {
            case HALF:
                gpu_sloppy_prec = QUDA_HALF_PRECISION;
                break;
            case SINGLE:
                gpu_sloppy_prec = QUDA_SINGLE_PRECISION;
                break;
            case DOUBLE:
                gpu_sloppy_prec = QUDA_DOUBLE_PRECISION;
                break;
            default:
                gpu_sloppy_prec = gpu_prec;
                break;
            }

            // Work out GPU Precondition Precision
            // Default: No Sloppy
            switch( invParam.cudaPreconditionPrecision ) {
            case HALF:
                gpu_precondition_prec = QUDA_HALF_PRECISION;
                break;
            case SINGLE:
                gpu_precondition_prec = QUDA_SINGLE_PRECISION;
                break;
            case DOUBLE:
                gpu_precondition_prec = QUDA_DOUBLE_PRECISION;
                break;
            default:
                gpu_precondition_prec = gpu_prec;
                break;
            }

            // 3) set lattice size
            const multi1d<int>& latdims = Layout::subgridLattSize();

            q_gauge_param.X[0] = latdims[0];
            q_gauge_param.X[1] = latdims[1];
            q_gauge_param.X[2] = latdims[2];
            q_gauge_param.X[3] = latdims[3];

            // 5) - set QUDA_WILSON_LINKS, QUDA_GAUGE_ORDER
            q_gauge_param.type = QUDA_WILSON_LINKS;
#ifndef BUILD_QUDA_DEVIFACE_GAUGE
            q_gauge_param.gauge_order = QUDA_QDP_GAUGE_ORDER; // gauge[mu], p
#else
            QDPIO::cout << "MDAGM Using QDP-JIT gauge order" << std::endl;
            q_gauge_param.location    = QUDA_CUDA_FIELD_LOCATION;
            q_gauge_param.gauge_order = QUDA_QDPJIT_GAUGE_ORDER;
#endif

            // 6) - set t_boundary
            // Convention: BC has to be applied already
            // This flag just tells QUDA that this is so,
            // so that QUDA can take care in the reconstruct
            if( invParam.AntiPeriodicT ) {
                q_gauge_param.t_boundary = QUDA_ANTI_PERIODIC_T;
            }
            else {
                q_gauge_param.t_boundary = QUDA_PERIODIC_T;
            }

            // Set cpu_prec, cuda_prec, reconstruct and sloppy versions
            q_gauge_param.cpu_prec  = cpu_prec;
            q_gauge_param.cuda_prec = gpu_prec;


            switch( invParam.cudaReconstruct ) {
            case RECONS_NONE:
                q_gauge_param.reconstruct = QUDA_RECONSTRUCT_NO;
                break;
            case RECONS_8:
                q_gauge_param.reconstruct = QUDA_RECONSTRUCT_8;
                break;
            case RECONS_12:
                q_gauge_param.reconstruct = QUDA_RECONSTRUCT_12;
                break;
            default:
                q_gauge_param.reconstruct = QUDA_RECONSTRUCT_12;
                break;
            };

            q_gauge_param.cuda_prec_sloppy = gpu_sloppy_prec;

            switch( invParam.cudaSloppyReconstruct ) {
            case RECONS_NONE:
                q_gauge_param.reconstruct_sloppy = QUDA_RECONSTRUCT_NO;
                break;
            case RECONS_8:
                q_gauge_param.reconstruct_sloppy = QUDA_RECONSTRUCT_8;
                break;
            case RECONS_12:
                q_gauge_param.reconstruct_sloppy = QUDA_RECONSTRUCT_12;
                break;
            default:
                q_gauge_param.reconstruct_sloppy = QUDA_RECONSTRUCT_12;
                break;
            };
            // Precondition Reconstruct
            q_gauge_param.cuda_prec_precondition = gpu_precondition_prec;

            switch( invParam.cudaPreconditionReconstruct ) {
            case RECONS_NONE:
                q_gauge_param.reconstruct_precondition = QUDA_RECONSTRUCT_NO;
                break;
            case RECONS_8:
                q_gauge_param.reconstruct_precondition = QUDA_RECONSTRUCT_8;
                break;
            case RECONS_12:
                q_gauge_param.reconstruct_precondition = QUDA_RECONSTRUCT_12;
                break;
            default:
                q_gauge_param.reconstruct_precondition = QUDA_RECONSTRUCT_12;
                break;
            };

            // Gauge fixing:

            // These are the links
            // They may be smeared and the BC's may be applied
            Q links_single(Nd);

            // Now downcast to single prec fields.
            for(int mu=0; mu < Nd; mu++) {
                links_single[mu] = (state_->getLinks())[mu];
            }

            // GaugeFix
            if( invParam.axialGaugeP ) {
                QDPIO::cout << "Fixing Temporal Gauge" << std::endl;
                temporalGauge(links_single, GFixMat, Nd-1);
                for(int mu=0; mu < Nd; mu++){
                    links_single[mu] = GFixMat*(state_->getLinks())[mu]*adj(shift(GFixMat, FORWARD, mu));
                }
                q_gauge_param.gauge_fix = QUDA_GAUGE_FIXED_YES;
            }
            else {
                // No GaugeFix
                q_gauge_param.gauge_fix = QUDA_GAUGE_FIXED_NO;  // No Gfix yet
            }

            // Don't support anisotorpy for Moebius
            q_gauge_param.anisotropy = 1.0;

            // Now onto the inv param:
            // Dslash type
            quda_inv_param.dslash_type = QUDA_MOBIUS_DWF_DSLASH;
            // Invert type:
            switch( invParam.solverType ) {
            case CG:
                quda_inv_param.inv_type = QUDA_CG_INVERTER;
                solver_string = "CG";
                break;
            case BICGSTAB:
                QDPIO::cerr << "Solver BICGSTAB not supported for MDWF" << std::endl;
                QDP_abort(1);
                break;
            case GCR:
                QDPIO::cerr << "Solver GCR not supported for MDWF" << std::endl;
                QDP_abort(1);
                break;
            default:
                QDPIO::cerr << "Unknown Solver type" << std::endl;
                QDP_abort(1);
                break;
            }

            // Mass
            quda_inv_param.mass = toDouble(invParam.NEFParams.Mass);
            // The sign convention of M5 for Chroma and QUDA are opposite
            quda_inv_param.m5=toDouble(-invParam.NEFParams.OverMass);
            quda_inv_param.Ls=invParam.NEFParams.N5;
            // M.A.C. made these static so no need to alloc them.
            if ( invParam.NEFParams.N5 >= QUDA_MAX_DWF_LS ) {
                QDPIO::cerr << "LS can be at most " << QUDA_MAX_DWF_LS << std::endl;
                QDP_abort(1);
            }

            // Copy b5 and c5 into static array
            QDPIO::cout << "Ls from matrix: " << A->size() << std::endl;
            QDPIO::cout << "Ls from params: " << invParam.NEFParams.N5 << std::endl;
            QDPIO::cout << "Ls from quda: " << quda_inv_param.Ls << std::endl;

            // Alias the b_5 complex array as an array
            // of doubles
            double* b5_arr = reinterpret_cast<double*>(quda_inv_param.b_5);
            double* c5_arr = reinterpret_cast<double*>(quda_inv_param.c_5);

            // Set the params
            for(int s=0; s < quda_inv_param.Ls; ++s) {
                // Real
                b5_arr[2*s] = toDouble(invParam.NEFParams.b5[s]);
                // Imag
                b5_arr[2*s+1] = 0;

                // Real
                c5_arr[2*s] = toDouble(invParam.NEFParams.c5[s]);
                // Imag
                c5_arr[2*s+1] = 0;
            }

            struct cmpx {
                double re;
                double im;
            };

            for(int s=0; s < quda_inv_param.Ls; ++s) {
                QDPIO::cout << "CHROMA_NEF_PARAM: b5[" <<s<<"] = " << toDouble(invParam.NEFParams.b5[s]) << "     c5[" << s << "] = " << toDouble(invParam.NEFParams.c5[s]);

                const cmpx& qb5 = reinterpret_cast<const cmpx&>(quda_inv_param.b_5[s]);
                const cmpx& qc5 = reinterpret_cast<const cmpx&>(quda_inv_param.c_5[s]);

                QDPIO::cout << "    QUDA_NEF_PARAM:   b5[" <<s<<"] = (" << qb5.re << ", " << qb5.im << ")"
                            << "  c5[" << s << "] = (" << qc5.re << ", " << qc5.im << ")" << std::endl;
            }

            quda_inv_param.tol            = toDouble(invParam.RsdTarget);
            quda_inv_param.maxiter        = invParam.MaxIter;
            quda_inv_param.reliable_delta = toDouble(invParam.Delta);

            // Solve type, only CG-types supported so far:
            switch( invParam.solverType ) {
            case CG:
                if( invParam.cgnrP ) {
                    QDPIO::cout << "Doing CGNR solve" << std::endl;
                    quda_inv_param.solve_type = QUDA_NORMOP_PC_SOLVE;
                }
                else {
                    QDPIO::cout << "Doing CGNE solve" << std::endl;
                    quda_inv_param.solve_type = QUDA_NORMERR_PC_SOLVE;
                }

                break;

            default:
                if( invParam.cgnrP ) {
                    QDPIO::cout << "Doing CGNR solve" << std::endl;
                    quda_inv_param.solve_type = QUDA_NORMOP_PC_SOLVE;
                }
                else {
                    QDPIO::cout << "Doing CGNE solve" << std::endl;
                    quda_inv_param.solve_type = QUDA_NORMERR_PC_SOLVE;
                }
                break;
            }

            quda_inv_param.dagger                 = QUDA_DAG_NO;
            if ( invParam.massNorm ) {
                quda_inv_param.mass_normalization = QUDA_MASS_NORMALIZATION;
                QDPIO::cout << "Using MASS Normalization" << std::endl;
            }
            else {
                quda_inv_param.mass_normalization = QUDA_KAPPA_NORMALIZATION;
                QDPIO::cout << "Using KAPPA Normalization" << std::endl;
            }
            quda_inv_param.cpu_prec               = cpu_prec;
            quda_inv_param.cuda_prec              = gpu_prec;
            quda_inv_param.cuda_prec_sloppy       = gpu_sloppy_prec;
            quda_inv_param.cuda_prec_precondition = gpu_precondition_prec;
            quda_inv_param.preserve_source        = QUDA_PRESERVE_SOURCE_NO;
            quda_inv_param.gamma_basis            = QUDA_DEGRAND_ROSSI_GAMMA_BASIS;

#ifndef BUILD_QUDA_DEVIFACE_SPINOR
            quda_inv_param.dirac_order            = QUDA_DIRAC_ORDER;
#else
            QDPIO::cout << "MDAGM Using QDP-JIT spinor order" << std::endl;
            quda_inv_param.dirac_order            = QUDA_QDPJIT_DIRAC_ORDER;
            quda_inv_param.input_location         = QUDA_CUDA_FIELD_LOCATION;
            quda_inv_param.output_location        = QUDA_CUDA_FIELD_LOCATION;
#endif

            // Autotuning
            if( invParam.tuneDslashP ) {
                QDPIO::cout << "Enabling Dslash Autotuning" << std::endl;
                quda_inv_param.tune = QUDA_TUNE_YES;
            }
            else {
                QDPIO::cout << "Disabling Dslash Autotuning" << std::endl;
                quda_inv_param.tune = QUDA_TUNE_NO;
            }

            // Setup padding
            multi1d<int> face_size(4);
            face_size[0] = latdims[1]*latdims[2]*latdims[3]/2;
            face_size[1] = latdims[0]*latdims[2]*latdims[3]/2;
            face_size[2] = latdims[0]*latdims[1]*latdims[3]/2;
            face_size[3] = latdims[0]*latdims[1]*latdims[2]/2;

            int max_face = face_size[0];
            for(int i=1; i < 4; i++) {
                if ( face_size[i] > max_face ) {
                    max_face = face_size[i];
                }
            }

            q_gauge_param.ga_pad = max_face;
            // PADDING
            quda_inv_param.sp_pad = 0;
            quda_inv_param.cl_pad = 0;

            QDPIO::cout << "Setting Precondition stuff to defaults for not using" << std::endl;
            quda_inv_param.inv_type_precondition  = QUDA_INVALID_INVERTER;
            quda_inv_param.tol_precondition       = 1.0e-1;
            quda_inv_param.maxiter_precondition   = 1000;
            quda_inv_param.verbosity_precondition = QUDA_SILENT;
            quda_inv_param.gcrNkrylov             = 1;

            if ( invParam.InvDeflate ) {
                // Lanczos Eigenvector Info
                QDPIO::cout << "Setting Lanczos Eigenvector Parameters" << std::endl;
                QDPIO::cout << "Telling QUDA to preserve deflation" << std::endl;
                quda_eig_param.preserve_deflation = QUDA_BOOLEAN_TRUE;

                switch( invParam.NEFLanczosParams.EigSpectrum ) {
                case SR:
                    quda_eig_param.spectrum            = QUDA_SPECTRUM_SR_EIG;
                    break;
                default:
                    QDPIO::cerr << "ONLY SR EigSpectrum type accepted at this time" << std::endl;
                }
                quda_eig_param.nConv               = invParam.NEFLanczosParams.EigNConv;
                quda_eig_param.nEv                 = invParam.NEFLanczosParams.EigNEv;
                quda_eig_param.nKr                 = invParam.NEFLanczosParams.EigNKr;
                quda_eig_param.tol                 = toDouble(invParam.NEFLanczosParams.EigTol);
                quda_eig_param.batched_rotate      = invParam.NEFLanczosParams.EigBatchedRotate;
                quda_eig_param.require_convergence = invParam.NEFLanczosParams.EigRequireConvergence == true
                                                     ? QUDA_BOOLEAN_TRUE : QUDA_BOOLEAN_FALSE;
                quda_eig_param.max_restarts        = invParam.NEFLanczosParams.EigMaxRestarts;
                quda_eig_param.cuda_prec_ritz      = gpu_precondition_prec;
                quda_eig_param.use_norm_op         = QUDA_BOOLEAN_TRUE;
                quda_eig_param.use_dagger          = QUDA_BOOLEAN_FALSE;
                quda_eig_param.compute_svd         = QUDA_BOOLEAN_FALSE;
                quda_eig_param.use_poly_acc        = invParam.NEFLanczosParams.EigUsePolyAcc == true
                                                     ? QUDA_BOOLEAN_TRUE : QUDA_BOOLEAN_FALSE;
                quda_eig_param.poly_deg            = invParam.NEFLanczosParams.EigPolyDeg;
                quda_eig_param.a_min               = toDouble(invParam.NEFLanczosParams.EigAmin);
                quda_eig_param.a_max               = toDouble(-1.0); // Instruct QUDA to estimate
                quda_eig_param.arpack_check        = QUDA_BOOLEAN_FALSE;
                //strcpy(quda_eig_param.vec_infile,  (char*)(invParam.NEFLanczosParams.EigLoadPath).c_str() );
                //strcpy(quda_eig_param.vec_outfile, (char*)(invParam.NEFLanczosParams.EigSavePath).c_str() );
                strcpy(quda_eig_param.vec_infile,  (invParam.NEFLanczosParams.EigLoadPath).c_str() );
                strcpy(quda_eig_param.vec_outfile, (invParam.NEFLanczosParams.EigSavePath).c_str() );

                quda_inv_param.eig_param           = &quda_eig_param;
            }
            else {
                quda_inv_param.eig_param           = nullptr;
            }

            if( invParam.verboseP ) {
                quda_inv_param.verbosity = QUDA_VERBOSE;
            }
            else {
                quda_inv_param.verbosity = QUDA_SUMMARIZE;
            }

            // Set up the links
            void* gauge[4];

            for(int mu=0; mu < Nd; mu++) {
#ifndef BUILD_QUDA_DEVIFACE_GAUGE
                gauge[mu] = (void *)&(links_single[mu].elem(all.start()).elem().elem(0,0).real());
#else
                gauge[mu] = GetMemoryPtr( links_single[mu].getId() );
                std::cout << "MDAGM CUDA gauge[" << mu << "] in = " << gauge[mu] << "\n";
#endif
            }

            loadGaugeQuda((void *)gauge, &q_gauge_param);

            END_CODE();
        }

        //! Destructor is automatic
        ~LinOpSysSolverQUDANEF()
        {
            QDPIO::cout << "Destructing" << std::endl;
            freeGaugeQuda();
        }

        //! Return the subset on which the operator acts
        const Subset& subset() const {return A->subset();}

        //return the size:
        int size() const {return invParam.NEFParams.N5;}

        //! Solver the linear system
        /*!
         * \param psi      solution ( Modify )
         * \param chi      source ( Read )
         * \return syssolver results
         */
        SystemSolverResults_t operator() (multi1d<T>& psi, const multi1d<T>& chi) const
        {
            SystemSolverResults_t res;

            START_CODE();
            StopWatch swatch;
            swatch.start();

            // Declare sub_domain which is either all for MAT or rb[1] for MATPC
            const Subset& sub_domain = A->subset();
            /*  If the user has asked QUDA to return a full solution, but used an UNPREF solver (NEF), then
                we should abort immediately to not waste anyone's time or compute time
             */
            if ( quda_returns_mat && sub_domain.numSiteTable() != Layout::sitesOnNode() ) {
                QDPIO::cerr << "You have requested a full MAT solution from QUDA, but you are using "
                            << "a preconditioned solver from Chroma.  In order for Chroma to "
                            << "accept a full MAT solution, you must pick an unpreconditioned Chroma "
                            << "solver for memory and solution checking" << std::endl;
                QDPIO::cerr << "NEF ---> UNPRECONDITIONED_NEF" << std::endl;
                QDP_abort(1);
            }
            if ( ! quda_returns_mat && sub_domain.numSiteTable() == Layout::sitesOnNode() ) {
                QDPIO::cerr << "You have requested a MATPC solution from QUDA, but you are using "
                            << "an UNPRECONDITIONED solver from Chroma.  You most likely did not intend "
                            << "to do this.  Add support for this option to the code if you really want it.  "
                            << "Otherwise, change" << std::endl;
                QDPIO::cerr << "UNPRECONDITIONED_NEF --> NEF" << std::endl;
                QDP_abort(1);
            }
            
            // 1/( 2 kappa_b ) =  b5 * ( Nd - M5 ) + 1
            Real invTwoKappaB = invParam.NEFParams.b5[0]*( Nd - invParam.NEFParams.OverMass) + Real(1);
            Real twoKappaB = Real(1)/invTwoKappaB;
            // QUDA defines kappa_b without b5.  We need to in the test code to assign a kappa_b value
            // 1/( 2 kappa_b ) =       ( Nd - M5 ) +1  NOTE: M5_QUDA = -M5_Chroma
            Double invTwoKappaBQuda = Nd - invParam.NEFParams.OverMass + 1.;
            Double twoKappaBQuda    = 1./invTwoKappaBQuda;
            
            const multi1d<int>& latdims = Layout::subgridLattSize();
            int fermsize = 2 * Nc * Ns * latdims[0]*latdims[1]*latdims[2]*latdims[3];
            int half_fermsize;
            if ( ! quda_returns_mat){
                fermsize = fermsize / 2.;
                half_fermsize = fermsize;
            }
            else {
                half_fermsize = fermsize / 2;
            }

#if 0
            {
                // Set Kappa since QUDA MatQuda used in test code does not call massRescale to set this parameter
                quda_inv_param.kappa = 0.5 * toDouble(twoKappaBQuda);
                bool VERBOSE_TEST=true;
                bool POINT_SOURCE=false;
                
                QDPIO::cout << "TEST the Application of MatChroma and MatQUDA to the source" << std::endl;

                // Set operator as MAT or MATPC based on size of sub_domain
                if ( sub_domain.numSiteTable() == Layout::sitesOnNode() ) {
                    QDPIO::cout << "Using NORMOP / NORMERR for comparison" << std::endl;
                    if( invParam.cgnrP ) { quda_inv_param.solve_type = QUDA_NORMOP_SOLVE;}
                    else { quda_inv_param.solve_type = QUDA_NORMERR_SOLVE;}
                }
                else{
                    QDPIO::cout << "Using NORMOP_PC / NORMERR_PC for comparison" << std::endl;
                    // no need to set this as it is the solve_type requested
                }

                // In_Quda is input to QUDA, Out_Quda is result
                multi1d<T> in_quda( this->size() );
                multi1d<T> out_quda(this->size() );            
                // In_Chroma is input to Chroma, Out_Chroma is result
                multi1d<T> in_chroma( this->size() );
                multi1d<T> out_chroma(this->size() );
            
                for(int s=0; s < this->size(); s++ ) {
                    // POINT SOURCE
                    if ( POINT_SOURCE ) {
                        in_quda[s] = zero;
                        if ( sub_domain.numSiteTable() == Layout::sitesOnNode() ) {
                            in_quda[s].elem(even.start()).elem(0).elem(0).real() = 1.;
                        }
                        else {
                            in_quda[s] = chi[s];
                        }
                    } // GAUSSIAN SOURCE
                    else {
                        gaussian(in_quda[s]);  // Gaussian into in_quda
                    }
                    in_chroma[s] = in_quda[s];   // copy to in_chroma
                }
                // used for memcpy to transform data
                int start_cb=0;
                int num_cb=2;
                int dest_cb;
                if ( A->subset().start() == rb[1].start() ) {
                    // A is on ODD Checker Board, so only 1 parity and start at 1
                    start_cb = 1;
                    num_cb = 1;
                }

                for(int d=0; d < 2; d++) {
                    for(int s=0; s < this->size(); s++ ) {
                        out_quda[s]   = zero;   // zero both out_quda and out_chroma
                        out_chroma[s] = zero;
                    }
                    // Apply Chroma A
                    if ( d==0 ) {
                        QDPIO::cout << "Doing Mat" << std::endl;
                        (*A)(out_chroma, in_chroma, PLUS);
                    }
                    else {
                        QDPIO::cout << "Doing MatDag" << std::endl;
                        (*A)(out_chroma, in_chroma, MINUS);
                    }                    
                    
                    REAL* spinorIn  = new REAL[quda_inv_param.Ls*fermsize];
                    REAL* spinorOut = new REAL[quda_inv_param.Ls*fermsize];
                    memset((spinorIn),  0, fermsize*quda_inv_param.Ls*sizeof(REAL));
                    memset((spinorOut), 0, fermsize*quda_inv_param.Ls*sizeof(REAL));

                    // Now compare out_quda and out_chroma
                        
                    // test the source to make sure the parity-re-assignment was handled properly
                    QDPIO::cout << "compare nomr2 of operators BEFORE handing to QUDA" << std::endl;
                    for(int s=0; s < this->size();s++) {
                        QDPIO::cout << "QUDA[" <<s << ", rb[0]] = " << norm2(in_quda[s], rb[0])
                                    << "  QUDA[" <<s << ", rb[1]] = " << norm2(in_quda[s], rb[1])
                                    << "  Chroma[" <<s << ", rb[0]] = " << norm2(in_chroma[s], rb[0]) 
                                    << "  Chroma[" <<s << ". rb[1]] = " << norm2(in_chroma[s], rb[1])
                                    << std::endl;
                    }
                    // See .cc file for explanation of this mapping
                    QDPIO::cout << "Source parity map QUDA <-- Chroma" << std::endl;
                    dest_cb=0;
                    for ( int cb=start_cb; cb <= num_cb; cb++) {
                        if ( cb < 2 ) {// we only want cb=0 and 1 for PREC solution
                            for ( int s=0; s < quda_inv_param.Ls; s++ ) {
                                QDPIO::cout << "    QUDA[p=" << dest_cb <<", s=" << s << ", i_sp= " << s+ quda_inv_param.Ls*cb << "]"
                                            << "   <--   "
                                            << "Chroma[s=" << s << "][p=" << cb << "]"
                                            << std::endl;
                                memcpy(reinterpret_cast<char*>(&spinorIn[ half_fermsize*(s + quda_inv_param.Ls*(dest_cb))]),
                                       reinterpret_cast<char*>(const_cast<REAL*>(&(in_quda[s].elem(rb[cb].start()).elem(0).elem(0).real()))),
                                       half_fermsize*sizeof(REAL));
                            }
                        }
                        dest_cb++;
                    }
                    // Apply QUDA
                    if( d==0 ) {
                        quda_inv_param.dagger = QUDA_DAG_NO;
                        MatQuda((void *)spinorOut, (void *)spinorIn, (QudaInvertParam*)&quda_inv_param);
                    }
                    else {
                        quda_inv_param.dagger = QUDA_DAG_YES;
                        MatQuda((void *)spinorOut, (void *)spinorIn, (QudaInvertParam*)&quda_inv_param);
                    }

                    // see .cc file for explanation of parity mapping
                    QDPIO::cout << "Solution parity map Chroma <-- QUDA" << std::endl;
                    //if ( d == 0 ) {
                    dest_cb = start_cb;
                    for( int cb=0; cb < num_cb; cb++) {
                        for( int s=0; s < quda_inv_param.Ls; s++){
                            QDPIO::cout << "    Chroma[s=" << s << "][p=" << dest_cb << "]"
                                        << "   <--   "
                                        << "QUDA[p=" << dest_cb <<", s=" << s << ", i_sp= " << s+ quda_inv_param.Ls*cb << "]"
                                        << std::endl;
                            memcpy(reinterpret_cast<char*>(const_cast<REAL*>(&(out_quda[s].elem(rb[dest_cb].start()).elem(0).elem(0).real()))),
                                   reinterpret_cast<char*>(&spinorOut[half_fermsize*(s + quda_inv_param.Ls*cb)]),
                                   half_fermsize*sizeof(REAL));
                        }
                        dest_cb++;
                    }
                    
                    // declare these in case we do a verbose test
                    multi1d<int> coord_0(4);
                    int Nt = QDP::Layout::lattSize()[3];
                    int mid_t = (int)Nt/2;
                    coord_0(0) = 0; coord_0(1) = 0; coord_0(2) = 0; coord_0(3) = 0;
                    multi1d<int> node_0 = Layout::nodeCoord(coord_0);
                    if ( VERBOSE_TEST ) {
                        QDPIO::cout << "compare nomr2 of operators and color=spin=0 of a few time slices" << std::endl;
                    }
                    // Now compare out_quda (QUDA) and out_chroma (Chroma)
                    // QUDA MatQuda does NOT use massRescale which sets kappa for us
                    //      therefore, we have to 
                    for(int s=0; s < this->size();s++) {
                        // QUDA to Chroma normalization
                        out_quda[s] *= invTwoKappaB;
                        // MASS normalization
                        if ( invParam.massNorm ) {
                            if ( ! quda_returns_mat ) {
                                out_quda[s] *= twoKappaBQuda * twoKappaBQuda;
                            }
                            else {
                                out_quda[s] *= twoKappaBQuda;
                            }
                        }

                        QDPIO::cout <<   "QUDA[" <<s << ", rb[0]] = "   << norm2(out_quda[s], rb[0])
                                    << "  QUDA[" <<s << ", rb[1]] = "   << norm2(out_quda[s], rb[1])
                                    << "  Chroma[" <<s << ", rb[0]] = " << norm2(out_chroma[s], rb[0])
                                    << "  Chroma[" <<s << ", rb[1]] = " << norm2(out_chroma[s], rb[1])
                                    << std::endl;
                        QDPIO::cout << "QUDA[" <<s << "] = " << norm2(out_quda[s])
                                    << "  Chroma[" <<s << "] = " << norm2(out_chroma[s])
                                    << "  QUDA/Chroma[" <<s << "] = " << norm2(out_quda[s]) / norm2(out_chroma[s])
                                    <<std::endl;

                        if ( VERBOSE_TEST ) {
                            QDPIO::cout << "L5 X Y Z T s c "<<std::endl;
                            int linearInd;
                            if (Layout::nodeCoord() == node_0) {
                                for (int xi=0; xi<QDP::Layout::lattSize()[0]; xi++){
                                    coord_0(0) = xi;
                                    for (int yi=0; yi<QDP::Layout::lattSize()[1]; yi++) {
                                        coord_0(1) = yi;
                                        for (int zi=0; zi<QDP::Layout::lattSize()[2]; zi++) {
                                            coord_0(2) = zi;
                                            for (int ti=0; ti<Nt; ti++){
                                                coord_0(3) = ti;
                                                linearInd = Layout::linearSiteIndex(coord_0);
                                                for (int si=0; si<Ns; si++){
                                                    for (int ci=0; ci<Nc; ci++){
                                                        double quda_re   = toDouble(QDP::real(out_quda[s].elem(linearInd).elem(si).elem(ci)));
                                                        double quda_im   = toDouble(QDP::imag(out_quda[s].elem(linearInd).elem(si).elem(ci)));
                                                        double chroma_re = toDouble(QDP::real(out_chroma[s].elem(linearInd).elem(si).elem(ci)));
                                                        double chroma_im = toDouble(QDP::imag(out_chroma[s].elem(linearInd).elem(si).elem(ci)));
                                                        if ( quda_re > 1.e-6 || chroma_re > 1.e-6){
                                                            std::cout << s << "  " << xi << " " << yi << " " << zi << " " << ti << " " << si << " " << ci
                                                                      <<"  QUDA = " << std::setprecision(5) << quda_re << " + "
                                                                      << std::setprecision(5)<< quda_im << " I "
                                                                      << "  Chroma = "<< std::setprecision(5)<< chroma_re << " + "
                                                                      << std::setprecision(5)<< chroma_im << " I "
                                                                      << std::endl;
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        } // end VERBOSE_TEST
                    }

                    delete [] spinorIn;
                    delete [] spinorOut;
                }// D D^dagger loop
                //QDPIO::cout << "Exiting after one solve" << std::endl;
                //QDP_abort(1);
                // Reset quda_inv_param params
                quda_inv_param.dagger = QUDA_DAG_NO;
                if ( sub_domain.numSiteTable() == Layout::sitesOnNode() ) {
                    QDPIO::cout << "Restoring NORMOP_PC / NORMERR_PC for solve" << std::endl;
                    if( invParam.cgnrP ) { quda_inv_param.solve_type = QUDA_NORMOP_PC_SOLVE;}
                    else { quda_inv_param.solve_type = QUDA_NORMERR_PC_SOLVE;}
                }// if NORMOP_PC is used in test, no need to restore as that is what we want
            }
#endif

            QDPIO::cout << "Norm of chi = " << norm2(chi,sub_domain) << std::endl;
            QDPIO::cout << "TwoKappa = " << twoKappaB << "   invTwoKappa_b = " << invTwoKappaB << std::endl;
            
            if ( invParam.axialGaugeP ) {
                multi1d<T> g_chi( this->size());
                multi1d<T> g_psi( this->size());

                // Gauge Fix source and initial guess
                QDPIO::cout << "Gauge Fixing source and initial guess" << std::endl;
                for(unsigned int s=0; s<invParam.NEFParams.N5; s++){
                    g_chi[s][ sub_domain ]  = GFixMat * chi[s];
                    g_psi[s][ sub_domain ]  = GFixMat * psi[s];
                }
                QDPIO::cout << "Solving" << std::endl;
                res = qudaInvert(g_chi,g_psi);
                for(int s=0; s < this->size(); s++) {
                    g_psi[s][sub_domain] *= twoKappaB;
                }

                QDPIO::cout << "Untransforming solution." << std::endl;
                for(unsigned int s=0; s<invParam.NEFParams.N5; s++){
                    psi[s][sub_domain]  = adj(GFixMat)*g_psi[s];
                }

            }
            else {
                QDPIO::cout << "Calling qudaInvert" << std::endl;
                res = qudaInvert(chi, psi);
                // convert the normalization from QUDA back to Chroma
                for(int s=0; s < this->size(); s++) {
                    // overall normalization relevant for KAPPA and MASS normalization
                    psi[s][sub_domain] *= twoKappaB;

                    // MASS normalization
                    if ( invParam.massNorm) {
                        // MATPC solve
                        if ( ! quda_returns_mat ) {
                            psi[s][sub_domain] *= invTwoKappaBQuda * invTwoKappaBQuda;
                        }
                        // MAT solve
                        else {
                            psi[s][sub_domain] *= invTwoKappaBQuda;
                        }
                    }
                }
            }

            swatch.stop();
            double time = swatch.getTimeInSeconds();

            // Check Solution
            if ( invParam.checkSolution ) {
                multi1d<T> r(A->size());
                multi1d<T> Ax(A->size());
                r=zero;
                Ax=zero;
                Double r_norm(zero);
                Double b_norm(zero);
                Double r_norm_s(zero);
                Double b_norm_s(zero);
                (*A)(Ax, psi, PLUS);
                if ( quda_returns_mat ){
                    QDPIO::cout << "===============================================================" << std::endl;
                    QDPIO::cout << "NOTE: with MAT solution, sources are ZERO in 5th dimension bulk" << std::endl;
                }
                for(int s=0; s < A->size(); s++) {
                    r[s][sub_domain] = chi[s] - Ax[s];
                    r_norm_s = norm2(r[s],   sub_domain);
                    b_norm_s = norm2(chi[s], sub_domain);
                    r_norm += r_norm_s;
                    b_norm += b_norm_s;
                    if( quda_returns_mat ){
                        if( s == 0 ) {
                            QDPIO::cout << " | r[" << s <<"] | = " << sqrt(r_norm_s)
                                        << " | b[" << s <<"] | = " << sqrt(b_norm_s)
                                        << " |r|/|b|[" <<s<< "] = " << sqrt(r_norm_s)/sqrt(b_norm_s)
                                        << std::endl;
                        }
                        else {
                            QDPIO::cout << " | r[" << s <<"] | = " << sqrt(r_norm_s)
                                        << " | b[" << s <<"] | = " << sqrt(b_norm_s)
                                        << std::endl;                            
                        }
                    }
                    else{
                        QDPIO::cout << " | r[" << s <<"] | = " << sqrt(r_norm_s)
                                    << " | b[" << s <<"] | = " << sqrt(b_norm_s)
                                    << " |r|/|b|[" <<s<< "] = " << sqrt(r_norm_s)/sqrt(b_norm_s)
                                    << std::endl;
                    }
                }

                Double resid = sqrt(r_norm);
                Double rel_resid = sqrt(r_norm/b_norm);

                res.resid = resid;
                QDPIO::cout << "QUDA_"<< solver_string <<"_NEF_SOLVER: "
                            << res.n_count << " iterations. Max. Rsd = " << res.resid
                            << " Max. Relative Rsd = " << rel_resid << std::endl;

                // Convergence Check/Blow Up
                if ( ! invParam.SilentFailP ) {
                    if (  toBool( rel_resid >  invParam.RsdToleranceFactor*invParam.RsdTarget) )
                        {
                            QDPIO::cerr << "ERROR: QUDA Solver residuum is outside tolerance: QUDA resid=" << rel_resid
                                        << " Desired =" << invParam.RsdTarget
                                        << " Max Tolerated = " << invParam.RsdToleranceFactor*invParam.RsdTarget << std::endl;
                            QDP_abort(1);
                        }
                }
            }
            else {
                QDPIO::cout << "QUDA_"<< solver_string <<"_NEF_SOLVER: " << res.n_count << " iterations" << std::endl;
            }
            END_CODE();
            return res;
        }

    private:
        // Hide default constructor
        LinOpSysSolverQUDANEF() {}

#if 1
        Q links_orig;
#endif

        U GFixMat;
        QudaPrecision_s cpu_prec;
        QudaPrecision_s gpu_prec;
        QudaPrecision_s gpu_sloppy_prec;
        QudaPrecision_s gpu_precondition_prec;

        Handle< LinearOperatorArray<T> > A;
        const SysSolverQUDANEFParams invParam;
        QudaGaugeParam q_gauge_param;
        mutable QudaInvertParam quda_inv_param;
        mutable QudaEigParam quda_eig_param;

        SystemSolverResults_t qudaInvert(const multi1d<T>& chi_s, multi1d<T>& psi_s)const;

        std::string solver_string;
    };


} // End namespace

#endif // BUILD_QUDA
#endif
