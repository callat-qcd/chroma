/*! \file
*  \brief Solve a MdagM*psi=chi linear system by CG2
*/

#include "actions/ferm/invert/syssolver_linop_factory.h"
#include "actions/ferm/invert/syssolver_linop_aggregate.h"
#include "actions/ferm/invert/quda_solvers/syssolver_quda_nef_params.h"
#include "actions/ferm/invert/quda_solvers/syssolver_linop_nef_quda_w.h"
#include "io/aniso_io.h"


#include "handle.h"
#include "actions/ferm/fermstates/periodic_fermstate.h"
#include "actions/ferm/linop/eoprec_nef_linop_array_w.h"
#include "meas/glue/mesplq.h"
// QUDA Headers
#include <quda.h>
// #include <util_quda.h>

namespace Chroma
{
    namespace LinOpSysSolverQUDANEFEnv
    {

        //! Anonymous namespace
        namespace
        {
            //! Name to be used
            const std::string name("QUDA_NEF_INVERTER");

            //! Local registration flag
            bool registered = false;
        }

        LinOpSystemSolverArray<LatticeFermion>* 
            createFerm(XMLReader& xml_in,
                       const std::string& path,
                       Handle< FermState< LatticeFermion, multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > > state,
                       Handle< LinearOperatorArray<LatticeFermion> > A)
        {
            return new LinOpSysSolverQUDANEF(A, state,SysSolverQUDANEFParams(xml_in, path));
        }

        //! Register all the factories
        bool registerAll()
        {
            bool success = true;
            if (! registered)
            {
                success &= Chroma::TheLinOpFermSystemSolverArrayFactory::Instance().registerObject(name, createFerm);
                registered = true;
            }
            return success;
        }
    }

    SystemSolverResults_t LinOpSysSolverQUDANEF::qudaInvert(const multi1d<T>& chi_s, multi1d<T>& psi_s) const
    {
        SystemSolverResults_t ret;

        /*  In order to support Chroma asking QUDA to return a full MAT solution, we need to adjust the
            data ordering so that ALL EVEN and ODD fields are grouped together vs in Chroma, where the
            EVEN and ODD are split for each s_5 in the 5th dimension.  To do this, we must
            - transform the source to QUDA's data ordering
            - transform the solution back to Chroma's data ordering
            
            This is a temporary fix in Chroma until QUDA supports the data re-ordering.  Therefore, we will
            not worry about the memory bloat of doubling the SOURCE and SOLUTION in memory unless the QUDA
            solution takes a long time to get sorted
        */

        // Declare sub_domain which is either all for MAT or rb[1] for MATPC
        const Subset& sub_domain = A->subset();

        //size of field to copy. If we do a MAT solution, need full size, else, half
        const multi1d<int>& latdims = Layout::subgridLattSize();
        //            complex * Nc * Ns * V4
        // fermsize starts as full 4D fermion
        int fermsize = 2 * Nc * Ns * latdims[0]*latdims[1]*latdims[2]*latdims[3];
        int half_fermsize;

        int start_cb=0;
        int num_cb=2;
        if ( A->subset().start() == rb[1].start() ) {
            // A is on ODD Checker Board, so only 1 parity and start at 1
            start_cb = 1;
            num_cb = 1;
            fermsize = fermsize / 2;  // for PREC operator, fermsize is 1/2 the size
            half_fermsize = fermsize; // and half_fermsize is the same
        }
        // fermsize = full if quda_returns_mat
        // fermsize = 1/2  if MATPC solution is requested
        else {
            half_fermsize = fermsize / 2; // for UNPREC operator, fermsize is full size
        }

        // Now - define the in and out spinor
        REAL* spinorIn  = new REAL[quda_inv_param.Ls*fermsize];
        REAL* spinorOut = new REAL[quda_inv_param.Ls*fermsize];
        // NOTE - we could move spinorOut to after the solution, and delete spinorIn before assigning spinorOut to save memory
        memset(reinterpret_cast<char*>(spinorIn),  0, fermsize*quda_inv_param.Ls*sizeof(REAL));
        memset(reinterpret_cast<char*>(spinorOut), 0, fermsize*quda_inv_param.Ls*sizeof(REAL));

        //copy the source into spinorIn
#ifndef BUILD_QUDA_DEVIFACE_SPINOR
        // If num_cb=2 (UNPREC solution) we will populate both parities into the spinor
        // from Chroma's layout [s][parity][fermion] --> QUDA's layout [parity][s][fermion]
        // if PREC solution - then there is only 1 parity
        QDPIO::cout << "Source parity map QUDA <-- Chroma" << std::endl;
        int dest_cb=0;
        for ( int cb=start_cb; cb <= num_cb; cb++) {
            if ( cb < 2 ) {// we only want cb=0 and 1 for PREC solution
                for ( int s=0; s < quda_inv_param.Ls; s++ ) {
                    QDPIO::cout << "    QUDA[p=" << dest_cb <<", s=" << s << ", i_sp= " << s+ quda_inv_param.Ls*cb << "]"
                                << "   <--   "
                                << "Chroma[s=" << s << "][p=" << cb << "]"
                                << std::endl;
                    memcpy(reinterpret_cast<char*>(&spinorIn[ half_fermsize*(s + quda_inv_param.Ls*(dest_cb))]),
                           reinterpret_cast<char*>(const_cast<REAL*>(&(chi_s[s].elem(rb[cb].start()).elem(0).elem(0).real()))),
                           half_fermsize*sizeof(REAL));
                }
            }
            dest_cb++;
        }
#else
        //not yet
        //for(unsigned int s=0; s<quda_inv_param.Ls; s++){
        //	spinorIn[s]=GetMemoryPtr( chi_s[s].getId() );
        //	spinorOut[s]=GetMemoryPtr( psi_s[s].getId() );
        //}
#endif

        // Do the solve here
        StopWatch swatch1;
        swatch1.reset();
        swatch1.start();
        invertQuda(reinterpret_cast<void*>(spinorOut), reinterpret_cast<void*>(spinorIn), (QudaInvertParam*)&quda_inv_param);
        swatch1.stop();

        //copy result
        psi_s.resize(quda_inv_param.Ls);
        psi_s = zero;

#ifndef BUILD_QUDA_DEVIFACE_SPINOR
        // If UNPREC solution (num_cb=2) we have to transform the solution from QUDA layout back to Chroma layout
        // QUDA : [parity][s][fermion] --> Chroma : [s][parity][fermion]
        QDPIO::cout << "Solution parity map Chroma <-- QUDA" << std::endl;
        dest_cb = start_cb;
        for( int cb=0; cb < num_cb; cb++) {
            for( int s=0; s < quda_inv_param.Ls; s++){
                QDPIO::cout << "    Chroma[s=" << s << "][p=" << dest_cb << "]"
                            << "   <--   "
                            << "QUDA[p=" << dest_cb <<", s=" << s << ", i_sp= " << s+ quda_inv_param.Ls*cb << "]"
                            << std::endl;
                memcpy(reinterpret_cast<char*>(const_cast<REAL*>(&(psi_s[s].elem(rb[dest_cb].start()).elem(0).elem(0).real()))),
                       reinterpret_cast<char*>(&spinorOut[half_fermsize*(s + quda_inv_param.Ls*cb)]),
                       half_fermsize*sizeof(REAL));
            }
            dest_cb++;
        }
#else
        //not yet implemented
        //for(unsigned int s=0; s<quda_inv_param.Ls; s++){
        //	spinorIn[s]=GetMemoryPtr( chi_s[s].getId() );
        //	spinorOut[s]=GetMemoryPtr( psi_s[s].getId() );
        //}
#endif

        //clean up
        delete [] spinorIn;
        delete [] spinorOut;

        // report Performance and iteration count
        if ( quda_inv_param.solution_type == QUDA_MATPC_SOLUTION ){
            if ( quda_inv_param.matpc_type == QUDA_MATPC_ODD_ODD_ASYMMETRIC || quda_inv_param.matpc_type == QUDA_MATPC_ODD_ODD ){
                QDPIO::cout << "Norm2 psi = " << norm2(psi_s, rb[1])  << std::endl;
            }
            else if (quda_inv_param.matpc_type == QUDA_MATPC_EVEN_EVEN_ASYMMETRIC || quda_inv_param.matpc_type == QUDA_MATPC_EVEN_EVEN){
                QDPIO::cout << "Norm2 psi = " << norm2(psi_s, rb[0])  << std::endl;
            }
        }
        else if ( quda_inv_param.solution_type == QUDA_MAT_SOLUTION ){
            QDPIO::cout << "Norm2 psi = " << norm2(psi_s)  << std::endl;
        }
        QDPIO::cout << "QUDA_" << solver_string << "_NEF_SOLVER: time=" << quda_inv_param.secs << " s" ;
        QDPIO::cout << "\tPerformance="<<  quda_inv_param.gflops/quda_inv_param.secs<<" GFLOPS" << std::endl;


        ret.n_count =quda_inv_param.iter;

        return ret;

    }


}
