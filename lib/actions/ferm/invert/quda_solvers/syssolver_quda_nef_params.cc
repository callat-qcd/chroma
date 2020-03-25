#include "actions/ferm/invert/quda_solvers/syssolver_quda_nef_params.h"
#include "chromabase.h"
#include "io/xml_group_reader.h"
#include "chroma_config.h"



using namespace QDP;

namespace Chroma {

    SysSolverQUDANEFParams::SysSolverQUDANEFParams(XMLReader& xml, 
                                                   const std::string& path)
    {
        XMLReader paramtop(xml, path);

        read(paramtop, "MaxIter", MaxIter);
        read(paramtop, "RsdTarget", RsdTarget);
        read(paramtop, "NEFParams", NEFParams);
        read(paramtop, "AntiPeriodicT", AntiPeriodicT);

        read(paramtop, "Delta", Delta);
  

        read(paramtop, "SolverType", solverType);

        if ( paramtop.count("Verbose") > 0 ) { 
            read(paramtop, "Verbose", verboseP);
        }
        else { 
            verboseP = false;
        }
        if ( paramtop.count("AsymmetricLinop") > 0 ) { 
            read(paramtop, "AsymmetricLinop", asymmetricP);
        }
        else { 
            asymmetricP = false; // Symmetric is default 
        }

        if( paramtop.count("CudaPrecision") > 0 ) {
            read(paramtop, "CudaPrecision", cudaPrecision);
        }
        else { 
            cudaPrecision = DEFAULT;
        }
        if( paramtop.count("CudaReconstruct") > 0 ) {
            read(paramtop, "CudaReconstruct", cudaReconstruct);
        }
        else { 
            cudaReconstruct = RECONS_12;
        }

        // Sloppy params
        if( paramtop.count("CudaSloppyPrecision") > 0 ) {
            read(paramtop, "CudaSloppyPrecision", cudaSloppyPrecision);
        }
        else { 
            cudaSloppyPrecision = DEFAULT;
        }

        if( paramtop.count("CudaSloppyReconstruct") > 0 ) {
            read(paramtop, "CudaSloppyReconstruct", cudaSloppyReconstruct);
        }
        else { 
            cudaSloppyReconstruct = RECONS_12;
        }
        // Precondition params
        if( paramtop.count("CudaPreconditionPrecision") > 0 ) {
            read(paramtop, "CudaPreconditionPrecision", cudaPreconditionPrecision);
        }
        else {
            cudaPreconditionPrecision = DEFAULT;
        }

        if( paramtop.count("CudaPreconditionReconstruct") > 0 ) {
            read(paramtop, "CudaPreconditionReconstruct", cudaPreconditionReconstruct);
        }
        else {
            cudaPreconditionReconstruct = RECONS_12;
        }

        if( paramtop.count("AxialGaugeFix") > 0 ) {
            read(paramtop, "AxialGaugeFix", axialGaugeP);
        }
        else { 
            axialGaugeP = false;
        }

        if( paramtop.count("SilentFail") > 0) {
            read(paramtop, "SilentFail", SilentFailP);
        }
        else { 
            SilentFailP = false;
        }

        if( paramtop.count("RsdToleranceFactor") > 0 ) { 
            read(paramtop, "RsdToleranceFactor", RsdToleranceFactor);
        }
        else { 
            RsdToleranceFactor = Real(10); // Tolerate an order of magnitude difference by default.
        }

        if( paramtop.count("AutotuneDslash") > 0 ) { 
            read(paramtop, "AutotuneDslash", tuneDslashP);
        }
        else { 
            tuneDslashP = false;
        }
        QDPIO::cout << "tuneDslasP = " << tuneDslashP << std::endl;


        if( paramtop.count("GCRInnerParams") > 0 ) {
            innerParams = new GCRInnerSolverParams(paramtop, "./GCRInnerParams");
            innerParamsP = true;
        }
        else { 
            innerParamsP = false;
        }

        if ( paramtop.count("BackupSolverParam") > 0 ) { 
            // If user specified a backup solver, let's read it
            backup_invP = true;
            backup_inv_param = readXMLGroup(paramtop, "./BackupSolverParam", "invType");
        }
        else { 
            // No backup
            backup_invP = false;
            backup_inv_param.xml = "";
        }

        if ( paramtop.count("DumpOnFail") > 0 ) { 
            read(paramtop, "DumpOnFail", dump_on_failP);
        }
        else { 
            dump_on_failP  = false;
        }
   
        if ( paramtop.count("DoCGNR")  > 0 ) { 
            read(paramtop, "DoCGNR", cgnrP);
        }
        else { 
            cgnrP = false ; // Do CGNE by default
        }
        if( paramtop.count("NEFLanczosParams") > 0 ) {
            read(paramtop, "NEFLanczosParams", NEFLanczosParams);
            //NEFLanczosParams = new LanczosParams(paramtop, "./NEFLanczosParams");
            InvDeflate = true;
        }
        else { 
            InvDeflate = false;
        }



    }

    void read(XMLReader& xml, const std::string& path, 
              SysSolverQUDANEFParams& p)
    {
        SysSolverQUDANEFParams tmp(xml, path);
        p = tmp;
    }

    void write(XMLWriter& xml, const std::string& path, 
               const SysSolverQUDANEFParams& p) {
        push(xml, path);
        write(xml, "MaxIter", p.MaxIter);
        write(xml, "RsdTarget", p.RsdTarget);
        write(xml, "NEFParams", p.NEFParams);
        write(xml, "AntiPeriodicT", p.AntiPeriodicT);
        write(xml, "Delta", p.Delta);
        write(xml, "SolverType", p.solverType);
        write(xml, "Verbose", p.verboseP);
        write(xml, "AsymmetricLinop", p.asymmetricP);
        write(xml, "CudaPrecision", p.cudaPrecision);
        write(xml, "CudaReconstruct", p.cudaReconstruct);
        write(xml, "CudaSloppyPrecision", p.cudaSloppyPrecision);
        write(xml, "CudaSloppyReconstruct", p.cudaSloppyReconstruct);
        write(xml, "CudaPreconditionPrecision", p.cudaPreconditionPrecision);
        write(xml, "CudaPreconditionReconstruct", p.cudaPreconditionReconstruct);
        write(xml, "AxialGaugeFix", p.axialGaugeP);
        write(xml, "SilentFail", p.SilentFailP);
        write(xml, "RsdToleranceFactor", p.RsdToleranceFactor);

        write(xml, "AutotuneDslash", p.tuneDslashP);
        if( p.innerParamsP ) { 
            write(xml, "GCRInnerParams", *(p.innerParams));
        }

        write(xml, "DumpOnFail", p.dump_on_failP);
        write(xml, "DoCGNR", p.cgnrP );

        // If we added LanczosParams, InvDeflate is set to true
        if ( p.InvDeflate ) {
            write(xml, "NEFLanczosParams", p.NEFLanczosParams);
        }

        if( p.backup_invP ) { 
            // Need to dump out the XML for the back up solver here...
            // Turn XML into an std::istringstream
            std::istringstream is( p.backup_inv_param.xml);
            XMLReader tmp_reader(is);
            // I wonder if this method works....
            write(xml, "BackupSolverParam", tmp_reader);
        }
        pop(xml);

    }



}
