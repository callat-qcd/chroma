// -*- C++ -*-
/*! \file
 *  \brief rmc gauge action monomial
 */

#include "chromabase.h"

#include "update/molecdyn/monomial/rmc_gauge_monomial.h"
#include "update/molecdyn/monomial/monomial_factory.h"
#include "actions/gauge/gaugeacts/gaugeacts_aggregate.h"
#include "actions/gauge/gaugeacts/gaugeact_factory.h"
#include "actions/gauge/gaugestates/gauge_createstate_factory.h"
#include "actions/gauge/gaugestates/gauge_createstate_aggregate.h"

namespace Chroma 
{ 
  namespace RMCGaugeMonomialEnv 
  {
    namespace
    {
      //! Callback function for the factory
      Monomial< multi1d<LatticeColorMatrix>,
		multi1d<LatticeColorMatrix> >*
      createMonomial(XMLReader& xml, const std::string& path) 
      {
	QDPIO::cout << "Create monomial: " << name << std::endl;

	return new RMCGaugeMonomial(RMCGaugeMonomialParams(xml, path));
      }

      //! Local registration flag
      bool registered = false;
    }

    const std::string name("RMC_GAUGE_MONOMIAL");
    
    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= GaugeActsEnv::registerAll();
	success &= TheMonomialFactory::Instance().registerObject(name, createMonomial);
	registered = true;
      }
      return success;
    }
    /*void refreshInternalFields(const AbsFieldState<P,Q>& s) 
    {
    }*/ 
  } //end namespace RMCGaugeMonomialEnv



  // Read the parameters
  RMCGaugeMonomialParams::RMCGaugeMonomialParams(XMLReader& xml_in, const std::string& path)
  {
    // Get the top of the parameter XML tree
    XMLReader paramtop(xml_in, path);
    
    try {
      // Read the inverter Parameters
      XMLReader xml_tmp(paramtop, "./GaugeAction");
      std::ostringstream os;
      xml_tmp.print(os);
      gauge_act = os.str();
   
    }
    catch(const std::string& s) {
      QDPIO::cerr << "Caught Exception while reading parameters: " << s <<std::endl;
      QDP_abort(1);
    }

    QDPIO::cout << "RMCGaugeMonomialParams: read \n" << gauge_act << std::endl;
  }

  //! Read Parameters
  void read(XMLReader& xml, const std::string& path,
	    RMCGaugeMonomialParams& params) 
  {
    RMCGaugeMonomialParams tmp(xml, path);
    params = tmp;
  }

  //! Write Parameters
  void write(XMLWriter& xml, const std::string& path,
	     const RMCGaugeMonomialParams& params) 
  {
    // Not implemented
    QDPIO::cerr << RMCGaugeMonomialEnv::name << ": write not implemented" << std::endl;
    QDP_abort(1);
  }


  // Constructor
  RMCGaugeMonomial::RMCGaugeMonomial(const RMCGaugeMonomialParams& param_) 
  {
    std::istringstream is(param_.gauge_act);
    XMLReader gaugeact_reader(is);

    // Get the name of the gauge act
    std::string gaugeact_string;
    try { 
      read(gaugeact_reader, "/GaugeAction/Name", gaugeact_string);
    }
    catch( const std::string& e) { 
      QDPIO::cerr << "Error grepping the gaugeact name: " << e<<  std::endl;
      QDP_abort(1);
    }

    // Throw an exception if not found
    /*gaugeact = TheGaugeActFactory::Instance().createObject(gaugeact_string, 
							   gaugeact_reader, 
							   "/GaugeAction");*/
    gaugeact = new RMCGaugeActEnv::RMCGaugeAct(CreateGaugeStateEnv::reader(gaugeact_reader, "/GaugeAction"),
	       RMCGaugeActEnv::Params(gaugeact_reader, "/GaugeAction"));
  }
  

} //end namespace Chroma


