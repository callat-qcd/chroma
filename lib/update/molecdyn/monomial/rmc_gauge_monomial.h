// -*- C++ -*-
/*!
 *  RMC Gauge Monomial
 */

#ifndef __rmc_gaugeact_monomial_h__
#define __rmc_gaugeact_monomial_h__

#include "chromabase.h"

#include "update/molecdyn/field_state.h"
#include "update/molecdyn/monomial/abs_monomial.h"
#include "update/molecdyn/monomial/force_monitors.h"
#include "io/xmllog_io.h"
#include "actions/gauge/gaugeacts/rmc_gaugeact.h"

namespace Chroma 
{
  /*! @ingroup monomial */
  namespace RMCGaugeMonomialEnv 
  {
    extern const std::string name;
    bool registerAll();
  }


  // Parameter structure
  /*! @ingroup monomial */
  struct RMCGaugeMonomialParams 
  {
    // Base Constructor
    RMCGaugeMonomialParams();

    // Read monomial from some root path
    RMCGaugeMonomialParams(XMLReader& in, const std::string& path);
    std::string gauge_act;
  };


  //! Wrapper class for  gauge monomials
  /*! @ingroup monomial
   *
   * Monomial is expected to be the same for these fermacts
   */
  class RMCGaugeMonomial :
    public ExactMonomial< multi1d<LatticeColorMatrix>,
                          multi1d<LatticeColorMatrix> >    
  {
  public: 
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    //! Construct out of a parameter struct. Check against the desired GaugeAct name
    RMCGaugeMonomial(const RMCGaugeMonomialParams& param_);

    //! Copy Constructor
    RMCGaugeMonomial(const RMCGaugeMonomial& m) : gaugeact((m.gaugeact)) {}

    //! Create a suitable state and compute F
    void dsdq(P& F, const AbsFieldState<P,Q>& s) 
    {
      START_CODE();

      XMLWriter& xml_out = TheXMLLogWriter::Instance();
      push(xml_out, "RMCGaugeMonomial");

      // Make a gauge connect state
      Handle< GaugeState<P,Q> > g_state(getGaugeAct().createState(s.getQ()));

      gaugeact->RMC_deriv(F,g_state);

      monitorForces(xml_out, "Forces", F);
      pop(xml_out);

      END_CODE();
    }


    //! Gauge action value
    Double S(const AbsFieldState<P,Q>& s)  
    {
      START_CODE();

      XMLWriter& xml_out = TheXMLLogWriter::Instance();
      push(xml_out, "RMCGaugeMonomial");

      Handle< GaugeState<P,Q> > g_state(getGaugeAct().createState(s.getQ()));
 
      Double action = gaugeact -> RMC_S(g_state);
   
      write(xml_out, "S", action);
      pop(xml_out);

      END_CODE();

      return action;
    }
	
	
    void refreshInternalFields(const AbsFieldState<P,Q>& s)
    {
      //No internal fields to refresh => Nop
    
      Handle< GaugeState<P,Q> > g_state(getGaugeAct().createState(s.getQ()));
      
      gaugeact->updateBeta(g_state);
      /*RMCGaugeActEnv::GaugeAct downcast=dynamic_cast<GaugeAction<P,Q>*>(&gaugeact);

      
      // Check success of the downcast 
      if( downcast == 0x0 ) 
      {
	QDPIO::cerr << __func__ << ": unable to downcast to RMC Gauge Act" << std::endl;
	QDP_abort(1);
      }*/

    } //Rmc update will go here.

    void setInternalFields(const Monomial<P,Q>& m) 
    {
      // No internal fields to refresh => Nop
    }

    protected:
      const GaugeAction<P,Q>& getGaugeAct(void) const { 
	return *gaugeact;
      }

    private:
      // Hide empty constructor and =
      RMCGaugeMonomial();
      void operator=(const RMCGaugeMonomial&);

    private:
      // A handle for the gaugeact
      //Handle< GaugeAction<P,Q> > gaugeact;
      Handle< RMCGaugeActEnv::RMCGaugeAct > gaugeact;
    };


}; //end namespace chroma

#endif
