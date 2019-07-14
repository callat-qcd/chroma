// -*- C++ -*-
/*!
 *  RMC Gauge Monomial
 */

#ifndef __rmc_coupled_gaugeact_monomial_h__
#define __rmc_coupled_gaugeact_monomial_h__

#include "chromabase.h"

#include "update/molecdyn/field_state.h"
#include "update/molecdyn/monomial/abs_monomial.h"
#include "update/molecdyn/monomial/force_monitors.h"
#include "io/xmllog_io.h"
#include "actions/gauge/gaugeacts/rmc_coupled_gaugeact.h"

namespace Chroma 
{
  /*! @ingroup monomial */
  namespace RMCCoupledGaugeMonomialEnv 
  {
    extern const std::string name;
    bool registerAll();
  }


  // Parameter structure
  /*! @ingroup monomial */
  struct RMCCoupledGaugeMonomialParams 
  {
    // Base Constructor
    RMCCoupledGaugeMonomialParams();

    // Read monomial from some root path
    RMCCoupledGaugeMonomialParams(XMLReader& in, const std::string& path);
    std::string gauge_act;
  };


  //! Wrapper class for  gauge monomials
  /*! @ingroup monomial
   *
   * Monomial is expected to be the same for these fermacts
   */
  class RMCCoupledGaugeMonomial :
    public ExactMonomial< multi1d<LatticeColorMatrix>,
                          multi1d<LatticeColorMatrix> >    
  {
  public: 
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    //! Construct out of a parameter struct. Check against the desired GaugeAct name
    RMCCoupledGaugeMonomial(const RMCCoupledGaugeMonomialParams& param_);

    //! Copy Constructor
    RMCCoupledGaugeMonomial(const RMCCoupledGaugeMonomial& m) : gaugeact((m.gaugeact)) {}

    //! Create a suitable state and compute F
    void dsdq(P& F, const AbsFieldState<P,Q>& s) 
    {
      START_CODE();

      XMLWriter& xml_out = TheXMLLogWriter::Instance();
      push(xml_out, "RMCCoupledGaugeMonomial");

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
      push(xml_out, "RMCCoupledGaugeMonomial");

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
      gaugeact->updateGamma(g_state);

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
      RMCCoupledGaugeMonomial();
      void operator=(const RMCCoupledGaugeMonomial&);

    private:
      // A handle for the gaugeact
      //Handle< GaugeAction<P,Q> > gaugeact;
      Handle< RMCCoupledGaugeActEnv::RMCCoupledGaugeAct > gaugeact;
    };


}; //end namespace chroma

#endif
