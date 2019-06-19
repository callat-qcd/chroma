// -*- C++ -*-
/*! \file
 *  \brief  gauge action as a sum of characters of the SU(3) irreps
 */

#ifndef __rmc_gaugeact_h__
#define __rmc_gaugeact_h__

#include "chromabase.h"
#include "gaugeact.h"
#include "gaugebc.h"

namespace Chroma
{

  /*! @ingroup gaugeacts */
  namespace RMCGaugeActEnv 
  { 
    bool registerAll();

    //! Parameter structure
    /*! @ingroup gaugeacts */
    struct Params 
    {
      //! Base Constructor
      Params() {}
    
      //! Read params from some root path
      Params(XMLReader& xml_in, const std::string& path);

      Real  beta;   // Original coupling constant
      Real  alpha;   // New coupling constant for RMC
      mutable multi2d<LatticeReal> RMCBeta{Nd,Nd}; //RMC coupling constant 
      mutable int curr_num; //Current Update Number
      int freeze_steps; //Number of Markovian steps that RMCBeta is frozen
      bool freeze_check; //Check if freezing is turned on or not

    };
  

    //! RMC gauge action
    /*! \ingroup gaugeacts
     *
     * The gauge action with Reverse Monte Carlo
     */
    class RMCGaugeAct : public LinearGaugeAction
    {
    public:
      // Typedefs to save typing
      typedef multi1d<LatticeColorMatrix>  P;
      typedef multi1d<LatticeColorMatrix>  Q;

      //! General CreateGaugeState<P,Q>
      //! Read coeff from a param struct
      RMCGaugeAct(Handle< CreateGaugeState<P,Q> > cgs_, const Params& p) :
	cgs(cgs_), param(p) {}

      //! Return the set on which the gauge action is defined
      /*! Defined on the even-off (red/black) set */
      const Set& getSet() const {return rb;}

      //! Compute staple
      /*! Default version. Derived class should override this if needed. */
      void staple(LatticeColorMatrix& result,
		  const Handle< GaugeState<P,Q> >& state,
		  int mu, int cb) const;

      //! Compute dS/dU
      void deriv(multi1d<LatticeColorMatrix>& result,
		 const Handle< GaugeState<P,Q> >& state) const;

      //! Compute the actions
      Double S(const Handle< GaugeState<P,Q> >& state) const;

      //! Produce a gauge create state object
      const CreateGaugeState<P,Q>& getCreateState() const {return *cgs;}

      //! Update coupling coefficient Beta using RMC
      void updateBeta(const Handle< GaugeState<P,Q> >& state) const;

      //! Restore/initialize coupling coefficient Beta to original values
      void restoreBeta() const;

      //! Compute the RMC Action
      Double RMC_S(/*multi2d<LatticeReal>& RMCBeta,*/
		      const Handle< GaugeState<P,Q> >& state) const;

      //! Compute dS_{RMC}/dU
      void RMC_deriv(multi1d<LatticeColorMatrix>& result,
		      const Handle< GaugeState<P,Q> >& state) const;


      //! Destructor is automatic
      ~RMCGaugeAct() {}

    protected:
      //! Hide assignment
      void operator=(const RMCGaugeAct& a) {}
      
      //! Compute the site-level action
      void siteAction(multi2d<LatticeColorMatrix>& site_act, 
		      const Handle< GaugeState<P,Q> >& state) const;

    private:
      Handle< CreateGaugeState<P,Q> >  cgs;  /*!< Create Gauge State */
      Params               param;            /*!< The parameters */
    };

  }

}


#endif
