// -*- C++ -*-
/*! \file
 *  \brief  gauge action as a sum of characters of the SU(3) irreps
 */

#ifndef __rmc_coupled_gaugeact_h__
#define __rmc_coupled_gaugeact_h__

#include "chromabase.h"
#include "gaugeact.h"
#include "gaugebc.h"

namespace Chroma
{

  /*! @ingroup gaugeacts */
  namespace RMCCoupledGaugeActEnv 
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

      Real beta;   // Coupling Constant for Lattice 1/Shadow Simulation
      Real alpha;   // Coupling Constant for Lattice 2/Driven Simulation
      Real gamma; // Coupling Constant for RMC Coupling bet. Lattices
      mutable multi2d<LatticeReal> RMCGamma{Nd,Nd}; //Coupling Lattice
      mutable multi1d<LatticeColorMatrix> shadow_gauge_field{Nd}; //Shadow Simulation Gauge Field
      mutable std::string shadow_gauge_id; //Shadow Simulation Gauge Field LIME ID
      QDP_serialparallel_t file_format; //File Format; either set as "SERIAL" or "PARALLEL" in xml
      mutable int curr_num; //Current Update Number
      int freeze_steps; //Number of Markovian steps that RMCGamma is frozen
      bool freeze_check; //Check if freezing is turned on or not

    };
  

    //! RMC gauge action
    /*! \ingroup gaugeacts
     *
     * The gauge action with Reverse Monte Carlo
     */
    class RMCCoupledGaugeAct : public LinearGaugeAction
    {
    public:
      // Typedefs to save typing
      typedef multi1d<LatticeColorMatrix>  P;
      typedef multi1d<LatticeColorMatrix>  Q;

      //! General CreateGaugeState<P,Q>
      //! Read coeff from a param struct
      RMCCoupledGaugeAct(Handle< CreateGaugeState<P,Q> > cgs_, const Params& p) :
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

      //! Update coupling coefficient Gamma using RMC
      void updateGamma(const Handle< GaugeState<P,Q> >& state) const;

      //! Restore/initialize coupling coefficient Gammaa to original values
      void restoreGamma() const;

      //! Compute the RMC Action
      Double RMC_S(const Handle< GaugeState<P,Q> >& state) const;

      //! Compute dS_{RMC}/dU
      void RMC_deriv(multi1d<LatticeColorMatrix>& result,
		      const Handle< GaugeState<P,Q> >& state) const;


      //! Destructor is automatic
      ~RMCCoupledGaugeAct() {}

    protected:
      //! Hide assignment
      void operator=(const RMCCoupledGaugeAct& a) {}
      
      //! Compute the site-level action
      void siteAction(multi2d<LatticeColorMatrix>& site_act, 
		      const Handle< GaugeState<P,Q> >& state) const;

      //! Compute the site-level action for the shadow field
      void siteAction_shadow(multi2d<LatticeColorMatrix>& site_act) const;

    private:
      Handle< CreateGaugeState<P,Q> >  cgs;  /*!< Create Gauge State */
      Params               param;            /*!< The parameters */
    };

  }

}


#endif
