#ifndef pg_leapfrog_h
#define pg_leapfrog_h

#include "chromabase.h"
#include "update/field_state.h"
#include "update/molecdyn/abs_symp_updates.h"
#include "update/molecdyn/abs_hyb_int.h"

using namespace QDP;
using namespace std;

//! A Concrete Leapfrog for Pure Gauge Systems.
//  Templated on 

class PureGaugePQPLeapFrog : public AbsLatColMatPQPLeapFrog<
  AbsFieldState<multi1d<LatticeColorMatrix>,
	        multi1d<LatticeColorMatrix> >,
  AbsHamiltonian< multi1d<LatticeColorMatrix>,
		   multi1d<LatticeColorMatrix> >
 >
{
 public: 
  // Virtual Destructor
  ~PureGaugePQPLeapFrog() {}

  PureGaugePQPLeapFrog(AbsPureGaugeSympUpdates& symp_updates_,
		       Real delta_tau_, 
		       Real tau_) : symp_updates(symp_updates_.clone()), delta_tau(delta_tau_), tau(tau_) {}
  
  // Get at the leap P and leap Q
  // through a Symplectic Updates step.
  // Use Base Class...
  const AbsPureGaugeSympUpdates& getSympUpdates(void) const {
    return *symp_updates;
  }

  // Get step size
  virtual Real getStepSize(void) const { return delta_tau; }

  // Get traj length
  //virtual Real getTrajLength(void) const { return tau; }
  virtual Real getTrajLength(void) const { return tau; }


  PureGaugePQPLeapFrog(const PureGaugePQPLeapFrog& f) : symp_updates(f.symp_updates->clone()), delta_tau(f.delta_tau), tau(f.tau) {}
  
  virtual PureGaugePQPLeapFrog* clone(void) const {
    return new PureGaugePQPLeapFrog(*this);
  }

protected:
  Handle<AbsPureGaugeSympUpdates> symp_updates;
  const Real delta_tau;
  const Real tau;
};
#endif
