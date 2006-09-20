// -*- C++ -*-
// $Id: prec_slic_fermact_w.h,v 3.2 2006-09-20 20:27:59 edwards Exp $
/*! \file
 *  \brief Even-odd preconditioned Clover fermion action
 */

#ifndef __prec_slic_fermact_w_h__
#define __prec_slic_fermact_w_h__

#include "fermact.h"
#include "actions/ferm/linop/lgherm_w.h"
#include "actions/ferm/fermacts/clover_fermact_params_w.h"

namespace Chroma 
{ 
  //! Name and registration
  /*! \ingroup fermacts */
  namespace EvenOddPrecSLICFermActEnv
  {
    extern const std::string name;
    bool registerAll();
  }
  

  //! Even-odd preconditioned Clover fermion action
  /*! \ingroup fermacts
   *
   * Even-odd preconditioned clover fermion action. 
   * Only defined on odd subset.
   */

  class EvenOddPrecSLICFermAct : public EvenOddPrecLogDetWilsonTypeFermAct<LatticeFermion, 
				   multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >
  {
  public:
    // Typedefs to save typing
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    //! Partial constructor
    EvenOddPrecSLICFermAct() {}

    //! General FermState
    EvenOddPrecSLICFermAct(Handle< CreateFermState<T,P,Q> > cfs_, 
			     const CloverFermActParams& param_) : 
      cfs(cfs_), param(param_) {}

    //! Copy constructor
    EvenOddPrecSLICFermAct(const EvenOddPrecSLICFermAct& a) : 
      cfs(a.cfs), param(a.param) {}

    //! Produce a linear operator for this action
    EvenOddPrecLogDetLinearOperator<T,P,Q>* linOp(Handle< FermState<T,P,Q> > state) const;

    //! Produce the gamma_5 hermitian operator H_w
    LinearOperator<LatticeFermion>* hermitianLinOp(Handle< FermState<T,P,Q> > state) const 
      { 
	return new lgherm<LatticeFermion>(linOp(state));
      }

    //! Destructor is automatic
    ~EvenOddPrecSLICFermAct() {}

  protected:
    //! Return the fermion BC object for this action
    const CreateFermState<T,P,Q>& getCreateState() const {return *cfs;}

    //! Assignment
    void operator=(const EvenOddPrecSLICFermAct& a) {}

  private:
    Handle< CreateFermState<T,P,Q> >  cfs;
    CloverFermActParams param;
  };

}; // End Namespace Chroma


#endif
