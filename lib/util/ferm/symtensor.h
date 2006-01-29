// -*- C++ -*-
// $Id: symtensor.h,v 1.1 2006-01-29 05:21:55 edwards Exp $
/*! \file
 *  \brief Compute symmetric tensors
 */

#ifndef __symtensor_h__
#define __symtensor_h__

namespace Chroma 
{

  //! Return 3d symmetric tensor
  /*!
   * \ingroup ferm
   *
   * \return  \f$s_{ijk} = |\epsilon_{ijk}|\f$
   */
  int symTensor3d(int i, int j, int k);
  
}  // end namespace Chroma

#endif
