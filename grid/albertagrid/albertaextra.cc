// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*  Non-inline Implementation--File for extra Albert Functions              */
/****************************************************************************/

#include <config.h>
#include <agrid.hh>

//#include <iostream>
//#include <fstream>
//#include <vector>
//#include <assert.h>

#ifdef __ALBERTApp__
namespace Albert {
#endif

namespace AlbertHelp
{
  MeshCallBack::callBackPointer_t * MeshCallBack::postRefinementPtr = 0;
  MeshCallBack::callBackPointer_t * MeshCallBack::preCoarseningPtr  = 0;
  MESH * MeshCallBack::lockMeshPtr = 0;
  void * MeshCallBack::dataHandler = 0;
}   // end namespace AlbertHelp

#ifdef __ALBERTApp__
} // end namespace Albert
#endif
