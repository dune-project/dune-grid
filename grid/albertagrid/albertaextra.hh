// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*  Header--File for extra Albert Functions                                 */
/****************************************************************************/
#ifndef DUNE_ALBERTAEXTRA_HH
#define DUNE_ALBERTAEXTRA_HH

#include <algorithm>
#include <cassert>
#include <cstring>

#include <dune/grid/albertagrid/albertaheader.hh>

#if HAVE_ALBERTA

#ifdef __ALBERTApp__
namespace Albert {
#endif

//*********************************************************************
//
//  Help Routines for the ALBERTA Mesh
//
//*********************************************************************
namespace AlbertHelp
{

  //*********************************************************************

  // Leaf Data for Albert, only the leaf elements have this data set
  template <int cdim, int vertices>
  struct AlbertLeafData
  {
    // type of stored data
    typedef struct
    {
      double determinant;
    } Data;

    // keep element numbers
    inline static void AlbertLeafRefine( EL *parent, EL *child[2] )
    {
      Data * ldata;
      int i;

      ldata = (Data *) parent->child[1];
      assert(ldata != 0);

      //std::cout << "Leaf refine for el = " << parent << "\n";

      double childDet = 0.5 * ldata->determinant;

      /* bisection ==> 2 children */
      for(i=0; i<2; i++)
      {
        Data *ldataChi = (Data *) child[i]->child[1];
        assert(ldataChi != 0);
        ldataChi->determinant = childDet;
      }
    }

    inline static void AlbertLeafCoarsen(EL *parent, EL *child[2])
    {
      Data *ldata;
      int i;

      ldata = (Data *) parent->child[1];
      assert(ldata != 0);
      double & det = ldata->determinant;
      det = 0.0;

      //std::cout << "Leaf coarsen for el = " << parent << "\n";

      /* bisection ==> 2 children */
      for(i=0; i<2; i++)
      {
        Data *ldataChi = (Data *) child[i]->child[1];
        assert(ldataChi != 0);
        det += ldataChi->determinant;
      }
    }

    // function for mesh_traverse, is called on every element
    inline static void setLeafData(const EL_INFO * elf)
    {
      assert( elf->el->child[0] == 0 );
      Data *ldata = (Data *) elf->el->child[1];
      assert(ldata != 0);

      ldata->determinant = ALBERTA el_det(elf);
    }

    // remember on which level an element realy lives
    inline static void initLeafDataValues( MESH * mesh, int proc )
    {
      // see ALBERTA Doc page 72, traverse over all hierarchical elements
      ALBERTA meshTraverse(mesh,-1, CALL_LEAF_EL|FILL_COORDS,setLeafData);
    }

  }; // end of AlbertLeafData

} // end namespace AlbertHelp

#ifdef __ALBERTApp__
} // end namespace Albert
#endif

#endif // HAVE_ALBERTA

#endif  /* !_ALBERTAEXTRA_H_ */
