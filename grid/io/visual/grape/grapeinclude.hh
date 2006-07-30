// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRIDDISPLAY_HH
#define DUNE_GRIDDISPLAY_HH

#include <stdarg.h>
#include <assert.h>
#include <string.h>

#include "grapewrapper.hh"
#include "grapecommon.hh"

enum { MAX_EL_DOF  = 8 };
enum { MAX_EL_FACE = 6 };

namespace GrapeInterface_two_two
{
#define GRAPE_DIM 2
#define GRAPE_DIMWORLD 2
#include "grapehmesh.hh"
}

namespace GrapeInterface_two_three
{
#define GRAPE_DIM 2
#define GRAPE_DIMWORLD 3
#include "grapehmesh.hh"
}

namespace GrapeInterface_three_three
{
#define GRAPE_DIM 3
#define GRAPE_DIMWORLD 3
#include "grapehmesh.hh"
}

namespace Dune
{

  static int __MaxPartition = 1;
  // the interface to dune
  template <int dim, int dimworld>
  struct GrapeInterface
  {
    static int called;
    typedef GrapeInterface_two_two::DUNE_ELEM DUNE_ELEM;
    typedef GrapeInterface_two_two::DUNE_FDATA DUNE_FDATA;
    typedef GrapeInterface_two_two::DUNE_DAT DUNE_DAT;
    typedef GrapeInterface_two_two::DUNE_FUNC DUNE_FUNC;

    inline static void init()
    {
      GrapeInterface_two_two::grape_add_remove_methods();
    }

    inline static void setThread(int t)
    {}

    inline static void handleMesh (void *hmesh, bool grdMode = false )
    {
      GrapeInterface_two_two::handleMesh(hmesh,grdMode);
    }
    inline static void addDataToHmesh(void  *hmesh, DUNE_FDATA * fe,
                                      void (* const func_real) (DUNE_ELEM *, DUNE_FDATA*, int ind, const
                                                                double *, double *)  )
    {
      GrapeInterface_two_two::addDataToHmesh(hmesh,fe,func_real);
    }

    inline static void *setupHmesh(
      void (* const func_real) (DUNE_ELEM *, DUNE_FDATA*, int ind, const double *coord,  double *),
      const int noe, const int nov, const int maxlev,
      DUNE_FDATA * fe, DUNE_DAT * dune )
    {
      return GrapeInterface_two_two::setupHmesh(
               func_real,noe,nov,maxlev,fe,dune);
    }

    inline static void deleteHmesh( void * hmesh )
    {
      GrapeInterface_two_two::deleteHmesh( hmesh );
    }

    inline static void addHmeshToTimeScene(void * timescene, double time, void  *hmesh , int proc)
    {
      GrapeInterface_two_two::addHmeshToTimeScene(timescene,time,hmesh,proc);
    }

    inline static void addHmeshToGlobalTimeScene(double time, void  *hmesh , int proc)
    {
      GrapeInterface_two_two::addHmeshToGlobalTimeScene(time,hmesh,proc);
    }

    inline static void colorBarMinMax(const double min, const double max)
    {
      GrapeInterface_two_two::colorBarMinMax(min,max);
    }
  };

  // the interface to dune for dim = dimworld = 3
  template <>
  struct GrapeInterface<3,3>
  {
    typedef GrapeInterface_three_three::DUNE_ELEM DUNE_ELEM;
    typedef GrapeInterface_three_three::DUNE_FDATA DUNE_FDATA;
    typedef GrapeInterface_three_three::DUNE_DAT DUNE_DAT;
    typedef GrapeInterface_three_three::DUNE_FUNC DUNE_FUNC;

    inline static void init()
    {
      GrapeInterface_three_three::initPartitionDisp(__MaxPartition);
      GrapeInterface_three_three::grape_add_remove_methods();
    }

    inline static void setThread(int t)
    {
      GrapeInterface_three_three::setThread(t);
    }

    inline static void handleMesh (void *hmesh, bool grdMode = false )
    {
      GrapeInterface_three_three::handleMesh(hmesh,grdMode);
    }
    inline static void addDataToHmesh(void  *hmesh, DUNE_FDATA * fe,
                                      void (* const func_real) (DUNE_ELEM *, DUNE_FDATA*, int ind, const
                                                                double *, double *)  )
    {
      GrapeInterface_three_three::addDataToHmesh(hmesh,fe,func_real);
    }

    inline static void *setupHmesh(
      void (* const func_real) (DUNE_ELEM *, DUNE_FDATA*, int ind, const double *coord,  double *),
      const int noe, const int nov, const int maxlev,
      DUNE_FDATA * fe, DUNE_DAT * dune )
    {
      return GrapeInterface_three_three::
             setupHmesh(func_real,noe,nov,maxlev,fe,dune);
    }

    inline static void deleteHmesh( void * hmesh )
    {
      GrapeInterface_three_three::deleteHmesh( hmesh );
    }

    inline static void addHmeshToTimeScene(void * timescene, double time, void  *hmesh , int proc)
    {
      GrapeInterface_three_three::addHmeshToTimeScene(timescene,time,hmesh,proc);
    }
    inline static void addHmeshToGlobalTimeScene(double time, void  *hmesh , int proc)
    {
      GrapeInterface_three_three::addHmeshToGlobalTimeScene(time,hmesh,proc);
    }

    inline static void colorBarMinMax(const double min, const double max)
    {
      GrapeInterface_three_three::colorBarMinMax(min,max);
    }
  };


} // end namespace Dune

#include "grapecommon.cc"

#endif
