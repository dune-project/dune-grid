// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DGFPARSERYASP_HH
#define DUNE_DGFPARSERYASP_HH
#include <dune/grid/yaspgrid.hh>
#include "dgfparser.hh"
namespace Dune {

  namespace dgf {

    /** \brief Grid parameters for YaspGrid
        \ingroup DGFGridParameter
        The YaspGridParameter class is in charge of passing YaspGrid specific
        parameters to grid construction. Current parameters are: \n \n
          1. \b overlap defining the overlap of the grid (default value is zero) \n
          2. \b periodic defining which dimension should have periodic
                boundaries, i.e. passing \b periodic 0 1 will set
                periodic boundaries for x and y direction. \n
        See the \b examplegrid5.dgf file for an example.
     */
    class YaspGridParameterBlock
      : public GridParameterBlock
    {
    protected:
#if 0
      std::set<int> _periodic; // periodic grid
#endif
      int _overlap; // overlap for YaspGrid

    public:
      //! constructor taking istream
      YaspGridParameterBlock( std::istream &in )
        : GridParameterBlock( in ),
#if 0
          _periodic(), // default is non periodic
#endif
          _overlap( 0 ) // default value
      {
        // check overlap
        if( findtoken( "overlap" ) )
        {
          int x;
          if( getnextentry(x) ) _overlap = x;
          else
          {
            dwarn << "GridParameterBlock: found keyword `overlap' but no value, defaulting to `" <<  _overlap  <<"' !\n";
          }

          if (_overlap < 0)
          {
            DUNE_THROW(DGFException,"Negative overlap specified!");
          }
        }
        else
        {
          dwarn << "YaspGridParameterBlock: Parameter 'overlap' not specified, "
                << "defaulting to '" << _overlap << "'." << std::endl;
        }

#if 0
        // check periodic grid
        if (findtoken("periodic"))
        {
          int x;
          while (getnextentry(x))
          {
            _periodic.insert(x);
          }
        }
        else
        {
          dwarn << "YaspGridParameterBlock: Parameter 'periodic' not specified, "
                << "defaulting to no periodic boundary." << std::endl;
        }
#endif
      }

      //! get dimension of world found in block
      int overlap () const
      {
        return _overlap;
      }

#if 0
      //! returns true if dimension is periodic
      bool isPeriodic ( const int dim ) const
      {
        return (_periodic.find(dim) != _periodic.end());
      }
#endif
    };

  }

  template <int dim>
  class MacroGrid::Impl<YaspGrid<dim> > {
    typedef MPIHelper::MPICommunicator MPICommunicatorType;
  public:
    static YaspGrid<dim>* generate(MacroGrid& mg,
                                   const char* filename, MPICommunicatorType MPICOMM = MPIHelper::getCommunicator() );
  };
  template <int dim>
  struct DGFGridInfo< YaspGrid<dim> > {
    static int refineStepsForHalf() {return 1;}
    static double refineWeight() {return pow(0.5,dim);}
  };
}
#include "dgfyasp.cc"
#endif
