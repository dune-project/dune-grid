// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DGFPARSERUG_HH
#define DUNE_DGFPARSERUG_HH

// only include if UG is used
#if defined ENABLE_UG
#include <dune/grid/uggrid.hh>
#include "dgfparser.hh"
namespace Dune
{

  namespace dgf {

    /** \brief Grid parameters for UGGrid
        \ingroup DGFGridParameter
        The UGGridParameter class is in charge of passing UGGrid specific
        parameters to grid construction. Current parameters are:
          1. \b closure (valid values are \b none or \b green,
             which is the default value) will set the closure
             type of the returned UGGrid. \n
          2. \b copies (valid values are \b yes or \b none,
             which is the default value) will enforce that non-refined
             element are copied to the next level on refinement of a UGGrid. \n
          3. \b heapsize set heap size for UGGrid (default value is 500 MB).
        See the \b examplegrid5.dgf file for an example.
     */
    class UGGridParameterBlock
      : public GridParameterBlock
    {
    public:
      static const Flags foundClosure  = 1 << 3;
      static const Flags foundCopies   = 1 << 4;
      static const Flags foundHeapSize = 1 << 5;

    protected:
      bool _noClosure; // no closure for UGGrid
      bool _noCopy; // no copies  for UGGrid
      size_t _heapsize; // heap size  for UGGrid

    public:
      //! constructor taking istream
      UGGridParameterBlock(std::istream &in)
        : GridParameterBlock(in, false),
          _noClosure( false ),  // default value
          _noCopy( true ),      // default value
          _heapsize( 500 )      // default value (see UGGrid constructor)
      {
        int val = foundToken("closure" , "NONE" );
        if( val > 0 )
        {
          if( val == 2 ) _noClosure = true;
          foundFlags_ |= foundClosure;
        }

        val = foundToken("copies" , "NONE" );
        if( val > 0 )
        {
          if( val == 2 ) _noCopy = true;
          else _noCopy = false;
          foundFlags_ |= foundCopies;
        }
      }

      //! returns true if no closure should be used for UGGrid
      bool noClosure () const
      {
        if( (foundFlags_ & foundClosure) == 0 )
        {
          dwarn << "UGGridParameterBlock: Parameter 'closure' not specified, "
                << "defaulting to 'GREEN'." << std::endl;
        }
        return _noClosure;
      }

      //! returns true if no closure should be used for UGGrid
      bool noCopy () const
      {
        if( (foundFlags_ & foundCopies) == 0 )
        {
          dwarn << "UGGridParameterBlock: Parameter 'copies' not specified, "
                << "no copies will be generated." << std::endl;
        }
        return _noCopy;
      }

      //! returns heap size used on construction of the grid
      size_t heapSize() const
      {
        if( (foundFlags_ & foundHeapSize) == 0 )
        {
          dwarn << "UGGridParameterBlock: Parameter 'heapsize' not specified, "
                << "defaulting to '500' MB." << std::endl;
        }
        return _heapsize;
      }

    protected:
      int foundToken(const char * key, const char * value)
      {
        std::string VALUE( value );
        makeupcase( VALUE );

        // check closure
        if (findtoken( key ))
        {
          std::string clo;
          if( getnextentry(clo) )
          {
            makeupcase(clo);
            if(clo == VALUE )
            {
              return 2;
            }
          }
          return 1;
        }
        return 0;
      }
    };

  } // end namespace dgf

  /** \cond */
  template <int dim>
  class MacroGrid::Impl<UGGrid<dim> > {
    typedef MPIHelper::MPICommunicator MPICommunicatorType;
  public:
    static UGGrid<dim>* generate(MacroGrid& mg,
                                 const char* filename, MPICommunicatorType MPICOMM = MPIHelper::getCommunicator() );
  };
  template <int dimw>
  struct DGFGridInfo< UGGrid<dimw> > {
    static int refineStepsForHalf() {return 1;}
    static double refineWeight() {return -1.;}
  };
  /** \endcond */

}
#include "dgfug.cc"
#endif

#endif
