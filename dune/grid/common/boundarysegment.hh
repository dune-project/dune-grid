// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_COMMON_BOUNDARY_SEGMENT_HH
#define DUNE_GRID_COMMON_BOUNDARY_SEGMENT_HH

#include <map>
#include <sstream>

#include <dune/common/singleton.hh>
#include <dune/common/parameterizedobject.hh>
#include <dune/common/fvector.hh>

/** \file
    \brief Base class for grid boundary segments of arbitrary geometry
 */

namespace Dune {

  /** \brief Base class for classes implementing geometries of boundary segments

     Some grid implementations, as for example UGGrid, allow to use boundary segments of
     arbitrary geometry.  That means that you can have grid boundaries approach smooth
     shapes when refining the grid.

     Such curved boundary segments are specified by giving classes that implement them.
     Each boundary segment is implemented by an object of a class derived from
     BoundarySegment.  The set of these objects is handed over to the grid upon grid
     construction.

     \tparam dim Dimension of the grid
     \tparam dimworld Dimension of the Euclidean space the grid is embedded in
     \tparam ctype type of coordinate storage (default is double)
   */
  template< int dim, int dimworld = dim, class ctype = double >
  struct BoundarySegment;

  template <class BndSeg>
  class BoundarySegmentBackupRestore
  {
  public:
    // type of object stream used for storing boundary segment information
    typedef std::stringstream ObjectStreamType ;

  protected:
    //! type of BoundarySegment interface class
    typedef BndSeg  BoundarySegment;

    //! type of factory creating a unique_ptr from an ObjectStreamType
    typedef Dune::ParameterizedObjectFactory< std::unique_ptr< BoundarySegment > ( ObjectStreamType& ), int > FactoryType;

    /** \brief create an object of BoundarySegment type from a previously
     *         registered factory linked to key.
     *  \param in stream buffer previously written with backup containing key and object data
     *
     *  \return Object derived from BoundarySegment.
     */
    static std::unique_ptr< BoundarySegment > restore( ObjectStreamType& in )
    {
      int key = -1;
      // read class key for restore
      in.read( (char *) &key, sizeof( int ) );

      // factory creates a unique_ptr which can be released later on
      return factory().create( key, in );
    }

    template <class DerivedType>
    static int registerFactory()
    {
      const int key = createKey();
      // create factory method that produces unique_ptr
      factory().template define< DerivedType >( key );
      // return key for storage in derived class
      return key;
    }

  private:
    static int createKey()
    {
      static int key = 0;
      return key++;
    }

    static FactoryType& factory()
    {
      return Dune::Singleton< FactoryType > :: instance();
    }
  };

  template< int dim, int dimworld, class ctype >
  struct BoundarySegment : public BoundarySegmentBackupRestore< BoundarySegment< dim, dimworld, ctype > >
  {
    typedef BoundarySegment< dim, dimworld, ctype > ThisType;
    typedef BoundarySegmentBackupRestore< BoundarySegment< dim, dimworld, ctype > > BaseType;

    typedef typename BaseType :: ObjectStreamType ObjectStreamType;

    using BaseType :: restore;
    using BaseType :: registerFactory;

    /** \brief Dummy virtual destructor */
    virtual ~BoundarySegment() {}

    /** \brief A function mapping local coordinates on a boundary segment to world coordinates
     */
    virtual FieldVector< ctype, dimworld >
    operator() ( const FieldVector< ctype, dim-1> &local ) const = 0;

    /** \brief write BoundarySegment's data to stream buffer
     *  \param buffer buffer to store data
     */
    virtual void backup( [[maybe_unused]] ObjectStreamType& buffer ) const
    {
      DUNE_THROW(NotImplemented,"BoundarySegment::backup needs to be overloaded!");
    }
  };


}  // end namespace Dune

#endif
