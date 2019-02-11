// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_COMMON_BOUNDARY_SEGMENT_HH
#define DUNE_GRID_COMMON_BOUNDARY_SEGMENT_HH

#include <map>
#include <sstream>

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
  protected:
    typedef BndSeg  BoundarySegment;
    static const unsigned int keyLength = 4 ;
  public:
    typedef BoundarySegment* BoundarySegmentFactoryType( std::stringstream& in );

    //typedef std::type_info  KeyType;
    typedef std::string KeyType;

    typedef std::map< KeyType, BoundarySegmentFactoryType* > FactoryStorage;

  protected:
    /** \brief create an object of BoundarySegment type from a previously
     *         registered factory linked to key.
     *  \param in stream buffer previously written with backup containing key and object data
     *
     *  \return Object derived from BoundarySegment.
     */
    static BoundarySegment* restore( std::stringstream& in )
    {
      FactoryStorage& fac = storage();
      char key[ keyLength ];
      in.read( key, keyLength );

      auto it = fac.find( std::string(key) );
      if( it == fac.end() )
      {
        DUNE_THROW(NotImplemented,"Requested BoundarySegment type " << std::string(key) << " not previously registered");
      }

      // call factory routine
      return it->second( in );
    }

    static void registerFactory( const KeyType& key, BoundarySegmentFactoryType* factory )
    {
      if( key.size() != keyLength )
      {
        DUNE_THROW(InvalidStateException,"Key for BoundarySegment type is expected to have length " << keyLength);
      }

      FactoryStorage& fac = storage();
      auto it = fac.find( key );
      if( it == fac.end() )
      {
        fac.insert( std::make_pair( key, factory ) );
      }
    }

  private:
    static FactoryStorage& storage()
    {
      static FactoryStorage s;
      return s;
    }
  };

  template< int dim, int dimworld, class ctype >
  struct BoundarySegment : public BoundarySegmentBackupRestore< BoundarySegment< dim, dimworld, ctype > >
  {
    typedef BoundarySegment< dim, dimworld, ctype > ThisType;
    typedef BoundarySegmentBackupRestore< BoundarySegment< dim, dimworld, ctype > > BaseType;

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
    virtual void backup( std::stringstream& buffer ) const
    {
      DUNE_THROW(NotImplemented,"BoundarySegment::backup needs to be overloaded!");
    }
  };


}  // end namespace Dune

#endif
