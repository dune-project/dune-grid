// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALUGRIDGEOMETRYSTORAGE_HH
#define DUNE_ALUGRIDGEOMETRYSTORAGE_HH

// Dune includes
#include <dune/grid/common/grid.hh>
#include <dune/grid/common/gridfactory.hh>

#if HAVE_ALUGRID

#include <dune/grid/alugrid/common/declaration.hh>
#include <dune/grid/alugrid/3d/alu3dinclude.hh>
#include <dune/grid/alugrid/2d/alu2dinclude.hh>

namespace Dune
{
  template<int dim, int dimw>
  class ALUCubeGrid;

  template<int dim, int dimw>
  class ALUSimplexGrid;

  template<int dim, int dimw>
  class ALUConformGrid;

  template< class GridImp, class GeometryImpl, int nChild >
  class ALULocalGeometryStorage
  {
    typedef ALULocalGeometryStorage< GridImp, GeometryImpl, nChild > ThisType;

    // array with pointers to the geometries
    Dune::array< GeometryImpl *, nChild > geoms_;

    // count local geometry creation
    int count_;

    // type of grid impl
    typedef typename GridImp :: ctype ctype;
    enum { dimension       = GridImp :: dimension };
    enum { dimensionworld  = GridImp :: dimensionworld };

    template <int dummy, int dim, int dimworld, int >
    struct CreateGeometries;

    template <int dummy, int dimworld>
    struct CreateGeometries<dummy, 2, dimworld, ALU2DSPACE triangle >
    {
      template <class Storage>
      static void createGeometries(Storage& storage,
                                   const GeometryType& type,
                                   const bool nonConform )
      {
        if( nonConform )
        {
          typedef ALUGrid< 2, dimworld, simplex, nonconforming, No_Comm > Grid;
          storage.template createGeometries< Grid > (type);
        }
        else
        {
          typedef ALUGrid< 2, dimworld, simplex, conforming, No_Comm > Grid;
          storage.template createGeometries< Grid > (type);
        }
      }
    };

    template <int dummy>
    struct CreateGeometries<dummy, 3, 3, ALU3DSPACE tetra >
    {
      template <class Storage>
      static void createGeometries(Storage& storage,
                                   const GeometryType& type,
                                   const bool nonConform )
      {
        assert ( nonConform ) ;
        {
          typedef ALUGrid< 3, 3, simplex, nonconforming, No_Comm > Grid;
          storage.template createGeometries< Grid > (type);
        }
        /*
           // TODO, implement this for refinement of all edges (conforming)
           else
           {
           typedef ALUGrid< 3, 3, simplex, conforming, No_Comm > Grid;
           storage.template createGeometries< Grid > (type);
           }
         */
      }
    };

    template <int dummy, int dimworld>
    struct CreateGeometries<dummy, 2, dimworld, ALU2DSPACE quadrilateral >
    {
      template <class Storage>
      static void createGeometries(Storage& storage,
                                   const GeometryType& type,
                                   const bool nonConform )
      {
        assert ( nonConform ) ;
        {
          typedef ALUGrid< 2, dimworld, cube, nonconforming, No_Comm > Grid;
          storage.template createGeometries< Grid > (type);
        }
      }
    };

    template <int dummy>
    struct CreateGeometries<dummy, 3, 3, ALU3DSPACE hexa >
    {
      template <class Storage>
      static void createGeometries(Storage& storage,
                                   const GeometryType& type,
                                   const bool nonConform )
      {
        assert( nonConform );
        {
          typedef ALUGrid< 3, 3, cube, nonconforming, No_Comm > Grid;
          storage.template createGeometries< Grid > (type);
        }
      }
    };

  public:
    // create empty storage
    ALULocalGeometryStorage ( const GeometryType type, const bool nonConform )
      : count_( 0 )
    {
      geoms_.fill( (GeometryImpl *) 0 );
      // the idea is to create a grid containing the reference element,
      // refine once and the store the father - child relations
      CreateGeometries<0, dimension, dimensionworld, GridImp :: elementType >
      ::createGeometries(*this, type, nonConform);
    }

    // check if geometry has been created
    bool geomCreated(int child) const { return geoms_[child] != 0; }

    // return reference to local geometry
    const GeometryImpl & operator [] (int child) const
    {
      assert( geomCreated(child) );
      return *(geoms_[child]);
    }

    template < class Grid >
    void createGeometries(const GeometryType& type)
    {
      // create factory without verbosity
      GridFactory< Grid > factory( false );

      const Dune::ReferenceElement< ctype, dimension > &refElem
        = Dune::ReferenceElements< ctype, dimension >::general( type );

      // insert vertices
      FieldVector<ctype, dimensionworld> pos( 0 );
      const int vxSize = refElem.size(dimension);
      for(int i=0; i<vxSize; ++i)
      {
        FieldVector<ctype, dimension> position = refElem.position(i, dimension );
        // copy position
        for(int d = 0; d<dimension; ++d )
          pos[ d ] = position[ d ];

        factory.insertVertex( pos );
      }

      std::vector< unsigned int > vertices( vxSize );
      // create grid with reference element
      for(size_t i=0; i<vertices.size(); ++i) vertices[ i ] = i;
      factory.insertElement(type, vertices);

      // save original sbuf
      std::streambuf* cerr_sbuf = std::cerr.rdbuf();
      std::stringstream tempout;
      // redirect 'cerr' to a 'fout' to avoid unnecessary output in constructors
      std::cerr.rdbuf(tempout.rdbuf());

      Grid* gridPtr = factory.createGrid();
      Grid& grid    = *gridPtr;

      // restore the original stream buffer
      std::cerr.rdbuf(cerr_sbuf);

      //std::cerr = savecerr;

      // refine once to get children
      const int level = 1;
      grid.globalRefine( level );

      {
        typedef typename Grid :: template Partition< All_Partition >:: LevelGridView MacroGridView;
        MacroGridView macroView = grid.template levelView< All_Partition > ( 0 );
        typedef typename MacroGridView :: template Codim< 0 > :: Iterator Iterator;

        Iterator it = macroView.template begin<0> ();

        if( it == macroView.template end<0>() )
          DUNE_THROW(InvalidStateException,"Empty Grid, should contain at least 1 element");

        typedef typename Iterator :: Entity EntityType;

        const EntityType& entity = *it;
        const typename EntityType :: Geometry& geo = entity.geometry();
        typedef typename EntityType :: HierarchicIterator HierarchicIteratorType;
        const HierarchicIteratorType end = entity.hend( level );

        int childNum = 0;
        for( HierarchicIteratorType child = entity.hbegin( level );
             child != end; ++child, ++childNum )
        {
          create( geo, child->geometry(), childNum );
        }
      }

      // delete grid
      delete gridPtr;
    }

    // create local geometry
    template< class Geometry >
    void create ( const Geometry &father,
                  const Geometry &son,
                  const int child )
    {
      assert( !geomCreated( child ) );
      assert( (child >= 0) && (child < nChild) );

      assert( (count_ < nChild) );
      ++count_;

      geoms_[ child ] = new GeometryImpl();
      geoms_[ child ]->buildGeomInFather( father, son );
    }

  public:
    // desctructor deleteing geometries
    ~ALULocalGeometryStorage ()
    {
      for(size_t i=0; i<geoms_.size(); ++i)
        if(geoms_[i]) delete geoms_[i];
    }

    //! access local geometries
    static const GeometryImpl& geom( const GeometryType type, const bool nonConforming, const int child )
    {
      // create static variable on heap
      static ThisType instance( type, nonConforming );
      // make sure the geometry type is the same
      assert( type == instance[ child ].type() );
      return instance[ child ];
    }
  };

} // namespace Dune

#endif // end HAVE_ALUGRID

#endif // #ifndef DUNE_ALUGRIDGEOMETRYSTORAGE_HH
