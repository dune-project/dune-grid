// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_GENERICREFERENCEELEMENTS_HH
#define DUNE_GENERICREFERENCEELEMENTS_HH

#include <dune/common/forloop.hh>
#include <dune/common/polyallocator.hh>
#include <dune/common/typetraits.hh>

#include <dune/grid/genericgeometry/subtopologies.hh>
#include <dune/grid/genericgeometry/referencedomain.hh>
#include <dune/grid/genericgeometry/conversion.hh>
#include <dune/grid/genericgeometry/hybridmapping.hh>
#include <dune/grid/genericgeometry/mappingprovider.hh>

namespace Dune
{

  // forward declaration, needed for friend decl
  template< class ctype, int dim >
  class GenericReferenceElementContainer;

  // GenericReferenceElement
  // -----------------------

  /** \class GenericReferenceElement
   *  \ingroup GridGenericReferenceElements
   *  \brief This class provides access to geometric and topological
   *  properties of a reference element. This includes its type,
   *  the number of subentities, the volume, and a method for checking
   *  if a point is inside.
   *  The embedding of each subentity into the reference element is also
   *  provided.
   *
   *  A singleton of this class for a given geometry type can be accessed
   *  through the GenericReferenceElements class.

   *  \tparam ctype  field type for coordinates
   *  \tparam dim    dimension of the reference element
   *
   */
  template< class ctype, int dim >
  class GenericReferenceElement
  {
    typedef GenericReferenceElement< ctype, dim > This;

    friend class GenericReferenceElementContainer< ctype, dim >;

    // make copy constructor private
    GenericReferenceElement(const GenericReferenceElement &);

    // make empty constructor
    GenericReferenceElement() {};

    ~GenericReferenceElement ()
    {
      ForLoop< Destroy, 0, dim >::apply( mappings_, allocator_ );
    }

    class SubEntityInfo;
    template< class Topology > class CornerStorage;
    template< class Topology > struct Initialize;
    template< int codim > struct Destroy;

    struct GeometryTraits
      : public GenericGeometry::DefaultGeometryTraits< ctype, dim, dim >
    {
      typedef GenericGeometry::DefaultGeometryTraits< ctype, dim, dim > Base;

      typedef typename Base::CoordTraits CoordTraits;

      template< class Topology >
      struct Mapping
      {
        typedef GenericGeometry::CornerMapping< CoordTraits, Topology, dim, CornerStorage< Topology >, true > type;
      };

      struct Caching
      {
        static const GenericGeometry::EvaluationType evaluateJacobianTransposed = GenericGeometry::PreCompute;
        static const GenericGeometry::EvaluationType evaluateJacobianInverseTransposed = GenericGeometry::PreCompute;
        static const GenericGeometry::EvaluationType evaluateIntegrationElement = GenericGeometry::PreCompute;
        static const GenericGeometry::EvaluationType evaluateNormal = GenericGeometry::PreCompute;
      };

      typedef PolyAllocator Allocator;
    };

  public:
    /** \brief Collection of types depending on the codimension */
    template< int codim >
    struct Codim
    {
      //! type of mapping embedding a subentity into the reference element
      typedef GenericGeometry::HybridMapping< dim-codim, GeometryTraits > Mapping;
    };

  private:

    /** \brief Stores all subentities of a given codimension */
    template< int codim >
    struct MappingArray
      : public std::vector< typename Codim< codim >::Mapping * >
    {};

    /** \brief Type to store all subentities of all codimensions */
    typedef GenericGeometry::CodimTable< MappingArray, dim > MappingsTable;

    typename GeometryTraits::Allocator allocator_;
    std::vector< SubEntityInfo > info_[ dim+1 ];

    /** \brief The reference element volume */
    double volume_;
    std::vector< FieldVector< ctype, dim > > volumeNormals_;

    /** \brief Stores all subentities of all codimensions */
    MappingsTable mappings_;

  public:
    /** \brief number of subentities of codimension c
     *
     *  \param[in]  c  codimension whose size is desired
     */
    int size ( int c ) const
    {
      assert( (c >= 0) && (c <= dim) );
      return info_[ c ].size();
    }

    /** \brief number of subentities of codimension cc of subentity (i,c)
     *
     *  Denote by E the i-th subentity of codimension c of the current
     *  reference element. This method returns the number of subentities
     *  of codimension cc of the current reference element, that are also
     *  a subentity of E.
     *
     *  \param[in]  i   number of subentity E (0 <= i < size( c ))
     *  \param[in]  c   codimension of subentity E
     *  \param[in]  cc  codimension whose size is desired (c <= cc <= dim)
     */
    int size ( int i, int c, int cc ) const
    {
      assert( (c >= 0) && (c <= dim) );
      return info_[ c ][ i ].size( cc );
    }

    /** \brief obtain number of ii-th subentity with codim cc of (i,c)
     *
     *  Denote by E the i-th subentity of codimension c of the current
     *  reference element. And denote by S the ii-th subentity of codimension
     *  (cc-c) of E. Then, S is a also a subentity of codimension c of the current
     *  reference element. This method returns the number of S with respect
     *  to the current reference element.
     *
     *  \param[in]  i   number of subentity E (0 <= i < size( c ))
     *  \param[in]  c   codimension of subentity E
     *  \param[in]  ii  number of subentity S (with respect to E)
     *  \param[in]  cc  codimension of subentity S (c <= cc <= dim)
     */
    int subEntity ( int i, int c, int ii, int cc ) const
    {
      assert( (c >= 0) && (c <= dim) );
      return info_[ c ][ i ].number( ii, cc );
    }

    /** \brief position of the barycenter of entity (i,c)
     *
     *  Denote by E the i-th subentity of codimension c of the current
     *  reference element. This method returns the coordinates of
     *  the center of gravity of E within the current reference element.
     *
     *  \param[in]  i   number of subentity E (0 <= i < size( c ))
     *  \param[in]  c   codimension of subentity E
     */
    const FieldVector< ctype, dim > &position( int i, int c ) const
    {
      assert( (c >= 0) && (c <= dim) );
      return info_[ c ][ i ].position();
    }

    /** \brief check if a coordinate is in the reference element
     *
     *  This method returns true if the given local coordinate is within this
     *  reference element.
     *
     *  \param[in]  local  coordinates of the point
     */
    bool checkInside ( const FieldVector< ctype, dim > &local ) const
    {
      return checkInside< 0 >( local, 0 );
    }

    /** \brief check if a local coordinate is in the reference element of
     *         the i-th subentity E with codimension c of the current
     *         reference element.
     *
     *  Denote by E the i-th subentity of codimension codim of the current
     *  reference element. This method return true, if the given local
     *  coordinate is within the reference element for the entity E.
     *
     *  \tparam     codim  codimension of subentity E
     *
     *  \param[in]  local  coordinates of the point with respect to the
     *                     reference element of E
     *  \param[in]  i      number of subentity E (0 <= i < size( c ))
     */
    template< int codim >
    bool checkInside ( const FieldVector< ctype, dim-codim > &local, int i ) const
    {
      return mapping< codim >( i ).checkInside( local );
    }

    /** \brief map a local coordinate on subentity (i,codim) into the reference
     *         element
     *
     *  Denote by E the i-th subentity of codimension codim of the current
     *  reference element. This method maps a point within the reference
     *  element of E into the current reference element.
     *
     *  \tparam     codim  codimension of subentity E
     *
     *  \param[in]  local  coordinates of the point with respect to the reference
     *                     element of E
     *  \param[in]  i      number of subentity E (0 <= i < size( c ))
     *  \param[in]  c      codimension of subentity E
     *
     *  \note The runtime argument c is redundant and must equal codim.
     *
     *  \note This method is just an alias for
     *  \code
     *  mapping< codim >( i ).global( local );
     *  \endcode
     */
    template< int codim >
    FieldVector< ctype, dim >
    global( const FieldVector< ctype, dim-codim > &local, int i, int c ) const
    {
      if( c != codim )
        DUNE_THROW( Exception, "Local Coordinate Type does not correspond to codimension c." );
      assert( c == codim );
      return mapping< codim >( i ).global( local );
    }

    /** \brief map a local coordinate on subentity (i,codim) into the reference
     *         element
     *
     *  Denote by E the i-th subentity of codimension codim of the current
     *  reference element. This method maps a point within the reference
     *  element of E into the current reference element.
     *
     *  \tparam     codim  codimension of subentity E
     *
     *  \param[in]  local  coordinates of the point with respect to the reference
     *                     element of E
     *  \param[in]  i      number of subentity E (0 <= i < size( codim ))
     *
     *  \note This method is just an alias for
     *  \code
     *  mapping< codim >( i ).global( local );
     *  \endcode
     */
    template< int codim >
    FieldVector< ctype, dim >
    global( const FieldVector< ctype, dim-codim > &local, int i ) const
    {
      return mapping< codim >( i ).global( local );
    }

    /** \brief obtain the embedding of subentity (i,codim) into the reference
     *         element
     *
     *  Denote by E the i-th subentity of codimension codim of the current
     *  reference element. This method returns a
     *  \ref Dune::GenericGeometry::HybridMapping HybridMapping that maps
     *  the reference element of E into the current reference element.
     *
     *  This method can be used in a GenericGeometry to represent subentities
     *  of the current reference element.
     *
     *  \tparam     codim  codimension of subentity E
     *
     *  \param[in]  i      number of subentity E (0 <= i < size( codim ))
     */
    template< int codim >
    typename Codim< codim >::Mapping &mapping( int i ) const
    {
      integral_constant< int, codim > codimVariable;
      return *(mappings_[ codimVariable ][ i ]);
    }

    /** \brief obtain the type of subentity (i,c)
     *
     *  Denote by E the i-th subentity of codimension c of the current
     *  reference element. This method returns the GeometryType of E.
     *
     *  \param[in]  i      number of subentity E (0 <= i < size( c ))
     *  \param[in]  c      codimension of subentity E
     */
    const GeometryType &type ( int i, int c ) const
    {
      assert( (c >= 0) && (c <= dim) );
      return info_[ c ][ i ].type();
    }

    unsigned int topologyId ( int i, int c ) const DUNE_DEPRECATED
    {
      assert( (c >= 0) && (c <= dim) );
      return info_[ c ][ i ].topologyId();
    }

    /** \brief obtain the volume of the reference element */
    double volume () const
    {
      return volume_;
    }

    /** \brief obtain the volume outer normal of the reference element
     *
     *  The volume outer normal is the outer normal whose length coincides
     *  with the face's volume.
     *
     *  \param[in]  face  index of the face, whose normal is desired
     */
    const FieldVector< ctype, dim > &volumeOuterNormal ( int face ) const
    {
      assert( (face >= 0) && (face < int( volumeNormals_.size())) );
      return volumeNormals_[ face ];
    }

    /** \brief initialize the reference element
     *
     *  \tparam  Topology  topology of the desired reference element
     *
     *  \note The dimension of the topology must match dim.
     */
    template< class Topology >
    void initializeTopology ()
    {
      dune_static_assert( (Topology::dimension == dim),
                          "Cannot initialize reference element for different dimension." );
      typedef Initialize< Topology > Init;
      typedef GenericGeometry::VirtualMapping< Topology, GeometryTraits > VirtualMapping;

      // set up subentities
      integral_constant< int, 0 > codim0Variable;
      mappings_[ codim0Variable ].resize( 1 );
      mappings_[ codim0Variable ][ 0 ]  = allocator_.create( VirtualMapping( codim0Variable ) );

      Dune::ForLoop< Init::template Codim, 0, dim >::apply( info_, mappings_, allocator_ );

      // compute reference element volume
      typedef GenericGeometry::ReferenceDomain< Topology > ReferenceDomain;
      volume_ = ReferenceDomain::template volume< double >();

      // compute normals
      volumeNormals_.resize( ReferenceDomain::numNormals );
      for( unsigned int i = 0; i < ReferenceDomain::numNormals; ++i )
        ReferenceDomain::integrationOuterNormal( i ,volumeNormals_[ i ] );
    }
  };


  /** \brief Topological and geometric information about the subentities
   *     of a reference element
   */
  template< class ctype, int dim >
  class GenericReferenceElement< ctype, dim >::SubEntityInfo
  {
    template< class Topology, int codim > struct Initialize
    {
      template< int subcodim > struct SubCodim;
    };

    int codim_;
    std::vector< int > numbering_[ dim+1 ];
    FieldVector< ctype, dim > baryCenter_;
    GeometryType type_;

  public:
    int size ( int cc ) const
    {
      assert( (cc >= codim_) && (cc <= dim) );
      return numbering_[ cc ].size();
    }

    int number ( int ii, int cc ) const
    {
      assert( (cc >= codim_) && (cc <= dim) );
      return numbering_[ cc ][ ii ];
    }

    const FieldVector< ctype, dim > &position () const
    {
      return baryCenter_;
    }

    const GeometryType &type () const
    {
      return type_;
    }

    unsigned int topologyId () const DUNE_DEPRECATED
    {
      return type_.id();
    }

    template< class Topology, unsigned int codim, unsigned int i >
    void initialize ()
    {
      typedef Initialize< Topology, codim > Init;
      typedef GenericGeometry::ReferenceDomain< Topology > RefDomain;

      codim_ = codim;

      const unsigned int iVariable = i;
      Dune::ForLoop< Init::template SubCodim, 0, dim-codim >::apply( iVariable, numbering_ );

      baryCenter_ = ctype( 0 );
      static const unsigned int numCorners = size( dim );
      for( unsigned int j = 0; j < numCorners; ++j )
      {
        FieldVector< ctype, dim > corner;
        RefDomain::corner( number( j, dim ), corner );
        baryCenter_ += corner;
      }
      baryCenter_ *= ctype( 1 ) / ctype( numCorners );

      typedef typename GenericGeometry::SubTopology< Topology, codim, i >::type SubTopology;
      type_ = GeometryType( SubTopology::id, SubTopology::dimension );
      // type_ = GenericGeometry::DuneGeometryType< SubTopology, GeometryType::simplex >::type();
    }
  };


  template< class ctype, int dim >
  template< class Topology >
  class GenericReferenceElement< ctype, dim >::CornerStorage
  {
    typedef GenericGeometry::ReferenceDomain< Topology > RefDomain;

  public:
    static const unsigned int size = Topology::numCorners;

    template< class SubTopology >
    struct SubStorage
    {
      typedef CornerStorage< SubTopology > type;
    };

    explicit CornerStorage ( const integral_constant< int, 0 > & )
    {
      for( unsigned int i = 0; i < size; ++i )
        RefDomain::corner( i, coords_[ i ] );
    }

    template< class Mapping, unsigned int codim >
    explicit
    CornerStorage ( const GenericGeometry::SubMappingCoords< Mapping, codim > &coords )
    {
      for( unsigned int i = 0; i < size; ++i )
        coords_[ i ] = coords[ i ];
    }

    const FieldVector< ctype, dim > &operator[] ( unsigned int i ) const
    {
      return coords_[ i ];
    }

  private:
    FieldVector< ctype, dim > coords_[ size ];
  };


  template< class ctype, int dim >
  template< class Topology, int codim >
  template< int subcodim >
  struct GenericReferenceElement< ctype, dim >::SubEntityInfo::Initialize< Topology, codim >::SubCodim
  {
    typedef GenericGeometry::SubTopologySize< Topology, codim, subcodim > SubSize;
    typedef GenericGeometry::GenericSubTopologyNumbering< Topology, codim, subcodim > SubNumbering;

    static void apply ( unsigned int i, std::vector< int > (&numbering)[ dim+1 ] )
    {
      const unsigned int size = SubSize::size( i );
      numbering[ codim+subcodim ].resize( size );
      for( unsigned int j = 0; j < size; ++j )
        numbering[ codim+subcodim ][ j ] = SubNumbering::number( i, j );
    }
  };


  template< class ctype, int dim >
  template< class Topology >
  struct GenericReferenceElement< ctype, dim >::Initialize
  {
    typedef Dune::GenericReferenceElement< ctype, dim > GenericReferenceElement;

    typedef typename GenericReferenceElement::template Codim< 0 >::Mapping ReferenceMapping;

    template< int codim >
    struct Codim
    {
      template< int i >
      struct SubTopology
      {
        static void apply ( std::vector< SubEntityInfo > &info )
        {
          info[ i ].template initialize< Topology, codim, i >();
        }
      };

      static void
      apply ( std::vector< SubEntityInfo > (&info)[ dim+1 ],
              MappingsTable &mappings, typename GeometryTraits::Allocator &allocator )
      {
        const unsigned int size = GenericGeometry::Size< Topology, codim >::value;
        info[ codim ].resize( size );
        Dune::ForLoop< SubTopology, 0, size-1 >::apply( info[ codim ] );

        if( codim > 0 )
        {
          integral_constant< int, 0 > codim0Variable;
          const ReferenceMapping &refMapping = *(mappings[ codim0Variable ][ 0 ]);

          typedef typename GenericGeometry::MappingProvider< ReferenceMapping, codim > MappingProvider;

          integral_constant< int, codim > codimVariable;
          mappings[ codimVariable ].resize( size );
          for( unsigned int i = 0; i < size; ++i ) {
            char* storage = new char[MappingProvider::maxMappingSize];
            mappings[ codimVariable ][ i ] = refMapping.template trace< codim >( i, storage );
          }
        }
      }
    };
  };



  template< class ctype, int dim >
  template< int codim >
  struct GenericReferenceElement< ctype, dim >::Destroy
  {
    static void apply ( MappingsTable &mappings, typename GeometryTraits::Allocator &allocator )
    {
      if (codim > 0 )
      {
        integral_constant< int, codim > codimVariable;
        for( size_t i = 0; i < mappings[ codimVariable ].size(); ++i ) {
          typedef typename Codim<codim>::Mapping Mapping;
          mappings[ codimVariable ][ i ]->~Mapping();
          char* storage = (char*)mappings[ codimVariable ][ i ];
          delete[](storage);
        }
      }
    }
  };



  // GenericReferenceElementContainer
  // --------------------------------

  template< class ctype, int dim >
  class GenericReferenceElementContainer
  {
    static const unsigned int numTopologies = (1u << dim);

  public:
    typedef GenericReferenceElement< ctype, dim > value_type;

    const value_type &operator() ( const unsigned int topologyId ) const DUNE_DEPRECATED
    {
      return values_[ topologyId ];
    }

    const value_type &operator() ( const GeometryType &type ) const
    {
      assert( type.dim() == dim );
      return values_[ type.id() ];
    }

    const value_type &simplex () const
    {
      return values_[ GenericGeometry::SimplexTopology< dim >::type::id ];
    }

    const value_type &cube () const
    {
      return values_[ GenericGeometry::CubeTopology< dim >::type::id ];
    }

    const value_type &pyramid () const
    {
      return values_[ GenericGeometry::PyramidTopology< dim >::type::id ];
    }

    const value_type &prism () const
    {
      return values_[ GenericGeometry::PrismTopology< dim >::type::id ];
    }

    static const GenericReferenceElementContainer &instance ()
    {
      static GenericReferenceElementContainer inst;
      return inst;
    }

  private:
    template< int topologyId >
    struct Builder
    {
      static void apply ( value_type (&values)[ numTopologies ] )
      {
        typedef typename GenericGeometry::Topology< topologyId, dim >::type Topology;
        values[ topologyId ].template initializeTopology< Topology >();
      }
    };

    GenericReferenceElementContainer ()
    {
      ForLoop< Builder, 0, numTopologies-1 >::apply( values_ );
    }

    value_type values_[ numTopologies ];
  };


  // GenericReferenceElements
  // ------------------------

  /** \brief Class providing access to the singletons of the generic
   *  reference elements. Special methods are available for
   *  simplex and cube elements of any dimension.
   *  The method general can be used to obtain the reference element
   *  for a given geometry type.
   *
   *  \ingroup GridGenericReferenceElements
   */
  template< class ctype, int dim >
  struct GenericReferenceElements
  {
    //! get general generic reference elements
    static const GenericReferenceElement< ctype, dim > &
    general ( const GeometryType &type )
    {
      return GenericReferenceElementContainer< ctype, dim >::instance() ( type );
    }

    //! get simplex generic reference elements
    static const GenericReferenceElement< ctype, dim > &simplex ()
    {
      return GenericReferenceElementContainer< ctype, dim >::instance().simplex();
    }

    //! get hypercube generic reference elements
    static const GenericReferenceElement< ctype, dim > &cube ()
    {
      return GenericReferenceElementContainer< ctype, dim >::instance().cube();
    }
  };

}

#endif // #ifndef DUNE_GENERICREFERENCEELEMENTS_HH
