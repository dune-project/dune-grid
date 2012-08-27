// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALU3DGRID_ENTITYSEED_HH
#define DUNE_ALU3DGRID_ENTITYSEED_HH

#include <dune/common/nullptr.hh>

namespace Dune
{

  // Internal Forward Declarations
  // -----------------------------

  template< int codim, class GridImp >
  class ALU3dGridEntitySeed;



  // ALU3dGridEntitySeedBase
  // -----------------------

  template< int codim, class GridImp >
  class ALU3dGridEntitySeedBase
  {
    typedef ALU3dGridEntitySeedBase< codim, GridImp > This;

  protected:
    enum { dim       = GridImp::dimension };
    enum { dimworld  = GridImp::dimensionworld };


    typedef typename GridImp::MPICommunicatorType Comm;

    friend class ALU3dGridEntity< codim, dim, GridImp >;
    friend class ALU3dGridEntity< 0,dim,GridImp>;
    friend class ALU3dGrid < GridImp::elementType, Comm >;

    typedef ALU3dImplTraits<GridImp::elementType, Comm > ImplTraits;
    typedef typename ImplTraits::template Codim<codim>::ImplementationType ImplementationType;
    typedef typename ImplTraits::template Codim<codim>::InterfaceType HElementType;
    typedef typename ImplTraits::template Codim<codim>::EntitySeedType KeyType ;

    typedef typename ImplTraits::BNDFaceType BNDFaceType;
    typedef typename ImplTraits::HBndSegType HBndSegType;

    template< int cd, class Key >
    struct Bnd
    {
      static Key *toKey ( const HBndSegType * ) { return nullptr; }

      static HElementType *getItem ( KeyType *key )
      {
        return static_cast< HElementType* >( key );
      }

      static bool isGhost ( KeyType * ) { return false; }

      static BNDFaceType *ghost ( KeyType*  ) { return nullptr; }
    };

    template< class Key >
    struct Bnd< 0, Key >
    {
      static Key *toKey ( const HBndSegType *ghostFace )
      {
        return static_cast< KeyType* >( const_cast< BNDFaceType * >( static_cast< const BNDFaceType * >( ghostFace ) ) );
      }

      static HElementType *getItem ( KeyType *key )
      {
        if( key )
        {
          if( key->isboundary() )
            return static_cast< BNDFaceType * >( key )->getGhost().first;
          else
          {
            // we cannot cast to HElement here, since only the implementation is derived
            // from hasFace
            return static_cast< HElementType * >( static_cast< ImplementationType * >( key ) );
          }
        }
        else
          return nullptr;
      }

      static bool isGhost( KeyType *key ) { assert( key ); return key->isboundary(); }

      static BNDFaceType *ghost( KeyType* key ) { assert( key ); return static_cast< BNDFaceType * >( key ); }
    };

  public:
    static const int defaultValue = -1 ;
    static const int defaultTwist = 0 ;

    enum { codimension = codim };

    //! type of Entity
    typedef typename GridImp::template Codim<codimension>::Entity Entity;
    //! underlying EntityImplementation
    typedef MakeableInterfaceObject<Entity> EntityObject;
    typedef typename EntityObject :: ImplementationType EntityImp;

    //! typedef of my type
    typedef This ALU3dGridEntitySeedType;

    //! make type of entity pointer implementation available in derived classes
    typedef ALU3dGridEntitySeed<codimension,GridImp> EntitySeedImp;

    //! Destructor
    ~ALU3dGridEntitySeedBase()
    {
#ifndef NDEBUG
      // clear pointer
      clear();
#endif
    }

    //! Constructor for EntitySeed that points to an element
    ALU3dGridEntitySeedBase();

    //! Constructor for EntitySeed that points to an element
    ALU3dGridEntitySeedBase(const HElementType& item);

    //! Constructor for EntitySeed that points to an element
    ALU3dGridEntitySeedBase(const HElementType* item, const HBndSegType* ghostFace );

    //! Constructor for EntitySeed that points to an ghost
    ALU3dGridEntitySeedBase(const HBndSegType& ghostFace );

    /////////////////////////////////////////////////////////////
    //
    //  interface methods
    //
    /////////////////////////////////////////////////////////////
    //! copy constructor
    ALU3dGridEntitySeedBase(const ALU3dGridEntitySeedType & org);

    //! equality operator
    bool operator== ( const This &other ) const { return equals( other ); }

    //! inequality operator
    bool operator!= ( const This &other ) const { return !equals( other ); }

    //! assignment operator
    const This &operator= ( const This &org );

    //////////////////////////////////////////////////////
    //
    //  non-interface methods
    //
    //////////////////////////////////////////////////////
    //! equality
    bool equals ( const This &other ) const;

    //! invalidate seed
    void clear () { item_ = 0; }

    //! get item from key
    HElementType *item () const { return Bnd<codim,KeyType>::getItem( item_ ); }

    //! return iterior item
    HElementType *interior () const
    {
      assert( !isGhost() );
      return static_cast< HElementType * >( static_cast< ImplementationType * >( item_ ) );
    }

    // methods for ghosts

    bool isGhost () const { return Bnd< codim, KeyType >::isGhost( item_ ); }

    BNDFaceType *ghost () const
    {
      assert( isGhost() );
      return Bnd< codim, KeyType >::ghost( item_ );
    }

    KeyType *toKey ( const HElementType *item )
    {
      return static_cast< KeyType * >( const_cast< ImplementationType * >( static_cast< const ImplementationType * >( item ) ) );
    }

    void set ( const HElementType &item ) { item_ = toKey( &item ); }

    KeyType *toKey ( const HBndSegType *ghostFace )
    {
      return Bnd< codim, KeyType >::toKey( ghostFace );
    }

    void set ( const HBndSegType &ghostFace ) { item_ = toKey( &ghostFace ); }

    int level () const { return defaultValue; }
    int twist () const { return defaultTwist; }
    int face () const { return defaultValue; }

  protected:
    // pointer to item
    KeyType *item_;
  };



  // ALU3dGridEntitySeed
  // -------------------

  template< int codim, class GridImp >
  class ALU3dGridEntitySeed
    : public ALU3dGridEntitySeedBase< codim, GridImp >
  {
    typedef ALU3dGridEntitySeed< codim, GridImp > This;
    typedef ALU3dGridEntitySeedBase< codim, GridImp > Base;

    enum { dim       = GridImp::dimension };
    enum { dimworld  = GridImp::dimensionworld };

    typedef typename GridImp::MPICommunicatorType Comm;

    friend class ALU3dGridEntity<codim,dim,GridImp>;
    friend class ALU3dGridEntity< 0,dim,GridImp>;
    friend class ALU3dGrid < GridImp::elementType, Comm >;

    typedef ALU3dImplTraits< GridImp::elementType, Comm > ImplTraits;
    typedef typename ImplTraits::template Codim<codim>::ImplementationType ImplementationType;
    typedef typename ImplTraits::template Codim<codim>::InterfaceType HElementType;

    typedef typename ImplTraits::BNDFaceType BNDFaceType;
    typedef ALU3dGridEntity<codim,dim,GridImp> ALU3dGridEntityType;

  public:
    using Base::defaultValue;
    using Base::defaultTwist;

    //! type of Entity
    typedef typename GridImp::template Codim<codim>::Entity Entity;

    //! typedef of my type
    typedef ALU3dGridEntitySeed<codim,GridImp> ALU3dGridEntitySeedType;

    //! Constructor for EntitySeed that points to an element
    ALU3dGridEntitySeed(const ImplementationType & item)
    {
      // this constructor should only be called by codim=0 entity keys
      assert( false );
      abort();
    }

    //! Constructor for EntitySeed that points to an element
    ALU3dGridEntitySeed(const HElementType & item,
                        const int level,
                        const int twist = defaultTwist,
                        const int duneFace = defaultValue
                        );

    //! Constructor for EntitySeed that points to an element
    ALU3dGridEntitySeed ()
      : Base(),
        level_( defaultValue ),
        twist_( defaultTwist ),
        face_( defaultValue )
    {}

    //! Constructor for EntitySeed that points to given entity
    ALU3dGridEntitySeed ( const ALU3dGridEntityType &entity )
      : Base( entity.getItem() ),
        level_( entity.level() ),
        twist_( defaultTwist ),
        face_( defaultValue )
    {}

    //! copy constructor
    ALU3dGridEntitySeed ( const This &org );

    //! assignment operator
    const This &operator= ( const This &org );

    //! clear the key data structure
    void clear ();

    //! set element and level
    void set ( const HElementType &item, const int level )
    {
      Base::set( item );
      level_ = level;
    }

    //! return level
    int level () const { return level_ ; }
    //! return twist
    int twist () const { return twist_ ; }
    //! return face
    int face  () const { return face_ ; }

    using Base::set;

    bool operator== ( const This &other ) const { return equals( other ); }
    bool operator!= ( const This &other ) const { return !equals( other ); }

    //! equality, calls BaseType equals
    bool equals ( const This &other ) const
    {
      return Base::equals( other ) && (level() == other.level());
    }

  protected:
    // level of entity
    int level_;
    // twist of face, for codim 1 only
    int twist_;
    // face number, for codim 1 only
    int face_;
  };



  // ALU3dGridEntitySeed (for codim = 0)
  // -----------------------------------

  //! ALUGridEntitySeed points to an entity
  //! this class is the specialisation for codim 0,
  //! it has exactly the same functionality as the ALU3dGridEntitySeedBase
  template<class GridImp>
  class ALU3dGridEntitySeed< 0, GridImp >
    : public ALU3dGridEntitySeedBase< 0, GridImp >
  {
    typedef ALU3dGridEntitySeed< 0, GridImp > This;
    typedef ALU3dGridEntitySeedBase< 0, GridImp > Base;

  protected:
    enum { codim = 0 };
    enum { dim       = GridImp::dimension };
    enum { dimworld  = GridImp::dimensionworld };

    typedef typename GridImp::MPICommunicatorType Comm;

    friend class ALU3dGridEntity<codim,dim,GridImp>;
    friend class ALU3dGridEntity< 0,dim,GridImp>;
    friend class ALU3dGrid < GridImp::elementType, Comm >;

    typedef ALU3dImplTraits<GridImp::elementType, Comm > ImplTraits;
    typedef typename ImplTraits::template Codim<codim>::ImplementationType ImplementationType;
    typedef typename ImplTraits::template Codim<codim>::InterfaceType HElementType;

    typedef typename ImplTraits::BNDFaceType BNDFaceType;
    typedef typename ImplTraits::HBndSegType HBndSegType;

    typedef ALU3dGridEntity< 0, dim, GridImp > ALU3dGridEntityType;

  public:
    using Base::defaultValue;
    using Base::defaultTwist;

    //! type of Entity
    typedef typename GridImp::template Codim< codim >::Entity Entity;

    //! typedef of my type
    typedef This ALU3dGridEntitySeedType;

    //! Constructor for EntitySeed that points to an element
    ALU3dGridEntitySeed ()
      : Base()
    {}

    //! Constructor for EntitySeed that points to an interior element
    ALU3dGridEntitySeed ( const HElementType &item )
      : Base( item )
    {}

    //! Constructor for EntitySeed that points to an interior element
    ALU3dGridEntitySeed ( const HElementType &item, int, int, int )
      : Base( item )
    {}

    //! Constructor for EntitySeed that points to an ghost
    ALU3dGridEntitySeed ( const HBndSegType &ghostFace )
      : Base( ghostFace )
    {}

    //! copy constructor
    ALU3dGridEntitySeed ( const This &org )
      : Base( org )
    {}
  };


  //! print alugrid entity key to std::stream
  template< int codim, class GridImp >
  inline std::ostream &
  operator<< ( std::ostream &out, const ALU3dGridEntitySeed< codim, GridImp > &seed )
  {
    out << seed.item() << " " << seed.level() << " " << seed.twist() << " " << seed.face();
    return out;
  }


  // Implementation of ALU3dGridEntitySeedBase
  // -----------------------------------------

  template<int codim, class GridImp >
  inline ALU3dGridEntitySeedBase< codim, GridImp >::ALU3dGridEntitySeedBase ()
    : item_( 0 )
  {}


  template<int codim, class GridImp >
  inline ALU3dGridEntitySeedBase< codim, GridImp >
  ::ALU3dGridEntitySeedBase ( const HElementType &item )
    : item_( toKey( &item ) )
  {}


  template<int codim, class GridImp >
  inline ALU3dGridEntitySeedBase< codim, GridImp >
  ::ALU3dGridEntitySeedBase ( const HBndSegType &ghostFace )
    : item_( toKey( &ghostFace ) )
  {}


  template<int codim, class GridImp >
  inline ALU3dGridEntitySeedBase< codim, GridImp >
  ::ALU3dGridEntitySeedBase ( const This &org )
    : item_( org.item_ )
  {}


  template<int codim, class GridImp >
  inline const typename ALU3dGridEntitySeedBase< codim, GridImp >::This &
  ALU3dGridEntitySeedBase< codim, GridImp >::operator= ( const This &org )
  {
    item_  = org.item_;
    return *this;
  }


  template<int codim, class GridImp >
  inline bool ALU3dGridEntitySeedBase< codim, GridImp >
  ::equals ( const ALU3dGridEntitySeedBase< codim, GridImp > &other ) const
  {
    // check equality of underlying items
    return (item_ == other.item_);
  }



  // Implementation of ALU3dGridEntity
  // ---------------------------------

  template< int codim, class GridImp >
  inline ALU3dGridEntitySeed< codim, GridImp >
  ::ALU3dGridEntitySeed ( const HElementType &item,
                          const int level,
                          const int twist,
                          const int duneFace )
    : Base( item ),
      level_( level ),
      twist_( twist ),
      face_( duneFace ) // duneFace can be -1 when face was created by face iterator
  {}


  template< int codim, class GridImp >
  inline ALU3dGridEntitySeed< codim, GridImp >
  ::ALU3dGridEntitySeed ( const This &org )
    : Base( org ),
      level_( org.level_ ),
      twist_( org.twist_ ),
      face_( org.face_ )
  {}


  template< int codim, class GridImp >
  inline const typename ALU3dGridEntitySeed< codim, GridImp >::This &
  ALU3dGridEntitySeed< codim, GridImp >::operator= ( const This &org )
  {
    Base::operator=( org );
    level_ = org.level_;
    twist_ = org.twist_;
    face_ = org.face_;
    return *this;
  }


  template< int codim, class GridImp >
  inline void ALU3dGridEntitySeed< codim, GridImp >::clear ()
  {
    Base::clear();
    level_ = defaultValue;
    twist_ = defaultTwist;
    face_ = defaultValue;
  }

} // namespace Dune

#endif // #ifndef DUNE_ALU3DGRID_ENTITYSEED_HH
