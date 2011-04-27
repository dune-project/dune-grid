// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef ALU2DGRID_ENTITYKEY_HH
#define ALU2DGRID_ENTITYKEY_HH

namespace Dune
{

  template<int cd, class GridImp>
  class ALU2dGridEntitySeed ;

  //**********************************************************************
  //
  // --ALU2dGridEntitySeed
  // --EntitySeed
  //**********************************************************************
  template< int codim, class GridImp >
  class ALU2dGridEntitySeedBase
  {
  protected:
    typedef ALU2dGridEntitySeedBase< codim, GridImp > ThisType;
    enum { dim       = GridImp::dimension };
    enum { dimworld  = GridImp::dimensionworld };


    friend class ALU2dGridEntity<codim,dim,GridImp>;
    friend class ALU2dGridEntity< 0,dim,GridImp>;
    friend class ALU2dGrid <GridImp :: dimension, GridImp::dimensionworld, GridImp::elementType >;

    typedef ALU2dImplTraits< GridImp::dimensionworld, GridImp::elementType > ImplTraits;
    typedef typename ImplTraits::template Codim<codim>::InterfaceType ImplementationType;
    typedef ImplementationType HElementType;
    typedef ImplementationType KeyType;

  public:
    static const int defaultValue = -1 ;

    enum { codimension = codim };

    //! type of Entity
    typedef typename GridImp::template Codim<codimension>::Entity Entity;
    //! underlying EntityImplementation
    typedef MakeableInterfaceObject<Entity> EntityObject;
    typedef typename EntityObject :: ImplementationType EntityImp;

    //! typedef of my type
    typedef ThisType ALU2dGridEntitySeedType;

    //! make type of entity pointer implementation available in derived classes
    typedef ALU2dGridEntitySeed<codimension,GridImp> EntitySeedImp;

    //! Destructor
    ~ALU2dGridEntitySeedBase()
    {
#ifndef NDEBUG
      // clear pointer
      clear();
#endif
    }

    //! Constructor for EntitySeed that points to an element
    ALU2dGridEntitySeedBase();

    //! Constructor for EntitySeed that points to an element
    ALU2dGridEntitySeedBase(const HElementType& item);

    /////////////////////////////////////////////////////////////
    //
    //  interface methods
    //
    /////////////////////////////////////////////////////////////
    //! copy constructor
    ALU2dGridEntitySeedBase(const ALU2dGridEntitySeedType & org);

    //! equality operator
    bool operator == (const ALU2dGridEntitySeedType& i) const
    {
      return equals( i );
    }

    //! inequality operator
    bool operator != (const ALU2dGridEntitySeedType& i) const
    {
      return ! equals( i );
    }

    //! assignment operator
    ThisType & operator = (const ThisType & org);

    //////////////////////////////////////////////////////
    //
    //  non-interface methods
    //
    //////////////////////////////////////////////////////
    //! equality
    bool equals (const ALU2dGridEntitySeedType& i) const;

    //! invalidate seed
    void clear()
    {
      item_ = 0;
    }

    //! get item from key
    HElementType* item() const { return item_; }

    KeyType* toKey(const HElementType* item)
    {
      return static_cast< KeyType* > (const_cast< ImplementationType* > (static_cast<const ImplementationType* > (item)));
    }

    void set(const HElementType& item, const int level = -1 , const int face = -1 )
    {
      item_ = toKey( &item );
    }

    int level () const { return ( item_ ) ? item_->level() : defaultValue; }
    int face  () const { return defaultValue; }

  protected:
    // pointer to item
    mutable KeyType* item_;
  };

  template<int cd, class GridImp>
  class ALU2dGridEntitySeed :
    public ALU2dGridEntitySeedBase<cd,GridImp>
  {
    typedef ALU2dGridEntitySeedBase<cd,GridImp> BaseType;

    typedef ALU2dGridEntitySeed <cd,GridImp> ThisType;
    enum { dim       = GridImp::dimension };
    enum { dimworld  = GridImp::dimensionworld };

    friend class ALU2dGridEntity<cd,dim,GridImp>;
    friend class ALU2dGridEntity< 0,dim,GridImp>;
    friend class ALU2dGrid <GridImp :: dimension, GridImp::dimensionworld, GridImp::elementType >;

    typedef ALU2dImplTraits< GridImp::dimensionworld, GridImp::elementType > ImplTraits;
    typedef typename ImplTraits::template Codim<cd>::InterfaceType ImplementationType;
    typedef ImplementationType HElementType;

    typedef ALU2dGridEntity<cd,dim,GridImp> ALU2dGridEntityType;

  public:
    using BaseType :: defaultValue ;

    //! type of Entity
    typedef typename GridImp::template Codim<cd>::Entity Entity;

    //! typedef of my type
    typedef ALU2dGridEntitySeed<cd,GridImp> ALU2dGridEntitySeedType;

    //! Constructor for EntitySeed that points to an element
    ALU2dGridEntitySeed(const ImplementationType & item)
    {
      // this constructor should only be called by codim=0 entity keys
      assert( false );
      abort();
    }

    //! Constructor for EntitySeed that points to an element
    ALU2dGridEntitySeed(const HElementType & item,
                        const int level,
                        const int duneFace = defaultValue
                        );

    //! Constructor for EntitySeed that points to an element
    ALU2dGridEntitySeed()
      : BaseType(), level_(defaultValue), face_(defaultValue) {}

    //! Constructor for EntitySeed that points to given entity
    ALU2dGridEntitySeed(const ALU2dGridEntityType& entity)
      : ALU2dGridEntitySeedBase<cd,GridImp> (entity.getItem()),
        level_(entity.level()), face_(defaultValue)
    {}

    //! copy constructor
    ALU2dGridEntitySeed(const ALU2dGridEntitySeedType & org);

    //! assignment operator
    ThisType & operator = (const ThisType & org);

    //! clear the key data structure
    void clear();

    //! set element and level
    void set(const HElementType & item, const int level, const int duneFace )
    {
      BaseType :: set( item );
      level_ = level ;
      face_  = duneFace ;
    }

    //! return level
    int level () const { return level_ ; }
    //! return face
    int face  () const { return face_ ; }

    using BaseType :: set ;

    bool operator == (const ALU2dGridEntitySeedType& i) const
    {
      return equals( i );
    }

    bool operator != (const ALU2dGridEntitySeedType& i) const
    {
      return ! equals( i );
    }

    //! equality, calls BaseType equals
    bool equals (const ALU2dGridEntitySeedType& key) const
    {
      // only compare the item pointer, this is the real key
      return BaseType :: equals( key ) && (level() == key.level());
    }

  protected:
    // level of entity
    int level_;
    // face number, for codim 1 only
    int face_;
  };

  //! ALUGridEntitySeed points to an entity
  //! this class is the specialisation for codim 0,
  //! it has exactly the same functionality as the ALU2dGridEntitySeedBase
  template<class GridImp>
  class ALU2dGridEntitySeed<0,GridImp> :
    public ALU2dGridEntitySeedBase<0,GridImp>
  {
  protected:
    typedef ALU2dGridEntitySeedBase<0,GridImp> BaseType;

    enum { cd = 0 };
    typedef ALU2dGridEntitySeed <cd,GridImp> ThisType;
    enum { dim       = GridImp::dimension };
    enum { dimworld  = GridImp::dimensionworld };

    friend class ALU2dGridEntity<cd,dim,GridImp>;
    friend class ALU2dGridEntity< 0,dim,GridImp>;
    friend class ALU2dGrid <GridImp :: dimension, GridImp::dimensionworld, GridImp::elementType >;

    typedef ALU2dImplTraits< GridImp::dimensionworld, GridImp::elementType > ImplTraits;
    typedef typename ImplTraits::template Codim<cd>::InterfaceType ImplementationType;
    typedef ImplementationType HElementType;

    typedef ALU2dGridEntity< 0,dim,GridImp> ALU2dGridEntityType ;

  public:
    using BaseType :: defaultValue ;

    //! type of Entity
    typedef typename GridImp::template Codim<cd>::Entity Entity;

    //! typedef of my type
    typedef ThisType ALU2dGridEntitySeedType;

    //! Constructor for EntitySeed that points to an element
    ALU2dGridEntitySeed() : BaseType() {}

    //! Constructor for EntitySeed that points to an interior element
    ALU2dGridEntitySeed(const HElementType& item)
      : ALU2dGridEntitySeedBase<cd,GridImp> (item) {}

    //! Constructor for EntitySeed that points to an interior element
    ALU2dGridEntitySeed(const HElementType& item, int , int )
      : ALU2dGridEntitySeedBase<cd,GridImp> (item) {}

    //! copy constructor
    ALU2dGridEntitySeed(const ALU2dGridEntitySeedType & org)
      : ALU2dGridEntitySeedBase<cd,GridImp> (org)
    {}
  };


  //! print alugrid entity key to std::stream
  template <int cd, class GridImp>
  inline std :: ostream &operator<< ( std :: ostream &out,
                                      const ALU2dGridEntitySeed<cd,GridImp>& key)
  {
    out << key.item() << " " << key.level() << " " << key.face();
    return out;
  }


  //*******************************************************************
  //
  //  Implementation
  //
  //*******************************************************************
  template<int codim, class GridImp >
  inline ALU2dGridEntitySeedBase<codim,GridImp> ::
  ALU2dGridEntitySeedBase()
    : item_( 0 )
  {}

  template<int codim, class GridImp >
  inline ALU2dGridEntitySeedBase<codim,GridImp> ::
  ALU2dGridEntitySeedBase(const HElementType &item)
    : item_( toKey(&item) )
  {}

  template<int codim, class GridImp >
  inline ALU2dGridEntitySeedBase<codim,GridImp> ::
  ALU2dGridEntitySeedBase(const ALU2dGridEntitySeedType & org)
    : item_(org.item_)
  {}

  template<int codim, class GridImp >
  inline ALU2dGridEntitySeedBase<codim,GridImp> &
  ALU2dGridEntitySeedBase<codim,GridImp> ::
  operator = (const ALU2dGridEntitySeedType & org)
  {
    item_  = org.item_;
    return *this;
  }

  template<int codim, class GridImp >
  inline bool ALU2dGridEntitySeedBase<codim,GridImp>::
  equals (const ALU2dGridEntitySeedBase<codim,GridImp>& i) const
  {
    // check equality of underlying items
    return (item_ == i.item_);
  }

  ///////////////////////////////////////////////////////////////////
  //
  //  specialisation for higher codims
  //
  ///////////////////////////////////////////////////////////////////

  template<int codim, class GridImp >
  inline ALU2dGridEntitySeed<codim,GridImp> ::
  ALU2dGridEntitySeed(const HElementType &item,
                      const int level,
                      const int duneFace )
    : ALU2dGridEntitySeedBase<codim,GridImp> (item)
      , level_(level)
      , face_(duneFace)
  {
    assert( (codim == 1) ? (face_ >= 0) : 1 );
  }

  template<int codim, class GridImp >
  inline ALU2dGridEntitySeed<codim,GridImp> ::
  ALU2dGridEntitySeed(const ALU2dGridEntitySeedType & org)
    : ALU2dGridEntitySeedBase<codim,GridImp>(org)
      , level_(org.level_)
      , face_(org.face_)
  {}

  template<int codim, class GridImp >
  inline ALU2dGridEntitySeed<codim,GridImp> &
  ALU2dGridEntitySeed<codim,GridImp>::
  operator = (const ALU2dGridEntitySeedType & org)
  {
    // docu and cleanup
    BaseType :: operator = ( org );

    // clone other stuff
    level_ = org.level_;
    face_  = org.face_;
    return *this;
  }

  template<int codim, class GridImp >
  inline void
  ALU2dGridEntitySeed<codim,GridImp>::clear ()
  {
    BaseType :: clear();
    level_ = defaultValue ;
    face_  = defaultValue ;
  }

} // end namespace Dune
#endif
