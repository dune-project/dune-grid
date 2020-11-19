#ifndef DUNE_GRID_SINGLETON_HH
#define DUNE_GRID_SINGLETON_HH

//- System includes
#include <cassert>
#include <iostream>
#include <memory>
#include <typeindex>
#include <utility>
#include <unordered_map>

namespace Dune
{
  class SingletonStorage
  {
    struct Item
    {
      virtual ~Item() {}
    };

  public:
    template <class Object>
    struct DUNE_PRIVATE ItemWrapper : public Item
    {
      Object obj_;
      template <class... Args>
      ItemWrapper(Args &&... args) : obj_(std::forward< Args >( args )...)
      {}
    };

    typedef std::shared_ptr< Item > PointerType;
    typedef std::type_index         KeyType;

    typedef std::unordered_map< KeyType, PointerType > StorageType;

  private:
    StorageType storage_;

  public:
    ~SingletonStorage() {}
    /** \brief return singleton instance of given Object type.
     */
    template <class Object, class... Args>
    Object& instance(Args &&... args)
    {
      {
        typedef ItemWrapper< Object > ItemWrapperType;
        PointerType& ptr = storage_[ std::type_index(typeid(Object)) ];
        if( ! ptr )
        {
          ptr.reset( new ItemWrapperType(std::forward< Args >( args )...) );
        }
        assert( dynamic_cast< ItemWrapperType* > (ptr.operator->()) );
        return static_cast< ItemWrapperType& > (*ptr).obj_;
      }
    }
  };
} // namespace Dune

#endif //  #ifndef DUNE_GRID_SINGLETONLIST_HH
