#ifndef DUNE_GRID_SINGLETON_HH
#define DUNE_GRID_SINGLETON_HH

//- System includes
#include <cassert>
#include <iostream>
#include <memory>
#include <typeindex>
#include <utility>
#include <map>
#include <unordered_map>
#include <vector>

#include <dune/common/visibility.hh>

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
    struct ItemWrapper : public Item
    {
      Object obj_;
      template <class... Args>
      ItemWrapper(Args &&... args) : obj_(std::forward< Args >( args )...)
      {}
    };

    struct NullDeleter
    {
      void operator()(Item *p) {}
    };

    typedef std::shared_ptr< Item > PointerType; // TODO fix issue with unique_ptr
    typedef std::type_index         KeyType;

    typedef std::unordered_map< KeyType, PointerType > StorageType;

  private:
    StorageType storage_;

    static const bool placeStaticVariableInline = false ;

  public:
    ~SingletonStorage() {std::cout << "SingletonStorage\n";}
    /** \brief return singleton instance of given Object type.
     */
    template <class Object, class... Args>
    Object& instance(Args &&... args)
    {
      // this way of creating variables only works with gcc, not with clang
      if constexpr ( placeStaticVariableInline )
      {
        static Object obj( std::forward< Args >( args )...);
        return obj;
      }
      else
      {
        //std::cout << "Accessing Object " << typeid(Object).name() << std::endl;
        //std::cout << "typeindex = " << std::type_index(typeid(Object)).hash_code() << std::endl;

        typedef ItemWrapper< Object > ItemWrapperType;
        PointerType& ptr = storage_[ std::type_index(typeid(Object)) ];
        if( ! ptr )
        {
          ptr.reset( new ItemWrapperType(std::forward< Args >( args )...) );
          //std::cout << "Create Object " << typeid(Object).name() << std::endl;
          //std::cout << "typeindex = " << std::type_index(typeid(Object)).hash_code() << std::endl;
        }
        assert( dynamic_cast< ItemWrapperType* > (ptr.operator->()) );
        return static_cast< ItemWrapperType& > (*ptr).obj_;
      }
    }
  };
} // namespace Dune

#endif //  #ifndef DUNE_GRID_SINGLETONLIST_HH
