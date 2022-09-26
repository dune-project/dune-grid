// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_COMMON_DATAHANDLEIF_HH
#define DUNE_GRID_COMMON_DATAHANDLEIF_HH

/** @file
   @author Robert Kloefkorn
   @brief Describes the parallel communication interface class for
   MessageBuffers and DataHandles
 */

#include <dune/common/bartonnackmanifcheck.hh>

namespace Dune
{

  /** @brief
     Communication message buffer interface. This class describes the
     interface for reading and writing data to the communication message
     buffer. As message buffers might be deeply implemented in various
     packages the message buffers implementations cannot be derived from
     this interface class. Therefore we just apply the engine concept to
     wrap the message buffer call and make sure that the interface is
     fulfilled.

     \tparam MessageBufferImp Implementation of message buffer used by the grids' communication method
     \ingroup GICommunication
   */
  template <class MessageBufferImp>
  class MessageBufferIF
  {
    MessageBufferImp & buff_;
  public:
    //! stores reference to original buffer \c buff
    MessageBufferIF(MessageBufferImp & buff) : buff_(buff) {}

    /** @brief just wraps the call of the internal buffer method write
       which writes the data of type T from the buffer by using the
       assigment operator of T
       @param val reference to object that is written
     */
    template <class T>
    void write(const T & val)
    {
      buff_.write(val);
    }

    /** @brief just wraps the call of the internal buffer method read
       which reads the data of type T from the buffer by using the
       assigment operator of T
       @param val reference to object that is read

       The method is not const, because calling it advances the iterator
       to the current data of the MessageBufferImp member.
     */
    template <class T>
    void read(T & val)
    {
      buff_.read(val);
    }
  }; // end class MessageBufferIF


  /** @brief CommDataHandleIF describes the features of a data handle for
     communication in parallel runs using the Grid::communicate methods.
     Here the Barton-Nackman trick is used to interprete data handle objects
     as its interface. Therefore usable data handle classes need to be
     derived from this class.

     \tparam DataHandleImp implementation of the users data handle
     \tparam DataTypeImp type of data that are going to be communicated which is exported as \c DataType (for example double)
     \ingroup GICommunication
   */
  template <class DataHandleImp, class DataTypeImp>
  class CommDataHandleIF
  {
    template <class M>
    class CheckFixedSizeMethod
    {
      // check for old signature of deprecated fixedsize method.
      template <class T>
      static std::true_type testSignature(bool (T::*)(int, int) const);

      template <class T>
      static decltype(testSignature(&T::fixedsize)) test(std::nullptr_t);

      template <class T>
      static std::false_type test(...);

      using type = decltype(test<M>(nullptr));
    public:
      static const bool value = type::value;
    };


    template <class DH, bool>
    struct CallFixedSize
    {
      static bool fixedSize( const DH& dh, int dim, int codim )
      {
        return dh.fixedSize( dim, codim );
      }
    };

    // old, deprecated implementation
    template <class DH>
    struct CallFixedSize< DH, true >
    {
      static bool fixedSize( const DH& dh, int dim, int codim )
      {
        return dh.overloaded_deprecated_fixedsize( dim, codim );
      }
    };

  public:
    //! data type of data to communicate
    typedef DataTypeImp DataType;

  protected:
    // one should not create an explicit instance of this interface object
    CommDataHandleIF() {}

  public:
    /** @brief
       returns true if data for given valid codim should be communicated
       @param dim valid dimension (i.e. the grids dimension)
       @param codim valid codimension of the entity set for which data should be communicated
     */
    bool contains (int dim, int codim) const
    {
      CHECK_INTERFACE_IMPLEMENTATION((asImp().contains(dim,codim)));
      return asImp().contains(dim,codim);
    }

    /** @brief
       returns true if size of data per entity of given dim and codim is a constant
       @param dim valid dimension (i.e. the grids dimension)
       @param codim valid codimension of the entity set for which data should be communicated

       This method calls 'fixedSize' (with a capital S) of the derived class,
       if it exists in the derived class.  Otherwise, it calls 'fixedsize'.

       @deprecated This method (with the lower-case 's') is deprecated.  Use 'fixedSize' instead!
     */
    [[deprecated("fixedsize (lower case s) will be removed after release 2.8. Implement and call fixedSize (camelCase) instead!")]]
    int fixedsize (int dim, int codim) const
    {
      return int(fixedSize( dim, codim ));
    }

    // if this deprecation appears then the DataHandle implementation
    // is overloaded in the old 'fixedsize' method but not the new 'fixedSize'
    // method.
    [[deprecated("fixedsize (lower case s) will be removed after release 2.8. Implement and call fixedSize (camelCase) instead!")]]
    bool overloaded_deprecated_fixedsize( int dim, int codim ) const
    {
      return asImp().fixedsize( dim, codim );
    }

    /** @brief
       returns true if size of data per entity of given dim and codim is a constant
       @param dim valid dimension (i.e. the grids dimension)
       @param codim valid codimension of the entity set for which data should be communicated

       This method calls 'fixedSize' of the derived class.
     */
    bool fixedSize (int dim, int codim) const
    {
      // this should be enabled once the old fixedsize is removed
      //CHECK_INTERFACE_IMPLEMENTATION((asImp().fixedSize(dim,codim)));
      return CallFixedSize< DataHandleImp,
                            CheckFixedSizeMethod< DataHandleImp >::value >::fixedSize( asImp(), dim, codim );
    }

    /** @brief how many objects of type DataType have to be sent for a given entity
        @note Only the sender side needs to know this size.
        @param e entity for which the size should be determined
     */
    template<class EntityType>
    size_t size (const EntityType& e) const
    {
      CHECK_INTERFACE_IMPLEMENTATION((asImp().size(e)));
      return asImp().size(e);
    }

    /** @brief pack data from user to message buffer
        @param buff message buffer provided by the grid
        @param e entity for which date should be packed to buffer
     */
    template<class MessageBufferImp, class EntityType>
    void gather (MessageBufferImp& buff, const EntityType& e) const
    {
      MessageBufferIF<MessageBufferImp> buffIF(buff);
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION((asImp().gather(buffIF,e)));
    }

    /*! \brief unpack data from message buffer to user.
        @param buff message buffer provided by the grid.  This is not const,
          because the buffer has an internal iterator that gets advanced
          when reading from the buffer.
        @param e entity for which date should be unpacked from buffer
        @param n number of data written to buffer for this entity before
     */
    template<class MessageBufferImp, class EntityType>
    void scatter (MessageBufferImp& buff, const EntityType& e, size_t n)
    {
      MessageBufferIF<MessageBufferImp> buffIF(buff);
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION((asImp().scatter(buffIF,e,n)));
    }

  private:
    //!  Barton-Nackman trick
    DataHandleImp& asImp () {return static_cast<DataHandleImp &> (*this);}
    //!  Barton-Nackman trick
    const DataHandleImp& asImp () const
    {
      return static_cast<const DataHandleImp &>(*this);
    }
  }; // end class CommDataHandleIF

#undef CHECK_INTERFACE_IMPLEMENTATION
#undef CHECK_AND_CALL_INTERFACE_IMPLEMENTATION

} // end namespace Dune
#endif
