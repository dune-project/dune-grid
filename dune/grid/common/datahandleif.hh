// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DATAHANDLEIF_HH
#define DUNE_DATAHANDLEIF_HH

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
     \ingroup GICollectiveCommunication
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
     */
    template <class T>
    void read(T & val) const
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
     \ingroup GICollectiveCommunication
   */
  template <class DataHandleImp, class DataTypeImp>
  class CommDataHandleIF
  {
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
    bool fixedsize (int dim, int codim) const
    {
      auto basePtr = &CommDataHandleIF<DataHandleImp,DataTypeImp>::fixedSize;
      auto derPtr = &DataHandleImp::fixedSize;
      bool hasOverwrittenFixedSize = basePtr != derPtr;
      if (hasOverwrittenFixedSize)
        return asImp().fixedSize(dim,codim);
      else
        return asImp().fixedsize(dim,codim);
    }

    /** @brief
       returns true if size of data per entity of given dim and codim is a constant
       @param dim valid dimension (i.e. the grids dimension)
       @param codim valid codimension of the entity set for which data should be communicated

       This method calls 'fixedSize' (with a capital S) of the derived class,
       if it exists in the derived class.  Otherwise, it calls 'fixedsize'.
     */
    bool fixedSize (int dim, int codim) const
    {
      auto basePtr = &CommDataHandleIF<DataHandleImp,DataTypeImp>::fixedSize;
      auto derPtr = &DataHandleImp::fixedSize;
      bool hasOverwrittenFixedSize = basePtr != derPtr;
      if (hasOverwrittenFixedSize)
        return asImp().fixedSize(dim,codim);
      else
        return asImp().fixedsize(dim,codim);
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
        @param buff message buffer provided by the grid
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
