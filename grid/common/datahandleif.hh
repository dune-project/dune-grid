// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DATAHANDLE_HH
#define DUNE_DATAHANDLE_HH

/** @file
   @author FR-dune
   @brief Describes the parallel communication interface class for
   MessageBuffers and DataHandles
 */

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
     \ingroup GICollectiveCommunication
   */
  template <class MessageBufferImp>
  class MessageBufferIF
  {
    MessageBufferImp & buff_;
  public:
    //! stores reference to original buffer
    MessageBufferIF(MessageBufferImp & buff) : buff_(buff) {}

    //! just wraps the call of the internal buffer method write
    //! which writes the data of type T from the buffer by using the
    //! assigment operator of T
    template <class T>
    void write(const T & val)
    {
      buff_.write(val);
    }

    //! just wraps the call of the internal buffer method read
    //! which reads the data of type T from the buffer by using the
    //! assigment operator of T
    template <class T>
    void read(T & val) const
    {
      buff_.read(val);
    }
  }; // end class MessageBufferIF


  /** @brief CommDataHandleIF describes the features of a data handle for
     communication in parallel runs using the Grid::communication methods.
     Here the Barton-Nackman trick is used to interprete data handle objects
     as it's interface. Therefore usable data handle classes need to be
     derived from this class.
     \ingroup GICollectiveCommunication
   */
  template <class DataHandleImp, class DataTypeImp>
  class CommDataHandleIF
  {
  public:
    //! data type of data to communicate
    typedef DataTypeImp DataType;

  protected:
    // one should not create an explicit instance of this inteface object
    CommDataHandleIF() {};

  public:
    //! returns true if data for this codim should be communicated
    bool contains (int dim, int codim) const
    {
      return asImp().contains(dim,codim);
    }

    //! returns true if size of data per entity of given dim and codim is a constant
    bool fixedsize (int dim, int codim) const
    {
      return asImp().fixedsize(dim,codim);
    }

    /** @brief how many objects of type DataType have to be sent for a given entity
        Note: Only the sender side needs to know this size.
        @param e for which the size should be dertermined
     */
    template<class EntityType>
    size_t size (const EntityType& e) const
    {
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
      asImp().gather(buffIF,e);
    }

    /*! unpack data from message buffer to user
        n is the number of objects sent by the sender
        @param buff message buffer provided by the grid
        @param e entity for which date should be unpacked from buffer
        @param n number of data written to buffer for this entity before
     */
    template<class MessageBufferImp, class EntityType>
    void scatter (MessageBufferImp& buff, const EntityType& e, size_t n)
    {
      MessageBufferIF<MessageBufferImp> buffIF(buff);
      asImp().scatter(buffIF,e,n);
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

} // end namespace Dune
#endif
