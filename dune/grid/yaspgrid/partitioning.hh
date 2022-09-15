// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_GRID_YASPGRID_PARTITIONING_HH
#define DUNE_GRID_YASPGRID_PARTITIONING_HH

/** \file
 *  \brief This file provides tools to partition YaspGrids.
 *  If you want to write your own partitioner, inherit from YLoadBalance
 *  and implement the loadbalance() method. You can also browse this file
 *  for already available useful partitioners, like YaspFixedSizePartitioner.
 */

#include<array>

#include<dune/common/math.hh>

namespace Dune
{

  namespace Yasp
  {

    /** \brief a base class for the yaspgrid partitioning strategy
     * The name might be irritating. It will probably change to YaspPartitionerBase
     * in a 3.0 release.
     */
    template<int d>
    class Partitioning
    {
    public:
      typedef std::array<int, d> iTupel;
      virtual ~Partitioning() {}
      virtual void partition(const iTupel&, int, iTupel&, int) const = 0;
    };

  }

  /** \brief a base class for the yaspgrid partitioning strategy
   * The name might be irritating. It will probably change to YaspPartitionerBase
   * in a 3.0 release.
   */
  template<int d>
  class YLoadBalance : public Yasp::Partitioning<d>
  {
  public:
    typedef std::array<int, d> iTupel;
    virtual ~YLoadBalance() {}
    void partition (const iTupel& size, int P, iTupel& dims, int overlap) const {
      this->loadbalance(size,P,dims);
    }
    virtual void loadbalance(const iTupel&, int, iTupel&) const = 0;
  };

  /** \brief Implement the default load balance strategy of yaspgrid
   */
  template<int d>
  class YLoadBalanceDefault : public YLoadBalance<d>
  {
  public:
    typedef std::array<int, d> iTupel;
    virtual ~YLoadBalanceDefault() {}

    /** \brief Distribute a structured grid across a set of processors
     *
     * \param [in] size Number of elements in each coordinate direction, for the entire grid
     * \param [in] P Number of processors
     */
    void loadbalance (const iTupel& size, int P, iTupel& dims) const final
    {
      int overlap = 1;
      this->partition(size,P,dims,overlap);
    }
    void partition (const iTupel& size, int P, iTupel& dims, int overlap) const final
    {
      double opt=1E100;
      iTupel trydims;

      trydims.fill(-1);
      dims.fill(-1);

      optimize_dims(d-1,size,P,dims,trydims,opt,overlap);
      if (dims[0] == -1)
        DUNE_THROW(Dune::GridError, "Failed to find a suitable partition");
    }
  private:
    void optimize_dims (int i, const iTupel& size, int P, iTupel& dims, iTupel& trydims, double &opt, int overlap ) const
    {
      if (i>0) // test all subdivisions recursively
      {
        for (int k=1; k<=P; k++)
          if (
            P%k==0 // k devides P
            and (
              k == 1 // no neighbors
              or
              size[i] / k >= 2*overlap // size sufficient for overlap
              )
            )
          {
            // P divisible by k
            trydims[i] = k;
            optimize_dims(i-1,size,P/k,dims,trydims,opt,overlap);
          }
      }
      else
      {
        // found a possible combination
        if (
          P == 1 // no neighbors
          or
          size[0] / P >= 2*overlap // size sufficient for overlap
          )
          trydims[0] = P;
        else
          return;

        // check for optimality
        double m = -1.0;

        for (int k=0; k<d; k++)
        {
          double mm=((double)size[k])/((double)trydims[k]);
          if (fmod((double)size[k],(double)trydims[k])>0.0001) mm*=3;
          if ( mm > m ) m = mm;
        }
        //if (_rank==0) std::cout << "optimize_dims: " << size << " | " << trydims << " norm=" << m << std::endl;
        if (m<opt)
        {
          opt = m;
          dims = trydims;
        }
      }
    }
  };

  /** \brief Implement yaspgrid load balance strategy for P=x^{dim} processors
   */
  template<int d>
  class YLoadBalancePowerD : public YLoadBalance<d>
  {
  public:
    typedef std::array<int, d> iTupel;
    virtual ~YLoadBalancePowerD() {}

    virtual void loadbalance (const iTupel& size, int P, iTupel& dims) const
    {
      for(int i=1; i<=P; ++i)
        if(Dune::power(i, d) == P) {
          std::fill(dims.begin(), dims.end(),i);
          return;
        }

      DUNE_THROW(GridError, "Loadbalancing failed: your number of processes needs to be a " << d << "-th power.");
    }
  };

  /** \brief Implement partitioner that gets a fixed partitioning from an array
   *  If the given partitioning doesn't match the number of processors, the grid should
   *  be distributed to, an exception is thrown.
   */
  template<int d>
  class YaspFixedSizePartitioner : public YLoadBalance<d>
  {
  public:
    YaspFixedSizePartitioner(const std::array<int,d>& dims) : _dims(dims) {}

    virtual ~YaspFixedSizePartitioner() {}

    virtual void loadbalance(const std::array<int,d>&, int P, std::array<int,d>& dims) const
    {
      int prod = 1;
      for (int i=0; i<d; i++)
        prod *= _dims[i];
      if (P != prod)
        DUNE_THROW(Dune::Exception,"Your processor number doesn't match your partitioning information");
      dims = _dims;
    }

  private:
    std::array<int,d> _dims;
  };

}

#endif
