// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ONE_D_GRID_LEVELITERATOR_HH
#define DUNE_ONE_D_GRID_LEVELITERATOR_HH

/** \file
 * \brief The OneDGridLevelIterator class
 */

#include <dune/grid/common/gridenums.hh>

#include <dune/grid/onedgrid/onedgridentitypointer.hh>

namespace Dune {



  //**********************************************************************
  //
  // --OneDGridLevelIterator
  // --LevelIterator
  /** \brief Iterator over all entities of a given codimension and level of a grid.
   * \ingroup OneDGrid
   */
  template<int codim, PartitionIteratorType pitype, class GridImp>
  class OneDGridLevelIterator :
    public OneDGridEntityPointer <codim, GridImp>
  {
  public:
    enum {dim=GridImp::dimension};
    friend class OneDGrid;
    friend class OneDGridEntity<codim,dim,GridImp>;
    friend class OneDGridEntity<0,dim,GridImp>;

    typedef typename GridImp::template Codim<codim>::Entity Entity;

  protected:

    /** \brief Constructor from a given iterator */
    OneDGridLevelIterator<codim,pitype, GridImp>(OneDEntityImp<dim-codim>* it)
      : OneDGridEntityPointer<codim, GridImp>(it)
    {}

  public:

    //! prefix increment
    void increment() {
      GridImp::getRealImplementation(this->virtualEntity_).setToTarget(GridImp::getRealImplementation(this->virtualEntity_).target_->succ_);
    }
  };

}  // namespace Dune

#endif
