// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRAPE_DATA_DISPLAY_HH
#define DUNE_GRAPE_DATA_DISPLAY_HH

//- system includes
#include <vector>
#include <limits>

//- local includes
#include "grapegriddisplay.hh"

/** @file
   @author Robert Kloefkorn
   @brief Provides a DataDisplay class using the GridDisplay and
   dune-fem module for class DiscreteFubctions support.
 */

namespace Dune
{

  // Forward Declarations
  // --------------------

  template< class ctype, int dim, int dimworld, int polOrd >
  class GrapeLagrangePoints;



  // GrapeFunction
  // -------------

  template< class GV, int dimR, int polOrd >
  struct GrapeFunction
  {
    typedef GV GridView;

    static const int dimDomain = GridView::Grid::dimension;
    static const int dimRange = dimR;

    typedef FieldVector< typename GridView::Grid::ctype, dimDomain > DomainVector;
    typedef FieldVector< typename GridView::Grid::ctype, dimRange > RangeVector;

    typedef typename GridView::template Codim< 0 >::Entity Entity;

    virtual ~GrapeFunction ()
    {}

    virtual void evaluate ( const Entity &entity, const DomainVector &x, RangeVector &y ) const = 0;

    virtual const GridView &gridView () const = 0;

    virtual std::string name () const = 0;
  };


#if HAVE_GRAPE

  // EvalFunctionData
  // ----------------

  template <class EvalImpTraits>
  struct EvalFunctionData
  {
    typedef typename EvalImpTraits :: GridType GridType;
    typedef typename EvalImpTraits :: EvalImp EvalImp;

    typedef typename GridType :: template Codim<0> :: Entity EntityType;
    enum { dim = GridType::dimension };
    enum { dimworld = GridType::dimensionworld };

    typedef typename GridType :: ctype ctype;

    typedef typename GrapeInterface<dim,dimworld>::DUNE_ELEM DUNE_ELEM;
    typedef typename GrapeInterface<dim,dimworld>::DUNE_FDATA DUNE_FDATA;

    // for the data visualization, call implementations evalCoordNow
    inline static void evalCoordNow (const EntityType &en, DUNE_FDATA *fdata, const double *coord, double * val)
    {
      EvalImp::evalCoordNow(en,fdata,coord,val);
    }

    // for the data visualization, call implementations evalDofNow
    inline static void evalDofNow (const EntityType &en, int geomType, DUNE_FDATA *fdata , int localNum, double * val)
    {
      EvalImp::evalDofNow(en,geomType,fdata,localNum,val);
    }

    // evaluate at given local coord
    inline static void evalCoord (DUNE_ELEM *he, DUNE_FDATA *df,
                                  const double *coord, double * val);

    // evaluate at dof
    inline static void evalDof (DUNE_ELEM *he, DUNE_FDATA *df, int localNum, double * val);

    // get min and max value for colorbar
    inline static void getMinMaxValues(DUNE_FDATA *df, double* min, double* max );
  };

  template <class GridImp, class DiscreteFunctionType>
  struct EvalDiscreteFunctions;

  template <class GridImp, class DiscreteFunctionType>
  struct EvalDiscreteFunctionsTraits
  {
    typedef GridImp GridType;
    typedef EvalDiscreteFunctions <GridImp, DiscreteFunctionType > EvalImp;
  };

  template <class GridImp, class DiscreteFunctionType>
  struct EvalDiscreteFunctions
    : public EvalFunctionData< EvalDiscreteFunctionsTraits <GridImp, DiscreteFunctionType > >
  {
    typedef GridImp GridType;
    typedef typename GridType :: template Codim<0> :: Entity EntityType;
    enum { dim = GridType::dimension };
    enum { dimworld = GridType::dimensionworld };

    typedef typename GridType :: ctype ctype;

    typedef typename DiscreteFunctionType :: LocalFunctionType LocalFunctionType;
    typedef typename DiscreteFunctionType :: DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

    typedef typename DiscreteFunctionSpaceType :: RangeType RangeType;
    typedef typename DiscreteFunctionSpaceType :: DomainType DomainType;

    typedef typename GrapeInterface<dim,dimworld>::DUNE_ELEM DUNE_ELEM;
    typedef typename GrapeInterface<dim,dimworld>::DUNE_FDATA DUNE_FDATA;

    // for the data visualization
    inline static void evalCoordNow (const EntityType &en, DUNE_FDATA *fdata, const double *coord, double * val);

    // for the data visualization
    inline static void evalDofNow (const EntityType &en, int geomType, DUNE_FDATA *fdata , int localNum, double * val);

    // for the data visualization
    inline static void evalScalar (const EntityType &en, int geomType,
                                   DiscreteFunctionType & func, LocalFunctionType &lf,
                                   const int * comp , int localNum, double * val);

    // for the data visualization
    inline static void evalVector (const EntityType &en, int geomType,
                                   DiscreteFunctionType & func, LocalFunctionType &lf,
                                   const int * comp, int vend, int localNum, double * val);

    // calculate min and max value of function
    inline static void calcMinMax(DUNE_FDATA * df);
  };



  // EvalGrapeFunction
  // -----------------

  template< class GV, int dimR, int polOrd >
  struct EvalGrapeFunction;

  template< class GV, int dimR, int polOrd >
  struct EvalGrapeFunctionTraits
  {
    typedef typename GV::Grid GridType;
    typedef EvalGrapeFunction< GV, dimR, polOrd > EvalImp;
  };

  template< class GV, int dimR, int polOrd >
  struct EvalGrapeFunction
    : public EvalFunctionData< EvalGrapeFunctionTraits< GV, dimR, polOrd > >
  {
    typedef GV GridView;

    typedef Dune::GrapeFunction< GV, dimR, polOrd > GrapeFunction;

    static const int dimDomain = GrapeFunction::dimDomain;
    static const int dimRange = GrapeFunction::dimRange;
    static const int dimWorld = GridView::Grid::dimensionworld;

    typedef typename GrapeFunction::DomainVector DomainVector;
    typedef typename GrapeFunction::RangeVector RangeVector;

    typedef typename GridView::template Codim< 0 >::Entity Entity;

    typedef typename GrapeInterface< dimDomain, dimWorld >::DUNE_ELEM DUNE_ELEM;
    typedef typename GrapeInterface< dimDomain, dimWorld >::DUNE_FDATA DUNE_FDATA;

    // for the data visualization
    static void evalCoordNow ( const Entity &entity, DUNE_FDATA *fdata, const double *coord, double *val );

    // for the data visualization
    static void evalDofNow ( const Entity &entity, int geomType, DUNE_FDATA *fdata, int localNum, double *val );

    // calculate min and max value of function
    static void calcMinMax ( DUNE_FDATA *fdata );
  };



  // EvalVectorData
  // --------------

  template <class GridImp, class VectorType, class IndexSetImp >
  struct EvalVectorData;

  template <class GridImp, class VectorType , class IndexSetImp >
  struct EvalVectorDataTraits
  {
    typedef GridImp GridType;
    typedef EvalVectorData <GridImp, VectorType, IndexSetImp > EvalImp;
  };

  template <class GridImp, class VectorType, class IndexSetImp >
  struct EvalVectorData
    : public EvalFunctionData< EvalVectorDataTraits <GridImp, VectorType, IndexSetImp > >
  {
    typedef GridImp GridType;
    typedef typename GridType :: template Codim<0> :: Entity EntityType;
    enum { dim = GridType::dimension };
    enum { dimworld = GridType::dimensionworld };

    typedef typename GridType :: ctype ctype;

    typedef typename GrapeInterface<dim,dimworld>::DUNE_ELEM DUNE_ELEM;
    typedef typename GrapeInterface<dim,dimworld>::DUNE_FDATA DUNE_FDATA;

    // for the data visualization
    inline static void evalCoordNow (const EntityType &en, DUNE_FDATA *fdata, const double *coord, double * val);

    // for the data visualization
    inline static void evalDofNow (const EntityType &en, int geomType, DUNE_FDATA *fdata , int localNum, double * val);

    // for the data visualization, evaluate linear funcs
    static void evalVectorLinear ( const EntityType &entity, int geomType,
                                   VectorType & func, const IndexSetImp &indexSet,
                                   const int *comp, int vend, int localNum, double *val );

    // for the data visualization, evaluate const funcs
    static void evalVectorConst ( const EntityType &entity, int geomType,
                                  VectorType & func, const IndexSetImp &indexSet,
                                  const int * comp, int vend, int localNum, double * val);

    // calculate min and max value of function
    inline static void calcMinMax(DUNE_FDATA * df);
  };
#endif

  /** \todo Please doc me!
      \ingroup Grape
   */
  template<class GridType>
  class GrapeDataDisplay : public GrapeGridDisplay < GridType >
  {
    typedef GrapeDataDisplay < GridType > MyDisplayType;

    typedef GrapeGridDisplay < GridType > BaseType;

    enum { dim = GridType::dimension };
    enum { dimworld = GridType::dimensionworld };

    typedef typename GridType :: ctype ctype;

#if HAVE_GRAPE
    typedef typename GrapeInterface<dim,dimworld>::DUNE_ELEM DUNE_ELEM;
    typedef typename GrapeInterface<dim,dimworld>::DUNE_FDATA DUNE_FDATA;
    typedef typename GrapeInterface<dim,dimworld>::DUNE_DAT DUNE_DAT;
    typedef typename GrapeInterface<dim,dimworld>::F_DATA F_DATA;
#endif

  public:
    typedef GridType MyGridType;

    //! Constructor, make a GrapeDataDisplay for given grid
    inline GrapeDataDisplay(const GridType &grid, const int myrank = -1);

    //! Constructor, make a GrapeDataDisplay for given grid
    template <class GridPartType>
    inline GrapeDataDisplay(const GridPartType & gridPart, const int myrank = -1);

    //! Desctructor
    inline ~GrapeDataDisplay();

    /*! display data stored in vector
       @param name Name of data (i.e. solution)
       @param data Data vector storing data to display
       @param indexSet The corresponding index set related to the data
       @param polOrd polynominal order of Lagrangespace, only 0 and 1 allowed
       at the momnent
       @param dimRange dimension of the result data (scalar: 1)
       @param continuous continuous or not (i.e polOrd = 0 ==> discontinuous) default is discontinuous
     */
    template< class VectorType, class IndexSetType >
    inline void displayVector ( const std::string name,
                                const VectorType &data,
                                const IndexSetType &indexSet,
                                const int polOrd,
                                const unsigned int dimRange,
                                bool continuous = false );

    //! Calls the display of the grid and draws the discrete function
    //! if discretefunction is NULL, then only the grid is displayed
    template <class DiscFuncType>
    inline void dataDisplay(const DiscFuncType &func, bool vector = false);

    //! display grid and data without grid mode
    inline void display();

    //! add discrete function to display
    template <class DiscFuncType>
    inline void addData(const DiscFuncType &func, double time = 0.0, bool vector = false );

    //! add discrete function to display
    template <class DiscFuncType>
    inline void addData(const DiscFuncType &func, std::string name , double time , bool vector = false );

    template< class GV, int dimR, int polOrd >
    void addData ( const GrapeFunction< GV, dimR, polOrd > &function );

#if HAVE_GRAPE
    //! add discrete function to display
    template <class DiscFuncType>
    inline void addData(const DiscFuncType &func, const DATAINFO * , double time );

    // retrun whether we have data or not
    bool hasData () { return (vecFdata_.size() > 0); }

    // return vector for copying in combined display
    std::vector < DUNE_FDATA * > & getFdataVec () { return vecFdata_; }

    /*! add vector to display
       @param data Data vector storing data to display
       @param indexSet The corresponding index set related to the data
       @param dinf GrapeDataDisplay internal data
       @param time simulation time of data
       @param polOrd polynominal order of Lagrangespace, only 0 and 1 allowed
       at the momnent
       @ param continuous continuous or not (i.e polOrd = 0 ==> discontinuous)
     */
    template<class VectorType, class IndexSetType >
    inline void addVector(const std::string name,
                          const VectorType & data, const IndexSetType & indexSet,
                          const double time , const int polOrd ,
                          const int dimRange, bool continuous );

    /*! add vector to display
       @param data Data vector storing data to display
       @param indexSet The corresponding index set related to the data
       @param dinf GrapeDataDisplay internal data
       @param time simulation time of data
       @param polOrd polynominal order of Lagrangespace, only 0 and 1 allowed
       at the momnent
       @ param continuous continuous or not (i.e polOrd = 0 ==> discontinuous)
     */

    template<class VectorType, class IndexSetType >
    inline void addVector(const VectorType & data, const IndexSetType & indexSet,
                          const DATAINFO * dinf, double time ,
                          const int polOrd , const int dimRange, bool continuous );

  private:
    //! hold the diffrent datas on this mesh
    std::vector < DUNE_FDATA * > vecFdata_;

    enum { polynomialOrder = 1 };
    // store lagrange points for evaluation
    GrapeLagrangePoints<ctype,dim,dimworld,polynomialOrder> lagrangePoints_;

    typedef typename GridType :: template Codim<0> :: Entity EntityCodim0Type;
    typedef void evalCoord_t (EntityCodim0Type &, DUNE_FDATA *, const double *, double * );
    typedef void evalDof_t (EntityCodim0Type &,int , DUNE_FDATA * , int , double * );

  public:
    // create object DUNE_FDATA
    static DUNE_FDATA * createDuneFunc ();
    // delete object DUNE_FDATA
    static void deleteDuneFunc (DUNE_FDATA *);
#endif
  };

  template <typename ctype, int dim, int dimworld, int polOrd>
  class GrapeLagrangePoints
  {
#if HAVE_GRAPE
    enum { maxPoints = 20 };
    enum { numberOfTypes = (dim == 2) ? 2 : 6 };

    std::vector < FieldMatrix<ctype,maxPoints,dim> > points_;
    GrapeLagrangePoints ( const GrapeLagrangePoints& );
  public:
    //! create lagrange points for given polyOrder and dim,dimworld
    GrapeLagrangePoints () : points_()
    {
      GrapeInterface_two_two::setupReferenceElements();
      GrapeInterface_two_three::setupReferenceElements();
      GrapeInterface_three_three::setupReferenceElements();

      for(int type=0; type<numberOfTypes; ++type)
      {
        FieldMatrix<ctype,maxPoints,dim> coords( ctype(0) );
        const int nvx = numberOfVertices(type);

        for(int i=0; i<nvx; ++i)
        {
          const double * p = getCoordinate(type,i);
          for(int j=0; j<dimworld; ++j)
          {
            assert( p );
            coords[i][j] = p[j];
          }
        }
        points_.push_back( coords );
      }
    }

    //! return lagrange point with localNum
    //! for given element type and polyOrder
    const FieldVector<ctype,dim> &
    getPoint (const int geomType, const int polyOrder , const int localNum ) const
    {
      assert( polOrd == polyOrder );
      assert( geomType >= 0 );
      assert( geomType < numberOfTypes );
      return points_[geomType][localNum];
    }

  private:
    static int numberOfVertices( const int type )
    {
      if(type < GrapeInterface_three_three::gr_tetrahedron)
        return GrapeInterface_two_two::getElementDescription(type)->number_of_vertices;
      else
        return GrapeInterface_three_three::getElementDescription(type)->number_of_vertices;
    }

    static const double * getCoordinate( const int type, const int i )
    {
      if(type < GrapeInterface_three_three::gr_tetrahedron)
      {
        return GrapeInterface_two_two::getElementDescription(type)->coord[i];
      }
      else
        return GrapeInterface_three_three::getElementDescription(type)->coord[i];
    }
#endif
  };

} // end namespace Dune

#include "grape/grapedatadisplay.cc"
#endif
