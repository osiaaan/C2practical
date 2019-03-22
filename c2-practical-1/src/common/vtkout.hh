#ifndef VTKOUT_HH_INCLUDED
#define VTKOUT_HH_INCLUDED

//- include vector implementation
#include <dune/grid/io/file/vtk/vtkwriter.hh>

// forward decleration
template <class> class DiscreteFunction;
template <class> class DiscreteVectorfield;


/** \brief output method to write a \b Paraview vtu file
    \param[in] discreteFunction   discrete function to visualize
    \param[in] step optional parameter to number discrete function output
    \param[in] name name of discrete function, default is "Uh".
*/
class Output
{
  private:
  //---------------------------------------------------------------------
  //  Output for visualization with Paraview
  //---------------------------------------------------------------------
  template <class DF>
  struct VTKWrapper : public Dune::VTKWriter<Dune::GridPartView>::VTKFunction
  {
    typedef typename Dune::VTKWriter<Dune::GridPartView>::VTKFunction Base;
    typedef typename Dune::GridPartView:: template Codim<0>::Entity Entity;

    VTKWrapper(const DF &df, const std::string name)
    : df_(df), name_(name)
    {
    }

    virtual int ncomps () const
    {
      return 1;
    }
    virtual double evaluate (int comp, const Entity& e,
                             const Dune::LocalCoordType &xi) const
    {
      ElementType element_(df_.grid(), e.template subEntity<0>(0));
      double val = df_.evaluate(element_,xi);
      if (std::abs(val)<1e-10) val = 0.;
      return val;
    }
    virtual std::string name () const
    {
      return name_;
    }
    private:
    const DF &df_;
    const std::string name_;
  };
  public:
  /** \brief constrctor
   * \param[in] grid the grid for which vtk output is to be generated
   * \param[in] step number used in file name
   * \param[in] name base name for output file
   */
  Output ( const GridType& grid,
           const int step = -1,
           const std::string name = "Uh" )
  : step_(step),
    name_(name),
    out_( grid.gridView() )
  {}
  /** \brief write file to disc
   */
  void write()
  {
    // create filename
    std::stringstream filename;
    filename << name_;
    if( step_ >= 0 )
      filename << "_" << step_;
    // write file
    out_.write( filename.str().c_str() );
  }
  /** \brief add a discretefunction for output
   * \param[in] discreteFunction the discrete function to output
   * \param[in] name name of function to use in vtk file
   */
  template <class BaseFunctionSet>
  void add(const DiscreteFunction<BaseFunctionSet>& discreteFunction,
           const std::string name = "Uh" )
  {
    VTKWrapper< DiscreteFunction<BaseFunctionSet> >* f =
      new VTKWrapper< DiscreteFunction<BaseFunctionSet> >(discreteFunction, name);
    out_.addVertexData( f );
  }
  /** \brief add a piecewise linear function given by a Dof vector
   * \param[in] v dof vector
   * \param[in] name name of function to use in vtk file
   */
  template <class V>
  void addVertexData( V &v, const std::string& name )
  {
    for (unsigned int i=0;i<v.size();++i)
      if (std::abs(v[i])<1e-10) v[i] = 0.;
    out_.addVertexData( v, name );
  }
  /** \brief add a piecewise constant function given by a Dof vector
   * \param[in] v dof vector
   * \param[in] name name of function to use in vtk file
   */
  template <class V>
  void addCellData( V &v, const std::string& name )
  {
    for (unsigned int i=0;i<v.size();++i)
      if (std::abs(v[i])<1e-10) v[i] = 0.;
    out_.addCellData( v, name );
  }
  private:
  int step_;
  std::string name_;
  Dune :: OutputType out_;
};


/**
    @}
*/
#endif
