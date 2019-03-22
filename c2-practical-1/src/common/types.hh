#ifndef TYPES_HH 
#define TYPES_HH

//- system includes 
#include <cassert> 
#include <iostream> 

// ****************************************************************************
// Central Includes from the dune-interface files
// ****************************************************************************
#include <dune/geometry/referenceelements.hh>
#include <dune/grid/common/gridview.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/io/file/dgfparser/dgfalu.hh>

#ifndef HAVE_DUNE_ALUGRID 
#error "ALUGrid not found, re-configure with ALUGrid support!"
#endif
#include <dune/alugrid/dgf.hh>

#ifndef NDEBUG 
#define EXERCISE(m) std::cout << m << " has not been \
implemented! See file: " << __FILE__ << "  line: " << __LINE__ << std::endl;
#else 
#define EXERCISE(m)
#endif


//-----------------------------------------------------------------------
// Some typedefs to be used in the programm:
// dimw=dimp=2:      we study grids in R^2
// LocalCoordType:   coords in reference element
// GlobalCoordType:  coord in world coordinats
// GridType:         implementation of a triangulare grid
// IteratorType: iterator over the leaf elements of the grid
// EntityPtrType:    pointer to one element of the grid
// GeomType:         geometry of one element of the grid
// IndexSetType:     a class for obtaining indicies for grid entities
//-----------------------------------------------------------------------
namespace Dune 
{
  static const int dimw = 2;
  static const int dimp = 2;
  // type of local cooordinate type 
  typedef FieldVector<double,dimp> LocalCoordType;

  // type of local coordinate type for faces
  typedef FieldVector<double,dimw-1> FaceCoordType;

  // type of global coordinate type 
  typedef FieldVector<double,dimw> GlobalCoordType;

  // type of global coordinate type 
  typedef FieldMatrix<double,dimw,dimp> JacobianType;

  // type of used grid 
  typedef ALUGrid< dimp, dimw, simplex, conforming > GridType;
  //typedef ALUGrid< dimp, dimw, simplex, nonconforming > GridType;

  // type of used grid part 
  typedef GridType::LeafGridView GridPartView;

  // type of used output 
  typedef Dune::VTKWriter<GridPartView> OutputType;

} // end namespace Dune 

// type of local coordinate 
typedef Dune :: LocalCoordType LocalCoordType;

// type of local coordinate on face
typedef Dune :: FaceCoordType FaceCoordType;

// type of global coordinate 
typedef Dune :: GlobalCoordType GlobalCoordType;

// type of jacobian
typedef Dune :: JacobianType JacobianType;

// not whether this should be removed 
using namespace std;

#endif
