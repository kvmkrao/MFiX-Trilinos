#include <iostream>
#include "stdio.h"

//#include <Teuchos_GlobalMPISession.hpp>
//#include <Teuchos_oblackholestream.hpp>
//#include <Tpetra_DefaultPlatform.hpp>
//#include <Tpetra_Version.hpp>

extern "C"
{
  double* CPPWrapper(double** A, double* b, double* x_vec, double** aloc, int* glob, int rows, int cols, int glbn, int itmax, double tole);
//  void* solve();
  void cwrapper_(double* br, double* xr, double* alc, double* A, int* globn,int *nAx, int *nAy, int *nn, int *maxit, double *tol, char Str[7]) // note extra underscore to keep linker happy
  {
int n = *nAx;
int m = *nAy;

double** Amat;
double** loc;
  
    Amat=new double*[m];
    for(int k=0;k<m;k++)
     {
       Amat[k]=new double[n];
     }

   for (int i=0; i<m; i++)
   {
   for(int j=0; j <n; j++) 
    {
    Amat[i][j] = A[i+j*m];
    }
   }

    loc=new double*[m];
    for(int k=0;k<m;k++)
     {
       loc[k]=new double[n-1];
     }



   for (int i=0; i<m; i++)
   {
   for(int j=0; j <n-1; j++)
    {
    loc[i][j] = alc[i+j*m];
    }
   }

  CPPWrapper(Amat, br, xr,loc, globn, m, n, *nn, *maxit, *tol);

}
}

//int* CPPWrapper(double** Ao, double* bo, double* x_vec ,int rows, int cols, int nit, int flag, int nn) // A[n*m]x[m] = b[n]
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_Version.hpp>

#include <Teuchos_TimeMonitor.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <iostream>


//#ifndef __IFPACK2_TEST_ADDITIVESCHWARZ_RILUK_TYPEDEFS_AND_INCLUDES_HPP
//#define __IFPACK2_TEST_ADDITIVESCHWARZ_RILUK_TYPEDEFS_AND_INCLUDES_HPP

#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>

#include <MueLu.hpp>

//#include <MueLu_EpetraOperator.hpp>
#include <MueLu_Exceptions.hpp>
#include <MueLu_Hierarchy.hpp>
#include <MueLu_MasterList.hpp>
#include <MueLu_MLParameterListInterpreter.hpp>
#include <MueLu_ParameterListInterpreter.hpp>
#include <MueLu_Utilities.hpp>
//#include <MueLu_HierarchyHelpers.hpp>

#include "MueLu_FactoryManagerBase.hpp"
#include "MueLu_Hierarchy.hpp"
#include "MueLu_PFactory.hpp"
#include "MueLu_SaPFactory.hpp"
#include "MueLu_TransPFactory.hpp"
#include "MueLu_RAPFactory.hpp"
#include "MueLu_AmesosSmoother.hpp"
#include "MueLu_TrilinosSmoother.hpp"
#include "MueLu_SmootherFactory.hpp"
#include "MueLu_CoupledAggregationFactory.hpp"
#include "MueLu_TentativePFactory.hpp"
#include "MueLu_AmesosSmoother.hpp"

//#include <MueLu.hpp>
#include <MueLu_Level.hpp>
#include <MueLu_MLParameterListInterpreter.hpp>
#include <MueLu_ParameterListInterpreter.hpp>
//#ifdef HAVE_MUELU_TPETRA
#include <Tpetra_Operator.hpp>
#include <MueLu_TpetraOperator.hpp>
#include <Xpetra_TpetraVector.hpp>
#include <MueLu_CreateTpetraPreconditioner.hpp>
//#endif
//#ifdef HAVE_MUELU_EPETRA
//#include <MueLu_EpetraOperator.hpp>
//#include <Xpetra_EpetraVector.hpp>
//#include <MueLu_CreateEpetraPreconditioner.hpp>
//#endif

// These files must be included last
#include <MueLu_UseDefaultTypes.hpp>

#include <BelosSolverFactory.hpp>
#include <BelosTpetraAdapter.hpp>

#include <iostream> // library that contain basic input/output functions
#include <fstream>  // library that contains file input/output functions

#include <iomanip>
#include <chrono>
#include <ctime>
#include <thread>


/*
typedef Tpetra::MultiVector<>::scalar_type scalar_type;
typedef scalar_type ST;
typedef Tpetra::MultiVector<>::local_ordinal_type LO;
typedef Tpetra::MultiVector<>::global_ordinal_type GO;
typedef Tpetra::MultiVector<>::node_type node_type;

typedef Teuchos::ScalarTraits<scalar_type> STS;
typedef STS::magnitudeType magnitude_type;
typedef Teuchos::ScalarTraits<magnitude_type> STM;

typedef Tpetra::Map<LO, GO, node_type> map_type;
typedef Tpetra::MultiVector<scalar_type, LO, GO, node_type> multivector_type;
typedef Tpetra::CrsMatrix<scalar_type, LO, GO, node_type> sparse_mat_type;

typedef Tpetra::Vector<>::scalar_type scalar_type;
typedef Tpetra::Vector<>::global_ordinal_type global_ordinal_type;
typedef Tpetra::Vector<>::local_ordinal_type local_ordinal_type;
typedef MueLu::TpetraOperator<scalar_type,LO,GO,node_type> mtoperator; 

*/
//#include "Solve.hpp"

//int main (int argc, char **argv)
//{
double* CPPWrapper(double** Ao, double* bo, double* x_vec, double** aloc, int* glob, int rows, int cols, int nn, int itmax, double tole) // A[n*m]x[m] = b[n]
{

typedef Tpetra::MultiVector<>::scalar_type scalar_type;
typedef scalar_type ST;
typedef Tpetra::MultiVector<>::local_ordinal_type LO;
typedef Tpetra::MultiVector<>::global_ordinal_type GO;
//typedef Tpetra::MultiVector<>::node_type node_type;

typedef KokkosClassic::DefaultNode::DefaultNodeType node_type;

typedef Teuchos::ScalarTraits<scalar_type> STS;
typedef STS::magnitudeType magnitude_type;
typedef Teuchos::ScalarTraits<magnitude_type> STM;

typedef Tpetra::Map<LO, GO, node_type> map_type;
typedef Tpetra::MultiVector<scalar_type, LO, GO, node_type> multivector_type;
typedef Tpetra::CrsMatrix<scalar_type, LO, GO, node_type> sparse_mat_type;

typedef Tpetra::Vector<>::scalar_type scalar_type;
typedef Tpetra::Vector<>::global_ordinal_type global_ordinal_type;
typedef Tpetra::Vector<>::local_ordinal_type local_ordinal_type;
typedef MueLu::TpetraOperator<scalar_type,local_ordinal_type,global_ordinal_type,node_type> mtoperator;

typedef Tpetra::Operator<scalar_type,local_ordinal_type,global_ordinal_type,node_type>    operator_type;
typedef Belos::LinearProblem<scalar_type, multivector_type, operator_type> linear_problem_type;
typedef Belos::SolverManager<scalar_type, multivector_type, operator_type> belos_solver_manager_type;
typedef Belos::PseudoBlockCGSolMgr<scalar_type, multivector_type, operator_type> belos_pseudocg_manager_type;
typedef Belos::BlockGmresSolMgr<scalar_type, multivector_type, operator_type> belos_gmres_manager_type;
typedef Belos::BiCGStabSolMgr<scalar_type, multivector_type, operator_type> belos_bicgstab_manager_type;

  typedef MueLu::TpetraOperator<scalar_type, local_ordinal_type, global_ordinal_type, node_type> muelu_tpetra_operator_type;
//  typedef MueLu::Utilities<scalar_type,local_ordinal_type,global_ordinal_type,node_type> MueLuUtilities;

//#include "MueLu_UseShortNames.hpp"


//  typedef Ifpack2::Test::ST ST;
//  typedef Ifpack2::Test::GO GO;
//  typedef Ifpack2::Test::LO LO;
//  typedef Ifpack2::Test::STS STS;
//  typedef Ifpack2::Test::map_type map_type;
//  typedef Ifpack2::Test::multivector_type multivector_type;
//  typedef Ifpack2::Test::sparse_mat_type sparse_mat_type;

//  typedef Tpetra::Vector<>::scalar_type scalar_type;
//  typedef Tpetra::Vector<>::global_ordinal_type global_ordinal_type;
//  typedef Tpetra::Vector<>::local_ordinal_type local_ordinal_type;

  // global_size_t: Tpetra defines this unsigned integer type big
  // enough to hold any global dimension or amount of data.
  using Tpetra::global_size_t;
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::ArrayRCP;
  using Teuchos::arcp;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::tuple;
  using std::cerr;
  using std::cout;
  using std::endl;
  using Teuchos::parameterList;

  using Teuchos::updateParametersFromXmlFile;
  using Teuchos::updateParametersFromXmlString; 
 typedef Tpetra::CrsMatrix<> crs_matrix_type;

//  typedef typename CrsMatrixType::scalar_type scalar_type;
//  typedef typename CrsMatrixType::local_ordinal_type LO;
//  typedef typename CrsMatrixType::global_ordinal_type GO;

 std::ostream &out = std::cout; 

  Teuchos::oblackholestream blackHole;
//  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &blackHole);
//  Teuchos::oblackholestream blackHole;
//  Teuchos::GlobalMPISession mpiSession (int argc, char **argv, &blackHole);
//  Teuchos::GlobalMPISession mpiSession ( &NULL, &NULL, &blackHole);
// Teuchos::GlobalMPISession mpiSession (NULL,NULL, &blackHole);

  Teuchos::GlobalMPISession mpiSession();
  RCP<const Teuchos::Comm<int> > comm =
    Tpetra::DefaultPlatform::getDefaultPlatform ().getComm ();
  const int myRank = comm->getRank ();
  const int numProcs = comm->getSize ();


  std::clock_t c_start = std::clock();

  // The number of rows and columns in the matrix.
  const global_size_t numGlobalElements = nn;  // 50;
//  int numGlobalElements = nn;  // 50;

  const size_t numMyElements =  rows; //map->getNodeNumElements ();
//  int numMyElements =  rows; //map->getNodeNumElements ();
  ArrayRCP<size_t> glon = arcp<size_t> (numMyElements);

  // Create the matrix's row Map.
//  RCP<const map_type> map (new map_type (numGlobalElements, numMyElements, 0, comm));

  const global_ordinal_type indexBase = 0;
  RCP<const map_type> map =
    rcp (new map_type (numGlobalElements, numMyElements, indexBase, comm));
//  const global_ordinal_type indexBase = 0;
/*
  RCP<const map_type> map;
 {
    Array<GO>::size_type numEltsPerProc = rows;
    Array<GO> myGlobalElements (numEltsPerProc); 
      for(int i = 0; i < numMyElements; ++i) {
      myGlobalElements[i] = glob[i] -1;  //myRank + k*numProcs;
    }
  map = rcp (new map_type (numGlobalElements, myGlobalElements, 0, comm));
}
*/

  for(int i = 0; i < numMyElements; ++i) {
      glon[i] = glob[i] -1;  //myRank + k*numProcs;
    }

// RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));
// map->describe(*fos,Teuchos::VERB_EXTREME);

  RCP<crs_matrix_type> A (new crs_matrix_type (map, 5));
//  RCP<crs_matrix_type> A (new crs_matrix_type (map, 7));
//   sparse_mat_type A (map, 0);

//  const ST two    = static_cast<ST> ( 2.0);
//  const ST negOne = static_cast<ST> (-1.0);

    std::clock_t c_start1 = std::clock();


  double *Values = new double[5];
//  global_ordinal_type *Indices = new global_ordinal_type[7];
  int  *Indices = new int[5];
//  global_ordinal_type  NumEntries;
  int NumEntries;
  int RowLess3, RowLess2, RowLess1;
  int RowPlus1, RowPlus2, RowPlus3;


    for (local_ordinal_type lcR = 0; lcR < static_cast<local_ordinal_type> (numMyElements);  ++lcR) {
//    for(int lcR = 0; lcR < numMyElements; ++lcR) {
    int gblRow = glon[lcR]; //map->getGlobalElement (lcR);

    NumEntries = 0; 

    RowLess3 = gblRow + aloc[lcR][0];
    RowLess2 = gblRow + aloc[lcR][1];
    RowLess1 = gblRow + aloc[lcR][2];
    RowPlus1 = gblRow + aloc[lcR][3];
    RowPlus2 = gblRow + aloc[lcR][4];
    RowPlus3 = gblRow + aloc[lcR][5];

/*
      if ((RowLess3 >=0 )&&  abs(aloc[lcR][0]) != 0.0 )
      {
      Values[NumEntries]  = Ao[lcR][0]; 
      Indices[NumEntries] = RowLess3; 
      NumEntries          = NumEntries + 1 ; 
       }
*/

//      if (RowLess2 >=0 &&  RowLess2 !=GlobalRow) //aloc[i][1] != 0.0)
      if ((RowLess2 >=0 ) && abs(aloc[lcR][1]) != 0)
      { 
      Values[NumEntries]  = Ao[lcR][1];
      Indices[NumEntries] = RowLess2;
      NumEntries          = NumEntries + 1 ;
//       A.insertGlobalValues(gblRow, 1, &Ao[i][1], &RowLess2);
       }
  
       if ((RowLess1 >=0 ) && abs(aloc[lcR][2]) != 0)
       {
      Values[NumEntries]  = Ao[lcR][2];
      Indices[NumEntries] = RowLess1;
      NumEntries          = NumEntries + 1 ;
//       A.insertGlobalValues(gblRow, 1, &Ao[i][2], &RowLess1);
       }
  
//       A.insertGlobalValues(gblRow, 1, &Ao[i][3], &gblRow);

//       if (RowPlus1 < NumGlobalElements &&  RowPlus1 !=GlobalRow) //aloc[i][3] != 0.0)
       if ((RowPlus1 < numGlobalElements) &&  abs(aloc[lcR][3]) != 0)
       {
      Values[NumEntries]  = Ao[lcR][4];
      Indices[NumEntries] = RowPlus1;
      NumEntries          = NumEntries + 1 ;
//       A.insertGlobalValues(gblRow, 1, &Ao[i][4], &RowPlus1);
       }
    
       if ((RowPlus2 < numGlobalElements) &&  abs(aloc[lcR][4]) != 0)
       {
      Values[NumEntries]  = Ao[lcR][5];
      Indices[NumEntries] = RowPlus2;
      NumEntries          = NumEntries + 1 ;
//       A.insertGlobalValues(gblRow, 1, &Ao[i][5], &RowPlus2);
       }
   
//      if (RowPlus3 < NumGlobalElements && RowPlus3 !=GlobalRow ) //aloc[i][5] != 0.0)
/*
      if ((RowPlus3 < numGlobalElements) && abs(aloc[lcR][5]) != 0)
       {
       Values[NumEntries]  = Ao[lcR][6];
      Indices[NumEntries] = RowPlus3;
      NumEntries          = NumEntries + 1 ;
*/
      Values[NumEntries]  = Ao[lcR][3];
      Indices[NumEntries] = gblRow;  
      NumEntries          = NumEntries + 1 ;
//       A.insertGlobalValues(gblRow, 1, &Ao[i][3], &gblRow);
      A->insertGlobalValues (gblRow, NumEntries, Values, Indices); 

  }

  // Tell the sparse matrix that we are done adding entries to it.
//  A->fillComplete ();
  A->fillComplete ();

 std::clock_t c_start2 = std::clock();

//   typedef Tpetra::Vector<> vector_type;
//vector_type b (map);

//  multivector_type b (map,1);// ,false);

  RCP<multivector_type> x = rcp(new multivector_type(map,1));
  RCP<multivector_type> b = rcp(new multivector_type(map,1));

  // Create the right-hand side and initial guess of the linear system to solve.
//  b.putScalar (STS::one ());
//  x->randomize();

/*
   for (size_t i = 0; i < numMyElements; ++i) {
   if( b.getMap()->isNodeLocalElement(i) ){
      b.replaceLocalValue(i,0,bo[i]);
      }
    }
*/

  for (local_ordinal_type lclRow = 0;
       lclRow < static_cast<local_ordinal_type> (numMyElements);
       ++lclRow) {
    const GO gblRow = map->getGlobalElement (lclRow);
   b->sumIntoGlobalValue(gblRow, 0, bo[lclRow]);
   x->sumIntoGlobalValue(gblRow, 0, x_vec[lclRow]);
   }


 std::clock_t c_start3 = std::clock();

Teuchos::ParameterList paramList;
paramList.set("verbosity", "low");
paramList.set("max levels", 2);
paramList.set("coarse: max size", 10);
paramList.set("multigrid algorithm", "sa");
/*
paramList.set("verbosity", "low");
paramList.set("max levels", 3);
paramList.set("coarse: max size", 20);
paramList.set("multigrid algorithm", "sa");
paramList.set("sa: damping factor", 0.9);
paramList.set("smoother: type", "CHEBYSHEV");
paramList.set("smoother: overlap", 10);
paramList.set("aggregation: type", "uncoupled");
//paramList.set("aggregation: min agg size", 2);
*/

std::string xmlFileName = "test.xml";

//tested ok 
//Teuchos::RCP<muelu_tpetra_operator_type> mueLuPreconditioner = MueLu::CreateTpetraPreconditioner(A, paramList);


//tested ok
RCP< mtoperator > mueLuPreconditioner = MueLu::CreateTpetraPreconditioner<scalar_type,local_ordinal_type,global_ordinal_type,node_type>(RCP<operator_type>(A), xmlFileName);

//RCP<MueLu::TpetraOperator> mueLuPreconditioner = MueLu::CreateTpetraPreconditioner(A, paramList);
//RCP<muelu_tpetra_operator_type> mueLuPreconditioner =

//MueLu::CreateTpetraPreconditioner(A, paramList);

//RCP<muelu_tpetra_operator_type> M = MueLu::CreateTpetraPreconditioner((RCP<operator_type>)A, mueluParams);
//RCP<muelu_tpetra_operator_type> M = MueLu::CreateTpetraPreconditioner((RCP<operator_type>)A, mueluParams);

//  b.assign(*bo);

//  multivector_type x (map, 1); //, true);
//   vector_type x (map);

/*
  *fos << "MATRIX" << std::endl;
  A->describe(*fos,Teuchos::VERB_EXTREME);
  *fos << std::endl;
*/

//std::cout << itmax << " " << tole << std::endl;

//  typedef Belos::SolverFactory<scalar_type, multivector_type, operator_type> belos_factory_type;
//  typedef Belos::SolverManager<scalar_type, multivector_type, op_type> solver_type;

//  RCP<belos_solver_manager_type> solver = belos_factory_type ().create ("GMRES", solverParams);

  // Create a LinearProblem struct with the problem to solve.
  // A, X, B, and M are passed by (smart) pointer, not copied.
// ------------- set up  the pre-conditioner --------------------------
//   commented : Murali
// typedef Belos::LinearProblem<scalar_type, multivector_type, op_type> problem_type;
   // RCP<linear_problem_type> problem (new linear_problem_type (A, x, b));
//  problem->setRightPrec (additiveSchwarz);
//// ---------------------------------------------------------------------

  // Tell the solver what problem you want to solve.
//  solver->setProblem (problem);
// solver->reset (Belos::Problem);
//  Belos::ReturnType result = solver->solve ();

   RCP<linear_problem_type> Problem = rcp(new linear_problem_type(A, x, b));
   Problem->setRightPrec (mueLuPreconditioner);
//   Problem->setLeftPrec (mueLuPreconditioner);
   Problem->setProblem();
  // Set up Krylov solver and iterate.

  RCP<ParameterList> belosList = rcp(new ParameterList());
  belosList->set("Maximum Iterations",    itmax);     //maxIts); // Maximum number of iterations allowed
  belosList->set( "Num Blocks", 200);       // Maximum number of blocks in Krylov factorization
  belosList->set("Convergence Tolerance", tole); //tol);    // Relative convergence tolerance requested
  belosList->set("Block Size",          1);
//  belosList->set("Verbosity",             Belos::Errors + Belos::Warnings + Belos::StatusTestDetails + Belos::TimingDetails);
  belosList->set("Verbosity",             Belos::Errors + Belos::Warnings);
  belosList->set("Output Frequency",      20);
//  belosList->set("Output Style",          Belos::None);
  belosList->set("Output Style",          Belos::Brief);
//  belosList->set("Implicit Residual Scaling", "None");
  RCP<belos_solver_manager_type> solver;

//  if (krylovSolverType == "cg")
//    solver = rcp(new belos_pseudocg_manager_type(Problem, belosList));
//  else if (krylovSolverType == "gmres")

//    solver = rcp(new belos_gmres_manager_type(Problem, belosList));

  
  std::string inputFile="input_param.xml";
  std::string vname; 
  RCP<Teuchos::ParameterList> myParams = rcp(new Teuchos::ParameterList());
  Teuchos::updateParametersFromXmlFile(inputFile, myParams.ptr());
  vname = myParams->get<std::string>("Solver");

   std::clock_t c_start4 = std::clock();

    if(vname=="GMRES") {
    solver = rcp(new belos_gmres_manager_type(Problem, belosList));
    }
    else
    {
    solver = rcp(new belos_bicgstab_manager_type(Problem, belosList));
    }
//    throw std::invalid_argument("bad Krylov solver type");

  solver->solve();

  std::clock_t c_start5 = std::clock();

  std::cout << "time" << c_start1-c_start   << std::endl; 
  std::cout << "time1" << c_start2-c_start1 << std::endl; 
  std::cout << "time2" << c_start3-c_start2 << std::endl; 
  std::cout << "time3" << c_start4-c_start3 << std::endl; 
  
/*

  typedef Belos::SolverFactory<scalar_type, multivector_type, op_type> belos_factory_type;
  typedef Belos::SolverManager<scalar_type, multivector_type, op_type> solver_type;

  RCP<solver_type> solver = belos_factory_type ().create ("GMRES", solverParams);

  // Create a LinearProblem struct with the problem to solve.
  // A, X, B, and M are passed by (smart) pointer, not copied.
// ------------- set up  the pre-conditioner --------------------------
//   commented : Murali
 typedef Belos::LinearProblem<scalar_type, multivector_type, op_type> problem_type;
  RCP<problem_type> problem (new problem_type (A, x, b));
  problem->setRightPrec (mueLuPreconditioner);
//// ---------------------------------------------------------------------

  // Tell the solver what problem you want to solve.
  solver->setProblem (problem);
  solver->reset (Belos::Problem);
  Belos::ReturnType result = solver->solve ();

*/


/*
  using Teuchos::reduceAll;
  using Teuchos::REDUCE_MIN;
  using Teuchos::outArg;
  int lclSuccess = 1;
  int gblSuccess = 1;
  std::ostringstream errStrm;
*/

//  try {
//    Ifpack2::Test::solve (A, b, x, 1, 1, 100, 1000, 1.0e-8, false, "RILUK");
//  }
//  catch (std::exception& e) {
//    lclSuccess = 0;
//    errStrm << e.what ();
//  }
//  reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));

   ArrayRCP<ST> view;
   int size = x->getLocalLength ();
   Array<ST> copy1(numMyElements);
   view = x->get1dViewNonConst();
   x->get1dCopy(copy1(),numMyElements);

   for(int i=0; i < numMyElements; i++) {
      x_vec[i] = copy1[i]; 
//      std::cout << myRank  << "  " << glob[i] << " " << i <<"  " <<copy1[i] << std::endl;
   }


// *fos << "LHS :" << std::endl;
//  x->describe(*fos,Teuchos::VERB_EXTREME);





/* // -----------------------commented for the time being -------------------
  if (gblSuccess != 1) {
    // We assume that it's OK for MPI processes other than Proc 0 in
    // MPI_COMM_WORLD to print to stderr.  That's generally true when
    // we run tests.
    cerr << "Belos solve with Ifpack2 preconditioner threw an exception "
      "on one or more processes!" << endl;
    for (int r = 0; r < numProcs; ++r) {
      if (r == myRank) {
        std::ostringstream os;
        os << "Process " << myRank << ": " << errStrm.str () << endl;
        cerr << os.str ();
      }
      comm->barrier (); // wait for output to finish
      comm->barrier ();
      comm->barrier ();
    }
  }

// ----------------------------------------------------------------------- */

/*
  if (comm->getRank () == 0) {
    if (gblSuccess == 1) {
      cout << "End Result: TEST PASSED" << endl;
    }
    else {
      cout << "End Result: TEST FAILED" << endl;
    }
  }

 return EXIT_SUCCESS;
*/

}


/*
  std::ostream &out = std::cout;

  RCP<const sparse_mat_type> A_ptr = rcpFromRef (A);
  RCP<const multivector_type> b_ptr = rcpFromRef (b);
  RCP<multivector_type> x_ptr = rcpFromRef (x);
*/

/* ------------  // Set RILUK parameters \\ ------------------
  ParameterList innerPlist;
  innerPlist.set ("fact: drop tolerance", 0.0);
  innerPlist.set ("fact: iluk level-of-fill", ilukFillLevel);
  innerPlist.set ("fact: relax value", 0.0);

  //Teuchos::RCP<prec_type> innerPrec;
  //Ifpack2::Factory factory;
  //innerPrec = factory.create (innerPrecondType, A_ptr);
  //innerPrec->setParameters (innerPlist);
// -------------------------------------------------------------*/

   // Create (outer) Additive Schwarz preconditioner
//  RCP<outer_prec_type> additiveSchwarz (new outer_prec_type (A_ptr));

/* ------------- Set outer preconditioner parameters -------------
  ParameterList ASlist;
  ASlist.set ("inner preconditioner name", innerPrecondType);
  ASlist.set ("inner preconditioner parameters", innerPlist);
//  ASlist.set ("schwarz: combine mode", "ZERO");
  ASlist.set ("schwarz: combine mode", "ADD");
  ASlist.set ("schwarz: overlap level", overlapLevel);
  ASlist.set ("schwarz: use reordering", reorder);

  additiveSchwarz->setParameters (ASlist);

  // Compute (set up) the (outer) preconditioner
  additiveSchwarz->initialize ();
  additiveSchwarz->compute ();
// ----------------------------------------------------------------*/

/*
    std::string optionsFile = "mueluOptions.xml";

//  RCP<MueLu::TpetraOperator> mueLuPreconditioner;
//  Teuchos::RCP<MueLu::TpetraOperator<Scalar,LocalOrdinal,GlobalOrdinal,Node> > MueLu::CreateTpetraPreconditioner
//    Teuchos::RCP<MueLu::TpetraOperator<scalar_type,LO,GO,node_type> >  
       RCP<const mtoperator > mueLuPreconditioner ( A_ptr, optionsFile); 
//     MueLu::CreateTpetraPreconditioner (const Teuchos::RCP< Tpetra::CrsMatrix<scalar_type, LO,GO, node_type> >& A_ptr, optionsFile);
//  std::string optionsFile = "mueluOptions.xml";
//  mueLuPreconditioner = MueLu::CreateTpetraPreconditioner(A_ptr, optionsFile);



  // Set GMRES (iterative linear solver) parameters
  RCP<ParameterList> solverParams = parameterList ();
  solverParams->set ("Num Blocks", 30); //numBlocks); krylov
//  solverParams->set ("Block Size", 1);
  solverParams->set ("Maximum Iterations", 20); //maxIters);
//  solverParams->set ("Maximum Restarts", 10);
  solverParams->set ("Convergence Tolerance", 1.0e-4); //tol);
//  solverParams->set ("Implicit Residual Scaling", "None");
//  solverParams->set ("Explicit Residual Scaling", "None");

  // Create the GMRES solver using a "factory" and
  // the list of solver parameters created above.
  typedef Belos::SolverFactory<scalar_type, multivector_type, op_type> belos_factory_type;
  typedef Belos::SolverManager<scalar_type, multivector_type, op_type> solver_type;

  RCP<solver_type> solver = belos_factory_type ().create ("GMRES", solverParams);

  // Create a LinearProblem struct with the problem to solve.
  // A, X, B, and M are passed by (smart) pointer, not copied.
// ------------- set up  the pre-conditioner --------------------------
//   commented : Murali
 typedef Belos::LinearProblem<scalar_type, multivector_type, op_type> problem_type;
  RCP<problem_type> problem (new problem_type (A_ptr, x_ptr, b_ptr));
//  problem->setRightPrec (additiveSchwarz);
//// ---------------------------------------------------------------------

  // Tell the solver what problem you want to solve.
  solver->setProblem (problem);
  solver->reset (Belos::Problem);
  Belos::ReturnType result = solver->solve ();

*/

/*
  RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));

  *fos << "MATRIX" << std::endl;
  A_ptr->describe(*fos,Teuchos::VERB_EXTREME);
  *fos << std::endl;

  *fos << "RHS :" << std::endl;
  b_ptr->describe(*fos,Teuchos::VERB_EXTREME);
  *fos << std::endl;

*/

/*
  // Ask the solver how many iterations the last solve() took.
  const int numIters = solver->getNumIters ();

  if (myRank == 0) {
    if (solver->isLOADetected ()) {
      cout << "Detected a loss of accuracy!" << endl;
    }
    cout << "The Belos solve took " << numIters << " iteration"
         << (numIters != 1 ? "s" : "");
    if (result == Belos::Converged) {
      cout << " to converge." << endl;
    } else {
      cout << ", but did not converge." << endl;
    }
    cout << "It achieved a tolerance of: " << solver->achievedTol () << endl;
  }
*/
/*
}


} // namespace Test
} // namespace Ifpack2
*/



#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_config.h"

#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_oblackholestream.hpp>

extern "C" {

  void xmlstrc_(char *varname, int *varlen, char *Str1, int *len)  {
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::parameterList;
  Teuchos::oblackholestream blackHole;
  Teuchos::GlobalMPISession mpiSession();

  std::string vname;
  vname.assign(varname, *varlen);

  std::string name1;
  name1.assign(Str1, *len);

  std::string inputFile;
  inputFile="input_param.xml";
  RCP<Teuchos::ParameterList> myParams = rcp(new Teuchos::ParameterList());
  Teuchos::updateParametersFromXmlFile(inputFile, myParams.ptr());
  vname = myParams->get<std::string>(name1);
//  std::cout << vname << name1 << std::endl;
  strncpy(varname, vname.c_str(), *varlen);
  varname[*varlen-1] = 0 ;
//  std::cout << vname.c_str() << std::endl;
}
}


#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_config.h"

#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_oblackholestream.hpp>

extern "C" {
  void xmlintc_(int *intvar1, char *Str1, int *len)  {
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::parameterList;

  Teuchos::oblackholestream blackHole;
  Teuchos::GlobalMPISession mpiSession();

  std::string name1;
  name1.assign(Str1, *len);
  std::string inputFile;
  inputFile="input_param.xml";
  RCP<Teuchos::ParameterList> myParams = rcp(new Teuchos::ParameterList());
  Teuchos::updateParametersFromXmlFile(inputFile, myParams.ptr());
   *intvar1 = myParams->get<int>(name1);
}
}


#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_config.h"

#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_oblackholestream.hpp>

extern "C" {
  void xmlfltc_(double *intvar1, char *Str1, int *len)  {
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::parameterList;

  Teuchos::oblackholestream blackHole;
  Teuchos::GlobalMPISession mpiSession();

  std::string name1;
  name1.assign(Str1, *len);
  std::string inputFile;
  inputFile="input_param.xml";
  RCP<Teuchos::ParameterList> myParams = rcp(new Teuchos::ParameterList());
  Teuchos::updateParametersFromXmlFile(inputFile, myParams.ptr());
   *intvar1 = myParams->get<double>(name1);
//   std::cout << *intvar1 << name1 << std::endl;
}
}



