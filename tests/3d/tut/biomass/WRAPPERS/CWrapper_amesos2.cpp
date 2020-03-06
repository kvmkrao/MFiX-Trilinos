#include <iostream>
#include "stdio.h"

//#include <Teuchos_GlobalMPISession.hpp>
//#include <Teuchos_oblackholestream.hpp>
//#include <Tpetra_DefaultPlatform.hpp>
//#include <Tpetra_Version.hpp>

extern "C"
{
  double* CPPWrapper(double** A, double* b, double* x_vec, double** aloc, int* glob, int rows, int cols, int glbn);
//  void* solve();
  void cwrapper_(double* br, double* xr, double* alc, double* A, int* globn,int *nAx, int *nAy, int *nn, char Str[7]) // note extra underscore to keep linker happy
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

  CPPWrapper(Amat, br, xr,loc, globn, m, n, *nn);

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

#include "Amesos2.hpp"
#include "Amesos2_Version.hpp"

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
double* CPPWrapper(double** Ao, double* bo, double* x_vec, double** aloc, int* glob, int rows, int cols, int nn) // A[n*m]x[m] = b[n]
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

//  int argc;
//  char **argv;
//  char* argv[];

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

  // The number of rows and columns in the matrix.
  const global_size_t numGlobalElements = nn;  // 50;

  const size_t numMyElements =  rows; //map->getNodeNumElements ();
  ArrayRCP<size_t> myGlobalElements = arcp<size_t> (numMyElements);

  // Create the matrix's row Map.
//  RCP<const map_type> map (new map_type (numGlobalElements, 0, comm));

//  const global_ordinal_type indexBase = 0;
  RCP<const map_type> map;
 {
    Array<GO>::size_type numEltsPerProc = rows;
    Array<GO> myGlobalElements (numEltsPerProc); 
      for(int i = 0; i < numMyElements; ++i) {
      myGlobalElements[i] = glob[i] -1;  //myRank + k*numProcs;
    }
  map = rcp (new map_type (numGlobalElements, myGlobalElements, 0, comm));
}


// RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));
// map->describe(*fos,Teuchos::VERB_EXTREME);


//  RCP<crs_matrix_type> A (new crs_matrix_type (map, 5));
  RCP<crs_matrix_type> A (new crs_matrix_type (map, 7));

//   sparse_mat_type A (map, 0);

  // Fill the sparse matrix, one row at a time.
//  const ST two    = static_cast<ST> ( 2.0);
//  const ST negOne = static_cast<ST> (-1.0);

  double *Values = new double[7];
//  global_ordinal_type *Indices = new global_ordinal_type[7];
  int  *Indices = new int[7];
//  global_ordinal_type  NumEntries;
  int NumEntries;
  int RowLess3, RowLess2, RowLess1;
  int RowPlus1, RowPlus2, RowPlus3;


  for (local_ordinal_type lcR = 0;
       lcR < static_cast<local_ordinal_type> (numMyElements);
       ++lcR) {
    const GO gblRow = map->getGlobalElement (lcR);

    NumEntries = 0; 

    RowLess3 = gblRow + aloc[lcR][0];
    RowLess2 = gblRow + aloc[lcR][1];
    RowLess1 = gblRow + aloc[lcR][2];
    RowPlus1 = gblRow + aloc[lcR][3];
    RowPlus2 = gblRow + aloc[lcR][4];
    RowPlus3 = gblRow + aloc[lcR][5];

      if ((RowLess3 >=0 )&&  abs(aloc[lcR][0]) != 0.0 )
      {
      Values[NumEntries]  = Ao[lcR][0]; 
      Indices[NumEntries] = RowLess3; 
      NumEntries          = NumEntries + 1 ; 
//      A.insertGlobalValues(gblRow, 1, &Ao[i][0], &RowLess3); 
//      A.insertGlobalValues (gblRow, tuple<GO> (gblRow, gblRow + aloc[lcR][3],gblRow + aloc[lcR][4]),
//                            tuple<ST> ( Ao[lcR][2], Ao[lcR][3],Ao[lcR][4]));

//      std::cout << "Inserting node  RowLess3"  <<  RowLess3 << std::endl; 
       }
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
//       std::cout << "Inserting node  RowLess1"  <<  RowLess1 << std::endl;
       }
  
//       A.insertGlobalValues(gblRow, 1, &Ao[i][3], &gblRow);

//       if (RowPlus1 < NumGlobalElements &&  RowPlus1 !=GlobalRow) //aloc[i][3] != 0.0)
       if ((RowPlus1 < numGlobalElements) &&  abs(aloc[lcR][3]) != 0)
       {
      Values[NumEntries]  = Ao[lcR][4];
      Indices[NumEntries] = RowPlus1;
      NumEntries          = NumEntries + 1 ;
//       A.insertGlobalValues(gblRow, 1, &Ao[i][4], &RowPlus1);
//       std::cout << "Inserting node  RowPlus1"  <<  RowPlus1 << std::endl;
       }
    
       if ((RowPlus2 < numGlobalElements) &&  abs(aloc[lcR][4]) != 0)
       {
      Values[NumEntries]  = Ao[lcR][5];
      Indices[NumEntries] = RowPlus2;
      NumEntries          = NumEntries + 1 ;
//       A.insertGlobalValues(gblRow, 1, &Ao[i][5], &RowPlus2);
//       std::cout << "Inserting node  RowPlus2"  <<  RowPlus2 << std::endl;
       }
   
//      if (RowPlus3 < NumGlobalElements && RowPlus3 !=GlobalRow ) //aloc[i][5] != 0.0)
      if ((RowPlus3 < numGlobalElements) && abs(aloc[lcR][5]) != 0)
       {
       Values[NumEntries]  = Ao[lcR][6];
      Indices[NumEntries] = RowPlus3;
      NumEntries          = NumEntries + 1 ;
//      A.insertGlobalValues(gblRow, 1, &Ao[i][6], &RowPlus3);
//      std::cout << "Inserting node  RowPlus3"  <<  RowPlus3 << std::endl;
       }
  
      Values[NumEntries]  = Ao[lcR][3];
      Indices[NumEntries] = gblRow;
      NumEntries          = NumEntries + 1 ;
//       A.insertGlobalValues(gblRow, 1, &Ao[i][3], &gblRow);
      A->insertGlobalValues (gblRow, NumEntries, Values, Indices); 


/*
    if(gblRow == 0 ) {
         A->insertGlobalValues (gblRow,
                            tuple<GO> (gblRow, gblRow + aloc[lcR][3],gblRow + aloc[lcR][4]),
                            tuple<ST> ( Ao[lcR][2], Ao[lcR][3],Ao[lcR][4]));
    }
   
    else if(gblRow == 1 ) {
         A->insertGlobalValues (gblRow,
                            tuple<GO> (gblRow+aloc[lcR][1],gblRow, gblRow + aloc[lcR][3],gblRow + aloc[lcR][4]),
                            tuple<ST> (Ao[lcR][1], Ao[lcR][2], Ao[lcR][3],Ao[lcR][4]));
    }
   
    else if (gblRow == numGlobalElements - 1) {
        A->insertGlobalValues (gblRow,        //myGlobalElements[i],
                            tuple<GO> (gblRow+aloc[lcR][0], gblRow+aloc[lcR][1], gblRow),
                            tuple<ST> (Ao[lcR][0],Ao[lcR][1], Ao[lcR][2]));
    }
  
 
    else if (gblRow == numGlobalElements - 2) {
        A->insertGlobalValues (gblRow,        //myGlobalElements[i],
                            tuple<GO> (gblRow+aloc[lcR][0], gblRow+aloc[lcR][1], gblRow,gblRow+aloc[lcR][3]),
                            tuple<ST> (Ao[lcR][0],Ao[lcR][1], Ao[lcR][2],Ao[lcR][3]));
    }
  
     else {
        A->insertGlobalValues (gblRow,
                            tuple<GO> (gblRow+aloc[lcR][0],gblRow +aloc[lcR][1], gblRow, gblRow +aloc[lcR][3],gblRow +aloc[lcR][4]),
                            tuple<ST> (Ao[lcR][0],Ao[lcR][1], Ao[lcR][2], Ao[lcR][3],Ao[lcR][4]));
    }
*/
  }

  // Tell the sparse matrix that we are done adding entries to it.
//  A->fillComplete ();
  A->fillComplete ();


//   typedef Tpetra::Vector<> vector_type;
//vector_type b (map);

//  multivector_type b (map,1);// ,false);

  RCP<multivector_type> x = rcp(new multivector_type(map,1));
  RCP<multivector_type> b = rcp(new multivector_type(map,1));

  // Create the right-hand side and initial guess of the linear system to solve.
//  b.putScalar (STS::one ());

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
   }

  // Create solver interface to Superlu with Amesos2 factory method
  RCP<Amesos2::Solver<MAT,MV> > solver = Amesos2::create<MAT,MV>("Superlu", A, X, B);

  solver->symbolicFactorization().numericFactorization().solve();


  /* Print the solution
   *
   * Should be:
   *
   *  [[1]
   *   [2]
   *   [3]
   *   [4]
   *   [5]
   *   [6]]
   */
  RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));


/*
Teuchos::ParameterList paramList;
paramList.set("verbosity", "low");
paramList.set("max levels", 3);
paramList.set("coarse: max size", 20);
paramList.set("multigrid algorithm", "sa");
paramList.set("sa: damping factor", 0.9);
paramList.set("smoother: type", "CHEBYSHEV");
paramList.set("smoother: overlap", 10);
paramList.set("aggregation: type", "uncoupled");
//paramList.set("aggregation: min agg size", 2);

std::string xmlFileName = "test.xml";

//RCP<MueLu::TpetraOperator<SC,LO,GO,NO> > tH = MueLu::CreateTpetraPreconditioner<SC,LO,GO,NO>(RCP<operator_type>(tpA), xmlFileName);
RCP< mtoperator > mueLuPreconditioner = MueLu::CreateTpetraPreconditioner<scalar_type,local_ordinal_type,global_ordinal_type,node_type>(RCP<operator_type>(A), xmlFileName);

  RCP<ParameterList> solverParams = parameterList ();
  solverParams->set ("Num Blocks", 200); //numBlocks); krylov
//  solverParams->set ("Block Size", 1);
  solverParams->set ("Maximum Iterations", 50); //maxIters);
//  solverParams->set ("Maximum Restarts", 10);
  solverParams->set ("Convergence Tolerance", 1.0e-4); //tol);

   RCP<linear_problem_type> Problem = rcp(new linear_problem_type(A, x, b));
//   Problem->setRightPrec (mueLuPreconditioner);
   Problem->setLeftPrec (mueLuPreconditioner);
   Problem->setProblem();

  RCP<ParameterList> belosList = rcp(new ParameterList());
  belosList->set("Maximum Iterations",    20);     //maxIts); // Maximum number of iterations allowed
  belosList->set("Convergence Tolerance", 1.0e-4); //tol);    // Relative convergence tolerance requested
  belosList->set("Block Size",          200);
//  belosList->set("Verbosity",             Belos::Errors + Belos::Warnings + Belos::StatusTestDetails);
  belosList->set("Verbosity",             Belos::Errors + Belos::Warnings);
  belosList->set("Output Frequency",      20);
  belosList->set("Output Style",          Belos::Brief);
//  belosList->set("Implicit Residual Scaling", "None");
  RCP<belos_solver_manager_type> solver;

//  if (krylovSolverType == "cg")
//    solver = rcp(new belos_pseudocg_manager_type(Problem, belosList));
//  else if (krylovSolverType == "gmres")
    solver = rcp(new belos_gmres_manager_type(Problem, belosList));
//  else
//    throw std::invalid_argument("bad Krylov solver type");

  solver->solve();
*/

   ArrayRCP<ST> view;
   int size = x->getLocalLength ();
   Array<ST> copy1(numMyElements);
   view = x->get1dViewNonConst();
   x->get1dCopy(copy1(),numMyElements);

   for(int i=0; i < numMyElements; i++) {
      x_vec[i] = copy1[i]; 
//      std::cout << myRank  << "  " << glob[i] << " " << i <<"  " <<copy1[i] << std::endl;
   }




}





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



