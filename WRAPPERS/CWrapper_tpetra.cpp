#include <iostream>
#include "stdio.h"

extern "C"
{
  double* CPPWrapper(double** A, double* b, double* x_vec, double** aloc, int* glob, int rows, int cols, int glbn, int itmax, double tole);
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
std::cout << *maxit <<  *tol << std::endl; 
  delete[] Amat;
   delete [] loc;

}
}


#include <iostream>

#include <Teuchos_TypeNameTraits.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Tpetra_DefaultPlatform.hpp>
#include <Teuchos_TimeMonitor.hpp>

#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_Operator.hpp>

#include <MueLu_ParameterListInterpreter.hpp>
#include <MueLu_TpetraOperator.hpp>
#include <MueLu_CreateTpetraPreconditioner.hpp>
#include <MueLu_UseDefaultTypes.hpp>

#include <BelosSolverFactory.hpp>
#include <BelosTpetraAdapter.hpp>
#include "BelosConfigDefs.hpp"

#include <ctime>

using namespace Teuchos;
using namespace std;

double* CPPWrapper(double** Ao, double* bo, double* x_vec, double** aloc, int* glob, int rows, int cols, int nn, int itmax, double tole) // A[n*m]x[m] = b[n]
{

  typedef Tpetra::MultiVector<>::scalar_type ST;
  typedef Tpetra::MultiVector<>::local_ordinal_type  LO;
  typedef Tpetra::MultiVector<>::global_ordinal_type GO;

  typedef KokkosClassic::DefaultNode::DefaultNodeType node_type;

  typedef Tpetra::Map<LO, GO, node_type> map_type;
  typedef Tpetra::MultiVector<ST, LO, GO, node_type> multivector_type;
  typedef Tpetra::CrsMatrix<ST, LO, GO, node_type> sparse_mat_type;

  typedef MueLu::TpetraOperator<ST,LO,GO,node_type> mtoperator;

  typedef Tpetra::Operator<ST,LO,GO,node_type>    operator_type;
  typedef Belos::LinearProblem<ST, multivector_type, operator_type> linear_problem_type;
  typedef Belos::SolverManager<ST, multivector_type, operator_type> belos_solver_manager_type;

  typedef Belos::BlockCGSolMgr<ST, multivector_type, operator_type>    belos_blockcg_manager_type;
  typedef Belos::PseudoBlockCGSolMgr<ST, multivector_type, operator_type> belos_pseudocg_manager_type;

// flexible 
  typedef Belos::BlockGmresSolMgr<ST, multivector_type, operator_type>       belos_flexgmres_manager_type;
// standard 
  typedef Belos::PseudoBlockGmresSolMgr<ST, multivector_type, operator_type> belos_stdgmres_manager_type;

  typedef Belos::BiCGStabSolMgr<ST, multivector_type, operator_type> belos_bicgstab_manager_type;
  typedef Belos::TFQMRSolMgr<ST, multivector_type, operator_type>      belos_tfqmr_manager_type;

  typedef Belos::LSQRSolMgr<ST, multivector_type, operator_type>       belos_lsqr_manager_type;
  typedef Belos::GCRODRSolMgr<ST, multivector_type, operator_type>     belos_recgmres_type;
  typedef Belos::GmresPolySolMgr<ST, multivector_type, operator_type>  belos_hybridgmres_type;
  typedef Belos::PseudoBlockTFQMRSolMgr<ST, multivector_type, operator_type> belos_psedotfqmr_type;
//  typedef Belos::PseudoBlockCGSolMgr<ST, multivector_type, operator_type> belos_psedobcg_type;  

  std::cout << "Execution space name: ";
  typedef Tpetra::Map<>::device_type::execution_space default_execution_space;
  std::cout << Teuchos::TypeNameTraits<default_execution_space>::name () << std::endl;


  using Tpetra::global_size_t;
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::ArrayRCP;
  using Teuchos::arcp;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::tuple;
  using Teuchos::parameterList;

  using Teuchos::updateParametersFromXmlFile;
  using Teuchos::updateParametersFromXmlString; 



 std::ostream &out = std::cout; 

  Teuchos::oblackholestream blackHole;

  Teuchos::GlobalMPISession mpiSession();
  RCP<const Teuchos::Comm<int> > comm =
    Tpetra::DefaultPlatform::getDefaultPlatform ().getComm ();
  const int myRank = comm->getRank ();
  const int numProcs = comm->getSize ();


  RCP<Time> insertva     = TimeMonitor::getNewCounter ("InsertValues ");
  RCP<Time> FillTimer    = TimeMonitor::getNewCounter ("FillComplete A");


  RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));

  // The number of rows and columns in the matrix.
  const global_size_t numGlobalElements = nn;  // 50;
//  int numGlobalElements = nn;  // 50;

  const size_t numMyElements =  rows; //map->getNodeNumElements ();
//  int numMyElements =  rows; //map->getNodeNumElements ();
  ArrayRCP<size_t> glon = arcp<size_t> (numMyElements);

  // Create the matrix's row Map.
//  RCP<const map_type> map (new map_type (numGlobalElements, numMyElements, 0, comm));

  const GO indexBase = 0;
  RCP<const map_type> map =
    rcp (new map_type (numGlobalElements, numMyElements, indexBase, comm));
//  const GO indexBase = 0;
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

    RCP<sparse_mat_type> A (new sparse_mat_type (map, 7));

  double *Values = new double[7];
//  GO *Indices = new GO[7];
  int  *Indices = new int[7];
//  GO  NumEntries;
  int NumEntries;
  int RowLess3, RowLess2, RowLess1;
  int RowPlus1, RowPlus2, RowPlus3;

  {
  TimeMonitor monitor (*insertva);
    for (LO lcR = 0; lcR < static_cast<LO> (numMyElements);  ++lcR) {
//    for(int lcR = 0; lcR < numMyElements; ++lcR) {
    int gblRow = glon[lcR]; //map->getGlobalElement (lcR);

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
   
      if ((RowPlus3 < numGlobalElements) && abs(aloc[lcR][5]) != 0)
       {
       Values[NumEntries]  = Ao[lcR][6];
      Indices[NumEntries] = RowPlus3;
      NumEntries          = NumEntries + 1 ;
      }

      Values[NumEntries]  = Ao[lcR][3];
      Indices[NumEntries] = gblRow;  
      NumEntries          = NumEntries + 1 ;
//       A.insertGlobalValues(gblRow, 1, &Ao[i][3], &gblRow);
      A->insertGlobalValues (gblRow, NumEntries, Values, Indices); 
  }
}


{
    TimeMonitor monitor (*FillTimer);
    A->fillComplete (map,map);
  }


//   typedef Tpetra::Vector<> vector_type;
//vector_type b (map);

//  multivector_type b (map,1);// ,false);

  RCP<multivector_type> x = rcp(new multivector_type(map,1));
  RCP<multivector_type> b = rcp(new multivector_type(map,1));

  for (LO lclRow = 0;
       lclRow < static_cast<LO> (numMyElements);
       ++lclRow) {
    const GO gblRow = map->getGlobalElement (lclRow);
   b->sumIntoGlobalValue(gblRow, 0, bo[lclRow]);
   }


std::string xmlFileName = "test.xml";

RCP< mtoperator > mueLuPreconditioner = MueLu::CreateTpetraPreconditioner<ST,LO,GO,node_type>(RCP<operator_type>(A), xmlFileName);

   RCP<linear_problem_type> Problem = rcp(new linear_problem_type(A, x, b));
   Problem->setRightPrec (mueLuPreconditioner);
   Problem->setProblem();

  RCP<ParameterList> belosList = rcp(new ParameterList());
  belosList->set("Maximum Iterations",    itmax);     //maxIts); // Maximum number of iterations allowed
  belosList->set( "Num Blocks", 200);       // Maximum number of blocks in Krylov factorization
  belosList->set("Convergence Tolerance", tole); //tol);    // Relative convergence tolerance requested
  belosList->set("Block Size",          1);
  belosList->set("Verbosity",             Belos::Errors + Belos::Warnings + Belos::StatusTestDetails + Belos::TimingDetails);
//  belosList->set("Verbosity",  Belos::Errors + Belos::Warnings + Belos::StatusTestDetails + Belos::TimingDetails + Belos::OrthoDetails + Belos::IterationDetails + Belos::Debug );
//  belosList->set("Verbosity",             Belos::Errors + Belos::Warnings);
  belosList->set("Output Frequency",      20);
//  belosList->set("Output Style",          Belos::None);
  belosList->set("Output Style",          Belos::Brief);
//  belosList->set("Implicit Residual Scaling", "None");
  RCP<belos_solver_manager_type> solver;

  
  std::string inputFile="precond.xml";
  std::string vname; 
  RCP<Teuchos::ParameterList> myParams = rcp(new Teuchos::ParameterList());
  Teuchos::updateParametersFromXmlFile(inputFile, myParams.ptr());
  vname = myParams->get<std::string>("Solver");

    if(vname=="GMRES") {
    solver = rcp(new belos_stdgmres_manager_type (Problem, belosList));
    }
    else
    {
    solver = rcp(new belos_bicgstab_manager_type(Problem, belosList));
    }

  solver->solve();

   ArrayRCP<ST> view;
   int size = x->getLocalLength ();
   Array<ST> copy1(numMyElements);
   view = x->get1dViewNonConst();
   x->get1dCopy(copy1(),numMyElements);

   for(int i=0; i < numMyElements; i++) {
      x_vec[i] = copy1[i]; 
//      std::cout << myRank  << "  " << glob[i] << " " << i <<"  " <<copy1[i] << std::endl;
   }

  TimeMonitor::summarize();

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