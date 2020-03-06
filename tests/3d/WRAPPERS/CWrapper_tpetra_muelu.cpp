#include <iostream>
#include "stdio.h"
#include <ctime>

extern "C"
{
  double* CPPWrapper(double** A, double* b, double* x_vec, double** aloc, int* glob, int rows, int cols, int glbn, int itmax, double tole);
  void cwrapper_(double* br, double* xr, double* alc, double* A, int* globn,int *nAx, int *nAy, int *nn, int *maxit, double *tol, char Str[7])
  {

   
   std::clock_t c_start1 = std::clock();

   int n = *nAx;
   int m = *nAy;

   double** Amat;
   double** loc;
  
   Amat=new double*[m];
    for(int k=0;k<m;k++)
     {
       Amat[k]=new double[n];
     }

    for (int i=0; i<m; i++)   {
      for(int j=0; j <n; j++) {
    Amat[i][j] = A[i+j*m];
      }
    }

   loc=new double*[m];
    for(int k=0;k<m;k++)  {
       loc[k]=new double[n-1];
    }


    for (int i=0; i<m; i++) {
       for(int j=0; j <n-1; j++) {
        loc[i][j] = alc[i+j*m];
       }
    }

   std::clock_t c_start2 = std::clock();
   CPPWrapper(Amat, br, xr,loc, globn, m, n, *nn, *maxit, *tol);
   std::clock_t c_start3 = std::clock();
   delete[] Amat; 
   delete [] loc;
   std::clock_t c_start4 = std::clock();

  std::cout << "CWrapper_total_time  "      << 1000*(c_start4-c_start1)/(double)CLOCKS_PER_SEC   << std::endl;
  std::cout << "CWrapper_form_arrays  "     << 1000*(c_start2-c_start1)/(double)CLOCKS_PER_SEC << std::endl;
  std::cout << "CppWrapper_total time  "    << 1000*(c_start3-c_start2)/(double)CLOCKS_PER_SEC << std::endl;

  }
}

#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_Version.hpp>

#include <Teuchos_TimeMonitor.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <iostream>

#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>

#include <MueLu.hpp>

#include <MueLu_Exceptions.hpp>
#include <MueLu_Hierarchy.hpp>
#include <MueLu_MasterList.hpp>
#include <MueLu_MLParameterListInterpreter.hpp>
#include <MueLu_ParameterListInterpreter.hpp>
#include <MueLu_Utilities.hpp>

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

#include <MueLu_Level.hpp>
#include <MueLu_MLParameterListInterpreter.hpp>
#include <MueLu_ParameterListInterpreter.hpp>
//#ifdef HAVE_MUELU_TPETRA
#include <Tpetra_Operator.hpp>
#include <MueLu_TpetraOperator.hpp>
#include <Xpetra_TpetraVector.hpp>
#include <MueLu_CreateTpetraPreconditioner.hpp>
#include <MueLu_UseDefaultTypes.hpp>

#include <BelosSolverFactory.hpp>
#include <BelosTpetraAdapter.hpp>
#include "BelosConfigDefs.hpp"

#include <ctime>

/*
#include <Ifpack2_Factory.hpp>
#include <Ifpack2_Preconditioner.hpp>
#include <Ifpack2_Condest.hpp>
#include <Ifpack2_CondestType.hpp>
*/

double* CPPWrapper(double** Ao, double* bo, double* x_vec, double** aloc, int* glob, int rows, int cols, int nn, int itmax, double tole)
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

  typedef Belos::BlockCGSolMgr<scalar_type, multivector_type, operator_type>    belos_blockcg_manager_type;
  typedef Belos::PseudoBlockCGSolMgr<scalar_type, multivector_type, operator_type> belos_pseudocg_manager_type;
// flexible 
  typedef Belos::BlockGmresSolMgr<scalar_type, multivector_type, operator_type>       belos_flexgmres_manager_type;
// standard 
  typedef Belos::PseudoBlockGmresSolMgr<scalar_type, multivector_type, operator_type> belos_stdgmres_manager_type;

  typedef Belos::BiCGStabSolMgr<scalar_type, multivector_type, operator_type> belos_bicgstab_manager_type;
  typedef Belos::TFQMRSolMgr<scalar_type, multivector_type, operator_type>      belos_tfqmr_manager_type;
  
  typedef Belos::LSQRSolMgr<scalar_type, multivector_type, operator_type>       belos_lsqr_manager_type;
  typedef Belos::GCRODRSolMgr<scalar_type, multivector_type, operator_type>     belos_recgmres_type; 
  typedef Belos::GmresPolySolMgr<scalar_type, multivector_type, operator_type>  belos_hybridgmres_type; 
  typedef Belos::PseudoBlockTFQMRSolMgr<scalar_type, multivector_type, operator_type> belos_psedotfqmr_type;  
//  typedef Belos::PseudoBlockCGSolMgr<scalar_type, multivector_type, operator_type> belos_psedobcg_type;  

  typedef MueLu::TpetraOperator<scalar_type, local_ordinal_type, global_ordinal_type, node_type> muelu_tpetra_operator_type;
  
//  typedef Ifpack2::Preconditioner<scalar_type, local_ordinal_type, global_ordinal_type, node_type> prec_type;

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
  using Teuchos::Time;
  using Teuchos::TimeMonitor;

  using Teuchos::updateParametersFromXmlFile;
  using Teuchos::updateParametersFromXmlString; 
  typedef Tpetra::CrsMatrix<> crs_matrix_type;

  Teuchos::oblackholestream blackHole;
// Teuchos::GlobalMPISession mpiSession (NULL,NULL, &blackHole);

  std::clock_t c_start = std::clock();
  
  Teuchos::GlobalMPISession mpiSession();
  RCP<const Teuchos::Comm<int> > comm = Tpetra::DefaultPlatform::getDefaultPlatform ().getComm ();
  const int myRank = comm->getRank ();
  const int numProcs = comm->getSize ();
 
  RCP<Time> insertva     = TimeMonitor::getNewCounter ("InsertValues ");
  RCP<Time> FillTimer    = TimeMonitor::getNewCounter ("FillComplete A");
//  RCP<Time> makevec      = TimeMonitor::getNewCounter ("Make Vecs ");
//  RCP<Time> makeprb      = TimeMonitor::getNewCounter ("Problem setup ");
//  RCP<Time> SolsTimer    = TimeMonitor::getNewCounter ("SolverSetupTime ");

  std::ostream &out = std::cout;
  RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));

  const global_size_t numGlobalElements = nn;  // 50;

  const size_t numMyElements =  rows; //map->getNodeNumElements ();
//  int numMyElements =  rows; //map->getNodeNumElements ();
  ArrayRCP<size_t> glon = arcp<size_t> (numMyElements);

  const global_ordinal_type indexBase = 0;

  double Values[7];
  int   Indices[7];
  int NumEntries;
  int RowLess3, RowLess2, RowLess1;
  int RowPlus1, RowPlus2, RowPlus3;

  int i;
  LO lcR;
  int gblRow;

  std::clock_t c_start1 = std::clock(); 

  RCP<const map_type> map = rcp (new map_type (numGlobalElements, numMyElements, indexBase, comm));
  std::clock_t c_start2 = std::clock();  
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
// map->describe(*fos,Teuchos::VERB_EXTREME);

  for(i = 0; i < numMyElements; ++i) {
      glon[i] = glob[i] -1;  //myRank + k*numProcs;
  }

  RCP<crs_matrix_type> A (new crs_matrix_type (map, 7));

//  std::clock_t c_start1 = std::clock();

  {
  TimeMonitor monitor (*insertva);
  for(lcR = 0; lcR < static_cast<LO> (numMyElements);  ++lcR) {
    gblRow = glon[lcR]; //map->getGlobalElement (lcR);
    NumEntries = 0; 
    RowLess3 = gblRow + aloc[lcR][0];
    RowLess2 = gblRow + aloc[lcR][1];
    RowLess1 = gblRow + aloc[lcR][2];
    RowPlus1 = gblRow + aloc[lcR][3];
    RowPlus2 = gblRow + aloc[lcR][4];
    RowPlus3 = gblRow + aloc[lcR][5];

    if ((RowLess3 >=0 )&&  abs(aloc[lcR][0]) != 0.0 ){
      Values[NumEntries]  = Ao[lcR][0]; 
      Indices[NumEntries] = RowLess3; 
      NumEntries          = NumEntries + 1 ; 
    }

//      if (RowLess2 >=0 &&  RowLess2 !=GlobalRow) //aloc[i][1] != 0.0)
    if ((RowLess2 >=0 ) && abs(aloc[lcR][1]) != 0) { 
      Values[NumEntries]  = Ao[lcR][1];
      Indices[NumEntries] = RowLess2;
      NumEntries          = NumEntries + 1 ;
    }
  
    if ((RowLess1 >=0 ) && abs(aloc[lcR][2]) != 0)  {
      Values[NumEntries]  = Ao[lcR][2];
      Indices[NumEntries] = RowLess1;
      NumEntries          = NumEntries + 1 ;
    }
  
//       if (RowPlus1 < NumGlobalElements &&  RowPlus1 !=GlobalRow) //aloc[i][3] != 0.0)
    if ((RowPlus1 < numGlobalElements) &&  abs(aloc[lcR][3]) != 0)   {
       Values[NumEntries]  = Ao[lcR][4];
       Indices[NumEntries] = RowPlus1;
       NumEntries          = NumEntries + 1 ;
    }
    
   if ((RowPlus2 < numGlobalElements) &&  abs(aloc[lcR][4]) != 0)       {
      Values[NumEntries]  = Ao[lcR][5];
      Indices[NumEntries] = RowPlus2;
      NumEntries          = NumEntries + 1 ;
   }
   
   if ((RowPlus3 < numGlobalElements) && abs(aloc[lcR][5]) != 0)   {
      Values[NumEntries]  = Ao[lcR][6];
      Indices[NumEntries] = RowPlus3;
      NumEntries          = NumEntries + 1 ;
      }

      Values[NumEntries]  = Ao[lcR][3];
      Indices[NumEntries] = gblRow;  
      NumEntries          = NumEntries + 1 ;
      A->insertGlobalValues (gblRow, NumEntries, Values, Indices); 
//    std::cout << gblRow << " "  << NumEntries << " "<< Indices << std::endl; 
    }
  }

// return 0; 
  {
    TimeMonitor monitor (*FillTimer);
//  A->fillComplete ();
    A->fillComplete (map,map);
  }


// A->describe(*fos,Teuchos::VERB_EXTREME);

  std::clock_t c_start3 = std::clock();


//   TimeMonitor monitor(*makevec);
   RCP<multivector_type> x = rcp(new multivector_type(map,1));
   RCP<multivector_type> b = rcp(new multivector_type(map,1));

//  RCP<multivector_type> x(map);
//  RCP<multivector_type> b(map);

//  b.putScalar (STS::one ());
//  x->randomize();

   for (LO lclRow = 0; lclRow < static_cast<LO> (numMyElements);++lclRow) {
    const GO gblRow = map->getGlobalElement (lclRow);
    b->sumIntoGlobalValue(gblRow, 0, bo[lclRow]);
//   x->sumIntoGlobalValue(gblRow, 0, x_vec[lclRow]);
   }

   std::clock_t c_start4 = std::clock();

/*
  Teuchos::ParameterList paramList;
  paramList.set("verbosity", "low");
  paramList.set("max levels", 2);
  paramList.set("coarse: max size", 10);
  paramList.set("multigrid algorithm", "sa");
*/
  
  std::string xmlFileName = "test.xml";
  
   std::clock_t c_start5 = std::clock();
//   TimeMonitor monitor(*makeprb); 
   RCP<linear_problem_type> Problem = rcp(new linear_problem_type(A, x, b));
   Problem->setProblem();

  std::clock_t c_start6 = std::clock();

 //  TimeMonitor monitor(*SolsTimer) ;   
  RCP<ParameterList> belosList = rcp(new ParameterList());
//  belosList->set("Maximum Iterations",    itmax);    // Maximum number of iterations allowed
  belosList->set("Maximum Iterations",    20);    // Maximum number of iterations allowed
  belosList->set( "Num Blocks", 300);                // Maximum number of blocks in Krylov factorization
//  belosList->set("Maximum Restarts", 1000); 
//  belosList->set("Convergence Tolerance", tole);     // Relative convergence tolerance requested
  belosList->set("Convergence Tolerance", 1.0e-4);     // Relative convergence tolerance requested
//  "Implicit" means "the left-preconditioned approximate a.k.a. 'recursive' residual as computed by the Krylov method."
//  belosList->set ("Implicit Residual Scaling", "Norm of RHS"); 
//  "Explicit" means ||B - A*X||, the unpreconditioned, "exact"     // residual.
//  belosList->set ("Explicit Residual Scaling", "Norm of RHS");
// belosList->set ("Explicit Residual Scaling", "Norm of Initial Residual");
  belosList->set("Block Size",          1);
  belosList->set("Verbosity",             Belos::Errors + Belos::Warnings + Belos::StatusTestDetails + Belos::TimingDetails);
//  belosList->set("Verbosity",             Belos::Errors + Belos::Warnings);
  belosList->set("Output Frequency",      20);
  //  belosList->set("Output Style",          Belos::None);
  belosList->set("Output Style",          Belos::Brief);
//  belosList->set("Output Style",         Belos::General);
   belosList->set("Estimate Condition Number"  , true) ; 
  //  belosList->set("Implicit Residual Scaling", "None");



  RCP<belos_solver_manager_type> solver;
    std::string inputFile="input_param.xml";
  std::string vname, ortho; 
  bool prec;
  RCP<Teuchos::ParameterList> myParams = rcp(new Teuchos::ParameterList());
  Teuchos::updateParametersFromXmlFile(inputFile, myParams.ptr());
  vname = myParams->get<std::string>("Solver");
  prec  = myParams->get<bool>("Precond");
  
//  std::clock_t c_start4 = std::clock();
  if(vname=="GMRES") {
//    ortho = myParams->get<std::string>("Orthogonal"); 
    if(prec) {
    RCP< mtoperator > mueLuPreconditioner = MueLu::CreateTpetraPreconditioner<scalar_type,local_ordinal_type,global_ordinal_type,node_type>(RCP<operator_type>(A), xmlFileName);
    Problem->setRightPrec (mueLuPreconditioner);
    }

/*    if(ortho=="IMGS") {
    belosList->set("Output Style", "IMGS");
    }
    else if(ortho=="ICGS") {
    belosList->set("Output Style", "ICGS");
    }
    else {
   belosList->set("Output Style", "DGKS");
    }
*/
    solver = rcp(new belos_stdgmres_manager_type(Problem, belosList));
    }

  else if(vname=="FlexGMRES") {
    if(prec) {
    RCP< mtoperator > mueLuPreconditioner = MueLu::CreateTpetraPreconditioner<scalar_type,local_ordinal_type,global_ordinal_type,node_type>(RCP<operator_type>(A), xmlFileName);
    Problem->setRightPrec (mueLuPreconditioner);
    }
    solver = rcp(new belos_flexgmres_manager_type(Problem, belosList));
    }
   
    else if(vname=="CG")  {
    if(prec) {
    RCP< mtoperator > mueLuPreconditioner = MueLu::CreateTpetraPreconditioner<scalar_type,local_ordinal_type,global_ordinal_type,node_type>(RCP<operator_type>(A), xmlFileName);
    Problem->setRightPrec (mueLuPreconditioner);
             }
    solver = rcp(new belos_blockcg_manager_type(Problem, belosList));
    }
    else if(vname=="TFQMR")     {
     solver = rcp(new belos_tfqmr_manager_type(Problem, belosList));
   if(prec) {
    RCP< mtoperator > mueLuPreconditioner = MueLu::CreateTpetraPreconditioner<scalar_type,local_ordinal_type,global_ordinal_type,node_type>(RCP<operator_type>(A), xmlFileName);
    Problem->setRightPrec (mueLuPreconditioner);
           }
    }
    else if(vname=="LSQR")     {
     solver = rcp(new  belos_lsqr_manager_type(Problem, belosList));
    }

    else if(vname=="RecyleGMRES") { 
      solver = rcp(new belos_recgmres_type (Problem, belosList)); 
    if(prec) {
    RCP< mtoperator > mueLuPreconditioner = MueLu::CreateTpetraPreconditioner<scalar_type,local_ordinal_type,global_ordinal_type,node_type>(RCP<operator_type>(A), xmlFileName);
    Problem->setRightPrec (mueLuPreconditioner);
           }
    }
    
    else if(vname=="HybridGMRES") { 
       solver = rcp(new belos_hybridgmres_type (Problem, belosList)); 
      if(prec) {
    RCP< mtoperator > mueLuPreconditioner = MueLu::CreateTpetraPreconditioner<scalar_type,local_ordinal_type,global_ordinal_type,node_type>(RCP<operator_type>(A), xmlFileName);
    Problem->setLeftPrec (mueLuPreconditioner);
    }

    } 
    else if (vname=="PsedoTFQMR") { 
       solver = rcp (new belos_psedotfqmr_type (Problem, belosList)); 
    if(prec) {
    RCP< mtoperator > mueLuPreconditioner = MueLu::CreateTpetraPreconditioner<scalar_type,local_ordinal_type,global_ordinal_type,node_type>(RCP<operator_type>(A), xmlFileName);
    Problem->setRightPrec (mueLuPreconditioner);
    }

    }
    else    {
    if(prec) {
    RCP< mtoperator > mueLuPreconditioner = MueLu::CreateTpetraPreconditioner<scalar_type,local_ordinal_type,global_ordinal_type,node_type>(RCP<operator_type>(A), xmlFileName);
//    Problem->setLeftPrec (mueLuPreconditioner);
    Problem->setRightPrec (mueLuPreconditioner);
    }
   solver = rcp(new belos_bicgstab_manager_type(Problem, belosList));
//    solver = rcp(new belos_psedobcg_type(Problem, belosList));
  }
 

  solver->solve();
//  solver->isLOADetected(); 
  
//  Belos::ReturnType(solver->solve());

// Perform solve
//      Belos::ReturnType ret = solver->solve();
//      TEST_EQUALITY(ret, Belos::Converged);
 
  std::clock_t c_start7 = std::clock();

  if(myRank == 0 ) {

  std::cout << "Initializing arrays "  << 1000*(c_start1-c_start)/CLOCKS_PER_SEC   << std::endl; 
  std::cout << "Creating Map  "        << 1000*(c_start2-c_start1)/CLOCKS_PER_SEC   << std::endl; 
  std::cout << "Make vecs  "           << 1000*(c_start4-c_start3)/CLOCKS_PER_SEC   << std::endl; 
  std::cout << "IntiMapFormMatVec "    << 1000*(c_start5-c_start1)/CLOCKS_PER_SEC   << std::endl;
  std::cout << "Problem Setup  "       << 1000*(c_start6-c_start5)/CLOCKS_PER_SEC   << std::endl; 
  std::cout << "Solver  "              << 1000*(c_start7-c_start6)/CLOCKS_PER_SEC   << std::endl; 
  
/*  std::cout << "before fill "   << 1000*(c_start1-c_start)/CLOCKS_PER_SEC   << std::endl; 
  std::cout << " Fill "         << 1000*(c_start2-c_start1)/CLOCKS_PER_SEC << std::endl; 
  std::cout << " Fill b "       << 1000*(c_start3-c_start2)/CLOCKS_PER_SEC << std::endl; 
  std::cout << " before Solve " << 1000*(c_start4-c_start3)/CLOCKS_PER_SEC << std::endl; 
  std::cout << "Solve "         << 1000*(c_start5-c_start4)/CLOCKS_PER_SEC << std::endl; 
  std::cout << "total"          << 1000*(c_start5-c_start)/CLOCKS_PER_SEC << std::endl;
*/
  int numIterations = solver->getNumIters();
  std::cout << " number of iterations" <<  numIterations << std::endl;
  double tolach = solver->achievedTol( );
  std::cout << " tolerance achieved" <<  tolach << std::endl;

  }


/*
  RCP<prec_type> prec;
  Ifpack2::Factory factory;
  // Set up the preconditioner of the given type.
  std::string  precondType="ILUT"; 
  prec = factory.create (precondType, A);

  const double fillLevel = 2.0;
  const double dropTol = 0.0;
  const double absThreshold = 0.1;

  RCP<Teuchos::ParameterList> plist  = rcp(new Teuchos::ParameterList());

  plist->set ("fact: ilut level-of-fill", fillLevel);
  plist->set ("fact: drop tolerance", dropTol);
  plist->set ("fact: absolute threshold", absThreshold);

  prec.setParameters (plist);

  prec.initialize();
  prec.compute();

  magnitude_type condest = STM::one();
  condest = prec.computeCondEst (Ifpack2::Cheap);
  out << endl << "Ifpack2 preconditioner's estimated condition number: " << condest << endl;
*/

// std::cout << "condition number" << solver->getConditionEstimate() << endl; 

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

}  // end of cpp wrapper 

// *fos << "LHS :" << std::endl;
//  x->describe(*fos,Teuchos::VERB_EXTREME);

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

