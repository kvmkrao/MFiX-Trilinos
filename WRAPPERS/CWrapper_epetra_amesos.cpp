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

  std::cout << "CWrapper_total_time in sec "      << double(c_start4-c_start1)/CLOCKS_PER_SEC   << std::endl;
  std::cout << "CWrapper_form_arrays  in sec "    << double(c_start2-c_start1)/CLOCKS_PER_SEC << std::endl;
  std::cout << "CppWrapper_total time in sec "    << double(c_start3-c_start2)/CLOCKS_PER_SEC << std::endl;

  }
}

#include "AztecOO.h"

#ifdef EPETRA_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "Epetra_Map.h"
//#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_CrsMatrix.h"
//#include "Epetra_FECrsMatrix.h"
#include "Epetra_Time.h"
#include "Teuchos_ParameterList.hpp"
#include<iostream>
#include<fstream>
#include<string>

#include "EpetraExt_RowMatrixOut.h"
#include "EpetraExt_OperatorOut.h"

#include "Epetra_FEVector.h"
#include "Epetra_Util.h"

#include "Epetra_Export.h"
#include "Epetra_Import.h"

#include "ml_include.h"
#include "ml_MultiLevelPreconditioner.h"

using namespace std;
using namespace Teuchos;

#include "Trilinos_Util.h"
#include "Teuchos_RCP.hpp"
//#include "Epetra_Export.h"


// Amesos
#include "Amesos_ConfigDefs.h"
#include "Amesos.h"
#include "Amesos_Klu.h"

#include "Amesos_Klu.h"
#include "Amesos_TestRowMatrix.h"

  double* CPPWrapper(double** Ao,double* bo,double* x_vec,double** aloc,int* globn,int rows,int cols,int tnrows,int itmax,double tole) 
  {

#ifdef HAVE_MPI
  Epetra_MpiComm Comm (MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

#ifdef EPETRA_NO_32BIT_GLOBAL_INDICES
  // Epetra was compiled only with 64-bit global index support, so use 64-bit global indices.
  typedef long long global_ordinal_type;
#else
  // Epetra was compiled with 32-bit global index support.  If
  // EPETRA_NO_64BIT_GLOBAL_INDICES is defined, it does not also
  // support 64-bit indices.
  typedef int global_ordinal_type;
  typedef int local_ordinal_type;
#endif // EPETRA_NO_32BIT_GLOBAL_INDICES

std::ostream &out = std::cout;

// std::cout << NumMyElements  << NumGlobalElements << std::endl; 

//  Epetra_Map Map(-1,MyGlobalElements.Length(), MyGlobalElements.Values(),0, Comm);
//   const global_ordinal_type indexBase = 0;
    const int indexBase = 0;

//   const global_ordinal_type  NumGlobalElements = tnrows;
    const int NumGlobalElements = tnrows;

    const int  NumMyElements = rows;
//  const global_ordinal_type   NumMyElements = rows;

//  global_ordinal_type* glon = new global_ordinal_type [NumMyElements];
   int * glon = new int[NumMyElements];
 
   for (int k = 0; k < NumMyElements; ++k) {
    glon[k] = globn[k] - 1;
//    std::cout << Comm.MyPID() << " " << k << " " << glon[k] <<  std::endl; 
  }

//     Epetra_Map Map(-1, NumMyElements, glon, 0 ,Comm);

      Epetra_Map Map(NumGlobalElements, NumMyElements, 0 ,Comm);


  int * NumNz = new int[NumMyElements];
//  global_ordinal_type* NumNz = new global_ordinal_type [NumMyElements];

  int glbrow, RowLess2, RowLess1, RowLess3;
  int RowPlus2, RowPlus1, RowPlus3;
  int gid ;
  int NumEntries;


     Epetra_CrsMatrix A(Copy,Map,7);
//     Epetra_CrsMatrix Mat(Copy,Map,7);

//   Epetra_FECrsMatrix A(Copy,Map,NumNz);
//   Epetra_FECrsMatrix A(Copy,Map,7);

//std::cout << "Created crs_matrix" << std::endl; 

 int lclerr = 0; 
// Teuchos::RefCountPtr<Epetra_CrsMatrix> A = Teuchos::rcp( new Epetra_CrsMatrix(Copy,Map,NumNz) );


  double Values[7];
  int Indices[7];
//  double *Values = new double[7];
//  int *Indices = new int[7];

  double two = -1.0;
  int nend = NumMyElements; 
  int nstart= 0; 

//   if(Comm.MyPID() != Comm.NumProc() -1) nend = nend - 22;

/*
 if(Comm.MyPID()==0) { 
  ofstream myfile ("fort800.txt");
  if (myfile.is_open())
  {
  for( int i=0 ; i<NumMyElements; ++i ) {
      myfile << glon[i]+1 << " " << bo[i] << std::endl;
  }
  }
  myfile.close();
  }
*/


  std::clock_t c_start1 = std::clock(); 

  for( int i=0 ; i<NumMyElements ; ++i) {
//  std::cout << " " << aloc[i][0] << " " << aloc[i][1] << " "<< aloc[i][3]<<" "<<aloc[i][4] << std::endl;
//    int gid = glon[i]; //Map.GID(i);
      gid = glon[i] ; //Map.GID(i);
        RowLess3 = gid + aloc[i][0];
        RowLess2 = gid + aloc[i][1];
        RowLess1 = gid + aloc[i][2];
        RowPlus1 = gid + aloc[i][3];
        RowPlus2 = gid + aloc[i][4];
        RowPlus3 = gid + aloc[i][5];

      NumEntries = 0 ;
      if (RowLess3 >=0 &&  abs(aloc[i][0]) != 0 )
      {
      Values[NumEntries]  = Ao[i][0];
      Indices[NumEntries] = RowLess3;
      NumEntries          = NumEntries + 1 ;
       }

      if (RowLess2 >=0 &&  abs(aloc[i][1]) != 0 )
      {
      Values[NumEntries]  = Ao[i][1];
      Indices[NumEntries] = RowLess2;
      NumEntries          = NumEntries + 1 ;
       }
       if (RowLess1 >=0  &&  abs(aloc[i][2]) != 0)
       {
       Values[NumEntries]  = Ao[i][2];
       Indices[NumEntries] = RowLess1;
       NumEntries          = NumEntries + 1 ;
       }
       if (RowPlus1 < NumGlobalElements && abs(aloc[i][3]) !=0)
       {
      Values[NumEntries]  = Ao[i][4];
      Indices[NumEntries] = RowPlus1;
      NumEntries          = NumEntries + 1 ;
       }

       if (RowPlus2 < NumGlobalElements  && abs(aloc[i][4]) !=0)
       {
       Values[NumEntries]  = Ao[i][5];
       Indices[NumEntries] = RowPlus2;
       NumEntries          = NumEntries + 1 ;
       }

      if (RowPlus3 < NumGlobalElements  && abs(aloc[i][5]) !=0)
       {
       Values[NumEntries]  = Ao[i][6];
       Indices[NumEntries] = RowPlus3;
       NumEntries          = NumEntries + 1 ;
       }

      Values[NumEntries]  = Ao[i][3];
      Indices[NumEntries] = gid;
      NumEntries          = NumEntries + 1 ;
      A.InsertGlobalValues (gid, NumEntries, Values, Indices);
  }

   std::clock_t c_start2 = std::clock();
   A.FillComplete(Map,Map);
   std::clock_t c_start3 = std::clock();

   A.FillComplete(Map,Map);

  int NumVectors = 1;
  Epetra_MultiVector x(Map,NumVectors);
  Epetra_MultiVector b(Map,NumVectors);

//    Epetra_MultiVector RHS(Map,1,false);
//    Epetra_MultiVector LHS(Map,1,true);
//  RHS.PutScalar(1.0);  

    for( int i=0 ; i< NumMyElements; ++i ) {
      glbrow = glon[i];
//      RHS->ReplaceGlobalValue(glbrow, 0, 1.0);
      b.SumIntoGlobalValue(glbrow, 0, bo[i]);
  }

  x.Random();
//  A.Multiply(false,x_exact,b);

  Epetra_LinearProblem Problem(&A, &x, &b);

  Amesos_BaseSolver* Solver;

  // Initializes the Factory. Factory is a function class (a
  // class that contains methods only, no data). Factory
  // will be used to create Amesos_BaseSolver derived objects.
  //
  Amesos Factory;

  // Specifies the solver. String ``SolverType'' can assume one 
  // of the following values:
  // - Lapack
  // - Klu
  // - Umfpack
  // - Pardiso
  // - Taucs
  // - Superlu
  // - Superludist
  // - Mumps
  // - Dscpack
  // 
  std::string SolverType = "Klu";
  Solver = Factory.Create(SolverType, Problem);

  Teuchos::ParameterList List;
  List.set("PrintTiming", true);
  List.set("PrintStatus", true);
  Solver->SetParameters(List);

  std::clock_t c_start4 = std::clock();
  Solver->SymbolicFactorization();
  Solver->NumericFactorization();
  Solver->Solve();

 
  std::clock_t c_start5 = std::clock();
  if (Comm.MyPID() == 0) std::cout << "Solver time in sec  "  << double(c_start5-c_start4)/CLOCKS_PER_SEC   << std::endl;
  if (Comm.MyPID() == 0) std::cout << "InsertValues in sec "  << double(c_start2-c_start1)/CLOCKS_PER_SEC   << std::endl;
  if (Comm.MyPID() == 0) std::cout << "FillComplete in sec "  << double(c_start3-c_start2)/CLOCKS_PER_SEC   << std::endl;

 // you can get the timings here
  Teuchos::ParameterList TimingsList;
  Solver->GetTiming( TimingsList );

  // you can find out how much time was spent in ...
  double sfact_time, nfact_time, solve_time;
  double mtx_conv_time, mtx_redist_time, vec_redist_time;

  // 1) The symbolic factorization 
  //    (parameter doesn't always exist)
  sfact_time = TimingsList.get( "Total symbolic factorization time", 0.0 );
  if (Comm.MyPID() == 0) std::cout << " Total symbolic factorization time" << sfact_time << endl;
  // 2) The numeric factorization 
  //    (always exists if NumericFactorization() is called)

//  nfact_time = Teuchos::getParameter<double>( TimingsList, "Total numeric factorization time" );
//  if (Comm.MyPID() == 0) std::cout << "Total numeric factorization time " << nfact_time << endl;

  // 3) Solving the linear system 
  //    (always exists if Solve() is called)
//   solve_time = Teuchos::getParameter<double>( TimingsList, "Total solve time" );
//   if (Comm.MyPID() == 0) std::cout << "Total solve time" << solve_time << std::endl; 

  // 4) Converting the matrix to the accepted format for the solver
  //    (always exists if SymbolicFactorization() is called)
//   mtx_conv_time = Teuchos::getParameter<double>( TimingsList, "Total solve time" );
//   if (Comm.MyPID() == 0) std::cout << "Total solve time" <<  mtx_conv_time << std::endl; 

  // 5) Redistributing the matrix for each solve to the accepted format for the solver

//   mtx_redist_time = TimingsList.get( "Total matrix redistribution time", 0.0 );
//   if (Comm.MyPID() == 0) std::cout << "Total matrix redistribution time" << mtx_redist_time << std::endl;

  // 6) Redistributing the vector for each solve to the accepted format for the solver

//  vec_redist_time = TimingsList.get( "Total vector redistribution time", 0.0 );
//  if (Comm.MyPID() == 0) std::cout << "Total vector redistribution time" << vec_redist_time  << endl; 


//  double norm = ComputeNorm(Matrix, &x_exact, &b);
//  if (Comm.MyPID() == 0)
//    std::cout << "norm = " << norm << std::endl;


/*
   Epetra_LinearProblem Problem;
   Problem.SetOperator(&A);

  Amesos_BaseSolver * Solver;
  Amesos Factory;
      LHS.PutScalar(1.0);
      // 2.- create corresponding rhs
      Mat->Multiply(false, LHS, RHS);

      // 3.- randomize solution vector
      LHS.Random();

      // 4.- create the amesos solver object
//      Amesos_BaseSolver* Solver = Factory.Create("Amesos_KLU", Problem);
//      Solver->SetParameters(List);
//      Solver->SetUseTranspose( true) ;


  Solver = Factory.Create("Amesos_Klu", Problem);
  Solver->SymbolicFactorization();
  Solver->SetUseTranspose( true) ;
*/

/*
    ParameterList List;
    Epetra_LinearProblem Problem(&A, &LHS, &RHS);

      Teuchos::ParameterList ParamList ;
      Amesos_Klu klu( Problem ) ;
      ParamList.set( "MaxProcs", -3 );
      klu.SetParameters( ParamList );
      klu.SetUseTranspose( &A ) ;

      klu.SymbolicFactorization(  ) ;
      klu.NumericFactorization(  ) ;
      klu.Solve( ) ;
*/

/*
 ParameterList List;
 Epetra_LinearProblem problem(&A, &LHS, &RHS);

// ----------  ML Preconditioner ------------------
  ParameterList MLList;


 ML_Epetra::SetDefaults("SA",MLList);
// overwrite with user’s defined parameters
 MLList.set("max levels",6);
 MLList.set("increasing or decreasing","decreasing");
 MLList.set("aggregation: type", "Uncoupled");
// MLList.set("aggregation: type", "MIS");
 MLList.set("coarse: type","Amesos-KLU");


  ML_Epetra::MultiLevelPreconditioner* MLPrec =
    new ML_Epetra::MultiLevelPreconditioner(A, MLList, true);

  AztecOO Solver(problem);
  Solver.SetPrecOperator(MLPrec);
  std::cout << " Build Preconditioner" << std::endl; 
  
  Solver.SetAztecOption(AZ_solver, AZ_bicgstab);
//  Solver.SetAztecOption(AZ_solver, AZ_gmres);
//            Solver.SetAztecOption(AZ_precond, AZ_dom_decomp);
//  Solver.SetAztecOption(AZ_subdomain_solve, AZ_bilu);
//          Solver.SetAztecOption(AZ_subdomain_solve, AZ_ilu);
//  Solver.SetAztecOption(AZ_subdomain_solve, AZ_ilut);
//  Solver.SetAztecOption(AZ_subdomain_solve, AZ_rilu);
//   Solver.SetAztecOption(AZ_omega, 0.7); 
//  Solver.SetAztecOption(AZ_subdomain_solve, AZ_rilu);
  Solver.SetAztecOption(AZ_kspace, 200);
//  Solver.SetAztecOption(AZ_output, tolr);
//  Solver.Iterate(nitr, tolr);
   Solver.Iterate(itmax, tole);

*/

  double xnew[NumMyElements];
  x.ExtractCopy(xnew , 1);

  for(int i=0; i < NumMyElements ; i++)
  {
  x_vec[i] = xnew[i];
 }
return 0;
  }  // end of cpp wrapper 


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

