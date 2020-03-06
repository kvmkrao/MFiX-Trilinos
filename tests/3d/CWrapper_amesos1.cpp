#include <iostream>
#include "stdio.h"


extern "C"
{
  double* CPPWrapper(double** A, double* b, double* x_vec, double** aloc, int* globn, int rows, int cols, int tnrows);
  void cwrapper_(double* br, double* xr, double* alc, double* A, int* globn, int *nAx, int *nAy, int *nrows, char Str[7]) // note extra underscore to keep linker happy
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
//  std::cout << globn[i] << std::endl; 
   }
 

/*  for (int i=0; i<m; i++)
   {
    globn[i] = globn[i] - 1; 
//  std::cout << globn[i] << std::endl; 
   }

*/


//  std::cout << alc[0] << alc[1]<< alc[2]  << alc[3] << alc[4] << alc[5] << std::endl; 
  CPPWrapper(Amat,br,xr,loc,globn,m,n,*nrows);

}
}
//##################################################################################

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
// The C++ interface of ML (more precisely,
// ML_Epetra::MultiLevelPreconditioner), requires Trilinos to be
// configured with --enable-epetra --enable-teuchos. This example
// required --enable-galeri (for the definition of the linear systems)

//#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_AZTECOO) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_GALERI) && defined(HAVE_ML_IFPACK)

//#include "Epetra_Vector.h"
//#include "Epetra_LinearProblem.h"
#include "ml_MultiLevelPreconditioner.h"

using namespace std;

using namespace Teuchos;
//using namespace Galeri;

#include "Trilinos_Util.h"
//#include "Trilinos_Util_ReadMatrixMarket2Epetra.h"
//#include "CrsMatrixTranspose.h"
#include "Teuchos_RCP.hpp"
//#include "Epetra_Export.h"

// Amesos
#include "Amesos_ConfigDefs.h"
#include "Amesos.h"

double* CPPWrapper(double** Ao, double* bo, double* x_vec ,double** aloc,int* globn, int rows, int cols,int tnrows) // A[n*m]x[m] = b[n]
  {

#ifdef HAVE_MPI
//  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm (MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

#ifdef EPETRA_NO_32BIT_GLOBAL_INDICES
  // Epetra was compiled only with 64-bit global index support, so use
  // 64-bit global indices.
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

//  Epetra_Map Map(-1,MyGlobalElements.Length(),
//         MyGlobalElements.Values(),0, Comm);

   const global_ordinal_type indexBase = 0;
//   const int indexBase = 0;

  const global_ordinal_type  NumGlobalElements = tnrows;
//   const int NumGlobalElements = tnrows;

//    const int  NumMyElements = rows;
  const global_ordinal_type   NumMyElements = rows;

  global_ordinal_type* glon = new global_ordinal_type [NumMyElements];
//   int * glon = new int[NumMyElements];
 
   for (int k = 0; k < NumMyElements; ++k) {
    glon[k] = globn[k] - 1;
//    std::cout << k << " " << glon[k] << std::endl; 
  }
  
   Epetra_Map Map(-1, NumMyElements, glon, 0 ,Comm);
//   Epetra_Map map1to1 =Epetra_Util::Create_OneToOne_Map(Map);

//   cout << Map ; 

// int * NumNz = new int[NumMyElements];
//  global_ordinal_type* NumNz = new global_ordinal_type [NumMyElements];


/*
  for ( int i=0; i<NumMyElements; i++){
//  std::cout << " proc " << Comm.MyPID() << "i "<< i << " " << glon[i] << std::endl; 
    if ( glon[i] == 0) {  //  || (glon[i] == NumGlobalElements-1)) {
      NumNz[i] = 2;
     }
    else if (glon[i] == NumGlobalElements-1) 
      NumNz[i] = 2;
    else {
      NumNz[i] = 3;
    }
  }
*/

//   Epetra_CrsMatrix A(Copy,Map,NumNz);
//   Epetra_FECrsMatrix A(Copy,Map,5);
//     Epetra_FECrsMatrix A(Copy,map1to1,5);

 int lclerr = 0; 
// Teuchos::RefCountPtr<Epetra_CrsMatrix> A = Teuchos::rcp( new Epetra_CrsMatrix(Copy,Map,NumNz) );
  Epetra_CrsMatrix A(Copy,Map,7);


//  double Values[5];
//  global_ordinal_type Indices[5];

  double *Values = new double[7];
  Values[0] = 0.0; Values[1] = 0.0;
  Values[2] = 0.0; Values[3] = 0.0;
  Values[4] = 0.0; Values[5] = 0.0;
  Values[6] = 0.0;
  int *Indices = new int[7];
  double two = 2.0;
  int NumEntries;
  int glbrow, RowLess2, RowLess1, RowLess3;
  int RowPlus2, RowPlus1, RowPlus3;

  for( int i=0 ; i< NumMyElements; ++i) {
//  std::cout << " " << aloc[i][0] << " " << aloc[i][1] << " "<< aloc[i][3]<<" "<<aloc[i][4] << std::endl;

    int gid = Map.GID(i);
//    std::cout << gid << " "<< glon[i] << std::endl; 
	RowLess3 = gid + aloc[i][0];
	RowLess2 = gid + aloc[i][1];
	RowLess1 = gid + aloc[i][2];
	RowPlus1 = gid + aloc[i][3];
	RowPlus2 = gid + aloc[i][4];
	RowPlus3 = gid + aloc[i][5];

	NumEntries = 0 ; 
      if (RowLess3 >=0 &&  abs(aloc[i][0]) != 0)
      {
//      Values[NumEntries]  = Ao[i][1];
      Indices[NumEntries] = RowLess3;
      NumEntries          = NumEntries + 1 ;
       }


      if (RowLess2 >=0 &&  abs(aloc[i][1]) != 0)
      {
//      Values[NumEntries]  = Ao[i][1];
      Indices[NumEntries] = RowLess2;
      NumEntries          = NumEntries + 1 ;
       }

       if (RowLess1 >=0  &&  abs(aloc[i][2]) != 0)
       {
//       Values[NumEntries]  = Ao[i][2];
       Indices[NumEntries] = RowLess1;
       NumEntries          = NumEntries + 1 ;
       }
//       A.insertGlobalValues(gblRow, 1, &Ao[i][3], &gblRow);
//       if (RowPlus1 < NumGlobalElements &&  RowPlus1 !=GlobalRow) //aloc[i][3] != 0.0)
       if (RowPlus1 < NumGlobalElements &&  abs(aloc[i][3]) != 0)
       {
//      Values[NumEntries]  = Ao[i][4];
      Indices[NumEntries] = RowPlus1;
      NumEntries          = NumEntries + 1 ;
       }

       if (RowPlus2 < NumGlobalElements &&  abs(aloc[i][4]) != 0)
       {
//       Values[NumEntries]  = Ao[i][5];
       Indices[NumEntries] = RowPlus2;
       NumEntries          = NumEntries + 1 ;
       }

     if (RowPlus3 < NumGlobalElements &&  abs(aloc[i][5]) != 0)
       {
//       Values[NumEntries]  = Ao[i][6];
       Indices[NumEntries] = RowPlus3;
       NumEntries          = NumEntries + 1 ;
       }

//      Values[NumEntries]  = Ao[i][3];
      Indices[NumEntries] = gid;
      NumEntries          = NumEntries + 1 ;
//       A.insertGlobalValues(gblRow, 1, &Ao[i][3], &gblRow);
//      A.SumIntoGlobalValues (gid, NumEntries, Values, Indices);
        A.InsertGlobalValues (gid, NumEntries, Values, Indices);
  }

//  	A.GlobalAssemble();
	A.FillComplete();

// Print the sparse matrix A
//  A.Print(std::cout); 
// EpetraExt::RowMatrixToMatrixMarketFile("strided.mm",*A);

//  Teuchos::RefCountPtr<Epetra_MultiVector> LHS = Teuchos::rcp( new Epetra_MultiVector(Map, 1) );
//  Teuchos::RefCountPtr<Epetra_MultiVector> RHS = Teuchos::rcp( new Epetra_MultiVector(Map, 1) );


 // ===================================================== //
  // B E G I N N I N G   O F  T H E   AM E S O S   P A R T //
  // ===================================================== //

  // For comments on the commands in this section, please
  // see file example_AmesosFactory.cpp.

  Epetra_LinearProblem Problem;

  Problem.SetOperator(&A);

  // Initializes Amesos solver. Here we solve with Amesos_Klu.

  Amesos_BaseSolver * Solver;
  Amesos Factory;
  Solver = Factory.Create("Amesos_Klu", Problem);

  Solver->SymbolicFactorization();

  for( int i=0 ; i< NumMyElements; ++i) {
//  std::cout << " " << aloc[i][0] << " " << aloc[i][1] << " "<< aloc[i][3]<<" "<<aloc[i][4] << std::endl;

    int gid = Map.GID(i);
//    std::cout << gid << " "<< glon[i] << std::endl; 
	RowLess3 = gid + aloc[i][0];
	RowLess2 = gid + aloc[i][1];
	RowLess1 = gid + aloc[i][2];
	RowPlus1 = gid + aloc[i][3];
	RowPlus2 = gid + aloc[i][4];
	RowPlus3 = gid + aloc[i][5];

	NumEntries = 0 ;

      if (RowLess3 >=0 &&  abs(aloc[i][0]) != 0)
      {
      Values[NumEntries]  = Ao[i][0];
      Indices[NumEntries] = RowLess3;
      NumEntries          = NumEntries + 1 ;
       }

      if (RowLess2 >=0 &&  abs(aloc[i][1]) != 0)
      {
      Values[NumEntries]  = Ao[i][1];
      Indices[NumEntries] = RowLess2;
      NumEntries          = NumEntries + 1 ;
//       A.insertGlobalValues(gblRow, 1, &Ao[i][1], &RowLess2);
       }
       if (RowLess1 >=0  &&  abs(aloc[i][2]) != 0)
       {
       Values[NumEntries]  = Ao[i][2];
       Indices[NumEntries] = RowLess1;
       NumEntries          = NumEntries + 1 ;
//       A.insertGlobalValues(gblRow, 1, &Ao[i][2], &RowLess1);
       }
//       A.insertGlobalValues(gblRow, 1, &Ao[i][3], &gblRow);
//       if (RowPlus1 < NumGlobalElements &&  RowPlus1 !=GlobalRow) //aloc[i][3] != 0.0)
       if (RowPlus1 < NumGlobalElements &&  abs(aloc[i][3]) != 0)
       {
      Values[NumEntries]  = Ao[i][4];
      Indices[NumEntries] = RowPlus1;
      NumEntries          = NumEntries + 1 ;
//       A.insertGlobalValues(gblRow, 1, &Ao[i][4], &RowPlus1);
       }

       if (RowPlus2 < NumGlobalElements &&  abs(aloc[i][4]) != 0)
       {
       Values[NumEntries]  = Ao[i][5];
       Indices[NumEntries] = RowPlus2;
       NumEntries          = NumEntries + 1 ;
//       A.insertGlobalValues(gblRow, 1, &Ao[i][5], &RowPlus2);
       }

       if (RowPlus3 < NumGlobalElements &&  abs(aloc[i][5]) != 0)
       {
       Values[NumEntries]  = Ao[i][6];
       Indices[NumEntries] = RowPlus3;
       NumEntries          = NumEntries + 1 ;
       }


      Values[NumEntries]  = -1.0; //Ao[i][3];
      Indices[NumEntries] = gid;
      NumEntries          = NumEntries + 1 ;
//       A.insertGlobalValues(gblRow, 1, &Ao[i][3], &gblRow);
//      A.SumIntoGlobalValues (gid, NumEntries, Values, Indices);
        A.ReplaceGlobalValues (gid, NumEntries, Values, Indices);
  }

   Solver->NumericFactorization();

    Epetra_MultiVector RHS(Map,1,false); 
    Epetra_MultiVector LHS(Map,1,true);

//    Epetra_FEVector RHS(Map,1,false); 
//      Epetra_FEVector RHS(map1to1,1,true);
//    Epetra_FEVector LHS(Map,1,true);

//     Epetra_FEVector LHS(map1to1,1,true);
//     Epetra_FEVector xv(Map);

//  Epetra_Vector x(Map);
//  Epetra_Vector b(Map);

//  RHS->PutScalar(0.0);  

  for( int i=0 ; i<NumMyElements; ++i ) {
      glbrow = glon[i];
//      RHS->ReplaceGlobalValue(glbrow, 0, 1.0);
        RHS.SumIntoGlobalValue(glbrow, 0, bo[i]);
//      RHS.ReplaceGlobalValue(glbrow, 0, bo[i]);
  }

  // Finally, we set up the LHS and the RHS vector (Random().
/*
  Epetra_Vector b(Map);
  b.Random();
  Epetra_Vector x(Map);
  x.PutScalar(0.0);
*/

  Problem.SetLHS(&LHS);
  Problem.SetRHS(&RHS);


// ParameterList List;
// Epetra_LinearProblem problem(&A, &LHS, &RHS);


  Solver->Solve();

  // Print out the timing information and get it from the solver
  Solver->PrintTiming();

/*
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
  Solver = Factory.Create(SolverType, problem);
*/

  // Factory.Create() returns 0 if the requested solver
  // is not available



// ----------  ML Preconditioner ------------------

/*
  ParameterList MLList;

 ML_Epetra::SetDefaults("SA",MLList);
// overwrite with user’s defined parameters
 MLList.set("max levels",6);
 MLList.set("increasing or decreasing","decreasing");
// MLList.set("aggregation: type", "Uncoupled");
 MLList.set("aggregation: type", "MIS");
 MLList.set("coarse: type","Amesos-KLU");
*/ 

/*  commented 
  ML_Epetra::SetDefaults("DD",MLList);
  MLList.set("smoother: pre or post", "both");
  MLList.set("PDE equations", 1);

  // fix the smoother to be IFPACK; can be set using (level X) syntax
  MLList.set("smoother: type","IFPACK");

  // now we have to specify which IFPACK preconditioner should be
  // built. Any value that is valid for the IFPACK factory. We also need
  // to define the overlap (>= 0).

  MLList.set("smoother: ifpack type", "ILU");
  MLList.set("smoother: ifpack overlap", 20);

  // Then, all parameters can will control the definition of the IFPACK
  // smoother are inserted in IFPACKList. In this case, we specify the fill-in
  // factor. For a list of supported parameters, please consult the IFPACK
  // documentation. For example, IFPACK preconditioner "Amesos" or
  // "Amesos stand-alone" can be used to solve with an LU
  // factorization on each domain.
  MLList.sublist("smoother: ifpack list").set("fact: level-of-fill", 5);
*/
  // we can now build the preconditioner...

/*
  ML_Epetra::MultiLevelPreconditioner* MLPrec =
    new ML_Epetra::MultiLevelPreconditioner(A, MLList, true);

  // ... and solve the linear system

  AztecOO Solver(problem);
  Solver.SetPrecOperator(MLPrec);
  Solver.SetAztecOption(AZ_solver, AZ_bicgstab);
//  Solver.SetAztecOption(AZ_solver, AZ_gmres);
//  Solver.SetAztecOption(AZ_precond, AZ_dom_decomp);
//  Solver.SetAztecOption(AZ_subdomain_solve, AZ_ilut);
  Solver.SetAztecOption(AZ_kspace, 200);
  Solver.SetAztecOption(AZ_output, 5);
  Solver.Iterate(20, 1e-4);

*/
//LHS.Print(std::cout);
//  Epetra_Import Importer(Map,map1to1);
//  xv.Import(LHS,Importer,Insert);

  double xnew[NumMyElements];
  LHS.ExtractCopy(xnew , 1);
//  xv.ExtractCopy(xnew , 1);

  for(int i=0; i < NumMyElements ; i++)
  {
  x_vec[i] = xnew[i];
 }

return 0;
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
    
