#include <iostream>
#include "stdio.h"

extern "C"
{
  double* CPPWrapper(double** A, double* b, double* x_vec, double** aloc, int* globn, int rows, int cols, int tnrows);
  void cwrapper_(double* br, double* xr, double* alc, double* A, int *globn, int *nAx, int *nAy, int *nrows, char Str[7]) // note extra underscore to keep linker happy
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
 

  for (int i=0; i<m; i++)
   {
    globn[i] = globn[i] - 1; 
//  std::cout << globn[i] << std::endl; 
   }


//  std::cout << m << " "  << n << " " << *nrows << std::endl; 
//std::cout << "time spent in Cwrapper" << " "<< clock() - time1  << std::endl; 

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
#include "Epetra_FECrsMatrix.h"
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

double* CPPWrapper(double** Ao, double* bo, double* x_vec ,double** aloc,int* glon, int rows, int cols,int tnrows) // A[n*m]x[m] = b[n]
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

//   const global_ordinal_type indexBase = 0;
   const int indexBase = 0;

//   const global_ordinal_type  NumGlobalElements = tnrows;
   const int NumGlobalElements = tnrows;

  const int  NumMyElements = rows;
//  const global_ordinal_type   NumMyElements = rows;

//  global_ordinal_type* glon = new global_ordinal_type [NumMyElements];
//    int * glon = new int[NumMyElements];

//    std::cout << "CPP Wrapper  "  << Comm.MyPID() << " " << NumGlobalElements << " " << rows << std::endl;
//     return 0 ; 


/*
   for (int k = 0; k < NumMyElements; k++) {
//    glon[k] = globn[k] ;
//    if(Comm.MyPID() == 0) 
   std::cout << Comm.MyPID() << "  " << k << " " << glon[k] << std::endl; 
  }
*/
  
   Epetra_Map Map(-1, NumMyElements, glon, 0 ,Comm);
   Epetra_Map map1to1 =Epetra_Util::Create_OneToOne_Map(Map);
 
//  cout << Map ; 
// cout << map1to1 ; 
  Map.Print(std::cout);

// return 0; 
//  void terminate();
//  std::terminate();    
   
// int * NumNz = new int[NumMyElements];
//  global_ordinal_type* NumNz = new global_ordinal_type [NumMyElements];

//   Epetra_CrsMatrix A(Copy,Map,NumNz);
//   Epetra_FECrsMatrix A(Copy,Map,5);
     Epetra_FECrsMatrix A(Copy,map1to1,5);

  // Create an Epetra_CrsMatrix using the Map, with dynamic allocation.
 int lclerr = 0; 
// Teuchos::RefCountPtr<Epetra_CrsMatrix> A = Teuchos::rcp( new Epetra_CrsMatrix(Copy,Map,NumNz) );
// Epetra_CrsMatrix A(Copy,Map,5);


//  double Values[5];
//  global_ordinal_type Indices[5];

  double *Values = new double[7];
//  Values[0] = -1.0; Values[1] = -1.0;
  int *Indices = new int[7];
  int NumEntries;
  int glbrow, RowLess2, RowLess1;
  int RowPlus2, RowPlus1;
  int RowPlus3, RowLess3;

  int istart, iend; 
  std::cout <<  "CPP wrapper" << Comm.MyPID() << Comm.NumProc() << std::endl; 

  for( int i=0 ; i< NumMyElements ; ++i) {

//    int gid = glon[i]; //Map.GID(i);
    int gid = Map.GID(i);
//    std::cout << gid << " "<< glon[i] << std::endl; 
	RowLess3 = gid + aloc[i][0];
	RowLess2 = gid + aloc[i][1];
	RowLess1 = gid + aloc[i][2];
	RowPlus1 = gid + aloc[i][3];
	RowPlus2 = gid + aloc[i][4];
	RowPlus3 = gid + aloc[i][5];

	NumEntries = 0 ;
//      if ((RowLess3 >=0 ) &&  abs(aloc[i][0]) != 0) {
//      Values[NumEntries]  = Ao[i][0];
//      Indices[NumEntries] = RowLess3;
//      NumEntries          = NumEntries + 1 ;
//       }

      if((RowLess2 >=0 ) &&  abs(aloc[i][1]) != 0)  {
      Values[NumEntries]  = Ao[i][1];
      Indices[NumEntries] = RowLess2;
      NumEntries          = NumEntries + 1 ;
       }

       if ((RowLess1 >=0 ) &&  abs(aloc[i][2]) != 0) {
       Values[NumEntries]  = Ao[i][2];
       Indices[NumEntries] = RowLess1;
       NumEntries          = NumEntries + 1 ;
//       A.insertGlobalValues(gblRow, 1, &Ao[i][2], &RowLess1);
       }
       if ((RowPlus1 < NumGlobalElements) &&  abs(aloc[i][3]) != 0) {
      Values[NumEntries]  = Ao[i][4];
      Indices[NumEntries] = RowPlus1;
      NumEntries          = NumEntries + 1 ;
//       A.insertGlobalValues(gblRow, 1, &Ao[i][4], &RowPlus1);
       }

       if ((RowPlus2 < NumGlobalElements) &&  abs(aloc[i][4]) != 0) {
       Values[NumEntries]  = Ao[i][5];
       Indices[NumEntries] = RowPlus2;
       NumEntries          = NumEntries + 1 ;
//       A.insertGlobalValues(gblRow, 1, &Ao[i][5], &RowPlus2);
       }

//       if ((RowPlus3 < NumGlobalElements) &&  abs(aloc[i][5]) != 0) {
//       Values[NumEntries]  = Ao[i][6];
//       Indices[NumEntries] = RowPlus3;
//       NumEntries          = NumEntries + 1 ;
////       A.insertGlobalValues(gblRow, 1, &Ao[i][5], &RowPlus2);
//       }
 

      Values[NumEntries]  = Ao[i][3];
      Indices[NumEntries] = gid;
      NumEntries          = NumEntries + 1 ;
//       A.insertGlobalValues(gblRow, 1, &Ao[i][3], &gblRow);
//      A.InsertGlobalValues (gid, NumEntries, Values, Indices);
      A.SumIntoGlobalValues (gid, NumEntries, Values, Indices);
      std::cout << "global row" << gid << " " << std::endl; 
  }
  A.GlobalAssemble();
  A.FillComplete();


// Print the sparse matrix A
//   A.Print(std::cout); 
//   return 0; 

// EpetraExt::RowMatrixToMatrixMarketFile("Sparsemat.mm",A);

//  Teuchos::RefCountPtr<Epetra_MultiVector> LHS = Teuchos::rcp( new Epetra_MultiVector(Map, 1) );
//  Teuchos::RefCountPtr<Epetra_MultiVector> RHS = Teuchos::rcp( new Epetra_MultiVector(Map, 1) );

//    Epetra_Vector RHS(Map);
//    Epetra_MultiVector RHS(Map,1,false); 
//    Epetra_MultiVector RHS(map1to1,1,false); 
      Epetra_FEVector RHS(map1to1,1,false); 
//    RHS->PutScalar(0.0);

//    Epetra_Vector LHS(Map);      
//    Epetra_MultiVector LHS(Map,1,true);
//    Epetra_MultiVector LHS(map1to1,1,true);
      Epetra_FEVector LHS(map1to1,1,true);
      Epetra_FEVector xv(Map);

//    Epetra_Export exporter(map1to1, Map);
//    xv.Export(LHS, exporter, Add);
//    xv.Import(LHS, exporter, Add);
//    std::cout << exporter ;

//    Epetra_Import Importer(Map,map1to1);
//    xv.Import(LHS,Importer,Insert);

//    xv.Print(std::cout);

      for( int i=0 ; i< NumMyElements ; i++ ) {
      glbrow = glon[i];
//      RHS.ReplaceGlobalValue(glbrow, 0, bo[i]);
//      LHS.ReplaceGlobalValue(glbrow, 0, x_vec[i]);
//       RHS.SumIntoGlobalValue(glbrow, 0, bo[i]);
        RHS.ReplaceGlobalValue(glbrow,  0, bo[i]);
//	LHS.SumIntoGlobalValue(glbrow, 0, x_vec[i]);
      }
      RHS.GlobalAssemble();

//     RHS.Print(std::cout);
//     LHS.Print(std::cout);
//     A.Print(std::cout); 
//     return 0;

// ----------  ML Preconditioner ------------------
  ParameterList MLList;

/*
 ML_Epetra::SetDefaults("SA",MLList);
// overwrite with user’s defined parameters
 MLList.set("max levels",6);
 MLList.set("increasing or decreasing","decreasing");
// MLList.set("aggregation: type", "Uncoupled");
 MLList.set("aggregation: type", "MIS");
 MLList.set("coarse: type","Amesos-KLU");
*/


/*
  ML_Epetra::SetDefaults("DD",MLList);
  MLList.set("smoother: pre or post", "both");
  MLList.set("PDE equations", 1);

  // fix the smoother to be IFPACK; can be set using (level X) syntax
  MLList.set("smoother: type","IFPACK");

  // now we have to specify which IFPACK preconditioner should be
  // built. Any value that is valid for the IFPACK factory. We also need
  // to define the overlap (>= 0).

  MLList.set("smoother: ifpack type", "ILU");
//  MLList.set("smoother: ifpack overlap", 20);

  // Then, all parameters can will control the definition of the IFPACK
  // smoother are inserted in IFPACKList. In this case, we specify the fill-in
  // factor. For a list of supported parameters, please consult the IFPACK
  // documentation. For example, IFPACK preconditioner "Amesos" or
  // "Amesos stand-alone" can be used to solve with an LU
  // factorization on each domain.
  MLList.sublist("smoother: ifpack list").set("fact: level-of-fill", 5);

  // we can now build the preconditioner...

  ML_Epetra::MultiLevelPreconditioner* MLPrec =
    new ML_Epetra::MultiLevelPreconditioner(A, MLList, true);

 // ... and solve the linear system
*/

  Epetra_LinearProblem problem(&A, &LHS, &RHS);
  AztecOO Solver(problem);
//  Solver.SetPrecOperator(MLPrec);
  Solver.SetAztecOption(AZ_solver, AZ_bicgstab);
//  Solver.SetAztecOption(AZ_solver, AZ_gmres);
  Solver.SetAztecOption(AZ_precond, AZ_dom_decomp);
  Solver.SetAztecOption(AZ_subdomain_solve, AZ_ilut);
//  Solver.SetAztecOption(AZ_overlap, 204);
  Solver.SetAztecOption(AZ_kspace, 100);
  Solver.SetAztecOption(AZ_output, 5);
  Solver.Iterate(20, 1e-5);

//  LHS.Print(std::cout);
  return 0; 

  Epetra_Import Importer(Map,map1to1);
  xv.Import(LHS,Importer,Insert);

  double xnew[NumMyElements];
  xv.ExtractCopy(xnew , 1);

  for(int i= 0; i < NumMyElements ; i++)
  {
   x_vec[i] = xnew[i];
//  std::cout << x_vec[i] << std::endl; 
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
    
