#include <iostream>
#include "stdio.h"


extern "C"
{
  double* CPPWrapper(double** A, double* b, double* x_vec, double** aloc, int* globn, int rows, int cols, int tnrows, int nitr, double tolr);
  void cwrapper_(double* br, double* xr, double* alc, double* A, int *globn, int *nAx, int *nAy, int *nrows, int *niter, double *tole, char Str[7]) 
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
 
//   std::cout << *niter << " " << *tole << std::endl; 
// std::cout << alc[0] << alc[1]<< alc[2]  << alc[3] << alc[4] << alc[5] << std::endl; 
  CPPWrapper(Amat,br,xr,loc,globn,m,n,*nrows,*niter,*tole);

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
#include "ml_MultiLevelPreconditioner.h"

using namespace std;

using namespace Teuchos;
//using namespace Galeri;

#include "Trilinos_Util.h"
//#include "Trilinos_Util_ReadMatrixMarket2Epetra.h"
//#include "CrsMatrixTranspose.h"
#include "Teuchos_RCP.hpp"
//#include "Epetra_Export.h"

double* CPPWrapper(double** Ao,double* bo,double* x_vec,double** aloc,int* globn,int rows,int cols,int tnrows,int nitr,double tolr) 
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
//     Epetra_Map Map(rows, NumMyElements, glon, 0 ,Comm);

//     std::cout << Map << std::endl; 
//     return 0; 

      Epetra_Map Map(NumGlobalElements, NumMyElements, 0 ,Comm);
//     Epetra_Map map1to1(NumGlobalElements, NumMyElements, 0 ,Comm);
//     Epetra_Map Map(-1, NumMyElements, 0, Comm);
//   Epetra_Map map1to1 =Epetra_Util::Create_OneToOne_Map(Map);

//    cout << Map ; 

//    Epetra_Map map1to1 =Epetra_Util::Create_OneToOne_Map(Map);

//    std::cout <<  NumGlobalElements << " " << NumMyElements << std::endl; 

  int * NumNz = new int[NumMyElements];
//  global_ordinal_type* NumNz = new global_ordinal_type [NumMyElements];

  int glbrow, RowLess2, RowLess1, RowLess3;
  int RowPlus2, RowPlus1, RowPlus3;
  int gid ;
  int NumEntries;

  for( int i=0 ; i<NumMyElements ; ++i) {
      gid = glon[i] ; //Map.GID(i);
      RowLess3 = gid + aloc[i][0];
      RowLess2 = gid + aloc[i][1];
      RowLess1 = gid + aloc[i][2];
      RowPlus1 = gid + aloc[i][3];
      RowPlus2 = gid + aloc[i][4];
      RowPlus3 = gid + aloc[i][5];

      NumEntries = 0 ;
      if (RowLess3 >=0 &&  abs(aloc[i][0]) != 0 )              NumEntries  = NumEntries + 1 ;
      if (RowLess2 >=0 &&  abs(aloc[i][1]) != 0 )               NumEntries = NumEntries + 1 ;
      if (RowLess1 >=0  &&  abs(aloc[i][2]) != 0)               NumEntries = NumEntries + 1 ;
      if (RowPlus1 <= NumGlobalElements && abs(aloc[i][3]) !=0) NumEntries = NumEntries + 1 ;
      if (RowPlus2 < NumGlobalElements  && abs(aloc[i][4]) !=0) NumEntries = NumEntries + 1 ;
      if (RowPlus3 < NumGlobalElements  && abs(aloc[i][5]) !=0) NumEntries = NumEntries + 1 ;
      NumEntries = NumEntries + 1 ;
      NumNz[i] = NumEntries;        
//      std::cout << gid << " " << NumEntries << std::endl; 
   }

// int * NumNz = new int[NumMyElements];
//  global_ordinal_type* NumNz = new global_ordinal_type [NumMyElements];

     Epetra_CrsMatrix A(Copy,Map,7);
//     Epetra_CrsMatrix A(Copy,Map,1);
//     Epetra_CrsMatrix A(Copy,Map,5);
//       Epetra_FECrsMatrix A(Copy,Map,NumNz);
//     Epetra_FECrsMatrix A(Copy,map1to1,5);
//       Epetra_FECrsMatrix A(Copy,Map,NumNz);
//      Epetra_FECrsMatrix A(Copy,Map,7);
//     Epetra_CrsMatrix A(Copy,map1to1,5);

//  std::cout << A << std::endl;
std::cout << "Created crs_matrix" << std::endl; 

 int lclerr = 0; 
// Teuchos::RefCountPtr<Epetra_CrsMatrix> A = Teuchos::rcp( new Epetra_CrsMatrix(Copy,Map,NumNz) );
// Epetra_CrsMatrix A(Copy,Map,5);


//  double Values[5];
//  global_ordinal_type Indices[5];

  double *Values = new double[7];
//  Values[0] = -1.0; Values[1] = -1.0;
  int *Indices = new int[7];
  double two = -1.0;
  int nend = NumMyElements; 
  int nstart= 0; 

//   if(Comm.MyPID() != Comm.NumProc() -1) nend = nend - 22;
//     if(Comm.MyPID() > 0) nstart = 22; 


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
//      A.SumIntoGlobalValues (gid, NumEntries, Values, Indices);
  }
   A.FillComplete();
  
  std::cout << " filled matrix A" << std::endl; 

//   A.GlobalAssemble();

// A.Print(std::cout); 
// return 0; 

// EpetraExt::RowMatrixToMatrixMarketFile("strided.mm",*A);

//  Teuchos::RefCountPtr<Epetra_MultiVector> LHS = Teuchos::rcp( new Epetra_MultiVector(Map, 1) );
//  Teuchos::RefCountPtr<Epetra_MultiVector> RHS = Teuchos::rcp( new Epetra_MultiVector(Map, 1) );


    Epetra_MultiVector RHS(Map,1,false); 
    Epetra_MultiVector LHS(Map,1,true);

//     Epetra_FEVector RHS(Map,1,false); 
//     Epetra_FEVector LHS(Map,1,true);
//     Epetra_FEVector RHS(map1to1,1,true);
//     Epetra_FEVector LHS(map1to1,1,true);
//      Epetra_FEVector xv(Map);

//    Epetra_MultiVector xv(map1to1,true);
//    Epetra_MultiVector RHS(Map,1,false); 
//    Epetra_MultiVector LHS(Map,1,true);


//  Epetra_Vector LHS(Map);
//  Epetra_Vector RHS(Map);

//  RHS.PutScalar(1.0);  
    for( int i=0 ; i< NumMyElements; ++i ) {
      glbrow = glon[i];
//      RHS->ReplaceGlobalValue(glbrow, 0, 1.0);
      RHS.SumIntoGlobalValue(glbrow, 0, bo[i]);
//      RHS.ReplaceGlobalValue(glbrow, 0, bo[i]);
  }


std::cout << " filled vector b"<< std::endl; 

//  RHS.SumIntoGlobalValues(NumMyElements, glon, bo, 0);
//  RHS.GlobalAssemble();

//  Epetra_FEVector LHS(RHS);

// RHS.Print(std::cout);
// return 0; 

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
  MLList.set("smoother: ifpack overlap", 50);

  // Then, all parameters can will control the definition of the IFPACK
  // smoother are inserted in IFPACKList. In this case, we specify the fill-in
  // factor. For a list of supported parameters, please consult the IFPACK
  // documentation. For example, IFPACK preconditioner "Amesos" or
  // "Amesos stand-alone" can be used to solve with an LU
  // factorization on each domain.
  //  MLList.sublist("smoother: ifpack list").set("fact: level-of-fill", 5);
*/

  // we can now build the preconditioner...

  ML_Epetra::MultiLevelPreconditioner* MLPrec =
    new ML_Epetra::MultiLevelPreconditioner(A, MLList, true);

  // ... and solve the linear system

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
  Solver.SetAztecOption(AZ_output, tolr);
  Solver.Iterate(nitr, tolr);

  std::cout << " Solved the problem" << std::endl; 

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
//  std::cout << nitr << " " << tolr << std::endl;
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
    
