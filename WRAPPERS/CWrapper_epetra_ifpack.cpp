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

#include "Ifpack_ConfigDefs.h"

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

#include "Trilinos_Util.h"
#include "Teuchos_RCP.hpp"
#include "Teuchos_RefCountPtr.hpp"
//#include "Epetra_Export.h"
#include "Ifpack.h"


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
//   Epetra_FECrsMatrix A(Copy,Map,NumNz);
//   Epetra_FECrsMatrix A(Copy,Map,7);

std::cout << "Created crs_matrix" << std::endl; 

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
   A.FillComplete(Map,Map);
  

  Teuchos::ParameterList List;

  // allocates an IFPACK factory. No data is associated 
  // to this object (only method Create()).
  Ifpack Factory;

  // create the preconditioner. For valid PrecType values,
  // please check the documentation
  std::string PrecType = "Amesos";
  int OverlapLevel = 2; // must be >= 0. If Comm.NumProc() == 1,
                        // it is ignored.

  Teuchos::RefCountPtr<Ifpack_Preconditioner> Prec = Teuchos::rcp( Factory.Create(PrecType, &A, OverlapLevel) );
//   Ifpack_Preconditioner Prec ( Factory.Create(PrecType, &A, OverlapLevel) );
//  assert(Prec != Teuchos::null);

  // specify the Amesos solver to be used. 
  // If the selected solver is not available,
  // IFPACK will try to use Amesos' KLU (which is usually always
  // compiled). Amesos' serial solvers are:
  // "Amesos_Klu", "Amesos_Umfpack", "Amesos_Superlu"
  List.set("amesos: solver type", "Amesos_Klu");

  // sets the parameters
  Prec->SetParameters(List);

  // initialize the preconditioner. At this point the matrix must
  // have been FillComplete()'d, but actual values are ignored.
  // At this call, Amesos will perform the symbolic factorization.
  Prec->Initialize();

  // Builds the preconditioners, by looking for the values of 
  // the matrix. At this call, Amesos will perform the
  // numeric factorization.
  Prec->Compute();

    Epetra_MultiVector RHS(Map,1,false); 
    Epetra_MultiVector LHS(Map,1,true);

   A.Apply(LHS,RHS); 
//  RHS.PutScalar(1.0);  

    for( int i=0 ; i< NumMyElements; ++i ) {
      glbrow = glon[i];
//      RHS->ReplaceGlobalValue(glbrow, 0, 1.0);
      RHS.SumIntoGlobalValue(glbrow, 0, bo[i]);
  }


// ParameterList List;
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

//  ML_Epetra::MultiLevelPreconditioner* MLPrec =
//    new ML_Epetra::MultiLevelPreconditioner(A, MLList, true);

  AztecOO Solver(problem);
//  Solver.SetPrecOperator(MLPrec);
  Solver.SetPrecOperator(&*Prec);
  std::cout << " Build Preconditioner" << std::endl; 
  
  Solver.SetAztecOption(AZ_solver, AZ_bicgstab);
//  Solver.SetAztecOption(AZ_solver, AZ_gmres);
  //          Solver.SetAztecOption(AZ_precond, AZ_dom_decomp);
//  Solver.SetAztecOption(AZ_subdomain_solve, AZ_bilu);
//          Solver.SetAztecOption(AZ_subdomain_solve, AZ_ilu);
//  Solver.SetAztecOption(AZ_subdomain_solve, AZ_ilut);
//  Solver.SetAztecOption(AZ_subdomain_solve, AZ_rilu);
//   Solver.SetAztecOption(AZ_omega, 0.7); 
 // Solver.SetAztecOption(AZ_subdomain_solve, AZ_rilu);
  Solver.SetAztecOption(AZ_kspace, 200);
//  Solver.SetAztecOption(AZ_output, tolr);
//  Solver.Iterate(nitr, tolr);
   Solver.Iterate(itmax, tole);

  std::cout << " Solved the problem" << std::endl; 

  double xnew[NumMyElements];
  LHS.ExtractCopy(xnew , 1);

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
