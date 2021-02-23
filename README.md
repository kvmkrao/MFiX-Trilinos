# MFiX-Trilinos
Integrated advanced linear solvers in Trilinos with MFiX 

#Commands to install TPL in ubuntu 
sudo apt-get install cmake 
sudo apt-get install mpi-default-dev
sudo apt-get install libblas-dev liblapack-dev
sudo apt-get install libboost-mpi-dev
sudo apt-get install libnetcdf-mpi-dev


#sample configuration and build file 

TRILSRCDIR=$(pwd)/src
TRILINSTALLDIR=$PATH/install
cmake \
 -D CMAKE_INSTALL_PREFIX:PATH=${TRILINSTALLDIR} \
 -D CMAKE_BUILD_TYPE:STRING=RELEASE \
 -D BUILD_SHARED_LIBS:BOOL=OFF \
 -D TPL_ENABLE_MPI:BOOL=ON \
 -D Trilinos_ENABLE_Fortran:BOOL=ON \
 -D CMAKE_VERBOSE_MAKEFILE:BOOL=OFF \
 -D Trilinos_ENABLE_ALL_PACKAGES:BOOL=OFF \
 -D Trilinos_WARNINGS_AS_ERRORS_FLAGS:STRING="" \
 -D Teuchos_ENABLE_LONG_LONG_INT:BOOL=ON \
\
 -D Trilinos_ENABLE_Teuchos:BOOL=ON \
 -D Trilinos_ENABLE_Belos:BOOL=ON \
\
 -D TPL_ENABLE_Boost:BOOL=ON \
 -D TPL_ENABLE_BoostLib:BOOL=ON \
\
 -D TPL_ENABLE_Netcdf:BOOL=ON \
\
 -D Trilinos_ENABLE_Tpetra:BOOL=ON \
 -D Trilinos_ENABLE_Kokkos:BOOL=ON \
 -D Trilinos_ENABLE_Ifpack2:BOOL=ON \
 -D Trilinos_ENABLE_Amesos2:BOOL=ON \
 -D Trilinos_ENABLE_Zoltan2:BOOL=ON \
 -D Trilinos_ENABLE_MueLu:BOOL=ON \
 -D Amesos2_ENABLE_KLU2:BOOL=ON \
  -D TPL_ENABLE_BLAS:BOOL=ON \
  -D TPL_ENABLE_LAPACK:BOOL=ON \
\
 -D Trilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON \
\
 -D Trilinos_ENABLE_Kokkos:BOOL=ON \
 -D Trilinos_ENABLE_KokkosCore:BOOL=ON \
 -D Kokkos_ENABLE_Serial:BOOL=ON \
 -D Kokkos_ENABLE_OpenMP:BOOL=OFF \
 -D Kokkos_ENABLE_Pthread:BOOL=OFF \
  -D TPL_ENABLE_BLAS:BOOL=ON \
  -D TPL_ENABLE_LAPACK:BOOL=ON \
  -D CMAKE_C_FLAGS:STRING="-O3 -fPIC" \
  -D CMAKE_CXX_FLAGS:STRING=" -std=c++14 -O3 -fPIC" \
  -D CMAKE_Fortran_FLAGS:STRING=" -O3 -fPIC" \
   -D Trilinos_ENABLE_Fortran:BOOL=ON \
${TRILSRCDIR}/src

make -j8
make install

