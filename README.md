# MFiX-Trilinos
Integrated advanced linear solvers in Trilinos with MFiX 

#Commands to install TPL in ubuntu 
```
sudo apt-get install cmake 
sudo apt-get install mpi-default-dev
sudo apt-get install libblas-dev liblapack-dev
sudo apt-get install libboost-mpi-dev
sudo apt-get install libnetcdf-mpi-dev
```

#sample configuration and build file 
```
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
 -D Trilinos_ENABLE_Teuchos:BOOL=ON \
 -D Trilinos_ENABLE_Belos:BOOL=ON \
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
 -D Trilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON \
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
```


#cat Dockerfile_mfix  (docker file to build and run a container)
```
# Download base image ubuntu 20.04
FROM ubuntu:20.04

# Install dependencies
RUN apt-get update \
        && apt-get install -y \
                vim \
                ssh \
                sudo \
                wget \
                cmake \
                gfortran \
                software-properties-common ;\
                rm -rf /var/lib/apt/lists/*

RUN useradd --user-group --create-home --shell /bin/bash mfix ;\
        echo "mfix ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers

RUN mkdir -p mfix && \
    wget https://mfix.netl.doe.gov/s3/ef9c88b3/3d6affc7bf2813475bf1d847f8cb8d02//source/mfix/mfix.t# Download base image ubuntu 20.04
FROM ubuntu:20.04

# Install dependencies
RUN apt-get update \
        && apt-get install -y \
                vim \
                ssh \
                sudo \
                wget \
                cmake \
                gfortran \
                software-properties-common ;\
                rm -rf /var/lib/apt/lists/*

RUN useradd --user-group --create-home --shell /bin/bash mfix ;\
        echo "mfix ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers

RUN mkdir -p mfix && \
    wget https://mfix.netl.doe.gov/ef9c88b3/3d6affc7bf2813475bf1d847f8cb8d02//source/mfix/mfix-21.4.tar.gz && \
    tar -xzvf mfix-21.4.tar.gz &&  \
    cd mfix-21.4/tutorials/tfm/fluid_bed_2d/  && \
    cmake ../../.. ; make;

# Configure apache
RUN echo 'installed mfix' && \
    echo 'check mfix-21.4/tutorials/tfm/fluid_bed_2d/' ;

RUN echo 'Feel free to make changes to the code and build mfix'

RUN chmod -R 755 /mfix-21.4
USER mfix
#RUN cp -r /mfix-21.4/tutorials/tfm/fluid_bed_2d/ /data/
#docker build -f dockerfile_mfix -t mfix:v21.4 .
#docker run -ti --rm -v $PWD:/data -w /data mfix:v21.4 ar.gz && \
    tar -xzvf mfix-21.4.tar.gz &&  \
    cd mfix-21.4/tutorials/tfm/fluid_bed_2d/  && \
    cmake ../../.. ; make;

# Configure apache
RUN echo 'installed mfix' && \
    echo 'check mfix-21.4/tutorials/tfm/fluid_bed_2d/' ;

RUN echo 'Feel free to make changes to the code and build mfix'

RUN chmod -R 755 /mfix-21.4
USER mfix
#RUN cp -r /mfix-21.4/tutorials/tfm/fluid_bed_2d/ /data/
#docker build -f dockerfile_mfix -t mfix:v21.4 .
#docker run -ti --rm -v $PWD:/data -w /data mfix:v21.4 
```



If you publish work which uses this code, please cite the following paper:
```
@article{kotteda2018performance,
  title={Performance of preconditioned iterative solvers in MFiX--Trilinos for fluidized beds},
  author={Kotteda, VM Krushnarao and Kumar, Vinod and Spotz, William},
  journal={The Journal of Supercomputing},
  volume={74},
  number={8},
  pages={4104--4126},
  year={2018},
  publisher={Springer}
}
```
