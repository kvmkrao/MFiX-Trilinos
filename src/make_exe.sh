
tripath=$HOME/Software/Trilinos/install
lapack=$HOME/Software/lapack/install
src=/home/vkotteda/junk/check5/mfix-trilinos/src
mfix=/home/vkotteda/junk/check5/mfix-trilinos/src
wrap=$(pwd)/WRAPPERS
parmetis=$HOME/Software/parmetis/install

FILE=./mfix-19.3.0
if [ -e $FILE ]; then
   echo "File $FILE exists."
else
   echo "File $FILE does not exist."
   wget https://mfix.netl.doe.gov/s3/ef9c88b3/2302e923d6d99aba81606df057407d94//source/mfix/mfix-19.3.1.tar.gz 
   tar -xzvf mfix-19.3.1.tar.gz ;
fi

SIZE=$(du -sb mfix-19.3.0 | cut -f1)
if [[ $SIZE -lt 1000000 ]]; then
   echo "Error in downloading the source file"
fi

#rsync -av   ../../../source/*  $(pwd)  --exclude=usr_rates.f --exclude=cases
MFIX_SRC=$(pwd)/mfix-19.3.1

#$mfix/configure --enable-dmp FCFLAGS=-O3 FFLAGS=-O3
cmake $MFIX_SRC -DENABLE_MPI=1 -DMPI_Fortran_COMPILER=mpif90 -DCMAKE_Fortran_FLAGS="-O2 -g"
make -j4
cp $wrap/solve_lin_eq.f   $src/model
cp $wrap/leq_bicgs.f      $src/model
cp $wrap/init_namelist.f  $src/model
cp $wrap/MFIXinterface.f  $src/model
cp $wrap/solve_pp_g.f     $src/model

mpif90 -I. -I$src    -I$src/model/include -I./model -J./model     -DMPI -O3 -ffree-form -ffree-line-length-none -fimplicit-none -cpp -c -o MFIXinterface.o  $wrap/MFIXinterface.f 
mpif90 -I. -I    -Iinclude -I./model -J./model     -DMPI -O3 -ffree-form -ffree-line-length-none -fimplicit-none -cpp -c -o model/solve_lin_eq.o solve_lin_eq.f

make 
#mpif90 -I. -I$src    -I$src/model/include -I./model -J./model     -DMPI -O3 -ffree-form -ffree-line-length-none -fimplicit-none -cpp -c -o  MFIXinterface.o  $src/model/MFIXinterface.f
mpif90 -DMPI -O3 -ffree-form -ffree-line-length-none -fimplicit-none -cpp -c -o FWrapper_ext.o $wrap/FWrapper_ext.f
mpicxx -c "-std=c++11" -o CWrapper_tpmu.o $wrap/CWrapper_tpetra_muelu.cpp  -I$tripath/include

ar cr libwrapper.a MFIXinterface.o CWrapper_tpmu.o FWrapper_ext.o

mpif90  -O3  -o mfixnew  model/mfix.o libmfix.a libwrapper.a -L$tripath/lib -I$tripath/include  -lmuelu-adapters -lmuelu-interface -lmuelu -llocathyra -llocaepetra -llocalapack -lloca  -lifpack2-adapters -lifpack2 -lanasazitpetra -lModeLaplace -lamesos2 -lbelostpetra -lbelosepetra -lbelos -lml -lifpack -lisorropia -lamesos -laztecoo -lxpetra-sup -lxpetra  -lepetraext -ltpetraext -ltpetrainout  -ltpetra -lkokkostsqr -ltpetrakernels -ltpetraclassiclinalg -ltpetraclassicnodeapi -ltpetraclassic -ltriutils -lglobipack  -lzoltan  -lzoltan2 -lepetra  -lteuchoskokkoscomm -lteuchoskokkoscompat -lteuchosremainder -lteuchosnumerics -lteuchoscomm -lteuchosparameterlist -lteuchoscore -lteuchoskokkoscomm -lteuchoskokkoscompat -lteuchosremainder -lteuchosnumerics -lteuchoscomm -lteuchosparameterlist -lteuchoscore -lkokkosalgorithms -lkokkoscontainers -lkokkoscore -lkokkosalgorithms -lkokkoscontainers -lkokkoscore $lapack/lib/liblapack.so.3 $lapack/lib/libblas.so.3  -I$parmetis/include $parmetis/lib/libparmetis.a $parmetis/lib/libmetis.a -lstdc++
