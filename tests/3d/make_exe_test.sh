
tripath=$HOME/Software/Trilinos/install
lapack=$HOME/Software/lapack/install
src=$(pwd)
mfix=$(pwd)
wrapper=$(pwd) #WRAPPERS
parmetis=$HOME/Software/parmetis/install

./configure  --disable-openmp   --enable-dmp  FCFLAGS=-O3 FFLAGS=-O3 
make
rm  *.o  mfixnew  *.a 
mpif90 -I. -I$src    -I$src/model/include -I$src/model -J$src/model -DMPI -O3 -ffree-form -ffree-line-length-none -fimplicit-none -cpp -c -o model/solve_lin_eq.o ./model/solve_lin_eq.f 
make
#mpif90 -I. -I$src    -I$src/model/include -I$src/model -J$src/model -DMPI -O3 -ffree-form -ffree-line-length-none -fimplicit-none -cpp -c -o MFIXinterface.o ./WRAPPERS/MFIXinterface5_bak.f
#mpif90 -I. -I$src    -I$src/model/include -I$src/model -J$src/model -DMPI -O3 -ffree-form -ffree-line-length-none -fimplicit-none -cpp -c -o MFIXinterface.o ./WRAPPERS/MFIXinterface6.f 
mpif90 -I. -I$src    -I$src/model/include -I$src/model -J$src/model -DMPI -O3 -ffree-form -ffree-line-length-none -fimplicit-none -cpp -c -o MFIXinterface.o ./WRAPPERS/MFIXinterface6_nodek.f 
mpif90  -I. -I$src    -I$src/model/include -I$src/model -J$src/model -DMPI -O3 -ffree-form -ffree-line-length-none -fimplicit-none -cpp -c -o FWrapper_ext.o ./WRAPPERS/FWrapper_ext.f
mpic++ -c "-std=c++11" -o CWrapper.o CWrapper_epetra_ifpack_flops.cpp -I$tripath/include
#mpic++ -c "-std=c++11" -o CWrapper.o CWrapper_tpetra.cpp  -I$tripath/include

ar cr libwrapper.a MFIXinterface.o CWrapper.o FWrapper_ext.o 

mpif90  -O3   -o mfixnew  $mfix/model/mfix.o $mfix/libmfix.a libwrapper.a -L$tripath/lib -I$tripath/include  -lpiro -lrol -lstokhos_muelu -lstokhos_ifpack2 -lstokhos_amesos2 -lstokhos_tpetra -lstokhos_sacado -lstokhos -lrythmos -lmuelu-adapters -lmuelu-interface -lmuelu -llocathyra -llocaepetra -llocalapack -lloca -lnoxepetra -lnoxlapack -lnox -lintrepid -lteko -lstratimikos -lstratimikosbelos -lstratimikosaztecoo -lstratimikosamesos -lstratimikosml -lstratimikosifpack -lifpack2-adapters -lifpack2 -lanasazitpetra -lModeLaplace -lanasaziepetra -lanasazi -lIonit -lIotr -lIohb -lIogn -lIopg -lIoexo_fac -lIofx -lIoex -lIoss -lexodus -lIonit -lIotr -lIohb -lIogn -lIopg -lIoexo_fac -lIofx -lIoex -lIoss -lexodus -lamesos2 -lbelostpetra -lbelosepetra -lbelos -lml -lifpack -lzoltan2 -lpamgen_extras -lpamgen -lamesos -lgaleri-xpetra -lgaleri-epetra -laztecoo -lisorropia -loptipack -lxpetra-sup -lxpetra -lthyratpetra -lthyraepetraext -lthyraepetra -lthyracore -lthyratpetra -lthyraepetraext -lthyraepetra -lthyracore -lepetraext -ltpetraext -ltpetrainout -ltpetra -lkokkostsqr -ltpetrakernels -ltpetraclassiclinalg -ltpetraclassicnodeapi -ltpetraclassic -ltpetraext -ltpetrainout -ltpetra -lkokkostsqr -ltpetrakernels -ltpetraclassiclinalg -ltpetraclassicnodeapi -ltpetraclassic -ltriutils -lglobipack -lshards -lzoltan -lepetra -lsacado -lrtop -lteuchoskokkoscomm -lteuchoskokkoscompat -lteuchosremainder -lteuchosnumerics -lteuchoscomm -lteuchosparameterlist -lteuchoscore -lteuchoskokkoscomm -lteuchoskokkoscompat -lteuchosremainder -lteuchosnumerics -lteuchoscomm -lteuchosparameterlist -lteuchoscore -lkokkosalgorithms -lkokkoscontainers -lkokkoscore -lkokkosalgorithms -lkokkoscontainers -lkokkoscore $lapack/lib/liblapack.so.3 $lapack/lib/libblas.so.3 -I$parmetis/include $parmetis/lib/libparmetis.a $parmetis/lib/libmetis.a -lstdc++

