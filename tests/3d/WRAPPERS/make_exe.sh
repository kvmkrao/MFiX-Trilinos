
tripath=$HOME/Software/Trilinos/install
lapack=$HOME/Software/lapack/install
#src=$HOME/Software/mfix-Trilinos/src
src=$HOME/Software/mfix/src
mfix=$HOME/Software/mfix/src
parmetis=$HOME/Software/parmetis/install
#mfix=$HOME/Software/mfix-Trilinos/src
#wrapper=$HOME/DOE/Trilinos-MFIX/SRC/WRAPPERS
wrapper=WRAPPERS

rm  *.o *.a mfixnew 
#mpif90 -I. -I$src    -I$src/model/include -I$src/model -J$src/model -DMPI -O3 -ffree-form -ffree-line-length-none -fimplicit-none -cpp -c -o MFIXinterface.o ./MFIXinterface4_bak.f
mpif90 -I. -I$src    -I$src/model/include -I$src/model -J$src/model -DMPI -O3 -ffree-form -ffree-line-length-none -fimplicit-none -cpp -c -o MFIXinterface.o ./MFIXinterface6.f

mpif90 -I. -I$src    -I$src/model/include -I$src/model -J$src/model -DMPI -O3 -ffree-form -ffree-line-length-none -fimplicit-none -cpp -c -o FWrapper.o      ./FWrapper.f

#------------------- git -----------------
#mpif90 -I. -I$src    -I$src/model/include -I$src/model -J$src/model -DMPI -O3 -ffree-form -ffree-line-length-none -fimplicit-none -cpp -c -o MFIXinterface.o $wrapper/MFIXinterface3.f
#mpif90 -I. -I$src    -I$src/model/include -I$src/model -J$src/model -DMPI -O3 -ffree-form -ffree-line-length-none -fimplicit-none -cpp -c -o FWrapper.o $wrapper/FWrapper.F90 
# ----------------------------------------

mpif90  -O3 -ffree-form -ffree-line-length-none -fimplicit-none -cpp -c -o FWrapper_ext.o FWrapper_ext.f
mpic++ -c "-std=c++11" CWrapper.cpp -I$tripath/include
#mpic++ -c "-std=c++11" -o CWrapper.o  CWrapper.cpp -L$tripath/lib -I$tripath/include  -lpiro -lrol -lstokhos_muelu -lstokhos_ifpack2 -lstokhos_amesos2 -lstokhos_tpetra -lstokhos_sacado -lstokhos -lrythmos -lmuelu-adapters -lmuelu-interface -lmuelu -llocathyra -llocaepetra -llocalapack -lloca -lnoxepetra -lnoxlapack -lnox -lintrepid -lteko -lstratimikos -lstratimikosbelos -lstratimikosaztecoo -lstratimikosamesos -lstratimikosml -lstratimikosifpack -lifpack2-adapters -lifpack2 -lanasazitpetra -lModeLaplace -lanasaziepetra -lanasazi -lIonit -lIotr -lIohb -lIogn -lIopg -lIoexo_fac -lIofx -lIoex -lIoss -lexodus -lIonit -lIotr -lIohb -lIogn -lIopg -lIoexo_fac -lIofx -lIoex -lIoss -lexodus -lamesos2 -lbelostpetra -lbelosepetra -lbelos -lml -lifpack -lzoltan2 -lpamgen_extras -lpamgen -lamesos -lgaleri-xpetra -lgaleri-epetra -laztecoo -lisorropia -loptipack -lxpetra-sup -lxpetra -lthyratpetra -lthyraepetraext -lthyraepetra -lthyracore -lthyratpetra -lthyraepetraext -lthyraepetra -lthyracore -lepetraext -ltpetraext -ltpetrainout -ltpetra -lkokkostsqr -ltpetrakernels -ltpetraclassiclinalg -ltpetraclassicnodeapi -ltpetraclassic -ltpetraext -ltpetrainout -ltpetra -lkokkostsqr -ltpetrakernels -ltpetraclassiclinalg -ltpetraclassicnodeapi -ltpetraclassic -ltriutils -lglobipack -lshards -lzoltan -lepetra -lsacado -lrtop -lteuchoskokkoscomm -lteuchoskokkoscompat -lteuchosremainder -lteuchosnumerics -lteuchoscomm -lteuchosparameterlist -lteuchoscore -lteuchoskokkoscomm -lteuchoskokkoscompat -lteuchosremainder -lteuchosnumerics -lteuchoscomm -lteuchosparameterlist -lteuchoscore -lkokkosalgorithms -lkokkoscontainers -lkokkoscore -lkokkosalgorithms -lkokkoscontainers -lkokkoscore $lapack/lib/liblapack.so.3 $lapack/lib/libblas.so.3

ar cr libwrapper.a MFIXinterface.o FWrapper.o FWrapper_ext.o CWrapper.o
#(cp $wrapper/solve_lin_eq.f  $wrapper/init_namelist.f  $mfix/model; cd $mfix; make)

mpif90  -O3   -o mfixnew  $mfix/model/mfix.o $mfix/libmfix.a libwrapper.a -L$tripath/lib -I$tripath/include  -lpiro -lrol -lstokhos_muelu -lstokhos_ifpack2 -lstokhos_amesos2 -lstokhos_tpetra -lstokhos_sacado -lstokhos -lrythmos -lmuelu-adapters -lmuelu-interface -lmuelu -llocathyra -llocaepetra -llocalapack -lloca -lnoxepetra -lnoxlapack -lnox -lintrepid -lteko -lstratimikos -lstratimikosbelos -lstratimikosaztecoo -lstratimikosamesos -lstratimikosml -lstratimikosifpack -lifpack2-adapters -lifpack2 -lanasazitpetra -lModeLaplace -lanasaziepetra -lanasazi -lIonit -lIotr -lIohb -lIogn -lIopg -lIoexo_fac -lIofx -lIoex -lIoss -lexodus -lIonit -lIotr -lIohb -lIogn -lIopg -lIoexo_fac -lIofx -lIoex -lIoss -lexodus -lamesos2 -lbelostpetra -lbelosepetra -lbelos -lml -lifpack -lzoltan2 -lpamgen_extras -lpamgen -lamesos -lgaleri-xpetra -lgaleri-epetra -laztecoo -lisorropia -loptipack -lxpetra-sup -lxpetra -lthyratpetra -lthyraepetraext -lthyraepetra -lthyracore -lthyratpetra -lthyraepetraext -lthyraepetra -lthyracore -lepetraext -ltpetraext -ltpetrainout -ltpetra -lkokkostsqr -ltpetrakernels -ltpetraclassiclinalg -ltpetraclassicnodeapi -ltpetraclassic -ltpetraext -ltpetrainout -ltpetra -lkokkostsqr -ltpetrakernels -ltpetraclassiclinalg -ltpetraclassicnodeapi -ltpetraclassic -ltriutils -lglobipack -lshards -lzoltan -lepetra -lsacado -lrtop -lteuchoskokkoscomm -lteuchoskokkoscompat -lteuchosremainder -lteuchosnumerics -lteuchoscomm -lteuchosparameterlist -lteuchoscore -lteuchoskokkoscomm -lteuchoskokkoscompat -lteuchosremainder -lteuchosnumerics -lteuchoscomm -lteuchosparameterlist -lteuchoscore -lkokkosalgorithms -lkokkoscontainers -lkokkoscore -lkokkosalgorithms -lkokkoscontainers -lkokkoscore $lapack/lib/liblapack.so.3 $lapack/lib/libblas.so.3 -I$parmetis/include $parmetis/lib/libparmetis.a $parmetis/lib/libmetis.a -lstdc++
#

cd LiquidFluidBed/3d
rm R.* fort.* 
mpirun -np 2 ../../mfixnew
