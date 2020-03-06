tripath=$HOME/Software/Trilinos/install
lapack=$HOME/Software/lapack/install
src=$(pwd)
mfix=$(pwd) 
wrapper=$(pwd) #WRAPPERS
parmetis=$HOME/Software/parmetis/install


#make distclean 
#make clean 
RUN_DIR=. config/rxn_preproc.sh
./configure  --disable-openmp   --enable-dmp  FCFLAGS=-O3 FFLAGS=-O3
make
SED="/bin/sed" GREP="/bin/grep" /bin/bash build-aux/fortran-depcomp model/usr_rates.f > model/usr_rates.d
SED="/bin/sed" GREP="/bin/grep" /bin/bash build-aux/fortran-depcomp model/des/usr_rates_des.f > model/des/usr_rates_des.d

mpif90 -I. -I$src    -I$src/model/include -I./model -J./model    -DMPI -g -O2 -ffree-form -ffree-line-length-none -fimplicit-none -cpp -c -o model/des/usr_rates_des.o model/des/usr_rates_des.f
mpif90 -I. -I$src   -I$src/model/include -I./model -J./model     -DMPI -g -O2 -ffree-form -ffree-line-length-none -fimplicit-none -cpp -c -o model/usr_rates.o  model/usr_rates.f

#rm  *.o  mfixnew  *.a 
mpif90 -I. -I$src    -I$src/model/include -I$src/model -J$src/model -DMPI -O3 -ffree-form -ffree-line-length-none -fimplicit-none -cpp -c -o model/solve_lin_eq.o ./model/solve_lin_eq.f 
make

mpif90 -I.  -I$src    -I$src/model/include -I$src/model -J$src/model -DMPI -O3 -ffree-form -ffree-line-length-none -fimplicit-none -cpp -c -o usr_rates.o usr_rates.f
mpif90 -I.  -I$src    -I$src/model/include -I$src/model -J$src/model -DMPI -O3 -ffree-form -ffree-line-length-none -fimplicit-none -cpp -c -o model/read_database.o model/read_database.f
mpif90 -I. -I$src    -I$src/model/include -I$src/model -J$src/model -DMPI -O3 -ffree-form -ffree-line-length-none -fimplicit-none -cpp -c -o model/read_namelist.o   model/read_namelist.f

#ar r $mfix/libmfix.a usr_rates.o   

mpif90 -I. -I$src    -I$src/model/include -I$src/model -J$src/model -DMPI -O3 -ffree-form -ffree-line-length-none -fimplicit-none -cpp -c -o MFIXinterface.o ./WRAPPERS/MFIXinterface6_nodek.f
mpif90  -I. -I$src    -I$src/model/include -I$src/model -J$src/model -DMPI -O3 -ffree-form -ffree-line-length-none -fimplicit-none -cpp -c -o FWrapper_ext.o ./WRAPPERS/FWrapper_ext.f
#mpic++ -c "-std=c++11" -o CWrapper.o CWrapper_aztec_map1bak.cpp  -I$tripath/include
mpic++  -c "-std=c++11" -o CWrapper.o CWrapper_tpetra.cpp   -I$tripath/include
#mpif90 -I.  -I$src    -I$src/model/include -I$src/model -J$src/model -DMPI -O3 -ffree-form -ffree-line-length-none -fimplicit-none -cpp -c -o usr_rates.o usr_rates.f

ar cr libwrapper.a MFIXinterface.o CWrapper.o FWrapper_ext.o 

cp $src/model/mfix.f build/mfix.f
#mpif90 -I. -I$src    -I$src/model/include -I$src/model  -DMPI   -g -O2 -ffree-form -ffree-line-length-none -fimplicit-none -cpp -c -o build/mfix.o build/mfix.f
mpif90 -I. -I$src    -I$src/model/include -I$src/model  -DMPI   -g -O2 -ffree-form -ffree-line-length-none -fimplicit-none -cpp -c -o usr_rates.o usr_rates.f
mpif90 -I. -I$src    -I$src/model/include -I$src/model  -DMPI   -g -O2 -ffree-form -ffree-line-length-none -fimplicit-none -cpp -DBURCAT_THR="'$src/model/thermochemical/BURCAT.THR'" -c -o  model/read_database.o $src/model/read_database.f
mpif90 -I. -I$src    -I$src/model/include -I$src/model  -DMPI   -g -O2 -ffree-form -ffree-line-length-none -fimplicit-none -cpp -c -o model/read_namelist.o $src/model/read_namelist.f

#mpif90  -g -O2  -DMPI  -o mfix  ./.build/mfix.o usr_rates.o ./.build/read_database.o ./.build/read_namelist.o $src/build/--enable-dmp_/libmfix.a


mpif90  -O3   -o mfixnew   usr_rates.o model/read_database.o  model/read_namelist.o $mfix/model/mfix.o $mfix/libmfix.a libwrapper.a -L$tripath/lib -I$tripath/include  -lpiro -lrol -lstokhos_muelu -lstokhos_ifpack2 -lstokhos_amesos2 -lstokhos_tpetra -lstokhos_sacado -lstokhos -lrythmos -lmuelu-adapters -lmuelu-interface -lmuelu -llocathyra -llocaepetra -llocalapack -lloca -lnoxepetra -lnoxlapack -lnox -lintrepid -lteko -lstratimikos -lstratimikosbelos -lstratimikosaztecoo -lstratimikosamesos -lstratimikosml -lstratimikosifpack -lifpack2-adapters -lifpack2 -lanasazitpetra -lModeLaplace -lanasaziepetra -lanasazi -lIonit -lIotr -lIohb -lIogn -lIopg -lIoexo_fac -lIofx -lIoex -lIoss -lexodus -lIonit -lIotr -lIohb -lIogn -lIopg -lIoexo_fac -lIofx -lIoex -lIoss -lexodus -lamesos2 -lbelostpetra -lbelosepetra -lbelos -lml -lifpack -lzoltan2 -lpamgen_extras -lpamgen -lamesos -lgaleri-xpetra -lgaleri-epetra -laztecoo -lisorropia -loptipack -lxpetra-sup -lxpetra -lthyratpetra -lthyraepetraext -lthyraepetra -lthyracore -lthyratpetra -lthyraepetraext -lthyraepetra -lthyracore -lepetraext -ltpetraext -ltpetrainout -ltpetra -lkokkostsqr -ltpetrakernels -ltpetraclassiclinalg -ltpetraclassicnodeapi -ltpetraclassic -ltpetraext -ltpetrainout -ltpetra -lkokkostsqr -ltpetrakernels -ltpetraclassiclinalg -ltpetraclassicnodeapi -ltpetraclassic -ltriutils -lglobipack -lshards -lzoltan -lepetra -lsacado -lrtop -lteuchoskokkoscomm -lteuchoskokkoscompat -lteuchosremainder -lteuchosnumerics -lteuchoscomm -lteuchosparameterlist -lteuchoscore -lteuchoskokkoscomm -lteuchoskokkoscompat -lteuchosremainder -lteuchosnumerics -lteuchoscomm -lteuchosparameterlist -lteuchoscore -lkokkosalgorithms -lkokkoscontainers -lkokkoscore -lkokkosalgorithms -lkokkoscontainers -lkokkoscore $lapack/lib/liblapack.so.3 $lapack/lib/libblas.so.3 -I$parmetis/include $parmetis/lib/libparmetis.a $parmetis/lib/libmetis.a -lstdc++

#cd tut/Liquidbed/TRI; rmmfixfiles; 
