InstallDir=$(pwd)
Soft=$InstallDir/Software;
(mkdir -p $Soft)
#(rm -rf $Soft/scripts) 
(cd  $Soft; git clone https://github.com/kvmkrao/scripts.git scripts; chmod +x $Soft/scripts/do-* $Soft/scripts/make*)

tarfiles=$InstallDir/Software/scripts;
cmake=$Soft/cmake
gcc=$Soft/gcc;
mpi=$Soft/mpi; 
lapack=$Soft/lapack;
netcdf=$Soft/netcdf; 
hdf=$Soft/hdf5; 
boost=$Soft/boost; 
zlib=$Soft/zlib; 
param=$Soft/parmetis; 
trilinos=$Soft/Trilinos; 
albany=$Soft/Albany; 
dakota=$Soft/Dakota; 
mfix=$Soft/mfix
mfixtri=$Soft/mfix-Trilinos

#export PATH="$gcc/install/bin:$mpi/install/bin:${PATH}"
#export LD_LIBRARY_PATH="$gcc/install/lib64:$mpi/install/lib:${LD_LIBRARY_PATH}"
##GCC5.4
#mkdir -p $gcc/build $gcc/install;
#if [ ! -d $gcc/install/lib ]; then
#  (cd $gcc; tar -xvf $tarfiles/gcc-5.4.0.tar.gz; ln -s $gcc/gcc-5.4.0 $gcc/src; cd $gcc/src; ./contrib/download_prerequisites)
#  (cd $gcc/build;  ../src/configure --prefix=$gcc/install --enable-lto --enable-libssp --enable-gold --with-arch=native --with-tune=native --enable-languages=c,c++,fortran --disable-multilib; make; make install);
#fi

# THe following commands are useful to install gcc and dependencies on Ubuntu 
#sudo apt-get install build-essential
#sudo apt-get install gcc
#sudo apt-get install gfortran
#sudo apt-get install zlib1g-dev

#MPI 
mkdir -p $mpi/build $mpi/install
if [ ! -d $mpi/install/lib ]; then
(cd $mpi; git clone git://git.mpich.org/mpich.git src; cd src; sh autogen.sh)
(cd $mpi/build;  ../src/configure --prefix=$mpi/install; make -j8; make install); 
fi

echo "export PATH=\$PATH:$mpi/install/bin" >> ~/.bashrc
echo "export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:$mpi/install/lib" >> ~/.bashrc

echo " ************ NOTICE ***************************"
echo "type the following command in the terminal" 
echo "source ~/.bashrc" 
echo "comment/remove the exit 1 command from the script and run the script again"
exit 1; 

mkdir -p $cmake
(cd $cmake; wget http://www.cmake.org/files/v3.2/cmake-3.2.2.tar.gz; tar xf cmake-3.2.2.tar.gz ;  mv cmake-3.2.2 src; cd src;  ./configure ; make ; sudo make install )

#Blas & Lapack
mkdir -p $lapack/build $lapack/install
if [ ! -d $lapack/install/lib ]; then
(cd $lapack; git clone https://github.com/kvmkrao/lapack.git src;)
(cd $lapack/build; $tarfiles/do-config_lapack $lapack/install $lapack/src; make -j6; make install );
fi

# ZLIB
mkdir -p $zlib/install 
if [ ! -d $zlib/install/lib ]; then
(cd $zlib; git clone https://github.com/kvmkrao/zlib.git src);
(cd $zlib/src; ./configure --prefix=$zlib/install ; make ; make install)
fi

#HDF5
mkdir -p $hdf/build $hdf/install
if [ ! -d $hdf/install/lib ]; then
(cd $hdf; wget https://www.hdfgroup.org/ftp/HDF5/current18/src/hdf5-1.8.19.tar.gz; tar -xvf hdf5-1.8.19.tar.gz; mv hdf5-1.8.19 src); 
(cd $hdf/build; $hdf/src/configure --with-zlib=$zlib/install --prefix=$hdf/install; make -j8; make install)
fi

# NETCDF
mkdir -p $netcdf/build $netcdf/install
if [ ! -d $netcdf/install/lib ]; then
(cd $netcdf; wget ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-4.4.0.tar.gz; tar -xvf netcdf-4.4.0.tar.gz; mv netcdf-4.4.0 src);    
(cd $netcdf/build; $netcdf/src/configure --prefix=$netcdf/install --disable-netcdf-4 --enable-parallel; make -j8; make install)
fi

export BOOST_BUILD_PATH=$boost/install
export BOOST_ROOT=$boost/src

### BOOST
mkdir -p  $boost/install
if [ ! -d $boost/install/lib ]; then
(cd $boost; wget https://sourceforge.net/projects/boost/files/boost/1.63.0/boost_1_63_0.tar.gz; tar -xzvf boost_1_63_0.tar.gz; mv boost_1_63_0 src);
(cd $boost/src; ./bootstrap.sh --prefix=$boost/install ac-library-zlib=$zlib/install)
(cd $boost/src; sed -i "s,using gcc ;,using gcc ; using mpi ;,g" project-config.jam; ./b2 -j8; ./b2 install)
fi
##

#mkdir -p $param/install;
#if [ ! -d $param/install/lib ]; then
#(cd $param ; tar -xvf $tarfiles/parmetis-4.0.3.tar.gz ; mv parmetis-4.0.3 $param/build);
#(cd $param/build; sed -i "s,define IDXTYPEWIDTH 32, define IDXTYPEWIDTH 64,g" metis/include/metis.h; sed -i "s,define REALTYPEWIDTH 32,define REALTYPEWIDTH 64,g" metis/include/metis.h )
#(cd $param/build/metis; make config prefix=$param/install cc=mpicc cxx=mpicxx; make -j8; make install)
#(cd $param/build; make config prefix=$param/install cc=mpicc cxx=mpicxx; make -j8; make install)
##(cd $boost/build; sed -i "s,using gcc ;,using gcc ; using mpi ;,g" project-config.jam; ./b2 install)
#fi

# TRILINOS
(mkdir -p  $trilinos/build $trilinos/install);
if [ ! -d $trilinos/install/lib ]; then 
(cd $trilinos; git clone https://github.com/trilinos/Trilinos.git src);
(cd $trilinos/build; $tarfiles/do-config_tri_taxila $InstallDir ; make -j8 ; make install)
fi

# ALBANY
#mkdir -p $albany/build
#if [ ! -d $albany/install ]; then
#(cd $albany; tar -xvf $tarfiles/Albanysep272016.tar.gz ); #ln -s Albanysep272016 src; ln -s build install);
#(cd $albany/build; $tarfiles/do-config_alb $Soft; make -j8;);
#fi
#echo 'exit.............'; exit 1;

#(mkdir -p  $mfix/install);
##if [ ! -d $mfix/install]; then
(mkdir -p $mfix )
(cd $mfix; tar -xzvf $tarfiles/2015.2.tar.gz ; mv mfix-2015.2 src);
(cd $mfix/src; ./configure -enable-dmp FFLAGS=-O3 FCFLAGS=-O3; make -j8)
##fi

#if [ ! -d $mfixtri/install/mfixnew ]; then
(mkdir -p $mfixtri) 
(cd $mfixtri; tar -xzvf $tarfiles/2015.2.tar.gz ;mv mfix-2015.2 src; cd src; git clone https://kvmkrao@bitbucket.org/kvmkrao/mfix-2016.git WRAPPERS; cp WRAPPERS/sol* WRAPPERS/leq* model/ ; cp WRAPPERS/make_exe.sh $mfixtri/src; cp WRAPPERS/MFIX* src/model);
(cd $mfixtri/src; ./configure -enable-dmp FFLAGS=-O3 FCFLAGS=-O3; make -j8; chmod +x make_exe.sh ; ./make_exe.sh $Soft)
#fi
 
# DAKOTA
#(mkdir -p $dakota/install;)
#if [ ! -d $dakota/install/lib ]; then
##(cd $dakota ; tar -xvf $tarfiles/dakota-6.4-public.src.tar.gz; cp $tarfiles/BuildDakotaCustom.cmake $dakota/dakota-6.4.0.src); 
#(cd $dakota/install; cmake -C $dakota/dakota-6.4.0.src/BuildDakotaCustom.cmake $dakota/dakota-6.4.0.src ; make -j4; make install)
#export PATH=$PATH:$dakota/install/bin:$dakota/install/test
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$dakota/install/bin:$dakota/install/lib
#fi
