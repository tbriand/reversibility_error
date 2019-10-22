# Minimal CMakeFile.txt for compiling nfft3
# Gabriele Facciolo <facciolo@cmla.ens-cachan.fr>
#
cmake_minimum_required(VERSION 2.6)
project(nfft3)

#NFFTW
#   CONFIGURE_COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/nfft-3.5.0/configure #--prefix=<INSTALL_DIR>
include(ExternalProject)
ExternalProject_Add(nfft-3.5.0
   SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
   CONFIGURE_COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/configure --prefix=<INSTALL_DIR> --disable-examples --disable-applications --disable-julia --enable-openmp --disable-nfct --disable-nfst -q
   BUILD_COMMAND ${MAKE})
ExternalProject_Get_Property(nfft-3.5.0 install_dir)
include_directories(${install_dir}/include)
link_directories(${install_dir}/lib)

# pass these variables to the parent 
set(NFFT3_INCPATH ${install_dir}/include PARENT_SCOPE)
set(NFFT3_LIBPATH ${install_dir}/lib PARENT_SCOPE)


