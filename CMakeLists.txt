cmake_minimum_required(VERSION 2.6)

# Enable C99
if(CMAKE_VERSION VERSION_LESS "3.1")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")
else()
  set(CMAKE_C_STANDARD 99)
endif()

# FindGSL.cmake not present before CMake 3.2, use provided one if necessary
if(CMAKE_VERSION VERSION_LESS "3.2")
  list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}")
endif()

set(GSL_FIND_QUIETLY TRUE)
find_package(GSL)
if(GSL_FOUND)
  include_directories("${GSL_INCLUDE_DIRS}")
  add_definitions(-DGSL_SUPPORT)
endif()

# Find OPENMP
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
endif()

# include source code directory
set(SRC src)
include_directories(${SRC})

# include external code directory
set(EXTERNAL ${SRC}/external)
include_directories(${EXTERNAL})

# include bspline code directory
set(BSPLINE ${EXTERNAL}/bspline)
include_directories(${BSPLINE})

# include NFFT
add_subdirectory(${EXTERNAL}/nfft-3.5.0)
include_directories(${NFFT3_INCPATH})
link_directories(   ${NFFT3_LIBPATH})

# set libraries
set(LIBS m jpeg png tiff)
set(LIBSFFT fftw3_threads fftw3)
if(GSL_FOUND)
set(LIBSINTERP nfft3_threads ${GSL_LIBRARIES})
else()
set(LIBSINTERP nfft3_threads)
endif()


# geometric transformation
add_executable(interpolation ${SRC}/main_interpolation.c ${SRC}/bicubic.c ${SRC}/fft_core.c ${SRC}/homography_core.c ${SRC}/tpi.c ${SRC}/periodic_plus_smooth.c ${SRC}/interpolation_core.c ${EXTERNAL}/iio.c ${BSPLINE}/splinter.c ${BSPLINE}/bspline.c)
add_dependencies(interpolation nfft-3.5.0)
target_link_libraries(interpolation ${LIBS} ${LIBSFFT} ${LIBSINTERP})

# create burst
add_executable(create_burst ${SRC}/main_create_burst.c ${SRC}/bicubic.c ${SRC}/fft_core.c ${SRC}/homography_core.c ${SRC}/tpi.c ${SRC}/periodic_plus_smooth.c ${SRC}/interpolation_core.c ${EXTERNAL}/iio.c ${BSPLINE}/splinter.c ${BSPLINE}/bspline.c)
add_dependencies(create_burst nfft-3.5.0)
target_link_libraries(create_burst ${LIBS} ${LIBSFFT} ${LIBSINTERP})

# spectrum clipping
add_executable(spectrum_clipping ${SRC}/main_spectrum_clipping.c ${SRC}/fft_core.c ${EXTERNAL}/iio.c)
target_link_libraries(spectrum_clipping ${LIBSFFT} ${LIBS})

# reversibility error
add_executable(reversibility_error ${SRC}/main_reversibility_error.c ${SRC}/fft_core.c ${EXTERNAL}/iio.c)
target_link_libraries(reversibility_error ${LIBSFFT} ${LIBS})

# crop
add_executable(crop ${SRC}/main_crop.c ${EXTERNAL}/iio.c)
target_link_libraries(crop ${LIBS})

