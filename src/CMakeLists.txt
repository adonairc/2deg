project(2DEG)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(GSL REQUIRED)
find_package(FFTW3 REQUIRED)
find_package(Boost 1.40 COMPONENTS program_options REQUIRED )




link_directories(/usr/local/lib)
include_directories(/usr/local/include)
include_directories( ${Boost_INCLUDE_DIR} )



add_executable(integration integration.cpp)
add_executable(current current.cpp)
add_executable(ldos ldos.cpp)
add_executable(magmoment magmoment.cpp)
add_executable(magmoment_dft magmoment_dft.cpp)
add_executable(two_impurities two_impurities.cpp)
add_executable(fringe_fields fringe_fields.cpp)
add_executable(fringe_fields_constant_mag fringe_fields_constant_mag.cpp)
add_executable(spinmag spinmag.cpp)
add_executable(poisson poisson.cpp)


target_link_libraries(magmoment_dft fftw3)
target_link_libraries(integration GSL::gsl GSL::gslcblas)
target_link_libraries(current GSL::gsl GSL::gslcblas)
target_link_libraries(ldos GSL::gsl GSL::gslcblas)
target_link_libraries(magmoment GSL::gsl GSL::gslcblas ${Boost_LIBRARIES})
target_link_libraries(spinmag GSL::gsl GSL::gslcblas)
target_link_libraries(fringe_fields GSL::gsl GSL::gslcblas fftw3)
target_link_libraries(two_impurities Eigen3::Eigen GSL::gsl GSL::gslcblas)
target_link_libraries(poisson fftw3)