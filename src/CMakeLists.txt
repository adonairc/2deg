project(2DEG)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(GSL REQUIRED)
#find_package(MPI REQUIRED)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})


find_package(Boost 1.40 COMPONENTS program_options mpi REQUIRED )




link_directories(/usr/local/lib)
include_directories(/usr/local/include)
include_directories( ${Boost_INCLUDE_DIR} )

#add_executable(magmoment_mpi magmoment_mpi.cpp)
#target_link_libraries(magmoment_energy GSL::gsl GSL::gslcblas Boost::mpi Boost::program_options)
#target_link_libraries(magmoment_mpi GSL::gsl GSL::gslcblas Boost::mpi Boost::program_options)
#target_link_libraries(magmoment_mpi ${Boost_LIBRARIES})
#target_link_libraries(magmoment_mpi ${MPI_LIBRARIES})

add_executable(ldos ldos.cpp)
add_executable(magmoment magmoment.cpp)
add_executable(fringe_fields_constant_mag fringe_fields_constant_mag.cpp)
add_executable(fringe_fields_envelope fringe_fields_envelope.cpp)
add_executable(fringe_fields_spin_envelope fringe_fields_spin_envelope.cpp)


target_link_libraries(ldos GSL::gsl GSL::gslcblas)
target_link_libraries(magmoment GSL::gsl GSL::gslcblas Boost::mpi Boost::program_options)
target_link_libraries(fringe_fields_constant_mag GSL::gsl GSL::gslcblas)
target_link_libraries(fringe_fields_envelope)
target_link_libraries(fringe_fields_spin_envelope)