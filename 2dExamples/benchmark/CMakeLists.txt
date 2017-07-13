#foreach(file
#        2d_benchmark_solvers
#        )
#
#    add_executable(${file} ${file}.cpp)
#    target_link_libraries(${file} NuToMechanics NuToMath NuToBase ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} ${ANN_LIBRARIES})
#    target_link_libraries(${file} NuToVisualize)
#    target_link_libraries(${file} ${MUMPS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
#endforeach()