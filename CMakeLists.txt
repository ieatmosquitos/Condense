cmake_minimum_required(VERSION 2.6)

project( Condense )

LIST(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules)


# find and include packages
include_directories( ${Condense_SOURCE_DIR}/include )

FIND_PACKAGE(CSparse)
include_directories( ${CSPARSE_INCLUDE_DIR} )

include_directories( $ENV{G2O_ROOT} )
link_directories( $ENV{G2O_ROOT}/lib )

FIND_PACKAGE(Eigen  REQUIRED)
include_directories(${EIGEN_INCLUDE_DIRS})

add_executable( condense src/condense.cpp src/edge_labeler.cpp)
#add_executable( multicondense src/multicondense.cpp src/edge_labeler.cpp)
#add_executable( visualize src/visualizer.cpp src/drawer.cpp )
#add_executable( listedges src/listedges)
add_executable( cond3nse src/cond3nse.cpp src/edge_labeler.cpp)
#add_executable( bestguess src/best_guess.cpp)
#add_executable( difference2d src/difference2d.cpp)
#add_executable( distance2d src/distance2d.cpp)
#add_executable( distance3d src/distance3d.cpp)
add_executable( compare src/compare.cpp)
add_executable( setzero src/setzero.cpp)
add_executable( analyze src/analyze.cpp )

target_link_libraries( condense ${CSPARSE_LIBRARY} ${Eigen_LIBS} g2o_core g2o_types_slam2d g2o_stuff g2o_solver_slam2d_linear g2o_solver_csparse g2o_csparse_extension)
#target_link_libraries( multicondense ${CSPARSE_LIBRARY} ${Eigen_LIBS} g2o_core g2o_types_slam2d g2o_stuff g2o_solver_slam2d_linear g2o_solver_csparse g2o_csparse_extension)
#target_link_libraries( visualize ${CSPARSE_LIBRARY} ${Eigen_LIBS} g2o_core g2o_types_slam2d g2o_stuff g2o_solver_slam2d_linear g2o_solver_csparse g2o_csparse_extension ${SFML_LIBRARIES})
#target_link_libraries( listedges ${CSPARSE_LIBRARY} ${Eigen_LIBS} g2o_core g2o_types_slam2d g2o_stuff g2o_solver_slam2d_linear g2o_solver_csparse g2o_csparse_extension g2o_types_slam3d)
target_link_libraries( cond3nse ${CSPARSE_LIBRARY} ${Eigen_LIBS} g2o_core g2o_types_slam3d g2o_stuff g2o_solver_csparse g2o_csparse_extension)
#target_link_libraries( bestguess ${CSPARSE_LIBRARY} ${Eigen_LIBS} g2o_core g2o_types_slam3d g2o_types_slam2d g2o_stuff g2o_solver_csparse g2o_csparse_extension)
#target_link_libraries( difference2d ${CSPARSE_LIBRARY} ${Eigen_LIBS} g2o_core g2o_types_slam3d g2o_types_slam2d g2o_stuff g2o_solver_csparse g2o_csparse_extension)
#target_link_libraries( distance2d ${CSPARSE_LIBRARY} ${Eigen_LIBS} g2o_core g2o_types_slam3d g2o_types_slam2d g2o_stuff g2o_solver_csparse g2o_csparse_extension)
#target_link_libraries( distance3d ${CSPARSE_LIBRARY} ${Eigen_LIBS} g2o_core g2o_types_slam3d g2o_types_slam2d g2o_stuff g2o_solver_csparse g2o_csparse_extension)
target_link_libraries( compare ${CSPARSE_LIBRARY} ${Eigen_LIBS} g2o_core g2o_types_slam3d g2o_types_slam2d g2o_stuff g2o_solver_csparse g2o_csparse_extension)
target_link_libraries( setzero ${CSPARSE_LIBRARY} ${Eigen_LIBS} g2o_core g2o_types_slam3d g2o_types_slam2d g2o_stuff g2o_solver_csparse g2o_csparse_extension)
target_link_libraries( analyze ${CSPARSE_LIBRARY} ${Eigen_LIBS} g2o_core g2o_types_slam3d g2o_types_slam2d g2o_stuff g2o_solver_csparse g2o_csparse_extension)


execute_process(COMMAND ${CMAKE_COMMAND} -E copy_directory
    ${CMAKE_SOURCE_DIR}/rd
    ${CMAKE_BINARY_DIR}/assets/rd
)

if(${CMAKE_BUILD_TYPE} MATCHES "Debug")
add_definitions( -DDEBUGMODE )
else()
endif()