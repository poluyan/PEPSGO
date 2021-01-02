cmake_minimum_required(VERSION 3.12)

set(PROJECT_NAME CXX)

project(${PROJECT_NAME})

OPTION (USE_OpenMP "Use OpenMP" ON)
IF(USE_OpenMP)
  FIND_PACKAGE(OpenMP)
  IF(OPENMP_FOUND)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  ENDIF()
ENDIF()

set( ROSETTAMAINPATH "/ssdwork/psp_builds/rosetta_src_2020.37.61417_bundle/main" )

if(IS_DIRECTORY ${ROSETTAMAINPATH})
	message("-- Found ROSETTAMAINPATH: ${ROSETTAMAINPATH}")
else()
	message( FATAL_ERROR "You must properly set ROSETTAMAINPATH varibale ${ROSETTAMAINPATH}" )
endif()

set( ROSETTALINUXVER "5.9" )
set( ROSETTACPPCOMP "clang" )
set( ROSETTACPPVER "10.0" )
set( MVEQFDIR "${CMAKE_SOURCE_DIR}" )

if(IS_DIRECTORY ${MVEQFDIR}/mveqf)
	message("-- Found MVEQFDIR: ${MVEQFDIR}/mveqf")
else()
	message( FATAL_ERROR "You must properly set MVEQFDIR varibale ${MVEQFDIR}/mveqf" )
endif()

set( RUNLIBS1 ${ROSETTAMAINPATH}/source/build/src/release/linux/${ROSETTALINUXVER}/64/x86/${ROSETTACPPCOMP}/${ROSETTACPPVER}/default )

if(IS_DIRECTORY ${RUNLIBS1})
	message("-- Found RUNLIBS1: ${RUNLIBS1}")
else()
	message( FATAL_ERROR "You must properly set RUNLIBS1 varibale ${RUNLIBS1}" )
endif()

set( RUNLIBS2 ${ROSETTAMAINPATH}/source/build/external/release/linux/${ROSETTALINUXVER}/64/x86/${ROSETTACPPCOMP}/${ROSETTACPPVER}/default )

if(IS_DIRECTORY ${RUNLIBS2})
	message("-- Found RUNLIBS2: ${RUNLIBS2}")
else()
	message( FATAL_ERROR "You must properly set RUNLIBS2 varibale ${RUNLIBS2}" )
endif()


LINK_DIRECTORIES(${RUNLIBS1})
LINK_DIRECTORIES(${RUNLIBS2})

#file(TOUCH ${CMAKE_BINARY_DIR}/bin/run.sh)
#file(WRITE ${CMAKE_BINARY_DIR}/bin/run.sh "#/bin/sh")


set(CMAKE_CXX_STANDARD 17)
set(CXX_STANDARD_REQUIRED ON)

add_compile_options(-pipe -pedantic -Wno-long-long -Wno-strict-aliasing -march=core2 -mtune=generic -O3 -ffast-math -funroll-loops -finline-functions -Wno-unused-variable -Wno-unused-parameter -DBOOST_ERROR_CODE_HEADER_ONLY -DBOOST_SYSTEM_NO_DEPRECATED -DBOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS -DPTR_STD -DNDEBUG)

add_compile_options(-L${ROSETTAMAINPATH}/source/src -L${ROSETTAMAINPATH}/external -L/usr/lib -L/usr/local/lib)

add_compile_options(-isystem ${ROSETTAMAINPATH}/source/external/boost_1_55_0/ -isystem ${ROSETTAMAINPATH}/source/external/ -isystem ${ROSETTAMAINPATH}/source/external/include/ -isystem ${ROSETTAMAINPATH}/source/external/dbio/)

add_compile_options(-Isrc -I${ROSETTAMAINPATH}/source/src -I${ROSETTAMAINPATH}/source/external/include -I${ROSETTAMAINPATH}/source/src/platform/linux/64/${ROSETTACPPCOMP}/${ROSETTACPPVER} -I${ROSETTAMAINPATH}/source/src/platform/linux/64/${ROSETTACPPCOMP} -I${ROSETTAMAINPATH}/source/src/platform/linux/64 -I${ROSETTAMAINPATH}/source/src/platform/linux -I${ROSETTAMAINPATH}/source/external/boost_1_55_0 -I${ROSETTAMAINPATH}/source/external/libxml2/include -I${ROSETTAMAINPATH}/source/external -I${ROSETTAMAINPATH}/source/external/dbio -I/usr/include -I/usr/local/include)

file(GLOB_RECURSE SOURCES "src/*")
add_library(pepsgo SHARED ${SOURCES})
set_property(TARGET pepsgo PROPERTY INCLUDE_DIRECTORIES ${CMAKE_SOURCE_DIR}/src)
set_property(TARGET pepsgo PROPERTY INCLUDE_DIRECTORIES ${MVEQFDIR}) # mveqf

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# test0
add_executable(test0 demos/test0/test0.cpp)
set_property(TARGET test0 PROPERTY INCLUDE_DIRECTORIES ${CMAKE_SOURCE_DIR}/src)
set_property(TARGET test0 PROPERTY INCLUDE_DIRECTORIES ${MVEQFDIR}) # mveqf

target_link_libraries(test0 pepsgo)

target_link_libraries(test0 -ldevel -lprotocols.8 -lprotocols.7 -lprotocols_e.6 -lprotocols_d.6 -lprotocols_c.6 -lprotocols_b.6 -lprotocols_a.6 -lprotocols_h.5 -lprotocols_g.5 -lprotocols_f.5 -lprotocols_e.5 -lprotocols_d.5 -lprotocols_c.5 -lprotocols_b.5 -lprotocols_a.5 -lprotocols.4 -lprotocols.3 -lprotocols_b.2 -lprotocols_a.2 -lprotocols.1 -lcore.6 -lcore.5 -lcore.4 -lcore.3 -lcore.2 -lcore.1 -lbasic -lnumeric -lutility -lObjexxFCL -lz -lcppdb -lsqlite3 -lcifparse -lxml2 -lzmq)

set_target_properties(test0 PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin/test0)


# test4
add_executable(test4 demos/test4/test4.cpp)
set_property(TARGET test4 PROPERTY INCLUDE_DIRECTORIES ${CMAKE_SOURCE_DIR}/src)
set_property(TARGET test4 PROPERTY INCLUDE_DIRECTORIES ${MVEQFDIR}) # mveqf

target_link_libraries(test4 pepsgo)

target_link_libraries(test4 -ldevel -lprotocols.8 -lprotocols.7 -lprotocols_e.6 -lprotocols_d.6 -lprotocols_c.6 -lprotocols_b.6 -lprotocols_a.6 -lprotocols_h.5 -lprotocols_g.5 -lprotocols_f.5 -lprotocols_e.5 -lprotocols_d.5 -lprotocols_c.5 -lprotocols_b.5 -lprotocols_a.5 -lprotocols.4 -lprotocols.3 -lprotocols_b.2 -lprotocols_a.2 -lprotocols.1 -lcore.6 -lcore.5 -lcore.4 -lcore.3 -lcore.2 -lcore.1 -lbasic -lnumeric -lutility -lObjexxFCL -lz -lcppdb -lsqlite3 -lcifparse -lxml2 -lzmq)

set_target_properties(test4 PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin/test4)

# test5
add_executable(test5 demos/test5/test5.cpp)
set_property(TARGET test5 PROPERTY INCLUDE_DIRECTORIES ${CMAKE_SOURCE_DIR}/src)
set_property(TARGET test5 PROPERTY INCLUDE_DIRECTORIES ${MVEQFDIR}) # mveqf

target_link_libraries(test5 pepsgo)

target_link_libraries(test5 -ldevel -lprotocols.8 -lprotocols.7 -lprotocols_e.6 -lprotocols_d.6 -lprotocols_c.6 -lprotocols_b.6 -lprotocols_a.6 -lprotocols_h.5 -lprotocols_g.5 -lprotocols_f.5 -lprotocols_e.5 -lprotocols_d.5 -lprotocols_c.5 -lprotocols_b.5 -lprotocols_a.5 -lprotocols.4 -lprotocols.3 -lprotocols_b.2 -lprotocols_a.2 -lprotocols.1 -lcore.6 -lcore.5 -lcore.4 -lcore.3 -lcore.2 -lcore.1 -lbasic -lnumeric -lutility -lObjexxFCL -lz -lcppdb -lsqlite3 -lcifparse -lxml2 -lzmq)

set_target_properties(test5 PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin/test5)


# test 0
file(WRITE "${CMAKE_SOURCE_DIR}/test0.cmake" "file(COPY \"${CMAKE_SOURCE_DIR}/runtemplate1.dat\" DESTINATION \"${CMAKE_BINARY_DIR}/bin/\" FILE_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ)\nfile(COPY \"${CMAKE_SOURCE_DIR}/demos/test0/\" DESTINATION \"${CMAKE_BINARY_DIR}/bin/test0/\")\nfile(REMOVE \"${CMAKE_BINARY_DIR}/bin/test0/test0.cpp\")\nfile(RENAME \"${CMAKE_BINARY_DIR}/bin/runtemplate1.dat\" \"${CMAKE_BINARY_DIR}/bin/test0/run_test0.sh\")\nfile(APPEND \"${CMAKE_BINARY_DIR}/bin/test0/run_test0.sh\" \"export LD_LIBRARY_PATH=${CMAKE_SOURCE_DIR}:${RUNLIBS1}:${RUNLIBS2}\n\")\nfile(APPEND \"${CMAKE_BINARY_DIR}/bin/test0/run_test0.sh\" \"./test0 -database ${ROSETTAMAINPATH}/database -in::file::fasta input/sequence.fasta -in::file::native $native -mute all\n\")")

add_custom_command(TARGET test0 POST_BUILD COMMAND ${CMAKE_COMMAND} -P "${CMAKE_BINARY_DIR}/test0.cmake")

# test 4
file(WRITE "${CMAKE_SOURCE_DIR}/test4.cmake" "file(COPY \"${CMAKE_SOURCE_DIR}/runtemplate1.dat\" DESTINATION \"${CMAKE_BINARY_DIR}/bin/\" FILE_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ)\nfile(COPY \"${CMAKE_SOURCE_DIR}/demos/test4/\" DESTINATION \"${CMAKE_BINARY_DIR}/bin/test4/\")\nfile(REMOVE \"${CMAKE_BINARY_DIR}/bin/test4/test4.cpp\")\nfile(RENAME \"${CMAKE_BINARY_DIR}/bin/runtemplate1.dat\" \"${CMAKE_BINARY_DIR}/bin/test4/run_test4.sh\")\nfile(APPEND \"${CMAKE_BINARY_DIR}/bin/test4/run_test4.sh\" \"export LD_LIBRARY_PATH=${CMAKE_SOURCE_DIR}:${RUNLIBS1}:${RUNLIBS2}\n\")\nfile(APPEND \"${CMAKE_BINARY_DIR}/bin/test4/run_test4.sh\" \"./test4 -database ${ROSETTAMAINPATH}/database -in::file::fasta input/sequence.fasta -in::file::native $native -mute all\n\")")

add_custom_command(TARGET test4 POST_BUILD COMMAND ${CMAKE_COMMAND} -P "${CMAKE_BINARY_DIR}/test4.cmake")

# test 5
file(WRITE "${CMAKE_SOURCE_DIR}/test5.cmake" "file(COPY \"${CMAKE_SOURCE_DIR}/runtemplate2.dat\" DESTINATION \"${CMAKE_BINARY_DIR}/bin/\" FILE_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ)\nfile(COPY \"${CMAKE_SOURCE_DIR}/demos/test5/\" DESTINATION \"${CMAKE_BINARY_DIR}/bin/test5/\")\nfile(REMOVE \"${CMAKE_BINARY_DIR}/bin/test5/test5.cpp\")\nfile(RENAME \"${CMAKE_BINARY_DIR}/bin/runtemplate2.dat\" \"${CMAKE_BINARY_DIR}/bin/test5/run_test5.sh\")\nfile(APPEND \"${CMAKE_BINARY_DIR}/bin/test5/run_test5.sh\" \"export LD_LIBRARY_PATH=${CMAKE_SOURCE_DIR}:${RUNLIBS1}:${RUNLIBS2}\n\")\nfile(APPEND \"${CMAKE_BINARY_DIR}/bin/test5/run_test5.sh\" \"./test5 -database ${ROSETTAMAINPATH}/database -in::file::fasta input/sequence.fasta -in::file::native $native -in::file::psipred_ss2 input/prediction.horiz -frags::frag_sizes $fragsize -in::file::frag_files input/fragments.Nmers -mute all\n\")")

add_custom_command(TARGET test5 POST_BUILD COMMAND ${CMAKE_COMMAND} -P "${CMAKE_BINARY_DIR}/test5.cmake")

