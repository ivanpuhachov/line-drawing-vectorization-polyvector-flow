cmake_minimum_required(VERSION 3.5)
set(CMAKE_VERBOSE_MAKEFILE ON)
project(polyvector_drawing_vectorization)
# Find includes in corresponding build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

# Instruct CMake to run moc automatically when needed.
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
add_executable(vectorize src/Resample.h src/FrameFieldFlowAux.h src/FrameFieldFlowAux.cpp src/distanceGraph.cpp src/distanceGraph.h)
#add_executable(${PROJECT_NAME} ${SOURCES})

set(CMAKE_FIND_LIBRARY_SUFFIXES ".so" ".a")

##### QT #####
find_package(Qt5 REQUIRED COMPONENTS Core Widgets Gui)

target_compile_definitions(vectorize PRIVATE WITH_GUI=1)
if (Qt5_FOUND)
 message(STATUS "Qt found")
 target_link_libraries(vectorize PUBLIC Qt5::Core Qt5::Widgets Qt5::Gui)
endif()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

target_include_directories(vectorize PUBLIC src)
#target_compile_options()

##### OpenMP #####
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(vectorize PUBLIC OpenMP::OpenMP_CXX)
endif()

##### OpenCV #####
find_package(OpenCV REQUIRED PATHS "/usr/lib/x86_64-linux-gnu")
target_link_libraries(vectorize PUBLIC ${OpenCV_LIBS})

##### Eigen #####
find_package (Eigen3 REQUIRED NO_MODULE)
target_link_libraries(vectorize PUBLIC Eigen3::Eigen)

##### Boost #####
find_package(Boost 1.48 REQUIRED)
message(STATUS "Boost_INCLUDE_DIRS = ${Boost_INCLUDE_DIRS}")
target_include_directories(vectorize PUBLIC ${Boost_INCLUDE_DIR})
target_link_libraries(vectorize PUBLIC ${Boost_LIBRARIES})

##### Gurobi #####
set(GUROBI_HOME "/opt/gurobi911/linux64/")
message(STATUS "GUROBI_HOME = ${GUROBI_HOME}")
find_package(GUROBI REQUIRED)
if(GUROBI_FOUND)
    message(STATUS "Gurobi found")
    target_include_directories(vectorize PUBLIC ${GUROBI_INCLUDE_DIRS})
    target_link_libraries(vectorize PUBLIC ${GUROBI_LIBRARIES})
else(GUROBI_FOUND)
    message(FATAL_ERROR "Gurobi not found, quitting")
endif(GUROBI_FOUND)

##### PAAL ##### http://paal.mimuw.edu.pl/
find_path(PAAL_INCLUDE_PATH lp/lp_base.hpp
        PATHS
            ../
            ../paal
            paal/
			D:\\libraries\\paal\\include\\)
message(STATUS "Paal = ${PAAL_INCLUDE_PATH}")
target_include_directories(vectorize PUBLIC ${PAAL_INCLUDE_PATH}/)

##### Sources #####
add_subdirectory(src)


