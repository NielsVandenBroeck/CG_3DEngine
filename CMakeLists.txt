cmake_minimum_required(VERSION 3.5)

project(engine)

############################################################
# Set compiler flags
############################################################
set(OWN_GXX_FLAGS "-std=c++14")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${OWN_GXX_FLAGS} -pg -g")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${OWN_GXX_FLAGS}")
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${OWN_GXX_FLAGS}")

############################################################
# List all sources
############################################################
set(engine_sources
        Libraries/easy_image.cc
        Libraries/vector3d.cc
        engine.cc
        Parsing/ini_configuration.cc
        Line2D.h
        Line2D.cpp
        L2D.h
        L2D.cpp
        Parsing/l_parser.h
        Parsing/l_parser.cc
        Figure3D.cpp
        Figure3D.h
        L3D.h
        L3D.cpp
        ZBuffer.h
        ZBuffer.cpp
        Light.cpp
        Light.h)

############################################################
# Create an executable
############################################################
set(exe_name "engine")
add_executable( ${exe_name} ${engine_sources} )
install( TARGETS ${exe_name} DESTINATION "../2D L-Systemen")
