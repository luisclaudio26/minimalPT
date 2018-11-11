cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
cmake_policy(SET CMP0015 NEW)
project(Raytracer)

#TODO:make less use of hardcoded GCC flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++11 -g -DDEBUG")

#Headers
include_directories(include)
include_directories(3rdparty)

#Source files
file(GLOB SOURCES #"3rdparty/tiny_obj_loader.cpp"
                  "src/*.cpp"
                  "src/core/camera.cpp"
                  "src/core/integrator.cpp"
                  "src/core/scene.cpp"
#                  "src/camera/*.cpp"
#                  "src/util/*.cpp"
#                  "src/shape/*.cpp"
#                  "src/texture/*.cpp"
#                  "src/integrator/*.cpp"
#                  "src/kdtree/*.cpp"
#                  "src/material/*.cpp"
                  "src/frontend/*.cpp")

#Link libraries
link_directories(3rdparty)
find_package(OpenGL REQUIRED)
set(LIBS freeimage m nanogui dl Xcursor X11 Xxf86vm Xinerama Xrandr Xi ${OPENGL_LIBRARIES} pthread)


add_executable(raytracer ${SOURCES})
target_link_libraries(raytracer ${LIBS})
