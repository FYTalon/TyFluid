set(OpenGL_GL_PREFERENCE GLVND)
find_package(OpenGL REQUIRED)
if(NOT OPENGL_FOUND)
    message(ERROR "OpenGL not found!")
endif()

find_package(OpenMP REQUIRED)
if(NOT OPENMP_FOUND)
    message(ERROR "OpenMP not found!")
endif()

find_package(MKL REQUIRED)

include(cmake/glad.cmake)
include(cmake/glfw.cmake)
include(cmake/imgui.cmake)
