set(GLAD_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/external/glad/include)
file(GLOB GLAD_SOURCES ${CMAKE_SOURCE_DIR}/external/glad/src/glad.c)
file(GLOB GLAD_HEADERS ${CMAKE_SOURCE_DIR}/external/glad/src/glad.h)
                 
add_library(glad STATIC ${GLAD_HEADERS} ${GLAD_SOURCES})

include_directories(${GLAD_INCLUDE_DIR})

set_target_properties(glad PROPERTIES LINKER_LANGUAGE C)
set_target_properties(glad PROPERTIES FOLDER external)

set(GLAD_LIBRARIES glad)

if(NOT WIN32)
    set(GLAD_LIBRARIES ${GLAD_LIBRARIES} dl)
endif()
