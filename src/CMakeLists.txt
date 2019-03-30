project(Lua CXX)

file(GLOB_RECURSE SOURCE_FILES *.c)
file(GLOB_RECURSE HEADER_FILES *.h)

list(REMOVE_ITEM SOURCE_FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/lua.c
        ${CMAKE_CURRENT_SOURCE_DIR}/luac.c
        )

set_source_files_properties(${SOURCE_FILES} PROPERTIES LANGUAGE CXX)

add_library(${PROJECT_NAME} ${SOURCE_FILES} ${HEADER_FILES})

target_compile_definitions(${PROJECT_NAME} PUBLIC
        -DLUA_COMPAT_ALL
        -DLUA_USE_POSIX
        -DLUA_USE_STRTODHEX
        -DLUA_USE_AFORMAT
        -DLUA_USE_LONGLONG
)

if (ENABLE_LUA_DEBUGGER)
    target_compile_definitions(${PROJECT_NAME} PUBLIC
            -DUSE_LUA_IO
            -DUSE_LUA_COROUTINE
    )
endif()

target_compile_options(${PROJECT_NAME} PRIVATE
        -include ${CMAKE_CURRENT_SOURCE_DIR}/override_printf.h
        -x c++
)

target_link_libraries(${PROJECT_NAME} trio)
