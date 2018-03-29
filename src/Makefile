EXCLUDE=lua.c luac.c
CFLAGS=-DLUA_COMPAT_ALL -DLUA_USE_POSIX -DLUA_USE_STRTODHEX -DLUA_USE_AFORMAT -DLUA_USE_LONGLONG
CFLAGS+=-I../trio -include override_printf.h
# Compiling C files as C++
CC:=${CXX}
CFLAGS+=-x c++
include ../common.make
