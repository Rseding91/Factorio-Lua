// lua.hpp
// Lua header files for C++
// <<extern "C">> not supplied automatically because Lua also compiles as C++

#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"
#include <Lua/LuaCPPUtilities.hpp>

#ifdef WIN32
#pragma comment(linker, "/export:luaL_addlstring")
#pragma comment(linker, "/export:luaL_addstring")
#pragma comment(linker, "/export:luaL_addvalue")
#pragma comment(linker, "/export:luaL_argerror")
#pragma comment(linker, "/export:luaL_buffinit")
#pragma comment(linker, "/export:luaL_buffinitsize")
#pragma comment(linker, "/export:luaL_callmeta")
#pragma comment(linker, "/export:luaL_checkany")
#pragma comment(linker, "/export:luaL_checkinteger")
#pragma comment(linker, "/export:luaL_checklstring")
#pragma comment(linker, "/export:luaL_checknumber")
#pragma comment(linker, "/export:luaL_checkoption")
#pragma comment(linker, "/export:luaL_checkstack")
#pragma comment(linker, "/export:luaL_checktype")
#pragma comment(linker, "/export:luaL_checkudata")
#pragma comment(linker, "/export:luaL_checkunsigned")
#pragma comment(linker, "/export:luaL_checkversion_")
#pragma comment(linker, "/export:luaL_error")
#pragma comment(linker, "/export:luaL_execresult")
#pragma comment(linker, "/export:luaL_fileresult")
#pragma comment(linker, "/export:luaL_getmetafield")
#pragma comment(linker, "/export:luaL_getsubtable")
#pragma comment(linker, "/export:luaL_gsub")
#pragma comment(linker, "/export:luaL_len")
#pragma comment(linker, "/export:luaL_loadbufferx")
#pragma comment(linker, "/export:luaL_loadfilex")
#pragma comment(linker, "/export:luaL_loadstring")
#pragma comment(linker, "/export:luaL_newmetatable")
#pragma comment(linker, "/export:luaL_newstate")
#pragma comment(linker, "/export:luaL_openlib")
#pragma comment(linker, "/export:luaL_openlibs")
#pragma comment(linker, "/export:luaL_optinteger")
#pragma comment(linker, "/export:luaL_optlstring")
#pragma comment(linker, "/export:luaL_optnumber")
#pragma comment(linker, "/export:luaL_optunsigned")
#pragma comment(linker, "/export:luaL_prepbuffsize")
#pragma comment(linker, "/export:luaL_pushmodule")
#pragma comment(linker, "/export:luaL_pushresult")
#pragma comment(linker, "/export:luaL_pushresultsize")
#pragma comment(linker, "/export:luaL_ref")
#pragma comment(linker, "/export:luaL_requiref")
#pragma comment(linker, "/export:luaL_setfuncs")
#pragma comment(linker, "/export:luaL_setmetatable")
#pragma comment(linker, "/export:luaL_testudata")
#pragma comment(linker, "/export:luaL_tolstring")
#pragma comment(linker, "/export:luaL_traceback")
#pragma comment(linker, "/export:luaL_unref")
#pragma comment(linker, "/export:luaL_where")
#pragma comment(linker, "/export:lua_absindex")
#pragma comment(linker, "/export:lua_arith")
#pragma comment(linker, "/export:lua_atpanic")
#pragma comment(linker, "/export:lua_callk")
#pragma comment(linker, "/export:lua_checkstack")
#pragma comment(linker, "/export:lua_close")
#pragma comment(linker, "/export:lua_compare")
#pragma comment(linker, "/export:lua_concat")
#pragma comment(linker, "/export:lua_copy")
#pragma comment(linker, "/export:lua_createtable")
#pragma comment(linker, "/export:lua_dump")
#pragma comment(linker, "/export:lua_error")
#pragma comment(linker, "/export:lua_gc")
#pragma comment(linker, "/export:lua_getallocf")
#pragma comment(linker, "/export:lua_getctx")
#pragma comment(linker, "/export:lua_getfield")
#pragma comment(linker, "/export:lua_getglobal")
#pragma comment(linker, "/export:lua_gethook")
#pragma comment(linker, "/export:lua_gethookcount")
#pragma comment(linker, "/export:lua_gethookmask")
#pragma comment(linker, "/export:lua_getinfo")
#pragma comment(linker, "/export:lua_getlocal")
#pragma comment(linker, "/export:lua_getmetatable")
#pragma comment(linker, "/export:lua_getstack")
#pragma comment(linker, "/export:lua_gettable")
#pragma comment(linker, "/export:lua_gettop")
#pragma comment(linker, "/export:lua_getupvalue")
#pragma comment(linker, "/export:lua_getuservalue")
#pragma comment(linker, "/export:lua_insert")
#pragma comment(linker, "/export:lua_iscfunction")
#pragma comment(linker, "/export:lua_isnumberorstringconvertabletonumber")
#pragma comment(linker, "/export:lua_isstringornumberconvertabletostring")
#pragma comment(linker, "/export:lua_isuserdata")
#pragma comment(linker, "/export:lua_len")
#pragma comment(linker, "/export:lua_load")
#pragma comment(linker, "/export:lua_newstate")
#pragma comment(linker, "/export:lua_newthread")
#pragma comment(linker, "/export:lua_newuserdata")
#pragma comment(linker, "/export:lua_next")
#pragma comment(linker, "/export:lua_pcallk")
#pragma comment(linker, "/export:lua_pushboolean")
#pragma comment(linker, "/export:lua_pushcclosure")
#pragma comment(linker, "/export:lua_pushfstring")
#pragma comment(linker, "/export:lua_pushinteger")
#pragma comment(linker, "/export:lua_pushlightuserdata")
#pragma comment(linker, "/export:lua_pushlstring")
#pragma comment(linker, "/export:lua_pushnil")
#pragma comment(linker, "/export:lua_pushnumber")
#pragma comment(linker, "/export:lua_pushstring")
#pragma comment(linker, "/export:lua_pushthread")
#pragma comment(linker, "/export:lua_pushunsigned")
#pragma comment(linker, "/export:lua_pushvalue")
#pragma comment(linker, "/export:lua_pushvfstring")
#pragma comment(linker, "/export:lua_rawequal")
#pragma comment(linker, "/export:lua_rawget")
#pragma comment(linker, "/export:lua_rawgeti")
#pragma comment(linker, "/export:lua_rawgetp")
#pragma comment(linker, "/export:lua_rawlen")
#pragma comment(linker, "/export:lua_rawset")
#pragma comment(linker, "/export:lua_rawseti")
#pragma comment(linker, "/export:lua_rawsetp")
#pragma comment(linker, "/export:lua_remove")
#pragma comment(linker, "/export:lua_replace")
#pragma comment(linker, "/export:lua_setallocf")
#pragma comment(linker, "/export:lua_setfield")
#pragma comment(linker, "/export:lua_setglobal")
#pragma comment(linker, "/export:lua_sethook")
#pragma comment(linker, "/export:lua_setlocal")
#pragma comment(linker, "/export:lua_setmetatable")
#pragma comment(linker, "/export:lua_settable")
#pragma comment(linker, "/export:lua_settop")
#pragma comment(linker, "/export:lua_setupvalue")
#pragma comment(linker, "/export:lua_setuservalue")
#pragma comment(linker, "/export:lua_status")
#pragma comment(linker, "/export:lua_toboolean")
#pragma comment(linker, "/export:lua_tocfunction")
#pragma comment(linker, "/export:lua_tointegerx")
#pragma comment(linker, "/export:lua_tolstring")
#pragma comment(linker, "/export:lua_tonumberx")
#pragma comment(linker, "/export:lua_topointer")
#pragma comment(linker, "/export:lua_tothread")
#pragma comment(linker, "/export:lua_tounsignedx")
#pragma comment(linker, "/export:lua_touserdata")
#pragma comment(linker, "/export:lua_type")
#pragma comment(linker, "/export:lua_typename")
#pragma comment(linker, "/export:lua_upvalueid")
#pragma comment(linker, "/export:lua_upvaluejoin")
#pragma comment(linker, "/export:lua_version")
#pragma comment(linker, "/export:lua_xmove")
#endif
