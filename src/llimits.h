/*
** $Id: llimits.h,v 1.99 2012/05/28 20:32:28 roberto Exp $
** Limits, basic types, and some other `installation-dependent' definitions
** See Copyright Notice in lua.h
*/

#ifndef llimits_h
#define llimits_h


#include <limits.h>
#include <stddef.h>


#include "lua.h"


typedef unsigned LUA_INT32 lu_int32;

typedef LUAI_UMEM lu_mem;

typedef LUAI_MEM l_mem;



/* chars used as small naturals (so that `char' is reserved for characters) */
typedef unsigned char lu_byte;


#define MAX_SIZET	((size_t)(~(size_t)0)-2)

#define MAX_LUMEM	((lu_mem)(~(lu_mem)0)-2)

#define MAX_LMEM	((l_mem) ((MAX_LUMEM >> 1) - 2))


#define MAX_INT (INT_MAX-2)  /* maximum value of an int (-2 for safety) */

/*
** conversion of pointer to integer
** this is for hashing only; there is no problem if the integer
** cannot hold the whole pointer value
*/
#define IntPoint(p)  ((unsigned int)(lu_mem)(p))



/* type to ensure maximum alignment */
#if !defined(LUAI_USER_ALIGNMENT_T)
#define LUAI_USER_ALIGNMENT_T	union { double u; void *s; long long l; }
#endif

typedef LUAI_USER_ALIGNMENT_T L_Umaxalign;


/* result of a `usual argument conversion' over lua_Number */
typedef LUAI_UACNUMBER l_uacNumber;


/* internal assertions for in-house debugging */
#if defined(lua_assert)
#define check_exp(c,e)		(lua_assert(c), (e))
/* to avoid problems with conditions too long */
#define lua_longassert(c)	{ if (!(c)) lua_assert(0); }
#else
#define lua_assert(c)		((void)0)
#define check_exp(c,e)		(e)
#define lua_longassert(c)	((void)0)
#endif

/*
** assertion for checking API calls
*/
#if !defined(luai_apicheck)

#if defined(LUA_USE_APICHECK)
  #ifndef NDEBUG
    #include <assert.h>
    #define luai_apicheck(L,e,msg)	assert((e) && msg)
  #else
    #ifndef LUA_MAKE_STRING
      #define LUA_MAKE_STRING(s) #s
    #endif

    #define luai_apicheck(L,e, msg)  (void)(!!(e) || lua_traceandabort(LUA_MAKE_STRING(e) " " msg))
  #endif
#else
  #define luai_apicheck(L,e,msg)	lua_assert(e)
#endif

#endif

#define api_check(l,e,msg)	luai_apicheck(l,e,msg)


#if !defined(UNUSED)
#define UNUSED(x)	((void)(x))	/* to avoid warnings */
#endif


#define lua_cast(t, exp)	((t)(exp))

#define cast_byte(i)	lua_cast(lu_byte, (i))
#define cast_num(i)	lua_cast(lua_Number, (i))
#define cast_int(i)	lua_cast(int, (i))
#define cast_uchar(i)	lua_cast(unsigned char, (i))


/*
** non-return type
*/
#if defined(__GNUC__)
#define l_noret		void __attribute__((noreturn))
#elif defined(_MSC_VER)
#define l_noret		void __declspec(noreturn)
#else
#define l_noret		void
#endif



/*
** maximum depth for nested C calls and syntactical nested non-terminals
** in a program. (Value must fit in an unsigned short int.)
*/
#if !defined(LUAI_MAXCCALLS)
#define LUAI_MAXCCALLS		200
#endif

/*
** maximum number of upvalues in a closure (both C and Lua). (Value
** must fit in an unsigned char.)
*/
#define MAXUPVAL	UCHAR_MAX


/*
** type for virtual-machine instructions
** must be an unsigned with (at least) 4 bytes (see details in lopcodes.h)
*/
typedef lu_int32 Instruction;



/* maximum stack for a Lua function */
#define MAXSTACK	250



/* minimum size for the string table (must be power of 2) */
#if !defined(MINSTRTABSIZE)
#define MINSTRTABSIZE	32
#endif


/* minimum size for string buffer */
#if !defined(LUA_MINBUFFER)
#define LUA_MINBUFFER	32
#endif


#if !defined(lua_lock)
#define lua_lock(L)     ((void) 0)
#define lua_unlock(L)   ((void) 0)
#endif

#if !defined(luai_threadyield)
#define luai_threadyield(L)     {lua_unlock(L); lua_lock(L);}
#endif


/*
** these macros allow user-specific actions on threads when you defined
** LUAI_EXTRASPACE and need to do something extra when a thread is
** created/deleted/resumed/yielded.
*/
#if !defined(luai_userstateopen)
#define luai_userstateopen(L)		((void)L)
#endif

#if !defined(luai_userstateclose)
#define luai_userstateclose(L)		((void)L)
#endif

#if !defined(luai_userstatethread)
#define luai_userstatethread(L,L1)	((void)L)
#endif

#if !defined(luai_userstatefree)
#define luai_userstatefree(L,L1)	((void)L)
#endif

#if !defined(luai_userstateresume)
#define luai_userstateresume(L,n)       ((void)L)
#endif

#if !defined(luai_userstateyield)
#define luai_userstateyield(L,n)        ((void)L)
#endif

/*
** lua_number2int is a macro to convert lua_Number to int.
** lua_number2integer is a macro to convert lua_Number to lua_Integer.
** lua_number2unsigned is a macro to convert a lua_Number to a lua_Unsigned.
** lua_unsigned2number is a macro to convert a lua_Unsigned to a lua_Number.
** luai_hashnum is a macro to hash a lua_Number value into an integer.
** The hash must be deterministic and give reasonable values for
** both small and large values (outside the range of integers).
*/

/* the following definitions always work, but may be slow */

#if !defined(lua_number2int)
#define lua_number2int(i,n)	((i)=(int)(n))
#endif

#if !defined(lua_number2integer)
#define lua_number2integer(i,n)	((i)=(lua_Integer)(n))
#endif

#if !defined(lua_number2unsigned)	/* { */
/* the following definition assures proper modulo behavior */
#if defined(LUA_NUMBER_DOUBLE)
#include <math.h>
#define SUPUNSIGNED	((lua_Number)(~(lua_Unsigned)0) + 1)
#define lua_number2unsigned(i,n)  \
	((i)=(lua_Unsigned)((n) - floor((n)/SUPUNSIGNED)*SUPUNSIGNED))
#else
#define lua_number2unsigned(i,n)	((i)=(lua_Unsigned)(n))
#endif
#endif				/* } */


#if !defined(lua_unsigned2number)
/* on several machines, coercion from unsigned to double is slow,
   so it may be worth to avoid */
#define lua_unsigned2number(u)  \
    (((u) <= (lua_Unsigned)INT_MAX) ? (lua_Number)(int)(u) : (lua_Number)(u))
#endif

#ifdef NINTENDO_SWITCH
//workaround for the undefined behaviors causing problems with the implementation below.
//This is taken from MS_ASMTRICK part of this original file
union luai_Cast { double l_d; LUA_INT32 l_p[2]; };
#define luai_hashnum(i,n)  \
  { volatile union luai_Cast u; u.l_d = (n) + 1.0;  /* avoid -0 */ \
    (i) = u.l_p[0]; (i) += u.l_p[1]; }  /* add double bits for his hash */
#endif

#if defined(ltable_c) && !defined(luai_hashnum)

#include <float.h>
#include <math.h>

#define luai_hashnum(i,n) { int e;  \
  n = frexp(n, &e) * (lua_Number)(INT_MAX - DBL_MAX_EXP);  \
  if (isinf(n)) i = INT_MAX; else lua_number2int(i, n); i += e; }

#endif



/*
** macro to control inclusion of some hard tests on stack reallocation
*/
#if !defined(HARDSTACKTESTS)
#define condmovestack(L)	((void)0)
#else
/* realloc stack keeping its size */
#define condmovestack(L)	luaD_reallocstack((L), (L)->stacksize)
#endif

#if !defined(HARDMEMTESTS)
#define condchangemem(L)	condmovestack(L)
#else
#define condchangemem(L)  \
	((void)(!(G(L)->gcrunning) || (G(L)->gcblocked) || (luaC_fullgc(L, 0), 1)))
#endif

#endif
