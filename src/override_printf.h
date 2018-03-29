#ifndef OVERRIDE_PRINTF
#define OVERRIDE_PRINTF

#include <trio.h>

#ifdef printf
#undef printf
#endif
#ifdef vprintf
#undef vprintf
#endif
#ifdef fprintf
#undef fprintf
#endif
#ifdef vfprintf
#undef vfprintf
#endif
#ifdef sprintf
#undef sprintf
#endif
#ifdef vsprintf
#undef vsprintf
#endif
#ifdef snprintf
#undef snprintf
#endif
#ifdef vsnprintf
#undef vsnprintf
#endif
#ifdef scanf
#undef scanf
#endif
#ifdef vscanf
#undef vscanf
#endif
#ifdef fscanf
#undef fscanf
#endif
#ifdef vfscanf
#undef vfscanf
#endif
#ifdef sscanf
#undef sscanf
#endif
#ifdef vsscanf
#undef vsscanf
#endif
#ifdef vasprintf
#undef vasprintf
#endif


#define printf trio_printf
#define vprintf trio_vprintf
#define fprintf trio_fprintf
#define vfprintf trio_vfprintf
#define sprintf trio_sprintf
#define vsprintf trio_vsprintf
#define snprintf trio_snprintf
#define vsnprintf trio_vsnprintf
#define scanf trio_scanf
#define vscanf trio_vscanf
#define fscanf trio_fscanf
#define vfscanf trio_vfscanf
#define sscanf trio_sscanf
#define vsscanf trio_vsscanf
#define vasprintf trio_vasprintf

#endif
