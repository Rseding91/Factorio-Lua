/*
** $Id: lmathlib.c,v 1.81 2012/05/18 17:47:53 roberto Exp $
** Standard mathematical library
** See Copyright Notice in lua.h
*/


#include <stdlib.h>
#include <math.h>
#include <stdint.h>

#define lmathlib_c
#define LUA_LIB

#include "lua.h"

#include "lauxlib.h"
#include "lualib.h"


/* macro 'l_tg' allows the addition of an 'l' or 'f' to all math operations */
#if !defined(l_tg)
#define l_tg(x)		(x)
#endif


#undef PI
#define PI (l_tg(3.1415926535897932384626433832795))
#define RADIANS_PER_DEGREE (PI/180.0)
#define M_1_2_PI 0.159154943091895335769 // 1 / (2 * pi)

#if defined(_MSC_VER)           /* Handle Microsoft VC++ compiler specifics. */
/* C4723: potential divide by zero. */
#pragma warning ( disable : 4723 )
#endif

static double custom_fabs(double v)
{
  union
  {
    double d;
    uint64_t i;
  } t;
  t.d = v;
  t.i &= 0x7FFF'FFFF'FFFF'FFFF;
  return t.d;
}

/** Evaluate atan on interval <0, 0.66> */
static double custom_xatan(double x)
{
  double P0 = -8.750608600031904122785e-01;
  double P1 = -1.615753718733365076637e+01;
  double P2 = -7.500855792314704667340e+01;
  double P3 = -1.228866684490136173410e+02;
  double P4 = -6.485021904942025371773e+01;
  double Q0 = +2.485846490142306297962e+01;
  double Q1 = +1.650270098316988542046e+02;
  double Q2 = +4.328810604912902668951e+02;
  double Q3 = +4.853903996359136964868e+02;
  double Q4 = +1.945506571482613964425e+02;

  double z = x * x;
  z = z * ((((P0 * z + P1) * z + P2) * z + P3) * z + P4) / (((((z + Q0) * z + Q1) * z + Q2) * z + Q3) * z + Q4);
  z = x * z + x;
  return z;
}

static double custom_satan(double x)
{
  double Morebits = 6.123233995736765886130e-17;
  double Tan3pio8 = 2.41421356237309504880;

  if (x <= 0.66)
    return custom_xatan(x);

  if (x > Tan3pio8)
    return M_PI_2 - custom_xatan(1 / x) + Morebits;

  return M_PI_4 + custom_xatan((x - 1) / (x + 1)) + 0.5 * Morebits;
}

static double custom_atan(double x)
{
  // Taken from http://golang.org/src/pkg/math/atan.go

  if (x > 0)
    return custom_satan(x);
  if (x < 0)
    return -custom_satan(-x);
  return x;
}

static double custom_atan2(double y, double x)
{
  if (isnan(y) || isnan(x))
    return NAN;
  if (x > 0)
    return custom_atan(y / x);
  else if (x < 0)
    return custom_atan(y / x) + copysign(M_PI, y);
  else
    return copysign(M_PI_2, y);
}

static double custom_asin(double x)
{
  // Taken from http://golang.org/src/pkg/math/asin.go

  if (x == 0)
    return x;

  int sign = 0;

  if (x < 0)
  {
    x = -x;
    sign = 1;
  }

  double temp = sqrt(1 - x * x);
  if (x > 0.7)
    temp = M_PI_2 - custom_satan(temp / x);
  else
    temp = custom_satan(x / temp);

  if (sign)
    temp = -temp;

  return temp;
}

static double custom_acos(double x)
{
  return M_PI_2 - custom_asin(x);
}

static double custom_xcos(double x)
{
  double y = x;
  // y = y - round(y); - it's not round to nearest even integer like the SIMD version, but it shouldn't matter due to the next line
  y -= y > 0 ? (int64_t)(y + 0.5) : (int64_t)(y - 0.5);
  x = 0.25 - custom_fabs(y);
  double a = 6.2831852696304121;
  double b = 4.1341675066657376e+1;
  double c = 8.1602015295955709e+1;
  double d = 7.6568876780232551e+1;
  double e = 3.965735524898863e+1;

  double x2 = x * x;
  double x4 = x2 * x2;
  double x8 = x4 * x4;
  return x * ((a - b * x2) + x4 * (c - d * x2) + e * x8);
}

static double custom_cos(double x)
{
  return custom_xcos(x * M_1_2_PI);
}

static double custom_sin(double x)
{
  return custom_xcos(x * M_1_2_PI - 0.25);
}

static double custom_tan(double x)
{
  return custom_sin(x) / custom_cos(x);
}

typedef union
{
  double value;
  struct
  {
    uint32_t lsw;
    uint32_t msw;
  } parts;
  struct
  {
    uint64_t w;
  } xparts;
} ieee_double_shape_type;

/* Get two 32 bit ints from a double.  */

#define EXTRACT_WORDS(ix0, ix1, d)                                \
    do {                                                            \
      ieee_double_shape_type ew_u;                                  \
      ew_u.value = (d);                                             \
      (ix0) = ew_u.parts.msw;                                       \
      (ix1) = ew_u.parts.lsw;                                       \
    } while (0)

  /* Get the more significant 32 bit int from a double.  */

#define GET_HIGH_WORD(i, d)                                      \
    do {                                                            \
      ieee_double_shape_type gh_u;                                  \
      gh_u.value = (d);                                             \
      (i) = gh_u.parts.msw;                                         \
    } while (0)

  /* Set the more significant 32 bits of a double from an int.  */

#define SET_HIGH_WORD(d, v)                                      \
    do {                                                            \
      ieee_double_shape_type sh_u;                                  \
      sh_u.value = (d);                                             \
      sh_u.parts.msw = (v);                                         \
      (d) = sh_u.value;                                             \
    } while (0)

/* Get the less significant 32 bit int from a double.  */

#define GET_LOW_WORD(i,d)                                          \
    do {                                                           \
      ieee_double_shape_type gl_u;                                 \
      gl_u.value = (d);                                            \
      (i) = gl_u.parts.lsw;                                        \
    } while (0)


  /* Set the less significant 32 bits of a double from an int.  */

#define SET_LOW_WORD(d, v)                                       \
    do {                                                            \
      ieee_double_shape_type sl_u;                                  \
      sl_u.value = (d);                                             \
      sl_u.parts.lsw = (v);                                         \
      (d) = sl_u.value;                                             \
    } while (0)

static double custom_log(double x)
{
  // Taken from SDL2/src/libm/e_log.c

  const double
    ln2_hi = 6.93147180369123816490e-01,	/* 3fe62e42 fee00000 */
    ln2_lo = 1.90821492927058770002e-10,	/* 3dea39ef 35793c76 */
    two54 = 1.80143985094819840000e+16,  /* 43500000 00000000 */
    Lg1 = 6.666666666666735130e-01,  /* 3FE55555 55555593 */
    Lg2 = 3.999999999940941908e-01,  /* 3FD99999 9997FA04 */
    Lg3 = 2.857142874366239149e-01,  /* 3FD24924 94229359 */
    Lg4 = 2.222219843214978396e-01,  /* 3FCC71C5 1D8E78AF */
    Lg5 = 1.818357216161805012e-01,  /* 3FC74664 96CB03DE */
    Lg6 = 1.531383769920937332e-01,  /* 3FC39A09 D078C69F */
    Lg7 = 1.479819860511658591e-01;  /* 3FC2F112 DF3E5244 */

  const double zero = 0.0;

  double hfsq, f, s, z, R, w, t1, t2, dk;
  int32_t k, hx, i, j;
  uint32_t lx;

  EXTRACT_WORDS(hx, lx, x);

  k = 0;
  if (hx < 0x00100000) {			/* x < 2**-1022  */
    if (((hx & 0x7fffffff) | lx) == 0)
      return -INFINITY;		/* log(+-0)=-inf */
    if (hx < 0) return NAN;	/* log(-#) = NaN */
    k -= 54; x *= two54; /* subnormal number, scale up x */
    GET_HIGH_WORD(hx, x);
  }
  if (hx >= 0x7ff00000) return x + x;
  k += (hx >> 20) - 1023;
  hx &= 0x000fffff;
  i = (hx + 0x95f64) & 0x100000;
  SET_HIGH_WORD(x, hx | (i ^ 0x3ff00000));	/* normalize x or x/2 */
  k += (i >> 20);
  f = x - 1.0;
  if ((0x000fffff & (2 + hx)) < 3) {	/* |f| < 2**-20 */
    if (f == zero) {
      if (k == 0) return zero;  else {
        dk = (double)k;
        return dk * ln2_hi + dk * ln2_lo;
      }
    }
    R = f * f*(0.5 - 0.33333333333333333*f);
    if (k == 0) return f - R; else {
      dk = (double)k;
      return dk * ln2_hi - ((R - dk * ln2_lo) - f);
    }
  }
  s = f / (2.0 + f);
  dk = (double)k;
  z = s * s;
  i = hx - 0x6147a;
  w = z * z;
  j = 0x6b851 - hx;
  t1 = w * (Lg2 + w * (Lg4 + w * Lg6));
  t2 = z * (Lg1 + w * (Lg3 + w * (Lg5 + w * Lg7)));
  i |= j;
  R = t2 + t1;
  if (i > 0) {
    hfsq = 0.5*f*f;
    if (k == 0) return f - (hfsq - s * (hfsq + R)); else
      return dk * ln2_hi - ((hfsq - (s*(hfsq + R) + dk * ln2_lo)) - f);
  }
  else {
    if (k == 0) return f - s * (f - R); else
      return dk * ln2_hi - ((s*(f - R) - dk * ln2_lo) - f);
  }
}

static double custom_log10(double x)
{
  // Taken from SDL2/src/libm/e_log10.c

  const double
    two54 = 1.80143985094819840000e+16, /* 0x43500000, 0x00000000 */
    ivln10 = 4.34294481903251816668e-01, /* 0x3FDBCB7B, 0x1526E50E */
    log10_2hi = 3.01029995663611771306e-01, /* 0x3FD34413, 0x509F6000 */
    log10_2lo = 3.69423907715893078616e-13; /* 0x3D59FEF3, 0x11F12B36 */

  static const double zero = 0.0;

  double y, z;
  int32_t i, k, hx;
  uint32_t lx;

  EXTRACT_WORDS(hx, lx, x);

  k = 0;
  if (hx < 0x00100000) {                  /* x < 2**-1022  */
    if (((hx & 0x7fffffff) | lx) == 0)
      return -INFINITY;             /* log(+-0)=-inf */
    if (hx < 0) return NAN;        /* log(-#) = NaN */
    k -= 54; x *= two54; /* subnormal number, scale up x */
    GET_HIGH_WORD(hx, x);
  }
  if (hx >= 0x7ff00000) return x + x;
  k += (hx >> 20) - 1023;
  i = ((uint32_t)k & 0x80000000) >> 31;
  hx = (hx & 0x000fffff) | ((0x3ff - i) << 20);
  y = (double)(k + i);
  SET_HIGH_WORD(x, hx);
  z = y * log10_2lo + ivln10 * custom_log(x);
  return  z + y * log10_2hi;
}

// not static because we use it also in luai_numpow macro
/* non-static */ double custom_pow(double x, double y)
{
  const double
    bp[] = { 1.0, 1.5, },
    dp_h[] = { 0.0, 5.84962487220764160156e-01, }, /* 0x3FE2B803, 0x40000000 */
    dp_l[] = { 0.0, 1.35003920212974897128e-08, }, /* 0x3E4CFDEB, 0x43CFD006 */
    zero = 0.0,
    one = 1.0,
    two = 2.0,
    two53 = 9007199254740992.0,  /* 0x43400000, 0x00000000 */
    huge = 1.0e300,
    tiny = 1.0e-300,
    /* poly coefs for (3/2)*(log(x)-2s-2/3*s**3 */
    L1 = 5.99999999999994648725e-01, /* 0x3FE33333, 0x33333303 */
    L2 = 4.28571428578550184252e-01, /* 0x3FDB6DB6, 0xDB6FABFF */
    L3 = 3.33333329818377432918e-01, /* 0x3FD55555, 0x518F264D */
    L4 = 2.72728123808534006489e-01, /* 0x3FD17460, 0xA91D4101 */
    L5 = 2.30660745775561754067e-01, /* 0x3FCD864A, 0x93C9DB65 */
    L6 = 2.06975017800338417784e-01, /* 0x3FCA7E28, 0x4A454EEF */
    P1 = 1.66666666666666019037e-01, /* 0x3FC55555, 0x5555553E */
    P2 = -2.77777777770155933842e-03, /* 0xBF66C16C, 0x16BEBD93 */
    P3 = 6.61375632143793436117e-05, /* 0x3F11566A, 0xAF25DE2C */
    P4 = -1.65339022054652515390e-06, /* 0xBEBBBD41, 0xC5D26BF1 */
    P5 = 4.13813679705723846039e-08, /* 0x3E663769, 0x72BEA4D0 */
    lg2 = 6.93147180559945286227e-01, /* 0x3FE62E42, 0xFEFA39EF */
    lg2_h = 6.93147182464599609375e-01, /* 0x3FE62E43, 0x00000000 */
    lg2_l = -1.90465429995776804525e-09, /* 0xBE205C61, 0x0CA86C39 */
    ovt = 8.0085662595372944372e-0017, /* -(1024-log2(ovfl+.5ulp)) */
    cp = 9.61796693925975554329e-01, /* 0x3FEEC709, 0xDC3A03FD =2/(3ln2) */
    cp_h = 9.61796700954437255859e-01, /* 0x3FEEC709, 0xE0000000 =(float)cp */
    cp_l = -7.02846165095275826516e-09, /* 0xBE3E2FE0, 0x145B01F5 =tail of cp_h*/
    ivln2 = 1.44269504088896338700e+00, /* 0x3FF71547, 0x652B82FE =1/ln2 */
    ivln2_h = 1.44269502162933349609e+00, /* 0x3FF71547, 0x60000000 =24b 1/ln2*/
    ivln2_l = 1.92596299112661746887e-08; /* 0x3E54AE0B, 0xF85DDF44 =1/ln2 tail*/

  double z, ax, z_h, z_l, p_h, p_l;
  double y1, t1, t2, r, s, t, u, v, w;
  int i, j, k, yisint, n;
  int hx, hy, ix, iy;
  unsigned lx, ly;

  EXTRACT_WORDS(hx, lx, x);
  EXTRACT_WORDS(hy, ly, y);
  ix = hx & 0x7fffffff; iy = hy & 0x7fffffff;

  /* y==zero: x**0 = 1 */
  if ((iy | ly) == 0)
    return one;

  /* +-NaN return x+y */
  if (ix > 0x7ff00000 || ((ix == 0x7ff00000) && (lx != 0)) ||
      iy > 0x7ff00000 || ((iy == 0x7ff00000) && (ly != 0)))
    return x + y;

  /* determine if y is an odd int when x < 0
   * yisint = 0 ... y is not an integer
   * yisint = 1 ... y is an odd int
   * yisint = 2 ... y is an even int
   */
  yisint = 0;
  if (hx < 0)
  {
    if (iy >= 0x43400000)
      yisint = 2; /* even integer y */
    else if (iy >= 0x3ff00000)
    {
      k = (iy >> 20) - 0x3ff;        /* exponent */
      if (k > 20)
      {
        j = ly >> (52 - k);
        if ((uint32_t)((j << (52 - k))) == ly)
          yisint = 2 - (j & 1);
      }
      else if (ly == 0)
      {
        j = iy >> (20 - k);
        if ((j << (20 - k)) == iy)
          yisint = 2 - (j & 1);
      }
    }
  }

  /* special value of y */
  if (ly == 0)
  {
    if (iy == 0x7ff00000)         /* y is +-inf */
    {
      if (((ix - 0x3ff00000) | lx) == 0)
        return y - y; /* inf**+-1 is NaN */
      else if (ix >= 0x3ff00000) /* (|x|>1)**+-inf = inf,0 */
        return (hy >= 0) ? y : zero;
      else /* (|x|<1)**-,+inf = inf,0 */
        return (hy < 0) ? -y : zero;
    }
    if (iy == 0x3ff00000)          /* y is  +-1 */
    {
      if (hy < 0)
        return one / x;
      else
        return x;
    }
    if (hy == 0x40000000)
      return x * x; /* y is  2 */
    if (hy == 0x3fe00000) /* y is  0.5 */
      if (hx >= 0) /* x >= +0 */
        return sqrt(x);
  }

  ax = fabs(x);
  /* special value of x */
  if (lx == 0)
    if (ix == 0x7ff00000 || ix == 0 || ix == 0x3ff00000)
    {
      z = ax;                 /*x is +-0,+-inf,+-1*/
      if (hy < 0)
        z = one / z; /* z = (1/|x|) */
      if (hx < 0)
      {
        if (((ix - 0x3ff00000) | yisint) == 0)
          z = (z - z) / (z - z); /* (-1)**non-int is NaN */
        else if (yisint == 1)
          z = -z; /* (x<0)**odd = -(|x|**odd) */
      }
      return z;
    }

  n = (hx >> 31) + 1;

  /* (x<0)**(non-int) is NaN */
  if ((n | yisint) == 0)
    return (x - x) / (x - x);

  s = one; /* s (sign of result -ve**odd) = -1 else = 1 */
  if ((n | (yisint - 1)) == 0)
    s = -one; /* (-ve)**(odd int) */

  /* |y| is huge */
  if (iy > 0x41e00000)   /* if |y| > 2**31 */
  {
    if (iy > 0x43f00000)   /* if |y| > 2**64, must o/uflow */
    {
      if (ix <= 0x3fefffff)
        return (hy < 0) ? huge * huge : tiny * tiny;
      if (ix >= 0x3ff00000)
        return (hy > 0) ? huge * huge : tiny * tiny;
    }
    /* over/underflow if x is not close to one */
    if (ix < 0x3fefffff)
      return (hy < 0) ? s * huge*huge : s * tiny*tiny;
    if (ix > 0x3ff00000)
      return (hy > 0) ? s * huge*huge : s * tiny*tiny;
    /* now |1-x| is tiny <= 2**-20, suffice to compute
       log(x) by x-x^2/2+x^3/3-x^4/4 */
    t = ax - one;         /* t has 20 trailing zeros */
    w = (t*t)*(0.5 - t * (0.3333333333333333333333 - t * 0.25));
    u = ivln2_h * t;      /* ivln2_h has 21 sig. bits */
    v = t * ivln2_l - w * ivln2;
    t1 = u + v;
    SET_LOW_WORD(t1, 0);
    t2 = v - (t1 - u);
  }
  else
  {
    double ss, s2, s_h, s_l, t_h, t_l;
    n = 0;
    /* take care subnormal number */
    if (ix < 0x00100000)
    {
      ax *= two53; n -= 53; GET_HIGH_WORD(ix, ax);
    }
    n += ((ix) >> 20) - 0x3ff;
    j = ix & 0x000fffff;
    /* determine interval */
    ix = j | 0x3ff00000;          /* normalize ix */
    if (j <= 0x3988E)
      k = 0; /* |x|<sqrt(3/2) */
    else if (j < 0xBB67A)
      k = 1; /* |x|<sqrt(3)   */
    else
    {
      k = 0; n += 1; ix -= 0x00100000;
    }
    SET_HIGH_WORD(ax, ix);

    /* compute ss = s_h+s_l = (x-1)/(x+1) or (x-1.5)/(x+1.5) */
    u = ax - bp[k];               /* bp[0]=1.0, bp[1]=1.5 */
    v = one / (ax + bp[k]);
    ss = u * v;
    s_h = ss;
    SET_LOW_WORD(s_h, 0);
    /* t_h=ax+bp[k] High */
    t_h = zero;
    SET_HIGH_WORD(t_h, ((ix >> 1) | 0x20000000) + 0x00080000 + (k << 18));
    t_l = ax - (t_h - bp[k]);
    s_l = v * ((u - s_h * t_h) - s_h * t_l);
    /* compute log(ax) */
    s2 = ss * ss;
    r = s2 * s2*(L1 + s2 * (L2 + s2 * (L3 + s2 * (L4 + s2 * (L5 + s2 * L6)))));
    r += s_l * (s_h + ss);
    s2 = s_h * s_h;
    t_h = 3.0 + s2 + r;
    SET_LOW_WORD(t_h, 0);
    t_l = r - ((t_h - 3.0) - s2);
    /* u+v = ss*(1+...) */
    u = s_h * t_h;
    v = s_l * t_h + t_l * ss;
    /* 2/(3log2)*(ss+...) */
    p_h = u + v;
    SET_LOW_WORD(p_h, 0);
    p_l = v - (p_h - u);
    z_h = cp_h * p_h;             /* cp_h+cp_l = 2/(3*log2) */
    z_l = cp_l * p_h + p_l * cp + dp_l[k];
    /* log2(ax) = (ss+..)*2/(3*log2) = n + dp_h + z_h + z_l */
    t = (double)n;
    t1 = (((z_h + z_l) + dp_h[k]) + t);
    SET_LOW_WORD(t1, 0);
    t2 = z_l - (((t1 - t) - dp_h[k]) - z_h);
  }

  /* split up y into y1+y2 and compute (y1+y2)*(t1+t2) */
  y1 = y;
  SET_LOW_WORD(y1, 0);
  p_l = (y - y1)*t1 + y * t2;
  p_h = y1 * t1;
  z = p_l + p_h;
  EXTRACT_WORDS(j, i, z);
  if (j >= 0x40900000)                              /* z >= 1024 */
  {
    if (((j - 0x40900000) | i) != 0) /* if z > 1024 */
      return s * huge*huge; /* overflow */
    else if (p_l + ovt > z - p_h)
      return s * huge*huge; /* overflow */
  }
  else if ((j & 0x7fffffff) >= 0x4090cc00)          /* z <= -1075 */
  {
    if (((j - 0xc090cc00) | i) != 0) /* z < -1075 */
      return s * tiny*tiny; /* underflow */
    else if (p_l <= z - p_h)
      return s * tiny*tiny; /* underflow */
  }
  /*
   * compute 2**(p_h+p_l)
   */
  i = j & 0x7fffffff;
  k = (i >> 20) - 0x3ff;
  n = 0;
  if (i > 0x3fe00000)                /* if |z| > 0.5, set n = [z+0.5] */
  {
    n = j + (0x00100000 >> (k + 1));
    k = ((n & 0x7fffffff) >> 20) - 0x3ff;     /* new k for n */
    t = zero;
    SET_HIGH_WORD(t, n&~(0x000fffff >> k));
    n = ((n & 0x000fffff) | 0x00100000) >> (20 - k);
    if (j < 0)
      n = -n;
    p_h -= t;
  }
  t = p_l + p_h;
  SET_LOW_WORD(t, 0);
  u = t * lg2_h;
  v = (p_l - (t - p_h))*lg2 + t * lg2_l;
  z = u + v;
  w = v - (z - u);
  t = z * z;
  t1 = z - t * (P1 + t * (P2 + t * (P3 + t * (P4 + t * P5))));
  r = (z*t1) / (t1 - two) - (w + z * w);
  z = one - (r - z);
  GET_HIGH_WORD(j, z);
  j += (n << 20);
  if ((j >> 20) <= 0)
    z = scalbn(z, n); /* subnormal output */
  else
    SET_HIGH_WORD(z, j);
  return s * z;
}

static double custom_exp(double x)	/* default IEEE double exp */
{
  // Taken from SDL2/src/libm/e_exp.c

  const double
    one = 1.0,
    halF[2] = { 0.5,-0.5, },
    huge = 1.0e+300,
    twom1000 = 9.33263618503218878990e-302,     /* 2**-1000=0x01700000,0*/
    o_threshold = 7.09782712893383973096e+02,  /* 0x40862E42, 0xFEFA39EF */
    u_threshold = -7.45133219101941108420e+02,  /* 0xc0874910, 0xD52D3051 */
    ln2HI[2] = { 6.93147180369123816490e-01,  /* 0x3fe62e42, 0xfee00000 */
           -6.93147180369123816490e-01, },/* 0xbfe62e42, 0xfee00000 */
    ln2LO[2] = { 1.90821492927058770002e-10,  /* 0x3dea39ef, 0x35793c76 */
           -1.90821492927058770002e-10, },/* 0xbdea39ef, 0x35793c76 */
    invln2 = 1.44269504088896338700e+00, /* 0x3ff71547, 0x652b82fe */
    P1 = 1.66666666666666019037e-01, /* 0x3FC55555, 0x5555553E */
    P2 = -2.77777777770155933842e-03, /* 0xBF66C16C, 0x16BEBD93 */
    P3 = 6.61375632143793436117e-05, /* 0x3F11566A, 0xAF25DE2C */
    P4 = -1.65339022054652515390e-06, /* 0xBEBBBD41, 0xC5D26BF1 */
    P5 = 4.13813679705723846039e-08; /* 0x3E663769, 0x72BEA4D0 */

  double y;
  double hi = 0.0;
  double lo = 0.0;
  double c;
  double t;
  int32_t k = 0;
  int32_t xsb;
  uint32_t hx;

  GET_HIGH_WORD(hx, x);
  xsb = (hx >> 31) & 1;		/* sign bit of x */
  hx &= 0x7fffffff;		/* high word of |x| */

    /* filter out non-finite argument */
  if (hx >= 0x40862E42) {			/* if |x|>=709.78... */
    if (hx >= 0x7ff00000) {
      uint32_t lx;
      GET_LOW_WORD(lx, x);
      if (((hx & 0xfffff) | lx) != 0)
        return x + x; 		/* NaN */
      else return (xsb == 0) ? x : 0.0;	/* exp(+-inf)={inf,0} */
    }
#if 1
    if (x > o_threshold) return huge * huge; /* overflow */
#else  /* !!! FIXME: check this: "huge * huge" is a compiler warning, maybe they wanted +Inf? */
    if (x > o_threshold) return INFINITY; /* overflow */
#endif

    if (x < u_threshold) return twom1000 * twom1000; /* underflow */
  }

    /* argument reduction */
  if (hx > 0x3fd62e42) {		/* if  |x| > 0.5 ln2 */
    if (hx < 0x3FF0A2B2) {	/* and |x| < 1.5 ln2 */
      hi = x - ln2HI[xsb]; lo = ln2LO[xsb]; k = 1 - xsb - xsb;
    }
    else {
      k = (int32_t)(invln2*x + halF[xsb]);
      t = k;
      hi = x - t * ln2HI[0];	/* t*ln2HI is exact here */
      lo = t * ln2LO[0];
    }
    x = hi - lo;
  }
  else if (hx < 0x3e300000) {	/* when |x|<2**-28 */
    if (huge + x > one) return one + x;/* trigger inexact */
  }
  else k = 0;

    /* x is now in primary range */
  t = x * x;
  c = x - t * (P1 + t * (P2 + t * (P3 + t * (P4 + t * P5))));
  if (k == 0) 	return one - ((x*c) / (c - 2.0) - x);
  else 		y = one - ((lo - (x*c) / (2.0 - c)) - hi);
  if (k >= -1021) {
    uint32_t hy;
    GET_HIGH_WORD(hy, y);
    SET_HIGH_WORD(y, hy + (k << 20));	/* add k to y's exponent */
    return y;
  }
  else {
    uint32_t hy;
    GET_HIGH_WORD(hy, y);
    SET_HIGH_WORD(y, hy + ((k + 1000) << 20));	/* add k to y's exponent */
    return y * twom1000;
  }
}

#undef EXTRACT_WORDS
#undef GET_HIGH_WORD
#undef SET_HIGH_WORD
#undef GET_LOW_WORD
#undef SET_LOW_WORD

static double custom_sinh(double x)
{
  double ex = custom_exp(x);
  return (ex * ex - 1) / (2 * ex);
}

static double custom_cosh(double x)
{
  double ex = custom_exp(x);
  return (ex * ex + 1) / (2 * ex);
}

static double custom_tanh(double x)
{
  double e2x = custom_exp(2 * x);
  return (e2x - 1) / (e2x + 1);
}

static int math_abs (lua_State *L) {
  lua_pushnumber(L, l_tg(custom_fabs)(luaL_checknumber(L, 1)));
  return 1;
}

static int math_sin (lua_State *L) {
  lua_pushnumber(L, l_tg(custom_sin)(luaL_checknumber(L, 1)));
  return 1;
}

static int math_sinh (lua_State *L) {
  lua_pushnumber(L, l_tg(custom_sinh)(luaL_checknumber(L, 1)));
  return 1;
}

static int math_cos (lua_State *L) {
  lua_pushnumber(L, l_tg(custom_cos)(luaL_checknumber(L, 1)));
  return 1;
}

static int math_cosh (lua_State *L) {
  lua_pushnumber(L, l_tg(custom_cosh)(luaL_checknumber(L, 1)));
  return 1;
}

static int math_tan (lua_State *L) {
  lua_pushnumber(L, l_tg(custom_tan)(luaL_checknumber(L, 1)));
  return 1;
}

static int math_tanh (lua_State *L) {
  lua_pushnumber(L, l_tg(custom_tanh)(luaL_checknumber(L, 1)));
  return 1;
}

static int math_asin (lua_State *L) {
  lua_pushnumber(L, l_tg(custom_asin)(luaL_checknumber(L, 1)));
  return 1;
}

static int math_acos (lua_State *L) {
  lua_pushnumber(L, l_tg(custom_acos)(luaL_checknumber(L, 1)));
  return 1;
}

static int math_atan (lua_State *L) {
  lua_pushnumber(L, l_tg(custom_atan)(luaL_checknumber(L, 1)));
  return 1;
}

static int math_atan2 (lua_State *L) {
  lua_pushnumber(L, l_tg(custom_atan2)(luaL_checknumber(L, 1),
                                       luaL_checknumber(L, 2)));
  return 1;
}

static int math_ceil (lua_State *L) {
  lua_pushnumber(L, l_tg(ceil)(luaL_checknumber(L, 1)));
  return 1;
}

static int math_floor (lua_State *L) {
  lua_pushnumber(L, l_tg(floor)(luaL_checknumber(L, 1)));
  return 1;
}

static int math_fmod (lua_State *L) {
  lua_pushnumber(L, l_tg(fmod)(luaL_checknumber(L, 1),
                               luaL_checknumber(L, 2)));
  return 1;
}

static int math_modf (lua_State *L) {
  lua_Number ip;
  lua_Number fp = l_tg(modf)(luaL_checknumber(L, 1), &ip);
  lua_pushnumber(L, ip);
  lua_pushnumber(L, fp);
  return 2;
}

static int math_sqrt (lua_State *L) {
  lua_pushnumber(L, l_tg(sqrt)(luaL_checknumber(L, 1)));
  return 1;
}

static int math_pow (lua_State *L) {
  lua_pushnumber(L, l_tg(custom_pow)(luaL_checknumber(L, 1),
                              luaL_checknumber(L, 2)));
  return 1;
}

static int math_log (lua_State *L) {
  lua_Number x = luaL_checknumber(L, 1);
  lua_Number res;
  if (lua_isnoneornil(L, 2))
    res = l_tg(custom_log)(x);
  else {
    lua_Number base = luaL_checknumber(L, 2);
    if (base == 10.0) res = l_tg(custom_log10)(x);
    else res = l_tg(custom_log)(x)/l_tg(custom_log)(base);
  }
  lua_pushnumber(L, res);
  return 1;
}

#if defined(LUA_COMPAT_LOG10)
static int math_log10 (lua_State *L) {
  lua_pushnumber(L, l_tg(custom_log10)(luaL_checknumber(L, 1)));
  return 1;
}
#endif

static int math_exp (lua_State *L) {
  lua_pushnumber(L, l_tg(custom_exp)(luaL_checknumber(L, 1)));
  return 1;
}

static int math_deg (lua_State *L) {
  lua_pushnumber(L, luaL_checknumber(L, 1)/RADIANS_PER_DEGREE);
  return 1;
}

static int math_rad (lua_State *L) {
  lua_pushnumber(L, luaL_checknumber(L, 1)*RADIANS_PER_DEGREE);
  return 1;
}

static int math_frexp (lua_State *L) {
  int e;
  lua_pushnumber(L, l_tg(frexp)(luaL_checknumber(L, 1), &e));
  lua_pushinteger(L, e);
  return 2;
}

static int math_ldexp (lua_State *L) {
  lua_pushnumber(L, l_tg(ldexp)(luaL_checknumber(L, 1),
                                luaL_checkint(L, 2)));
  return 1;
}



static int math_min (lua_State *L) {
  int n = lua_gettop(L);  /* number of arguments */
  lua_Number dmin = luaL_checknumber(L, 1);
  int i;
  for (i=2; i<=n; i++) {
    lua_Number d = luaL_checknumber(L, i);
    if (d < dmin)
      dmin = d;
  }
  lua_pushnumber(L, dmin);
  return 1;
}


static int math_max (lua_State *L) {
  int n = lua_gettop(L);  /* number of arguments */
  lua_Number dmax = luaL_checknumber(L, 1);
  int i;
  for (i=2; i<=n; i++) {
    lua_Number d = luaL_checknumber(L, i);
    if (d > dmax)
      dmax = d;
  }
  lua_pushnumber(L, dmax);
  return 1;
}


static int math_random (lua_State *L) {
  /* the `%' avoids the (rare) case of r==1, and is needed also because on
     some systems (SunOS!) `rand()' may return a value larger than RAND_MAX */
  lua_Number r = (lua_Number)(rand()%RAND_MAX) / (lua_Number)RAND_MAX;
  switch (lua_gettop(L)) {  /* check number of arguments */
    case 0: {  /* no arguments */
      lua_pushnumber(L, r);  /* Number between 0 and 1 */
      break;
    }
    case 1: {  /* only upper limit */
      lua_Number u = luaL_checknumber(L, 1);
      luaL_argcheck(L, 1.0 <= u, 1, "interval is empty");
      lua_pushnumber(L, l_tg(floor)(r*u) + 1.0);  /* int in [1, u] */
      break;
    }
    case 2: {  /* lower and upper limits */
      lua_Number l = luaL_checknumber(L, 1);
      lua_Number u = luaL_checknumber(L, 2);
      luaL_argcheck(L, l <= u, 2, "interval is empty");
      lua_pushnumber(L, l_tg(floor)(r*(u-l+1)) + l);  /* int in [l, u] */
      break;
    }
    default: return luaL_error(L, "wrong number of arguments");
  }
  return 1;
}


static int math_randomseed (lua_State *L) {
  srand(luaL_checkunsigned(L, 1));
  (void)rand(); /* discard first value to avoid undesirable correlations */
  return 0;
}


static const luaL_Reg mathlib[] = {
  {"abs",   math_abs},
  {"acos",  math_acos},
  {"asin",  math_asin},
  {"atan2", math_atan2},
  {"atan",  math_atan},
  {"ceil",  math_ceil},
  {"cosh",   math_cosh},
  {"cos",   math_cos},
  {"deg",   math_deg},
  {"exp",   math_exp},
  {"floor", math_floor},
  {"fmod",   math_fmod},
  {"frexp", math_frexp},
  {"ldexp", math_ldexp},
#if defined(LUA_COMPAT_LOG10)
  {"log10", math_log10},
#endif
  {"log",   math_log},
  {"max",   math_max},
  {"min",   math_min},
  {"modf",   math_modf},
  {"pow",   math_pow},
  {"rad",   math_rad},
  {"random",     math_random},
  {"randomseed", math_randomseed},
  {"sinh",   math_sinh},
  {"sin",   math_sin},
  {"sqrt",  math_sqrt},
  {"tanh",   math_tanh},
  {"tan",   math_tan},
  {NULL, NULL}
};


/*
** Open math library
*/
LUAMOD_API int luaopen_math (lua_State *L) {
  luaL_newlib(L, mathlib);
  lua_pushnumber(L, PI);
  lua_setfield(L, -2, "pi");
  lua_pushnumber(L, HUGE_VAL);
  lua_setfield(L, -2, "huge");
  return 1;
}

