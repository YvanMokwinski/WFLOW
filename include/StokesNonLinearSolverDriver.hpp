#pragma once

#include "nsGLOBAL.h"
#include "Monitor.h"
#include "NonLinearSolverNewton.hpp"
#include "INonLinearSolverDriver.hpp"
#include "SolverConstants.hpp"
#include "ns_automatic_constantes.h"
#include "libns_exact.h"



static const R L01[3] = {
  ((R)0.16666666666666666)
  ,((R)0.16666666666666666)
  ,((R)0.16666666666666666)
};

static const R L11[9] = {
  ((R)8.3333333333333329E-2)
  ,((R)4.1666666666666664E-2)
  ,((R)4.1666666666666664E-2)
  ,((R)4.1666666666666664E-2)
  ,((R)8.3333333333333329E-2)
  ,((R)4.1666666666666664E-2)
  ,((R)4.1666666666666664E-2)
  ,((R)4.1666666666666664E-2)
  ,((R)8.3333333333333329E-2)
};

static const R L21[18] = {
  ((R)1.6666666666666666E-2)
  ,((R)-8.3333333333333332E-3)
  ,((R)-8.3333333333333332E-3)
  ,((R)-8.3333333333333332E-3)
  ,((R)1.6666666666666666E-2)
  ,((R)-8.3333333333333332E-3)
  ,((R)-8.3333333333333332E-3)
  ,((R)-8.3333333333333332E-3)
  ,((R)1.6666666666666666E-2)
  ,((R)6.6666666666666666E-2)
  ,((R)6.6666666666666666E-2)
  ,((R)3.3333333333333333E-2)
  ,((R)3.3333333333333333E-2)
  ,((R)6.6666666666666666E-2)
  ,((R)6.6666666666666666E-2)
  ,((R)6.6666666666666666E-2)
  ,((R)3.3333333333333333E-2)
  ,((R)6.6666666666666666E-2)
};

static const R L2b1[21] = {
  ((R)1.6666666666666666E-2)
  ,((R)-8.3333333333333332E-3)
  ,((R)-8.3333333333333332E-3)
  ,((R)-8.3333333333333332E-3)
  ,((R)1.6666666666666666E-2)
  ,((R)-8.3333333333333332E-3)
  ,((R)-8.3333333333333332E-3)
  ,((R)-8.3333333333333332E-3)
  ,((R)1.6666666666666666E-2)
  ,((R)6.6666666666666666E-2)
  ,((R)6.6666666666666666E-2)
  ,((R)3.3333333333333333E-2)
  ,((R)3.3333333333333333E-2)
  ,((R)6.6666666666666666E-2)
  ,((R)6.6666666666666666E-2)
  ,((R)6.6666666666666666E-2)
  ,((R)3.3333333333333333E-2)
  ,((R)6.6666666666666666E-2)
  ,((R)7.4999999999999997E-2)
  ,((R)7.4999999999999997E-2)
  ,((R)7.4999999999999997E-2)
};

static const R L31[30] = {
  ((R)8.3333333333333332E-3)
  ,((R)4.1666666666666666E-3)
  ,((R)4.1666666666666666E-3)
  ,((R)4.1666666666666666E-3)
  ,((R)8.3333333333333332E-3)
  ,((R)4.1666666666666666E-3)
  ,((R)4.1666666666666666E-3)
  ,((R)4.1666666666666666E-3)
  ,((R)8.3333333333333332E-3)
  ,((R)3.7499999999999999E-2)
  ,((R)0.)
  ,((R)0.)
  ,((R)0.)
  ,((R)3.7499999999999999E-2)
  ,((R)0.)
  ,((R)0.)
  ,((R)3.7499999999999999E-2)
  ,((R)0.)
  ,((R)0.)
  ,((R)0.)
  ,((R)3.7499999999999999E-2)
  ,((R)0.)
  ,((R)0.)
  ,((R)3.7499999999999999E-2)
  ,((R)3.7499999999999999E-2)
  ,((R)0.)
  ,((R)0.)
  ,((R)7.4999999999999997E-2)
  ,((R)7.4999999999999997E-2)
  ,((R)7.4999999999999997E-2)
};

static const R L02[6] = {
  ((R)0.)
  ,((R)0.)
  ,((R)0.)
  ,((R)0.16666666666666666)
  ,((R)0.16666666666666666)
  ,((R)0.16666666666666666)
};

static const R L12[18] = {
  ((R)1.6666666666666666E-2)
  ,((R)-8.3333333333333332E-3)
  ,((R)-8.3333333333333332E-3)
  ,((R)6.6666666666666666E-2)
  ,((R)3.3333333333333333E-2)
  ,((R)6.6666666666666666E-2)
  ,((R)-8.3333333333333332E-3)
  ,((R)1.6666666666666666E-2)
  ,((R)-8.3333333333333332E-3)
  ,((R)6.6666666666666666E-2)
  ,((R)6.6666666666666666E-2)
  ,((R)3.3333333333333333E-2)
  ,((R)-8.3333333333333332E-3)
  ,((R)-8.3333333333333332E-3)
  ,((R)1.6666666666666666E-2)
  ,((R)3.3333333333333333E-2)
  ,((R)6.6666666666666666E-2)
  ,((R)6.6666666666666666E-2)
};

static const R L22[36] = {
  ((R)1.6666666666666666E-2)
  ,((R)-2.7777777777777779E-3)
  ,((R)-2.7777777777777779E-3)
  ,((R)0.)
  ,((R)-1.1111111111111112E-2)
  ,((R)0.)
  ,((R)-2.7777777777777779E-3)
  ,((R)1.6666666666666666E-2)
  ,((R)-2.7777777777777779E-3)
  ,((R)0.)
  ,((R)0.)
  ,((R)-1.1111111111111112E-2)
  ,((R)-2.7777777777777779E-3)
  ,((R)-2.7777777777777779E-3)
  ,((R)1.6666666666666666E-2)
  ,((R)-1.1111111111111112E-2)
  ,((R)0.)
  ,((R)0.)
  ,((R)0.)
  ,((R)0.)
  ,((R)-1.1111111111111112E-2)
  ,((R)8.8888888888888892E-2)
  ,((R)4.4444444444444446E-2)
  ,((R)4.4444444444444446E-2)
  ,((R)-1.1111111111111112E-2)
  ,((R)0.)
  ,((R)0.)
  ,((R)4.4444444444444446E-2)
  ,((R)8.8888888888888892E-2)
  ,((R)4.4444444444444446E-2)
  ,((R)0.)
  ,((R)-1.1111111111111112E-2)
  ,((R)0.)
  ,((R)4.4444444444444446E-2)
  ,((R)4.4444444444444446E-2)
  ,((R)8.8888888888888892E-2)
};

static const R L2b2[42] = {
  ((R)1.6666666666666666E-2)
  ,((R)-2.7777777777777779E-3)
  ,((R)-2.7777777777777779E-3)
  ,((R)0.)
  ,((R)-1.1111111111111112E-2)
  ,((R)0.)
  ,((R)-2.7777777777777779E-3)
  ,((R)1.6666666666666666E-2)
  ,((R)-2.7777777777777779E-3)
  ,((R)0.)
  ,((R)0.)
  ,((R)-1.1111111111111112E-2)
  ,((R)-2.7777777777777779E-3)
  ,((R)-2.7777777777777779E-3)
  ,((R)1.6666666666666666E-2)
  ,((R)-1.1111111111111112E-2)
  ,((R)0.)
  ,((R)0.)
  ,((R)0.)
  ,((R)0.)
  ,((R)-1.1111111111111112E-2)
  ,((R)8.8888888888888892E-2)
  ,((R)4.4444444444444446E-2)
  ,((R)4.4444444444444446E-2)
  ,((R)-1.1111111111111112E-2)
  ,((R)0.)
  ,((R)0.)
  ,((R)4.4444444444444446E-2)
  ,((R)8.8888888888888892E-2)
  ,((R)4.4444444444444446E-2)
  ,((R)0.)
  ,((R)-1.1111111111111112E-2)
  ,((R)0.)
  ,((R)4.4444444444444446E-2)
  ,((R)4.4444444444444446E-2)
  ,((R)8.8888888888888892E-2)
  ,((R)-1.0714285714285714E-2)
  ,((R)-1.0714285714285714E-2)
  ,((R)-1.0714285714285714E-2)
  ,((R)8.5714285714285715E-2)
  ,((R)8.5714285714285715E-2)
  ,((R)8.5714285714285715E-2)
};

static const R L32[60] = {
  ((R)5.9523809523809521E-3)
  ,((R)5.9523809523809529E-4)
  ,((R)5.9523809523809529E-4)
  ,((R)2.3809523809523812E-3)
  ,((R)4.7619047619047623E-3)
  ,((R)2.3809523809523812E-3)
  ,((R)5.9523809523809529E-4)
  ,((R)5.9523809523809521E-3)
  ,((R)5.9523809523809529E-4)
  ,((R)2.3809523809523812E-3)
  ,((R)2.3809523809523812E-3)
  ,((R)4.7619047619047623E-3)
  ,((R)5.9523809523809529E-4)
  ,((R)5.9523809523809529E-4)
  ,((R)5.9523809523809521E-3)
  ,((R)4.7619047619047623E-3)
  ,((R)2.3809523809523812E-3)
  ,((R)2.3809523809523812E-3)
  ,((R)1.607142857142857E-2)
  ,((R)-1.0714285714285714E-2)
  ,((R)-3.5714285714285713E-3)
  ,((R)2.8571428571428571E-2)
  ,((R)-7.1428571428571426E-3)
  ,((R)1.4285714285714285E-2)
  ,((R)-1.0714285714285714E-2)
  ,((R)1.607142857142857E-2)
  ,((R)-3.5714285714285713E-3)
  ,((R)2.8571428571428571E-2)
  ,((R)1.4285714285714285E-2)
  ,((R)-7.1428571428571426E-3)
  ,((R)-3.5714285714285713E-3)
  ,((R)1.607142857142857E-2)
  ,((R)-1.0714285714285714E-2)
  ,((R)1.4285714285714285E-2)
  ,((R)2.8571428571428571E-2)
  ,((R)-7.1428571428571426E-3)
  ,((R)-3.5714285714285713E-3)
  ,((R)-1.0714285714285714E-2)
  ,((R)1.607142857142857E-2)
  ,((R)-7.1428571428571426E-3)
  ,((R)2.8571428571428571E-2)
  ,((R)1.4285714285714285E-2)
  ,((R)-1.0714285714285714E-2)
  ,((R)-3.5714285714285713E-3)
  ,((R)1.607142857142857E-2)
  ,((R)-7.1428571428571426E-3)
  ,((R)1.4285714285714285E-2)
  ,((R)2.8571428571428571E-2)
  ,((R)1.607142857142857E-2)
  ,((R)-3.5714285714285713E-3)
  ,((R)-1.0714285714285714E-2)
  ,((R)1.4285714285714285E-2)
  ,((R)-7.1428571428571426E-3)
  ,((R)2.8571428571428571E-2)
  ,((R)-1.0714285714285714E-2)
  ,((R)-1.0714285714285714E-2)
  ,((R)-1.0714285714285714E-2)
  ,((R)8.5714285714285715E-2)
  ,((R)8.5714285714285715E-2)
  ,((R)8.5714285714285715E-2)
};

static const R L02b[7] = {
  ((R)0.)
  ,((R)0.)
  ,((R)0.)
  ,((R)0.16666666666666666)
  ,((R)0.16666666666666666)
  ,((R)0.16666666666666666)
  ,((R)0.22500000000000001)
};

static const R L12b[21] = {
  ((R)1.6666666666666666E-2)
  ,((R)-8.3333333333333332E-3)
  ,((R)-8.3333333333333332E-3)
  ,((R)6.6666666666666666E-2)
  ,((R)3.3333333333333333E-2)
  ,((R)6.6666666666666666E-2)
  ,((R)7.4999999999999997E-2)
  ,((R)-8.3333333333333332E-3)
  ,((R)1.6666666666666666E-2)
  ,((R)-8.3333333333333332E-3)
  ,((R)6.6666666666666666E-2)
  ,((R)6.6666666666666666E-2)
  ,((R)3.3333333333333333E-2)
  ,((R)7.4999999999999997E-2)
  ,((R)-8.3333333333333332E-3)
  ,((R)-8.3333333333333332E-3)
  ,((R)1.6666666666666666E-2)
  ,((R)3.3333333333333333E-2)
  ,((R)6.6666666666666666E-2)
  ,((R)6.6666666666666666E-2)
  ,((R)7.4999999999999997E-2)
};

static const R L22b[42] = {
  ((R)1.6666666666666666E-2)
  ,((R)-2.7777777777777779E-3)
  ,((R)-2.7777777777777779E-3)
  ,((R)0.)
  ,((R)-1.1111111111111112E-2)
  ,((R)0.)
  ,((R)-1.0714285714285714E-2)
  ,((R)-2.7777777777777779E-3)
  ,((R)1.6666666666666666E-2)
  ,((R)-2.7777777777777779E-3)
  ,((R)0.)
  ,((R)0.)
  ,((R)-1.1111111111111112E-2)
  ,((R)-1.0714285714285714E-2)
  ,((R)-2.7777777777777779E-3)
  ,((R)-2.7777777777777779E-3)
  ,((R)1.6666666666666666E-2)
  ,((R)-1.1111111111111112E-2)
  ,((R)0.)
  ,((R)0.)
  ,((R)-1.0714285714285714E-2)
  ,((R)0.)
  ,((R)0.)
  ,((R)-1.1111111111111112E-2)
  ,((R)8.8888888888888892E-2)
  ,((R)4.4444444444444446E-2)
  ,((R)4.4444444444444446E-2)
  ,((R)8.5714285714285715E-2)
  ,((R)-1.1111111111111112E-2)
  ,((R)0.)
  ,((R)0.)
  ,((R)4.4444444444444446E-2)
  ,((R)8.8888888888888892E-2)
  ,((R)4.4444444444444446E-2)
  ,((R)8.5714285714285715E-2)
  ,((R)0.)
  ,((R)-1.1111111111111112E-2)
  ,((R)0.)
  ,((R)4.4444444444444446E-2)
  ,((R)4.4444444444444446E-2)
  ,((R)8.8888888888888892E-2)
  ,((R)8.5714285714285715E-2)
};

static const R L2b2b[49] = {
  ((R)1.6666666666666666E-2)
  ,((R)-2.7777777777777779E-3)
  ,((R)-2.7777777777777779E-3)
  ,((R)0.)
  ,((R)-1.1111111111111112E-2)
  ,((R)0.)
  ,((R)-1.0714285714285714E-2)
  ,((R)-2.7777777777777779E-3)
  ,((R)1.6666666666666666E-2)
  ,((R)-2.7777777777777779E-3)
  ,((R)0.)
  ,((R)0.)
  ,((R)-1.1111111111111112E-2)
  ,((R)-1.0714285714285714E-2)
  ,((R)-2.7777777777777779E-3)
  ,((R)-2.7777777777777779E-3)
  ,((R)1.6666666666666666E-2)
  ,((R)-1.1111111111111112E-2)
  ,((R)0.)
  ,((R)0.)
  ,((R)-1.0714285714285714E-2)
  ,((R)0.)
  ,((R)0.)
  ,((R)-1.1111111111111112E-2)
  ,((R)8.8888888888888892E-2)
  ,((R)4.4444444444444446E-2)
  ,((R)4.4444444444444446E-2)
  ,((R)8.5714285714285715E-2)
  ,((R)-1.1111111111111112E-2)
  ,((R)0.)
  ,((R)0.)
  ,((R)4.4444444444444446E-2)
  ,((R)8.8888888888888892E-2)
  ,((R)4.4444444444444446E-2)
  ,((R)8.5714285714285715E-2)
  ,((R)0.)
  ,((R)-1.1111111111111112E-2)
  ,((R)0.)
  ,((R)4.4444444444444446E-2)
  ,((R)4.4444444444444446E-2)
  ,((R)8.8888888888888892E-2)
  ,((R)8.5714285714285715E-2)
  ,((R)-1.0714285714285714E-2)
  ,((R)-1.0714285714285714E-2)
  ,((R)-1.0714285714285714E-2)
  ,((R)8.5714285714285715E-2)
  ,((R)8.5714285714285715E-2)
  ,((R)8.5714285714285715E-2)
  ,((R)0.14464285714285716)
};

static const R L32b[70] = {
  ((R)5.9523809523809521E-3)
  ,((R)5.9523809523809529E-4)
  ,((R)5.9523809523809529E-4)
  ,((R)2.3809523809523812E-3)
  ,((R)4.7619047619047623E-3)
  ,((R)2.3809523809523812E-3)
  ,((R)2.6785714285714286E-3)
  ,((R)5.9523809523809529E-4)
  ,((R)5.9523809523809521E-3)
  ,((R)5.9523809523809529E-4)
  ,((R)2.3809523809523812E-3)
  ,((R)2.3809523809523812E-3)
  ,((R)4.7619047619047623E-3)
  ,((R)2.6785714285714286E-3)
  ,((R)5.9523809523809529E-4)
  ,((R)5.9523809523809529E-4)
  ,((R)5.9523809523809521E-3)
  ,((R)4.7619047619047623E-3)
  ,((R)2.3809523809523812E-3)
  ,((R)2.3809523809523812E-3)
  ,((R)2.6785714285714286E-3)
  ,((R)1.607142857142857E-2)
  ,((R)-1.0714285714285714E-2)
  ,((R)-3.5714285714285713E-3)
  ,((R)2.8571428571428571E-2)
  ,((R)-7.1428571428571426E-3)
  ,((R)1.4285714285714285E-2)
  ,((R)1.2053571428571429E-2)
  ,((R)-1.0714285714285714E-2)
  ,((R)1.607142857142857E-2)
  ,((R)-3.5714285714285713E-3)
  ,((R)2.8571428571428571E-2)
  ,((R)1.4285714285714285E-2)
  ,((R)-7.1428571428571426E-3)
  ,((R)1.2053571428571429E-2)
  ,((R)-3.5714285714285713E-3)
  ,((R)1.607142857142857E-2)
  ,((R)-1.0714285714285714E-2)
  ,((R)1.4285714285714285E-2)
  ,((R)2.8571428571428571E-2)
  ,((R)-7.1428571428571426E-3)
  ,((R)1.2053571428571429E-2)
  ,((R)-3.5714285714285713E-3)
  ,((R)-1.0714285714285714E-2)
  ,((R)1.607142857142857E-2)
  ,((R)-7.1428571428571426E-3)
  ,((R)2.8571428571428571E-2)
  ,((R)1.4285714285714285E-2)
  ,((R)1.2053571428571429E-2)
  ,((R)-1.0714285714285714E-2)
  ,((R)-3.5714285714285713E-3)
  ,((R)1.607142857142857E-2)
  ,((R)-7.1428571428571426E-3)
  ,((R)1.4285714285714285E-2)
  ,((R)2.8571428571428571E-2)
  ,((R)1.2053571428571429E-2)
  ,((R)1.607142857142857E-2)
  ,((R)-3.5714285714285713E-3)
  ,((R)-1.0714285714285714E-2)
  ,((R)1.4285714285714285E-2)
  ,((R)-7.1428571428571426E-3)
  ,((R)2.8571428571428571E-2)
  ,((R)1.2053571428571429E-2)
  ,((R)-1.0714285714285714E-2)
  ,((R)-1.0714285714285714E-2)
  ,((R)-1.0714285714285714E-2)
  ,((R)8.5714285714285715E-2)
  ,((R)8.5714285714285715E-2)
  ,((R)8.5714285714285715E-2)
  ,((R)0.14464285714285716)
};

static const R L03[10] = {
  ((R)1.6666666666666666E-2)
  ,((R)1.6666666666666666E-2)
  ,((R)1.6666666666666666E-2)
  ,((R)3.7499999999999999E-2)
  ,((R)3.7499999999999999E-2)
  ,((R)3.7499999999999999E-2)
  ,((R)3.7499999999999999E-2)
  ,((R)3.7499999999999999E-2)
  ,((R)3.7499999999999999E-2)
  ,((R)0.22500000000000001)
};

static const R L13[30] = {
  ((R)8.3333333333333332E-3)
  ,((R)4.1666666666666666E-3)
  ,((R)4.1666666666666666E-3)
  ,((R)3.7499999999999999E-2)
  ,((R)0.)
  ,((R)0.)
  ,((R)0.)
  ,((R)0.)
  ,((R)3.7499999999999999E-2)
  ,((R)7.4999999999999997E-2)
  ,((R)4.1666666666666666E-3)
  ,((R)8.3333333333333332E-3)
  ,((R)4.1666666666666666E-3)
  ,((R)0.)
  ,((R)3.7499999999999999E-2)
  ,((R)3.7499999999999999E-2)
  ,((R)0.)
  ,((R)0.)
  ,((R)0.)
  ,((R)7.4999999999999997E-2)
  ,((R)4.1666666666666666E-3)
  ,((R)4.1666666666666666E-3)
  ,((R)8.3333333333333332E-3)
  ,((R)0.)
  ,((R)0.)
  ,((R)0.)
  ,((R)3.7499999999999999E-2)
  ,((R)3.7499999999999999E-2)
  ,((R)0.)
  ,((R)7.4999999999999997E-2)
};

static const R L23[60] = {
  ((R)5.9523809523809521E-3)
  ,((R)5.9523809523809529E-4)
  ,((R)5.9523809523809529E-4)
  ,((R)1.607142857142857E-2)
  ,((R)-1.0714285714285714E-2)
  ,((R)-3.5714285714285713E-3)
  ,((R)-3.5714285714285713E-3)
  ,((R)-1.0714285714285714E-2)
  ,((R)1.607142857142857E-2)
  ,((R)-1.0714285714285714E-2)
  ,((R)5.9523809523809529E-4)
  ,((R)5.9523809523809521E-3)
  ,((R)5.9523809523809529E-4)
  ,((R)-1.0714285714285714E-2)
  ,((R)1.607142857142857E-2)
  ,((R)1.607142857142857E-2)
  ,((R)-1.0714285714285714E-2)
  ,((R)-3.5714285714285713E-3)
  ,((R)-3.5714285714285713E-3)
  ,((R)-1.0714285714285714E-2)
  ,((R)5.9523809523809529E-4)
  ,((R)5.9523809523809529E-4)
  ,((R)5.9523809523809521E-3)
  ,((R)-3.5714285714285713E-3)
  ,((R)-3.5714285714285713E-3)
  ,((R)-1.0714285714285714E-2)
  ,((R)1.607142857142857E-2)
  ,((R)1.607142857142857E-2)
  ,((R)-1.0714285714285714E-2)
  ,((R)-1.0714285714285714E-2)
  ,((R)2.3809523809523812E-3)
  ,((R)2.3809523809523812E-3)
  ,((R)4.7619047619047623E-3)
  ,((R)2.8571428571428571E-2)
  ,((R)2.8571428571428571E-2)
  ,((R)1.4285714285714285E-2)
  ,((R)-7.1428571428571426E-3)
  ,((R)-7.1428571428571426E-3)
  ,((R)1.4285714285714285E-2)
  ,((R)8.5714285714285715E-2)
  ,((R)4.7619047619047623E-3)
  ,((R)2.3809523809523812E-3)
  ,((R)2.3809523809523812E-3)
  ,((R)-7.1428571428571426E-3)
  ,((R)1.4285714285714285E-2)
  ,((R)2.8571428571428571E-2)
  ,((R)2.8571428571428571E-2)
  ,((R)1.4285714285714285E-2)
  ,((R)-7.1428571428571426E-3)
  ,((R)8.5714285714285715E-2)
  ,((R)2.3809523809523812E-3)
  ,((R)4.7619047619047623E-3)
  ,((R)2.3809523809523812E-3)
  ,((R)1.4285714285714285E-2)
  ,((R)-7.1428571428571426E-3)
  ,((R)-7.1428571428571426E-3)
  ,((R)1.4285714285714285E-2)
  ,((R)2.8571428571428571E-2)
  ,((R)2.8571428571428571E-2)
  ,((R)8.5714285714285715E-2)
};

static const R L2b3[70] = {
  ((R)5.9523809523809521E-3)
  ,((R)5.9523809523809529E-4)
  ,((R)5.9523809523809529E-4)
  ,((R)1.607142857142857E-2)
  ,((R)-1.0714285714285714E-2)
  ,((R)-3.5714285714285713E-3)
  ,((R)-3.5714285714285713E-3)
  ,((R)-1.0714285714285714E-2)
  ,((R)1.607142857142857E-2)
  ,((R)-1.0714285714285714E-2)
  ,((R)5.9523809523809529E-4)
  ,((R)5.9523809523809521E-3)
  ,((R)5.9523809523809529E-4)
  ,((R)-1.0714285714285714E-2)
  ,((R)1.607142857142857E-2)
  ,((R)1.607142857142857E-2)
  ,((R)-1.0714285714285714E-2)
  ,((R)-3.5714285714285713E-3)
  ,((R)-3.5714285714285713E-3)
  ,((R)-1.0714285714285714E-2)
  ,((R)5.9523809523809529E-4)
  ,((R)5.9523809523809529E-4)
  ,((R)5.9523809523809521E-3)
  ,((R)-3.5714285714285713E-3)
  ,((R)-3.5714285714285713E-3)
  ,((R)-1.0714285714285714E-2)
  ,((R)1.607142857142857E-2)
  ,((R)1.607142857142857E-2)
  ,((R)-1.0714285714285714E-2)
  ,((R)-1.0714285714285714E-2)
  ,((R)2.3809523809523812E-3)
  ,((R)2.3809523809523812E-3)
  ,((R)4.7619047619047623E-3)
  ,((R)2.8571428571428571E-2)
  ,((R)2.8571428571428571E-2)
  ,((R)1.4285714285714285E-2)
  ,((R)-7.1428571428571426E-3)
  ,((R)-7.1428571428571426E-3)
  ,((R)1.4285714285714285E-2)
  ,((R)8.5714285714285715E-2)
  ,((R)4.7619047619047623E-3)
  ,((R)2.3809523809523812E-3)
  ,((R)2.3809523809523812E-3)
  ,((R)-7.1428571428571426E-3)
  ,((R)1.4285714285714285E-2)
  ,((R)2.8571428571428571E-2)
  ,((R)2.8571428571428571E-2)
  ,((R)1.4285714285714285E-2)
  ,((R)-7.1428571428571426E-3)
  ,((R)8.5714285714285715E-2)
  ,((R)2.3809523809523812E-3)
  ,((R)4.7619047619047623E-3)
  ,((R)2.3809523809523812E-3)
  ,((R)1.4285714285714285E-2)
  ,((R)-7.1428571428571426E-3)
  ,((R)-7.1428571428571426E-3)
  ,((R)1.4285714285714285E-2)
  ,((R)2.8571428571428571E-2)
  ,((R)2.8571428571428571E-2)
  ,((R)8.5714285714285715E-2)
  ,((R)2.6785714285714286E-3)
  ,((R)2.6785714285714286E-3)
  ,((R)2.6785714285714286E-3)
  ,((R)1.2053571428571429E-2)
  ,((R)1.2053571428571429E-2)
  ,((R)1.2053571428571429E-2)
  ,((R)1.2053571428571429E-2)
  ,((R)1.2053571428571429E-2)
  ,((R)1.2053571428571429E-2)
  ,((R)0.14464285714285716)
};

static const R L33[100] = {
  ((R)5.6547619047619046E-3)
  ,((R)8.1845238095238097E-4)
  ,((R)8.1845238095238097E-4)
  ,((R)1.3392857142857143E-3)
  ,((R)0.)
  ,((R)2.0089285714285712E-3)
  ,((R)2.0089285714285712E-3)
  ,((R)0.)
  ,((R)1.3392857142857143E-3)
  ,((R)2.6785714285714286E-3)
  ,((R)8.1845238095238097E-4)
  ,((R)5.6547619047619046E-3)
  ,((R)8.1845238095238097E-4)
  ,((R)0.)
  ,((R)1.3392857142857143E-3)
  ,((R)1.3392857142857143E-3)
  ,((R)0.)
  ,((R)2.0089285714285712E-3)
  ,((R)2.0089285714285712E-3)
  ,((R)2.6785714285714286E-3)
  ,((R)8.1845238095238097E-4)
  ,((R)8.1845238095238097E-4)
  ,((R)5.6547619047619046E-3)
  ,((R)2.0089285714285712E-3)
  ,((R)2.0089285714285712E-3)
  ,((R)0.)
  ,((R)1.3392857142857143E-3)
  ,((R)1.3392857142857143E-3)
  ,((R)0.)
  ,((R)2.6785714285714286E-3)
  ,((R)1.3392857142857143E-3)
  ,((R)0.)
  ,((R)2.0089285714285712E-3)
  ,((R)4.0178571428571432E-2)
  ,((R)-1.40625E-2)
  ,((R)-1.0044642857142858E-2)
  ,((R)-4.0178571428571425E-3)
  ,((R)-1.0044642857142858E-2)
  ,((R)2.0089285714285716E-2)
  ,((R)1.2053571428571429E-2)
  ,((R)0.)
  ,((R)1.3392857142857143E-3)
  ,((R)2.0089285714285712E-3)
  ,((R)-1.40625E-2)
  ,((R)4.0178571428571432E-2)
  ,((R)2.0089285714285716E-2)
  ,((R)-1.0044642857142858E-2)
  ,((R)-4.0178571428571425E-3)
  ,((R)-1.0044642857142858E-2)
  ,((R)1.2053571428571429E-2)
  ,((R)2.0089285714285712E-3)
  ,((R)1.3392857142857143E-3)
  ,((R)0.)
  ,((R)-1.0044642857142858E-2)
  ,((R)2.0089285714285716E-2)
  ,((R)4.0178571428571432E-2)
  ,((R)-1.40625E-2)
  ,((R)-1.0044642857142858E-2)
  ,((R)-4.0178571428571425E-3)
  ,((R)1.2053571428571429E-2)
  ,((R)2.0089285714285712E-3)
  ,((R)0.)
  ,((R)1.3392857142857143E-3)
  ,((R)-4.0178571428571425E-3)
  ,((R)-1.0044642857142858E-2)
  ,((R)-1.40625E-2)
  ,((R)4.0178571428571432E-2)
  ,((R)2.0089285714285716E-2)
  ,((R)-1.0044642857142858E-2)
  ,((R)1.2053571428571429E-2)
  ,((R)0.)
  ,((R)2.0089285714285712E-3)
  ,((R)1.3392857142857143E-3)
  ,((R)-1.0044642857142858E-2)
  ,((R)-4.0178571428571425E-3)
  ,((R)-1.0044642857142858E-2)
  ,((R)2.0089285714285716E-2)
  ,((R)4.0178571428571432E-2)
  ,((R)-1.40625E-2)
  ,((R)1.2053571428571429E-2)
  ,((R)1.3392857142857143E-3)
  ,((R)2.0089285714285712E-3)
  ,((R)0.)
  ,((R)2.0089285714285716E-2)
  ,((R)-1.0044642857142858E-2)
  ,((R)-4.0178571428571425E-3)
  ,((R)-1.0044642857142858E-2)
  ,((R)-1.40625E-2)
  ,((R)4.0178571428571432E-2)
  ,((R)1.2053571428571429E-2)
  ,((R)2.6785714285714286E-3)
  ,((R)2.6785714285714286E-3)
  ,((R)2.6785714285714286E-3)
  ,((R)1.2053571428571429E-2)
  ,((R)1.2053571428571429E-2)
  ,((R)1.2053571428571429E-2)
  ,((R)1.2053571428571429E-2)
  ,((R)1.2053571428571429E-2)
  ,((R)1.2053571428571429E-2)
  ,((R)0.14464285714285716)
};
#include "../src/libns_exact_f3_u6_p3.c"

#define	clr_mem(_n,_x) { I _k; for (_k=0;_k<(_n);++_k) { (_x)[_k] = regal0;  }}
#define	cpy_mem(_n,_x,_y) { I _k; for (_k=0;_k<(_n);++_k) { (_y)[_k] = (_x)[_k];  }}
#define	cpy_mem_indirect(_n,_x,_b,_y) { I _k; for (_k=0;_k<(_n);++_k) { (_y)[_k] = (_x)[(_b)[_k]];  }}
#define mvass(_n,_m,_ax,_a,_x,_y)  dgemv("N",_n,_m,(_ax),(_a),_n,(_x),&negal1,&regal1,(_y),&negal1)





class StokesNonLinearSolverDriver : public INonLinearSolverDriver
{
protected:
  I			m_nelm;
  I 			required_rwork_n;
  I 			required_iwork_n;
  I 			iwork_n;
  I 			blank_n;
  
  pI 			iwork;
  pI 			blank;

  const Parameters * 	m_parameters;
  
  cst_pR		m_jacelm;
  Fem 			m_fem;
  
  SparseStokes*		m_sparseStokes;
  pSpaceReadOnly	m_space_u;
  pSpaceReadOnly	m_space_p;
  pSparse		m_sparseB;
  pSparse		m_sparseF;
  pSparse		m_sparseS;
  pSparse		m_sparseC;

  Workelm		m_workelm;
  Time			m_timeInfo;
  SolverConstants	m_solverConstants;

  /* \brief shift to access to u	*/	
  I 			m_dec_ddlu;
  /* \brief shift to access to p	*/	
  I 			m_dec_ddlp;
  I 			m_nddlu;
  I 			m_nedge;
  pR 			m_jacedge;
  pR 			m_normaledge;
  pI			m_slip_iperm;
  I 			m_total_nddl;
  void *		m_mesh_usrptr;  

public:  StokesNonLinearSolverDriver(pGlobal self);

public: void rhselm2rhs(const I 	jelm_,
			pR 		global_rhs_);
  
public:   virtual void * GetUsrPtr();

public:   inline void SetParams(const Parameters* parameters_);
  
public:   inline const Parameters* Params();
  
public:   virtual void BuildSystem(cst_pI 		nx_,
			   cst_pR 		x_,
			   cst_pI 		xoff_);
    

public:   void info();
  
public:   virtual void BuildResidu(pR 			global_rhs_,
				   cst_pR 		x_,
				   cst_pR 		xi_,
				   cst_pR 		xii_);

public:   void BuildResiduelmVelocity(const I 	jelm_,
				      cst_pR 	denselm_,
				      cst_pR 	viscelm_,
				      cst_pR 	uelm_,
				      cst_pR 	velm_,
				      cst_pR 	pelm_,
				      pR 	residu_u_,
				      pR 	residu_v_);  
  

public:   void ResiduExactDensityTimeDerivative(const I 		jelm_,
						cst_pR 			xa_,
						cst_pR 			denselm_,
						cst_pR 			uelm_,
						cst_pR 			uelmi_,
						cst_pR 			uelmii_,
						pR 			residu_,
						const KindTransientMethod::EnumType kindTransientMethod);
  
public:   void 	Setelm			(const I 		jelm_,
					 cst_pR 		x_,
					 cst_pR 		xi_,
					 cst_pR 		xii_);

  
public:   void 	get_ddlcnc_u		(const I 		ielm_,
					 pI 			ddl_);

public:   void 	get_ddlcnc_v		(const I 		ielm_,
					 pI 			ddl_);
  
public:   void 	get_ddlcnc_p		(const I 		ielm_,
					 pI 			ddl_);
  
  
public:   void get_ddlcnc(const I 		ielm_,
		  pI 			ddl_);


public:   void Linsys_clear();
public:   void Build_system_stokes_quadrature_free(cst_pR 	x_,
					   cst_pR 	xi_,
					   cst_pR 	xii_);

public:   void assmatelm_full(const I 	jelm_);
public:   void compute_viscosity_density_nodalelm2(cst_pR felm_);

};



void StokesNonLinearSolverDriver::rhselm2rhs(const I 	jelm_,
					     pR 	global_rhs_)
{
  auto local_ddlcnc = this->m_workelm.locnumer;
  this->get_ddlcnc(jelm_,
		   local_ddlcnc);  
  for (I i=0;i<_total_nddlelm;++i)
    {
      global_rhs_[local_ddlcnc[i]] += this->m_workelm.locresidu[i];
    } 
};


StokesNonLinearSolverDriver::StokesNonLinearSolverDriver(pGlobal self)
{
  //
  // Copy from self
  //
  this->m_sparseStokes = self->sparseStokes;
  this->m_dec_ddlu =  self->dec_ddlu;
  this->m_nddlu =  self->nddlu;
  this->m_dec_ddlp =  self->dec_ddlp;

  m_mesh_usrptr = self->mesh_usrptr;
  m_space_u = self->space_u;
  m_space_p = self->space_p;

  m_solverConstants[SolverConstants::iweber]=    self->rv[__ens_rv_iweber];
  m_solverConstants[SolverConstants::ireynold]=    self->rv[__ens_rv_ireynold];
  m_solverConstants[SolverConstants::ifroude]=    self->rv[__ens_rv_ifroude];
  m_solverConstants[SolverConstants::vmin]=    self->rv[__ens_rv_vmin];
  m_solverConstants[SolverConstants::vmax]=    self->rv[__ens_rv_vmax];
  m_solverConstants[SolverConstants::coeff_ratio_viscosity]=    self->rv[__ens_rv_coeff_ratio_viscosity];
  m_solverConstants[SolverConstants::coeff_ratio_density]=    self->rv[__ens_rv_coeff_ratio_density];
  m_solverConstants[SolverConstants::dt]=    self->rv[__ens_rv_dt];
  m_solverConstants[SolverConstants::dti]=    self->rv[__ens_rv_dti];
  m_solverConstants[SolverConstants::idt]=    self->rv[__ens_rv_idt];
  m_solverConstants[SolverConstants::midt]=    self->rv[__ens_rv_midt];
  m_solverConstants[SolverConstants::ratio_dti]=    self->rv[__ens_rv_ratio_dti];
  m_solverConstants[SolverConstants::iratio_dti]=    self->rv[__ens_rv_iratio_dti];
  m_solverConstants[SolverConstants::eps]=    self->rv[__ens_rv_eps];
  m_solverConstants[SolverConstants::ieps]=    self->rv[__ens_rv_ieps];
  
  m_total_nddl = self->total_nddl;
  m_normaledge = (pR)self->normaledge;
  m_jacedge = (pR)self->jacedge;
  blank = self->blank;
  blank_n = self->blank_n;
  iwork = self->iwork;
  iwork_n = self->iwork_n;
  m_nelm = self->nelm;
  m_jacelm = self->jacelm;
  m_nedge =   self->nedge;
  Workelm_init(&this->m_workelm);
  m_fem = self->m_fem;

  
};
  
void * StokesNonLinearSolverDriver::GetUsrPtr()
{
  return NULL;
};

inline void StokesNonLinearSolverDriver::SetParams(const Parameters* parameters_)
{
  this->m_parameters = parameters_;
}
  
inline const Parameters* StokesNonLinearSolverDriver::Params()
{
  return this->m_parameters;
}
  
 void StokesNonLinearSolverDriver::BuildSystem(cst_pI 		nx_,
						      cst_pR 		x_,
						      cst_pI 		xoff_)
{
#ifndef NDEBUG
  std::cout << "BUILD SYSTEM(nx = " << nx_ << ", x" << x_ << ", xoff " << xoff_[0] << ")"  <<   std::endl;
#endif
  const bool pressure_uncoupled		= this->Params()->GetInfoLogical(InfoLogical::pressure_uncoupled);
  const bool pressure_freematrix		= this->Params()->GetInfoLogical(InfoLogical::pressure_freematrix);
  const bool slip				= this->Params()->GetInfoLogical(InfoLogical::slip);
#ifndef NDEBUG
  std::cout << " - pressure_uncoupled  : " << pressure_uncoupled << std::endl;
  std::cout << " - pressure_freematrix : " << pressure_freematrix << std::endl;
  std::cout << " - slip                : " << slip << std::endl;
  std::cout << " - nddlu               : " << m_nddlu << std::endl;
  std::cout << " - dec_ddlu            : " << m_dec_ddlu << std::endl;
  std::cout << " - dec_ddlp            : " << m_dec_ddlp << std::endl;
  std::cout << " - nddlp_dirichelt     : " << this->m_sparseStokes->nddlp_dirichlet << std::endl;
  std::cout << " - nddlu_dirichlet     : " << this->m_sparseStokes->nddlu_dirichlet << std::endl;
#endif    
  //
  // Clean the linear system.
  //
#ifndef NDEBUG
  std::cout << "BUILD SYSTEM: CLEAR" << std::endl;
#endif
  this->Linsys_clear();

  //
  // build jacobian
  //
#ifndef NDEBUG
  std::cout << "BUILD SYSTEM QUADRATURE FREE " << std::endl;
#endif
  this->Build_system_stokes_quadrature_free(x_,
					    &x_[xoff_[0]],
					    &x_[2*xoff_[0]]);
#ifndef NDEBUG
  std::cout << "BUILD SYSTEM APPLY DIRICHLET CONDITION " << std::endl;
#endif
  //
  // dirichlet condition velocity
  //
  Sparse_dirichlet(this->m_sparseStokes->A,
		   this->m_sparseStokes->nddlu_dirichlet,
		   this->m_sparseStokes->ddlu_dirichlet,
		   this->m_dec_ddlu);
    
  Sparse_dirichlet(this->m_sparseStokes->A,
		   this->m_sparseStokes->nddlv_dirichlet,
		   this->m_sparseStokes->ddlv_dirichlet,
		   this->m_dec_ddlu+this->m_nddlu);
    
  /* 
     TANT QUE LES CONDITIONS DE DIICHLETS NE CHANGENT PAS SYMBOLIQUEMENT EN TEMPS
     CE QUI EST CI DESSOUS EST FAIT A CHAQUE PAS DE TEMPS POUR RIEN
     SI LES CONDITIONS CHANGENT IL FAUDRA RECALCULER LA MATRICE ECRASEE PRECEDENTE
     SUR LES MATRICES LINEAIRES
  */
    
  /* CONDITION DIRICHLET EN P */
  if (false == pressure_uncoupled)
    {
#ifndef NDEBUG
      std::cout << "BUILD SYSTEM APPLY DIRICHLET CONDITION ON P " << std::endl;
#endif
      Sparse_dirichlet(this->m_sparseStokes->A,
		       this->m_sparseStokes->nddlp_dirichlet,
		       this->m_sparseStokes->ddlp_dirichlet,this->m_dec_ddlp);	
    }
  else if ( (true == pressure_uncoupled)
	    &&
	    (false == pressure_freematrix) )
    {
      Sparse_dirichlet(this->m_sparseB,
		       this->m_sparseStokes->nddlp_dirichlet,
		       this->m_sparseStokes->ddlp_dirichlet,
		       this->m_dec_ddlp);
    }
    
  /* conditions de slip*/
  if (true == slip)
    {	    
      cst_pI 	Ab = Sparse_get_b(this->m_sparseStokes->A);
      cst_pI 	Ai = Sparse_get_i(this->m_sparseStokes->A);
      pR 	Ax = Sparse_get_ownx(this->m_sparseStokes->A);

      for (I i=0;i<this->m_sparseStokes->nddls_dirichlet;++i)
	{		      
	  const I k = this->m_slip_iperm[this->m_sparseStokes->ddls_dirichlet[i]];
	  if (0 == k)
	    {
	      std::cerr << "erreur slip iperm, iddl n est pas frontiere ..." << std::endl;
	      exit(1);
	    }
	  else
	    {
	      const I k2=this->m_total_nddl+k-1;
	      for (I j=Ab[k2];j<Ab[k2+1];++j)
		{
		  Ax[j-1] = ((R)0.0);
		}

	      bool found = false;
	      for (I j=Ab[k2];j<Ab[k2+1];++j)
		{
		  if (Ai[j-1]!=k2+1)
		    {
		      found = true;
		      Ax[j-1] =  ((R)1.0);
		      break;
		    }
		}
	      if (false == found)
		{
		  printf("wrong algo build system missing diagonal\n");
		}		
	    }
	} 
    }

#if 0
  {
    FILE * f = fopen("spy.txt","w");
    cst_pI 	Ab = Sparse_get_b(m_sparseStokes->A);
    cst_pI 	Ai = Sparse_get_i(m_sparseStokes->A);
    pR 	Ax = Sparse_get_ownx(m_sparseStokes->A);
    I M = Sparse_get_n(m_sparseStokes->A);
    for (I i=0;i<M;++i)
      {
	for (I at = Ab[i]-1;at < Ab[i+1]-1;++at)
	  {
	    I index = Ai[at];
	    R x = Ax[at];
	    fprintf(f,""ifmt" "ifmt" %e\n",M-i,index,x);
	  }
      }
    fclose(f);
    printf("see spy.txt\n");
    exit(1);
  }
#endif

    
};
    

void StokesNonLinearSolverDriver::info()
{  
  std::cout << "INFO ---" << std::endl;
}
 void StokesNonLinearSolverDriver::BuildResidu(pR 			global_rhs_,
						      cst_pR 		x_,
						      cst_pR 		xi_,
						      cst_pR 		xii_)
{
#ifndef NDEBUG
  std::cout << "BUILD RESIDU(global_rhs_ = " << global_rhs_ << ", x = " << x_ << ", xi_ = " << xi_ << ", xii_ = " << xii_ << ")"  <<   std::endl;
#endif
  pWorkelm const 		workelm		= &this->m_workelm;
  pFem const			fem 		= &this->m_fem;

  const KindTransientMethod::EnumType kindTransientMethod = this->Params()->GetKindTransientMethod(KindEquation::VELOCITY);
  const bool verbose 				= this->Params()->GetInfoLogical(InfoLogical::verbose);
#ifndef NDEBUG    
  std::cout << " - kindTransientMethod  : " << kindTransientMethod << std::endl;
  std::cout << " - verbose  : " << verbose << std::endl;
#endif
    
  if (verbose)
    {
      std::cout << "ns_build_residu:compute linear pressure residual ..."<< std::endl;
    }
    
#ifndef NDEBUG
  std::cout << "BUILD RESIDU - PRESSURE RHS" << std::endl;
#endif
  SparseStokes_pressure_rhs(this->m_sparseStokes,
			    x_,
			    global_rhs_);
#ifndef NDEBUG    
  std::cout << "BUILD RESIDU - PRESSURE RHS DONE." << std::endl;
#endif

  const R coeff_ratio_viscosity = m_solverConstants[SolverConstants::coeff_ratio_viscosity];
  const R coeff_ratio_density = m_solverConstants[SolverConstants::coeff_ratio_density];

  /*
    on calcule residu en fuv
  */
  for (I jelm=0;jelm<this->m_nelm;++jelm)
    {      
	
      this->Setelm(jelm,
		   x_,
		   xi_,
		   xii_);
#if 0
      printf("UELM line 1238 %p\n",this->m_workelm.uelm);
#endif	  
      Fem_setelm(fem,
		 &this->m_workelm);	
	  
#if 0
      integrator_compute_viscosity_density_nodalelm(self_->iproc,
						    elm->dgelm);
#endif
      /* printf("%e %e\n",Z->rv[__ens_rv_coeff_ratio_viscosity],Z->rv[__ens_rv_coeff_ratio_density]);*/
#if 0
      nsDG_transport_get_nodalvalues(self_->transport_usrptr,
				     &jelm,
				     self_->workelm.denselm);
      nsDG_transport_get_nodalvalues(self_->transport_usrptr,
				     &jelm,
				     self_->workelm.viscelm);
#endif
#if 0
      printf("UELM line 1238 %p\n",this->m_workelm.uelm);
#endif
      for (I k=0;k<_ndg;++k)						
	{
	  R hh = workelm->viscelm[k];
	  if (hh<((R)0.0))
	    {
	      hh=((R)0.0);
	    }
	  if (hh>((R)1.0))
	    {
	      hh=((R)1.0);
	    }
	      
	  workelm->viscelm[k] 	= 1.0+hh*(coeff_ratio_viscosity-1.0);
	  workelm->denselm[k] 	= 1.0+hh*(coeff_ratio_density-1.0);
	}

	  
      /*
	v1 + F(v2-v1) 
	r + F(1-r)
	printf("%e\n",Z->rv[__ens_rinfo_ratio_viscosity]); 
      */
	
      clr_mem(_total_nddlelm,
	      this->m_workelm.locresidu);

	
      this->BuildResiduelmVelocity		(jelm,
						 this->m_workelm.denselm,
						 this->m_workelm.viscelm,
						 this->m_workelm.uelm,
						 this->m_workelm.velm,
						 this->m_workelm.pelm,
						 &this->m_workelm.locresidu[offu],
						 &this->m_workelm.locresidu[offv]); 

      static const R requal1 = ((R)1.0);
	
      this->ResiduExactDensityTimeDerivative	(jelm,
						 &requal1,
						 this->m_workelm.denselm,
						 this->m_workelm.uelm,
						 this->m_workelm.uelmi,
						 this->m_workelm.uelmii,
						 &this->m_workelm.locresidu[offu],
						 kindTransientMethod);

      this->ResiduExactDensityTimeDerivative	(jelm,
						 &requal1,
						 this->m_workelm.denselm,
						 this->m_workelm.velm,
						 this->m_workelm.velmi,
						 this->m_workelm.velmii,
						 &this->m_workelm.locresidu[offv],
						 kindTransientMethod);
	
      this->rhselm2rhs			(jelm,
					 global_rhs_);
	
    }
    
  //  printf("putain de ta race\n");
  //  exit(1);

  /* RESIDU SLIP CONDITION */
    
  const bool hasSlipCondition = this->Params()->GetInfoLogical(InfoLogical::slip);
  if (hasSlipCondition)
    {

      if (verbose)
	{
	  std::cout << "ns_build_residu:compute slip boundary condition ..." << std::endl;
	}

      SparseStokes_slip_rhs(this->m_sparseStokes,
			    x_,
			    global_rhs_);

    }   
  else
    {
      if (verbose)
	{
	  std::cout << "ns_build_residu:no slip boundary condition to compute" << std::endl;
	}
    }
  
  /* IMPOSITION DIRICHLET */
  if (verbose)
    {
      std::cout << "ns_build_residu:compute boundary dirichlet condition" << std::endl;
    }
  
  /* DIRICHLET CONDITION */
  SparseStokes_dirichlet_rhs(this->m_sparseStokes,
			     global_rhs_);

#if 0
  for (I i=0;i<this->m_sparseStokes->A->n;++i)
    {
      std::cout << "residu["<< i << "] = "  << global_rhs_[i] << std::endl;
    }
  exit(1);
#endif
};


void StokesNonLinearSolverDriver::BuildResiduelmVelocity(const I 	jelm_,
							 cst_pR 	denselm_,
							 cst_pR 	viscelm_,
							 cst_pR 	uelm_,
							 cst_pR 	velm_,
							 cst_pR 	pelm_,
							 pR 	residu_u_,
							 pR 	residu_v_)
{

#if 0
#if ( (_ndg==1) AND (_nu==6) )
  cst_pR gravity_term = L02;
#elif ( (_ndg==3) AND (_nu==6) )
  cst_pR gravity_term = L12;
#elif ( (_ndg==6) AND (_nu==6) )
  cst_pR gravity_term = L22;
#elif ( (_ndg==10) AND (_nu==6) )
  cst_pR gravity_term = L32;
#elif ( (_ndg==1) AND (_nu==7) )
  cst_pR gravity_term = L02b;
#elif ( (_ndg==3) AND (_nu==7) )
  cst_pR gravity_term = L12b;
#elif ( (_ndg==6) AND (_nu==7) )
  cst_pR gravity_term = L22b;
#elif ( (_ndg==10) AND (_nu==7) ) 
  cst_pR gravity_term = L32b;
#endif
#endif 
  pFem const	fem 		= &this->m_fem;
    
  static const I nequal1 = 1;
  static const R requal0 = ((R)0.0);
  static const R requal1 = ((R)1.0);
  static const R requal2 = ((R)2.0);

#if 0
  const R 		ifroude 	= this->m_solverConstants[SolverConstants::ifroude];
#endif
  const R 		ireynold 	= this->m_solverConstants[SolverConstants::ireynold];

  nsblas_dgemv(transN,&nuxnu,&ndg,&requal1,fem->exact_visc_dxu_dxu,&nuxnu,viscelm_,&nequal1,&requal0,fem->visc_dxu_dxu,&nequal1);
  nsblas_dgemv(transN,&nuxnu,&ndg,&requal1,fem->exact_visc_dyu_dyu,&nuxnu,viscelm_,&nequal1,&requal0,fem->visc_dyu_dyu,&nequal1);
  nsblas_dgemv(transN,&nuxnu,&ndg,&requal1,fem->exact_visc_dyu_dxu,&nuxnu,viscelm_,&nequal1,&requal0,fem->visc_dyu_dxu,&nequal1);
  nsblas_dgemv(transN,&nuxnu,&ndg,&requal1,fem->exact_visc_dxu_dyu,&nuxnu,viscelm_,&nequal1,&requal0,fem->visc_dxu_dyu,&nequal1);
  
#if 0
  std::cout << "ifroude: "<< ifroude << std::endl;
  /* terme de gravite (0,-1) */
  { const R a = this->m_jacelm[jelm_]*ifroude;
    std::cout << "a: "<< a << std::endl;
    nsblas_dgemv(transN,
		 &nu,
		 &ndg,
		 &a,
		 gravity_term,
		 &nu,
		 denselm_,
		 &nequal1,
		 &requal1,
		 residu_v_,
		 &nequal1); }
#endif
    
  /* residu pressure -/- u */
  { const R a = ((R)-1.0);

    // RESIDU_U -= DXU_P * P_elm
    // RESIDU_V -= DYU_P * P_elm
    
    nsblas_dgemv(transT,
		 &np,
		 &nu,
		 &a,
		 fem->exact_dxu_p,
		 &np,
		 pelm_,
		 &nequal1,
		 &requal1,
		 residu_u_,
		 &nequal1);

    nsblas_dgemv(transT,
		 &np,
		 &nu,
		 &a,
		 fem->exact_dyu_p,
		 &np,
		 pelm_,
		 &nequal1,
		 &requal1,
		 residu_v_,
		 &nequal1); }

  /* BUILD STRESS stress */

  // RESIDU_U += 2.0 * ireynold * VISC_DXU_DXU * U_elm
  // RESIDU_V += ireynold * VISC_DXU_DXU * V_elm
  
  { const R a = ireynold * requal2;
    nsblas_dgemv(transN,&nu,&nu,&a,fem->visc_dxu_dxu,&nu,uelm_,&nequal1,&requal1,residu_u_,&nequal1); }
  { const R a = ireynold;
    nsblas_dgemv(transN,&nu,&nu,&a,fem->visc_dxu_dxu,&nu,velm_,&nequal1,&requal1,residu_v_,&nequal1); }
    
  // RESIDU_U += ireynold *       VISC_DYU_DYU * U_elm
  // RESIDU_V += 2.0 * ireynold * VISC_DYU_DYU * V_elm
    

  { const R a = ireynold * requal2;
    nsblas_dgemv(transN,&nu,&nu,&a,fem->visc_dyu_dyu,&nu,velm_,&nequal1,&requal1,residu_v_,&nequal1); }
  { const R a = ireynold;
    nsblas_dgemv(transN,&nu,&nu,&a,fem->visc_dyu_dyu,&nu,uelm_,&nequal1,&requal1,residu_u_,&nequal1); }
    
  // RESIDU_U += ireynold * VISC_DXU_DYU * U_elm
  // RESIDU_V += ireynold * VISC_DYU_DXU * V_elm
    
  { const R a = ireynold;
    nsblas_dgemv(transN,&nu,&nu,&a,fem->visc_dyu_dxu,&nu,uelm_,&nequal1,&requal1,residu_v_,&nequal1); 
    nsblas_dgemv(transN,&nu,&nu,&a,fem->visc_dxu_dyu,&nu,velm_,&nequal1,&requal1,residu_u_,&nequal1); }
    
  /* residu d'acceleration */

#if 0
  nsblas_dgemv(transN,&nuxnuxnu,&ndg,&requal1,fem->exact_dens_u_dxu_u,&nuxnuxnu,denselm_,&nequal1,&requal0,fem->dens_u_dxu_u,&nequal1);
  nsblas_dgemv(transN,&nuxnuxnu,&ndg,&requal1,fem->exact_dens_u_dyu_u,&nuxnuxnu,denselm_,&nequal1,&requal0,fem->dens_u_dyu_u,&nequal1);
  /**/
  nsblas_dgemv(transN,&nuxnu,&nu,&requal1,fem->dens_u_dxu_u,&nuxnu,uelm_,&nequal1,&requal0,fem->dens_dxu_u,&nequal1);
  nsblas_dgemv(transN,&nuxnu,&nu,&requal1,fem->dens_u_dyu_u,&nuxnu,velm_,&nequal1,&requal0,fem->dens_dyu_u,&nequal1);
  /**/	  
  nsblas_dgemv(transN,&nu,&nu,&requal1,fem->dens_dxu_u,&nu,uelm_,&nequal1,&requal1,residu_u_,&nequal1);
  nsblas_dgemv(transN,&nu,&nu,&requal1,fem->dens_dyu_u,&nu,uelm_,&nequal1,&requal1,residu_u_,&nequal1);
  nsblas_dgemv(transN,&nu,&nu,&requal1,fem->dens_dxu_u,&nu,velm_,&nequal1,&requal1,residu_v_,&nequal1);
  nsblas_dgemv(transN,&nu,&nu,&requal1,fem->dens_dyu_u,&nu,velm_,&nequal1,&requal1,residu_v_,&nequal1);
#endif
  
};



void StokesNonLinearSolverDriver::ResiduExactDensityTimeDerivative(const I 		jelm_,
								   cst_pR 		xa_,
								   cst_pR 		denselm_,
								   cst_pR 		uelm_,
								   cst_pR 		uelmi_,
								   cst_pR 		uelmii_,
								   pR 			residu_,
								   const KindTransientMethod::EnumType kindTransientMethod)
{    
  pTimeReadOnly const gTimeInfo 	= &this->m_timeInfo;
  pFem const	fem = &this->m_fem;      
  const R aa = xa_[0] * this->m_jacelm[jelm_];
  switch(kindTransientMethod)
    {
	
    case KindTransientMethod::UNDEFINED:
      {
	break;
      }
	
    case KindTransientMethod::EULER:
      {

	const R idt = this->m_solverConstants[SolverConstants::idt];
	const R a = idt * aa;
	  
	R tmp[_nu];
	{ I k;
	  for (k=0;k<_nu;++k)
	    {		
	      tmp[k] = (uelm_[k]-uelmi_[k]);
	    } }

	static const R requal0 = ((R)0.0);
	static const R requal1 = ((R)1.0);
	static const I nequal1 = 1;
	  
	nsblas_dgemv	(transN,
			 &nuxnu,
			 &ndg,
			 &requal1,
			 exactref_dens_u_u,
			 &nuxnu,
			 denselm_,
			 &nequal1,
			 &requal0,
			 Fem_get_dens_u_u(fem),
			 &nequal1);
	
	mvass(&nu,
	      &nu,
	      &a,
	      FemReadOnly_get_dens_u_u(fem),
	      tmp,
	      residu_);
	  
	break;
      }
	
    case KindTransientMethod::IMR:
      {
	const R idt = this->m_solverConstants[SolverConstants::idt];
	const R a = idt*((R)2.0);
	  
	R tmp[_nu];
	{ I k;
	  for (k=0;k<_nu;++k)
	    {		
	      tmp[k] = (uelm_[k]-uelmi_[k]);
	    } }

	static const R requal0 = ((R)0.0);
	static const R requal1 = ((R)1.0);
	static const I nequal1 = 1;

	  
	nsblas_dgemv	(transN,
			 &nuxnu,
			 &ndg,
			 &requal1,
			 exactref_dens_u_u,
			 &nuxnu,
			 denselm_,
			 &nequal1,
			 &requal0,
			 Fem_get_dens_u_u(fem),
			 &nequal1);

	mvass		(&nu,
			 &nu,
			 &a,
			 FemReadOnly_get_dens_u_u(fem),
			 tmp,
			 residu_);
	
	break;
      }
      
    case KindTransientMethod::GEAREULER:
      {

	static const R requal0 = ((R)0.0);
	static const R requal1 = ((R)1.0);
	static const R requal2 = ((R)2.0);
	static const I nequal1 = 1;
	
	if (gTimeInfo->itime>0)
	  {
	    const R a = aa;

	    const R ratio_dti = this->m_solverConstants[SolverConstants::ratio_dti];
	    const R dt = this->m_solverConstants[SolverConstants::dt];
	    const R dti = this->m_solverConstants[SolverConstants::dti];
	    const R iratio_dti = this->m_solverConstants[SolverConstants::iratio_dti];
	    
	    R tmp[_nu];
	    for (I k=0;k<_nu;++k)
	      {		
		tmp[k] = ( (requal2+ratio_dti)*uelm_[k]-(requal2 + ratio_dti + iratio_dti)*uelmi_[k] + iratio_dti*uelmii_[k] )/(dt+dti);
	      } 

	    nsblas_dgemv	(transN,
				 &nuxnu,
				 &ndg,
				 &requal1,
				 exactref_dens_u_u,
				 &nuxnu,
				 denselm_,
				 &nequal1,
				 &requal0,
				 Fem_get_dens_u_u(fem),
				 &nequal1);
	    
	    mvass		(&nu,
				 &nu,
				 &a,
				 FemReadOnly_get_dens_u_u(fem),
				 tmp,
				 residu_);
	    
	  }
	else
	  {
	    const R idt = this->m_solverConstants[SolverConstants::idt];
	    const R a = idt*aa;
	    R tmp[_nu];
	    for (I k=0;k<_nu;++k)
	      {		
		tmp[k] = (uelm_[k]-uelmi_[k]);
	      } 
	    nsblas_dgemv(transN,&nuxnu,&ndg,&requal1,exactref_dens_u_u,&nuxnu,denselm_,&nequal1,&requal0,Fem_get_dens_u_u(fem),&nequal1);
	    mvass(&nu,&nu,&a,FemReadOnly_get_dens_u_u(fem),tmp,residu_);
	  }
	break;
      }

    case KindTransientMethod::GEARIMR:
      {
	
	static const R requal0 = ((R)0.0);
	static const R requal1 = ((R)1.0);
	static const R requal2 = ((R)2.0);
	static const I nequal1 = 1;

	if (gTimeInfo->itime>0)
	  {
	    const R a = aa;

	    
	    const R ratio_dti = this->m_solverConstants[SolverConstants::ratio_dti];
	    const R dt = this->m_solverConstants[SolverConstants::dt];
	    const R dti = this->m_solverConstants[SolverConstants::dti];
	    const R iratio_dti = this->m_solverConstants[SolverConstants::iratio_dti];
	    
	    R tmp[_nu];
	    for (I k=0;k<_nu;++k)
	      {		
		tmp[k] = ( (requal2+ratio_dti)*uelm_[k]-(requal2 + ratio_dti + iratio_dti)*uelmi_[k] + iratio_dti*uelmii_[k] )/(dt+dti);
	      }
	    
	    nsblas_dgemv(transN,
			 &nuxnu,
			 &ndg,
			 &requal1,
			 exactref_dens_u_u,
			 &nuxnu,
			 denselm_,
			 &nequal1,
			 &requal0,
			 Fem_get_dens_u_u(fem),
			 &nequal1);

	    mvass(&nu,
		  &nu,
		  &a,
		  FemReadOnly_get_dens_u_u(fem),
		  tmp,
		  residu_);	    
	  }
	else
	  {
	    const R idt = this->m_solverConstants[SolverConstants::idt];
	    const R a = aa*idt*((R)2.0);
	    R tmp[_nu];
	    { I k;
	      for (k=0;k<_nu;++k)
		{		
		  tmp[k] = (uelm_[k]-uelmi_[k]);
		} }
	    
	    nsblas_dgemv(transN,
			 &nuxnu,
			 &ndg,
			 &requal1,
			 exactref_dens_u_u,
			 &nuxnu,
			 denselm_,
			 &nequal1,
			 &requal0,
			 Fem_get_dens_u_u(fem),
			 &nequal1);
	    
	    mvass(&nu,
		  &nu,
		  &a,
		  FemReadOnly_get_dens_u_u(fem),
		  tmp,
		  residu_);
	    
	  }
	break;
      }
      
    case KindTransientMethod::TRAPEZE:
      {
	std::cerr << "__eTransientMethod_TRAPEZE not yet" << std::endl;
	exit(1);
	break;
      }
      
    case KindTransientMethod::ERROR:
    case KindTransientMethod::ALL_VALUES:
      {
	std::cerr << "ns_residu_exact_dens_time_derivative:wrong case" << std::endl;
	exit(1);
	break;
      }
      
    };
};

  
void StokesNonLinearSolverDriver::Setelm	(const I 		jelm_,
						 cst_pR 		x_,
						 cst_pR 		xi_,
						 cst_pR 		xii_)
{
  I cellToEdges[16];
    
  const I itime 		= TimeReadOnly_get_itime(&this->m_timeInfo);

  ns_mesh_cooelm		((ns_mesh*)this->m_mesh_usrptr,
				 &jelm_,
				 this->m_workelm.cooelm);
    
  ns_mesh_get_cellToEdges	((ns_mesh*)this->m_mesh_usrptr,
				 &jelm_,
				 cellToEdges);
    
  ns_mesh_get_trelm		((ns_mesh*)this->m_mesh_usrptr,
				 &jelm_,
				 this->m_workelm.belm);
    
  this->m_workelm.alpha2 		= this->m_jacelm[jelm_];
  this->m_workelm.ialpha2 		= regal1/this->m_workelm.alpha2;
  this->m_workelm.edgelm[0]     	= this->m_workelm.cooelm[1]-this->m_workelm.cooelm[0];
  this->m_workelm.edgelm[1]     	= this->m_workelm.cooelm[4]-this->m_workelm.cooelm[3];
  this->m_workelm.edgelm[2]     	= this->m_workelm.cooelm[2]-this->m_workelm.cooelm[1];
  this->m_workelm.edgelm[3]     	= this->m_workelm.cooelm[5]-this->m_workelm.cooelm[4];
  this->m_workelm.edgelm[4]     	= this->m_workelm.cooelm[0]-this->m_workelm.cooelm[2];
  this->m_workelm.edgelm[5]     	= this->m_workelm.cooelm[3]-this->m_workelm.cooelm[5];
    
  this->m_workelm.sbelm[3]   	= this->m_workelm.belm[3]*this->m_workelm.alpha2;
  this->m_workelm.sbelm[1]   	= this->m_workelm.belm[1]*this->m_workelm.alpha2;
  this->m_workelm.sbelm[2]   	= this->m_workelm.belm[2]*this->m_workelm.alpha2;
  this->m_workelm.sbelm[0]   	= this->m_workelm.belm[0]*this->m_workelm.alpha2;
    
    
  {
#if 0
    I jedges[3];    
    jedges[0] 		= cncelm[3] - nvertex;
    jedges[1] 		= cncelm[4] - nvertex;
    jedges[2] 		= cncelm[5] - nvertex;
#endif
    this->m_workelm.normalelmx[0] 	= this->m_normaledge[_dim*cellToEdges[0]+0];
    this->m_workelm.normalelmy[0] 	= this->m_normaledge[_dim*cellToEdges[0]+1];
    this->m_workelm.normalelmx[1] 	= this->m_normaledge[_dim*cellToEdges[1]+0];
    this->m_workelm.normalelmy[1] 	= this->m_normaledge[_dim*cellToEdges[1]+1];
    this->m_workelm.normalelmx[2] 	= this->m_normaledge[_dim*cellToEdges[2]+0];
    this->m_workelm.normalelmy[2] 	= this->m_normaledge[_dim*cellToEdges[2]+1];
      
#if 0
#ifndef NDEBUG
    printf("-> "ifmt","ifmt","ifmt","ifmt"\n",this->m_nedge,cellToEdges[0],cellToEdges[1],cellToEdges[2]);
    { I i;
      for (i=0;i<3;++i)
	{
	  if (cellToEdges[i] >= this->m_nedge)
	    {
	      printf("problem cellToEdges["ifmt"] >= "ifmt"\n",
		     cellToEdges[i],
		     this->m_nedge);
	    }
	} }
#endif
#endif
      
      
    this->m_workelm.alpha3[0]     	= this->m_jacedge[cellToEdges[0]];
    this->m_workelm.alpha3[1]     	= this->m_jacedge[cellToEdges[1]];
    this->m_workelm.alpha3[2]     	= this->m_jacedge[cellToEdges[2]];
  }
    
  this->get_ddlcnc(jelm_,
		   this->m_workelm.locnumer);
    

  cpy_mem_indirect(_total_nddlelm,x_,
		   this->m_workelm.locnumer,
		   this->m_workelm.ddlelm);
    
  const I ntime = this->Params()->GetInfoInteger(InfoInteger::ntime);

  if (ntime>1)
    {
      cpy_mem_indirect	(_total_nddlelm,xi_,this->m_workelm.locnumer,this->m_workelm.ddlelmi);
    }
    
  if ( (ntime>1) AND (itime>0) )					
    {
      cpy_mem_indirect	(_total_nddlelm,xii_,this->m_workelm.locnumer,this->m_workelm.ddlelmii);
    }
    
};	

  
void 	StokesNonLinearSolverDriver::get_ddlcnc_u		(const I 		ielm_,
								 pI 			ddl_)
{
  nsSPACE_cncelm(this->m_space_u,
		 &ielm_,
		 ddl_);
};
  
void 	StokesNonLinearSolverDriver::get_ddlcnc_v		(const I 		ielm_,
								 pI 			ddl_)
{
  nsSPACE_cncelm(this->m_space_u,
		 &ielm_,
		 ddl_);

  for (I i=0;i<_nu;++i)
    {
      ddl_[i] += this->m_space_u->nddl;
    }
};

void 	StokesNonLinearSolverDriver::get_ddlcnc_p		(const I 		ielm_,
								 pI 			ddl_)
{
  nsSPACE_cncelm	(this->m_space_p,
			 &ielm_,
			 ddl_);  
  for (I i=0;i<_np;++i)
    {
      ddl_[i] += _dim*this->m_space_u->nddl;
    }
};
  
  
void StokesNonLinearSolverDriver::get_ddlcnc(const I 		ielm_,
					     pI 			ddl_)
{
  this->get_ddlcnc_u	(ielm_,ddl_ + _ju);
  this->get_ddlcnc_v	(ielm_,ddl_ + _jv);
  this->get_ddlcnc_p	(ielm_,ddl_ + _jp);
};


void StokesNonLinearSolverDriver::Linsys_clear()
{
  const bool pressure_uncoupled		= this->Params()->GetInfoLogical(InfoLogical::pressure_uncoupled);
  const bool pspg				= this->Params()->GetInfoLogical(InfoLogical::pspg);

  SparseBlock_clear(&this->m_sparseStokes->block_UU);  
  if (!pressure_uncoupled)
    {
      if (pspg)
	{
#ifndef NDEBUG
	  std::cout <<"debug:clear sparse_A_block_BB" << std::endl;
	  std::cout <<"debug:clear sparse_A_block_BU" << std::endl;
#endif
	  SparseBlock_clear(&this->m_sparseStokes->block_BB);
	  SparseBlock_clear(&this->m_sparseStokes->block_BU);
	}        
    }
  else
    {
      /* 
	 uncoupled pressure
	 we do nothing
      */
    }
};


void StokesNonLinearSolverDriver::Build_system_stokes_quadrature_free(cst_pR 	x_,
								      cst_pR 	xi_,
								      cst_pR 	xii_)
{ 
  pFem const	fem = &this->m_fem;
    
  pR J_UU = Fem_getJ_UU(fem);
  pR J_VV = Fem_getJ_VV(fem);
  pR J_UV = Fem_getJ_UV(fem);
  pR J_VU = Fem_getJ_VU(fem);

#if 0
#if ( (_ndg==1) AND (_nu==6) )
  cst_pR gravity_term = L02;
#endif
#if ( (_ndg==3) AND (_nu==6) )
  cst_pR gravity_term = L12;
#endif
#if ( (_ndg==6) AND (_nu==6) )
  cst_pR gravity_term = L22;
#endif
#if ( (_ndg==10) AND (_nu==6) )
  cst_pR gravity_term = L32;
#endif
#if ( (_ndg==1) AND (_nu==7) )
  cst_pR gravity_term = L02b;
#endif
#if ( (_ndg==3) AND (_nu==7) )
  cst_pR gravity_term = L12b;
#endif
#if ( (_ndg==6) AND (_nu==7) )
  cst_pR gravity_term = L22b;
#endif
#if ( (_ndg==10) AND (_nu==7) )
  cst_pR gravity_term = L32b;
#endif
#endif

  const bool 	hasTimeAdaptivity	= this->Params()->GetInfoLogical(InfoLogical::time_adaptivity);
  const R 	coeff_ratio_density 	= this->m_solverConstants[SolverConstants::coeff_ratio_density];
  const R 	coeff_ratio_viscosity 	= this->m_solverConstants[SolverConstants::coeff_ratio_viscosity];
  const R 	ireynold 		= this->m_solverConstants[SolverConstants::ireynold];
#if 0	  
  const bool color = this->Params()->GetInfoLogical(InfoLogical::color);  
#endif  
  /**OR (Z->newton_iter>0)*/
  if ( false == hasTimeAdaptivity )
    {
      pWorkelm elm = &this->m_workelm;
      
      for (I jelm=0;jelm<this->m_nelm;++jelm)
	{ 
	  clr_mem(total_nddlelmxtotal_nddlelm,
		  this->m_workelm.locmat);
	    
	  this->Setelm(jelm,
		       x_,
		       xi_,
		       xii_);
	    
	  Fem_setelm	(fem,
			 elm);
	    
#if 0
	  if (this->transport_usrptr)
	    {
	      transport_get_viscelm(this->transport_usrptr,jelm,this->m_workelm.viscelm);
	      transport_get_denselm(this->transport_usrptr,jelm,this->m_workelm.denselm);
	    }
	  else
	    {
	      this->m_workelm.viscelm[0]=((R)1.0)/nsSQRT(2.0);
	      for (i=1;i<_ndg;++i)
		this->m_workelm.viscelm[i]=((R)0.0);
	      this->m_workelm.denselm[0]=((R)1.0)/nsSQRT(2.0);
	      for (i=1;i<_ndg;++i)
		this->m_workelm.denselm[i]=((R)0.0);
	    }

	  if (color)							
	    {									
	      { I k;							
		for (k=0;k<_nf;++k)						
		  {
		    this->m_workelm.denselm[k] = changement_variable2(coeff_ratio_density,felm_[k]);	
		    this->m_workelm.viscelm[k] = changement_variable2(coeff_ratio_viscosity,felm_[k]);	
		  } }
	    }								
	  else									
	    {			
	      cpy_mem(_nf,felm_,this->m_workelm.viscelm);
	      vheaviside(iproc_,_nf,this->m_workelm.viscelm);
	      { I k;							
		for (k=0;k<_nf;++k)						
		  {								
		    const R aa = this->m_workelm.viscelm[k];
		    this->m_workelm.viscelm[k] = changement_variable2(coeff_ratio_viscosity,aa);	
		    this->m_workelm.denselm[k] = changement_variable2(coeff_ratio_density,aa);	
		  } }								
	    }
#endif

	  for (I k=0;k<_ndg;++k)						
	    {
	      elm->denselm[k] = 1.0;
	    }
	  
	  for (I k=0;k<_ndg;++k)						
	    {
	      elm->viscelm[k] = 1.0;
	    } 

#if 0
	  nsDG_transport_get_nodalvalues(self_->transport_usrptr,&jelm,elm->denselm);
	  nsDG_transport_get_nodalvalues(self_->transport_usrptr,&jelm,elm->viscelm);
#endif
	  

	  for (I k=0;k<_ndg;++k)						
	    {
	      R hh = elm->viscelm[k];
	      //
	      // Take care of overshoots.
	      //
	      if (hh<0.0)
		{
		  hh=0.0;
		}
	      if (hh>1.0)
		{
		  hh=1.0;
		}
	      elm->viscelm[k] 	= 1.0+hh*(coeff_ratio_viscosity-1.0);
	      elm->denselm[k] 	= 1.0+hh*(coeff_ratio_density-1.0);
#if 0
	      elm->viscelm[k] 	= coeff_ratio_viscosity+hh*(1.0-coeff_ratio_viscosity);
	      elm->denselm[k] 	= coeff_ratio_density+hh*(1.0-coeff_ratio_density);
#endif
	    } 
	  
	  
	  clr_mem(_nu*_nu,J_VU);
	  clr_mem(_nu*_nu,J_UV);
	  clr_mem(_nu*_nu,J_VV);
	  clr_mem(_nu*_nu,J_UU);


	  /* CONSTRUCTION DES MATRICES */
	  /* jacobien exact terme de stress -/- u */
	    
	  nsblas_dgemv(transN,&nuxnu,&ndg,&regal1,fem->exact_visc_dxu_dxu,&nuxnu,elm->viscelm,&negal1,&regal0,fem->visc_dxu_dxu,&negal1);
	  nsblas_dgemv(transN,&nuxnu,&ndg,&regal1,fem->exact_visc_dxu_dyu,&nuxnu,elm->viscelm,&negal1,&regal0,fem->visc_dxu_dyu,&negal1);
	  nsblas_dgemv(transN,&nuxnu,&ndg,&regal1,fem->exact_visc_dyu_dxu,&nuxnu,elm->viscelm,&negal1,&regal0,fem->visc_dyu_dxu,&negal1);
	  nsblas_dgemv(transN,&nuxnu,&ndg,&regal1,fem->exact_visc_dyu_dyu,&nuxnu,elm->viscelm,&negal1,&regal0,fem->visc_dyu_dyu,&negal1);
	    
	  { const R a = ireynold * regal2;
	    nsblas_daxpy(&nuxnu,
			 &a,
			 fem->visc_dxu_dxu,
			 &negal1,
			 J_UU,
			 &negal1); }
	  
	  { const R a = ireynold;
	    nsblas_daxpy(&nuxnu,&a,fem->visc_dyu_dyu,&negal1,J_UU,&negal1); }
	  { const R a = ireynold;
	    nsblas_daxpy(&nuxnu,&a,fem->visc_dxu_dyu,&negal1,J_UV,&negal1); }
	  { const R a = ireynold * regal2;
	    nsblas_daxpy(&nuxnu,&a,fem->visc_dyu_dyu,&negal1,J_VV,&negal1); }
	  { const R a = ireynold;
	    nsblas_daxpy(&nuxnu,&a,fem->visc_dxu_dxu,&negal1,J_VV,&negal1); }
	  { const R a = ireynold;
	    nsblas_daxpy(&nuxnu,&a,fem->visc_dyu_dxu,&negal1,J_VU,&negal1); }  

#if 0
	  ns_jacobian_exact_dens_time_derivative(iproc_,jelm);
	  nsblas_dgemv(transN,&nuxnuxnu,&ndg,&regal1,fem->exact_dens_u_dxu_u,&nuxnuxnu,elm->denselm,&negal1,&regal0,fem->dens_u_dxu_u,&negal1);
	  nsblas_dgemv(transN,&nuxnuxnu,&ndg,&regal1,fem->exact_dens_u_dyu_u,&nuxnuxnu,elm->denselm,&negal1,&regal0,fem->dens_u_dyu_u,&negal1);
	  nsblas_dgemv(transN,&nuxnuxnu,&ndg,&regal1,fem->exact_dens_dxu_u_u,&nuxnuxnu,elm->denselm,&negal1,&regal0,fem->dens_dxu_u_u,&negal1);
	  nsblas_dgemv(transN,&nuxnuxnu,&ndg,&regal1,fem->exact_dens_dyu_u_u,&nuxnuxnu,elm->denselm,&negal1,&regal0,fem->dens_dyu_u_u,&negal1);
	  { const R a = ireynold * regal2;
	    nsblas_daxpy(&nuxnu,&a,fem->visc_dxu_dxu,&negal1,J_UU,&negal1); }
	  { const R a = ireynold;
	    nsblas_daxpy(&nuxnu,&a,fem->visc_dyu_dyu,&negal1,J_UU,&negal1); }
	  { const R a = ireynold;
	    nsblas_daxpy(&nuxnu,&a,fem->visc_dxu_dyu,&negal1,J_UV,&negal1); }
	  { const R a = ireynold * regal2;
	    nsblas_daxpy(&nuxnu,&a,fem->visc_dyu_dyu,&negal1,J_VV,&negal1); }
	  { const R a = ireynold;
	    nsblas_daxpy(&nuxnu,&a,fem->visc_dxu_dxu,&negal1,J_VV,&negal1); }
	  { const R a = ireynold;
	    nsblas_daxpy(&nuxnu,&a,fem->visc_dyu_dxu,&negal1,J_VU,&negal1); }  
#endif

	  /* jacobien acceleration */
#if 0
	  nsblas_dgemv(transN,&nuxnu,&nu,&regal1,fem->dens_u_dxu_u,&nuxnu,&elm->ddlelm[_ju],&negal1,&regal1,J_UU,&negal1);
	  nsblas_dgemv(transN,&nuxnu,&nu,&regal1,fem->dens_u_dyu_u,&nuxnu,&elm->ddlelm[_jv],&negal1,&regal1,J_UV,&negal1);
	  nsblas_dgemv(transN,&nuxnu,&nu,&regal1,fem->dens_u_dxu_u,&nuxnu,&elm->ddlelm[_ju],&negal1,&regal1,J_VU,&negal1);
	  nsblas_dgemv(transN,&nuxnu,&nu,&regal1,fem->dens_u_dyu_u,&nuxnu,&elm->ddlelm[_jv],&negal1,&regal1,J_VV,&negal1);
	  nsblas_dgemv(transN,&nuxnu,&nu,&regal1,fem->dens_dxu_u_u,&nuxnu,&elm->ddlelm[_ju],&negal1,&regal1,J_UU,&negal1);
	  nsblas_dgemv(transN,&nuxnu,&nu,&regal1,fem->dens_dyu_u_u,&nuxnu,&elm->ddlelm[_jv],&negal1,&regal1,J_VV,&negal1);
#endif
	  /**/
	  /* jacobien terme de stress -/- visc  facultatif */
#if 0
	  
	  const R coeff_ratio_viscosity = this->m_solverConstants[SolverConstants::coeff_ratio_viscosity];

	  nsblas_dgemv(transN,&nuxndg,&nu,&regal1,fem->exact_dxu_visc_dxu,&nuxndg,uelm,&negal1,&regal0,dxu_visc_dxu,&negal1);
	  nsblas_dgemv(transN,&nuxndg,&nu,&regal1,fem->exact_dyu_visc_dyu,&nuxndg,uelm,&negal1,&regal0,dyu_visc_dyu,&negal1);
	  nsblas_dgemv(transN,&nuxndg,&nu,&regal1,fem->exact_dxu_visc_dyu,&nuxndg,velm,&negal1,&regal0,dxu_visc_dyu,&negal1);
	  { const R a = ireynold * regal2* (-coeff_ratio_viscosity);
	    nsblas_daxpy(&nuxndg,&a,dxu_visc_dxu,&negal1,J_UF,&negal1); }
	  { const R a = ireynold* (-coeff_ratio_viscosity);
	    nsblas_daxpy(&nuxndg,&a,dyu_visc_dyu,&negal1,J_UF,&negal1); }
	  { const R a = ireynold* (-coeff_ratio_viscosity);
	    nsblas_daxpy(&nuxndg,&a,dxu_visc_dyu,&negal1,J_UF,&negal1); }      
	  nsblas_dgemv(transN,&nuxndg,&nu,&regal1,fem->exact_dxu_visc_dxu,&nuxndg,velm,&negal1,&regal0,dxu_visc_dxu,&negal1);
	  nsblas_dgemv(transN,&nuxndg,&nu,&regal1,fem->exact_dyu_visc_dyu,&nuxndg,velm,&negal1,&regal0,dyu_visc_dyu,&negal1);
	  nsblas_dgemv(transN,&nuxndg,&nu,&regal1,fem->exact_dyu_visc_dxu,&nuxndg,uelm,&negal1,&regal0,dyu_visc_dxu,&negal1);
	  { const R a = ireynold * regal2* (-coeff_ratio_viscosity);
	    nsblas_daxpy(&nuxndg,&a,dyu_visc_dyu,&negal1,J_VF,&negal1); }
	  { const R a = ireynold* (-coeff_ratio_viscosity);
	    nsblas_daxpy(&nuxndg,&a,dxu_visc_dxu,&negal1,J_VF,&negal1); }
	  { const R a = ireynold* (-coeff_ratio_viscosity);
	    nsblas_daxpy(&nuxndg,&a,dyu_visc_dxu,&negal1,J_VF,&negal1); }
#endif
	  /* 
	     jacobien gravite -/- denselm
	  */
#if 0
	  const R ifroude = this->m_solverConstants[SolverConstants::ifroude];
	  const R coeff_ratio_density = this->m_solverConstants[SolverConstants::coeff_ratio_density];

	  { const R a = ifroude*self_->jacelm[jelm] * (-coeff_ratio_density);
	    nsblas_daxpy(&nuxndg,&a,gravity_term,&negal1,J_VF,&negal1); }
#endif

	  for (I j=0;j<_nu;++j)
	    {
	      for (I i=0;i<_nu;++i)
		{
		  elm->locmat[(offu+j)*off+(offu+i)] += J_UU[j*_nu+i];
		  elm->locmat[(offv+j)*off+(offu+i)] += J_UV[j*_nu+i];
		  elm->locmat[(offu+j)*off+(offv+i)] += J_VU[j*_nu+i];
		  elm->locmat[(offv+j)*off+(offv+i)] += J_VV[j*_nu+i];
		}
	    }

	  this->assmatelm_full(jelm); 
	} 
    }
  else
    {
      pWorkelm elm = &this->m_workelm;
      
      for (I jelm=0;jelm<this->m_nelm;++jelm)
	{      
	  clr_mem(total_nddlelmxtotal_nddlelm,
		  elm->locmat);
	  
	  this->Setelm(jelm,
		       x_,
		       xi_,
		       xii_);
	  
	  Fem_setelm(fem,
		     elm);
	  
	  compute_viscosity_density_nodalelm2(elm->ddlelm);

	  clr_mem(_nu*_nu,J_VU);
	  clr_mem(_nu*_nu,J_UV);
	  clr_mem(_nu*_nu,J_VV);
	  clr_mem(_nu*_nu,J_UU);

	  for (I j=0;j<_nu;++j)
	    {
	      for (I i=0;i<_nu;++i)
		{
		  elm->locmat[(offu+j)*off+(offu+i)] += J_UU[j*_nu+i];
		  elm->locmat[(offv+j)*off+(offu+i)] += J_UV[j*_nu+i];
		  elm->locmat[(offu+j)*off+(offv+i)] += J_VU[j*_nu+i];
		  elm->locmat[(offv+j)*off+(offv+i)] += J_VV[j*_nu+i];
		}
	    }
	  
	  this->assmatelm_full(jelm);
	} 
    }
};


void StokesNonLinearSolverDriver::assmatelm_full(const I 	jelm_)
{ 
  cst_pI Ab = Sparse_get_b(this->m_sparseStokes->A);
  cst_pI Ai = Sparse_get_i(this->m_sparseStokes->A);
  pR Ax = Sparse_get_ownx(this->m_sparseStokes->A);

  const bool hasPspg = this->Params()->GetInfoLogical(InfoLogical::pspg);  
  pWorkelm elm = &this->m_workelm;    
  this->get_ddlcnc(jelm_,
		   elm->locnumer);
  
  I N = 0;  
  for (I i=0;i<_total_nddlelm;++i)
    {
      if (0 == this->blank[elm->locnumer[i]])
	{
	  this->blank[elm->locnumer[i]] 	= ++N;
	  this->iwork[N-1] 	  		= elm->locnumer[i];
	}
    } 

  //
  // loop over local velocity.
  //
  for (I j=0;j<_total_nddlelm-_np;++j)
    {
      I jddl 	= elm->locnumer[j];
      const I n		= Ab[jddl+1]-Ab[jddl];
      cst_pI it 		= &Ai[Ab[jddl]-1];
      pR rt 			= &Ax[Ab[jddl]-1];
      for (I i=0;i<n;++i)
	{
	  I k;
	  if ( (k=this->blank[it[i]-1]) )
	    {
	      if (hasPspg)
		{
		  rt[i] += elm->locmat[j + (k-1)*off];
		}
	      else
		{
		  if (k-1<offp)
		    {
		      rt[i] += elm->locmat[j + (k-1)*off];
		    }
		}
#if 0
	      if (
		  ( (active_save_B_computation) AND (NOT save_B_computation) )
		  OR 
		  (NOT active_save_B_computation) )
		{
		  rt[i] += elm->locmat[j + (k-1)*off];
		}
	      else
		{
		  if (k-1<offp)
		    {
		      rt[i] += elm->locmat[j + (k-1)*off];
		    }
		}
#endif
		
	    }
	} 
    } 
  //
  // loop over the pressure
  //
  for (I j=_nuxdim;j<_total_nddlelm;++j)
    {
      I jddl 	= elm->locnumer[j];
      
      const I n		= Ab[jddl+1]-Ab[jddl];
      cst_pI it 		= &Ai[Ab[jddl]-1];
      pR rt 			= &Ax[Ab[jddl]-1];
      
#if 0
      const I n		= self_->sparse_stokes.A.b[jddl+1]-self_->sparse_stokes.A.b[jddl];
      cst_pI it 		= &self_->sparse_stokes.A.i[self_->sparse_stokes.A.b[jddl]-1];
      pR rt 			= &self_->sparse_stokes.A.own_x[self_->sparse_stokes.A.b[jddl]-1];
#endif
      
      for (I i=0;i<n;++i)
	{
	  I k;
	  if ( (k=this->blank[it[i]-1]) )
	    {
	      if (k-1<offp)
		{	
#if 0
		  if (
		      ( (active_save_B_computation) AND (NOT save_B_computation) )
		      OR 
		      (NOT active_save_B_computation) )
		    {}		      
#endif
		  
		  if (hasPspg)
		    {
		      rt[i] += elm->locmat[j + (k-1)*off];
		    }
		} 
	      else
		{
		  if (hasPspg)
		    {
		      rt[i] += elm->locmat[j + (k-1)*off];
		    }
		}
	    }
	}
    }    
  
  for (I i=0;i<N;++i) 
    {
      this->blank[this->iwork[i]] = 0;
    }
};

void StokesNonLinearSolverDriver::compute_viscosity_density_nodalelm2(cst_pR felm_)
{    

  for (I k=0;k<_ndg;++k)						
    {
      this->m_workelm.denselm[k] = 1.0;
    }

  for (I k=0;k<_ndg;++k)						
    {
      this->m_workelm.viscelm[k] = 1.0;
    }
  
#if 0

  const bool transport = this->Params()->GetInfoLogical(InfoLogical::transport);  
  if (true == transport)
    {
      const bool transport_uncoupled = this->Params()->GetInfoLogical(InfoLogical::transport_uncoupled);  
      if (false == transport_uncoupled)
	{
	  const bool color = this->Params()->GetInfoLogical(InfoLogical::color);  
	  const R coeff_ratio_density = this->m_solverConstants[SolverConstants::coeff_ratio_density];
	  const R coeff_ratio_viscosity = this->m_solverConstants[SolverConstants::coeff_ratio_viscosity];
	  if (true == color)							
	    {
 
	      for (I k=0;k<_nf;++k)						
		{
		  workelm->denselm[k] = changement_variable2(coeff_ratio_density,
							     felm_[k]);
		  workelm->viscelm[k] = changement_variable2(coeff_ratio_viscosity,
							     felm_[k]);
		}
	    }								
	  else									
	    {			
	      cpy_mem(_nf,felm_,workelm->viscelm);
	      vheaviside(self_,_nf,workelm->viscelm);
	      for (I k=0;k<_nf;++k)
		{
		  const R aa = workelm->viscelm[k];
		  workelm->viscelm[k] = changement_variable2(coeff_ratio_viscosity,
							     aa);
		  workelm->denselm[k] = changement_variable2(coeff_ratio_density,
							     aa);
		}
	    }
	}
      else
	{
	  
	}
    }
  else
    {

    }
  
#endif

};
