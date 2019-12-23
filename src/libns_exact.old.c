#include "ns_sys.h"
#include "ns_config.h"
#include "ns_constantes.h"
#include "nsGLOBAL.h"
#include "ns_extern.h"
#include "ns_automatic_constantes.h"

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



#if ( ( _ndg == 3 ) AND ( _nu == 6 ) AND ( _np == 3 ) )
#include "libns_exact_f3_u6_p3.c"
#else
#if ( ( _ndg == 6 ) AND ( _nu == 6 ) AND ( _np == 3 ) )
#include "libns_exact_f6_u6_p3.c"
#else
#error wrong mf wrong nu wrong np
#endif
#endif


static const I 	build_partialoff 	= ((I)_build_partialoff);
static const I 	build_twice_partialoff 	= ((I)_build_twice_partialoff);

void nsFEM_DATA_setelm_quadrature(nsFEM_DATA_ST*const 	fem_,
				  nsWORKELM_ST*const 	elm_)
{
}


Err nsFEM_DATA_def_from_quadrature(nsFEM_DATA_ST*const 	fem_,
				       cst_pI			qn_,
				       cst_pR 			qw_,
				       cst_pR 			qp_,
				       cst_pI			qoff_)
{
  Err err = __eErr_no;
#if 0
  extern void basis_	(cst_pI ,cst_pI ,cst_pR ,cst_pI,pR,cst_pI);
  extern void drbasis_	(cst_pI ,cst_pI ,cst_pR ,cst_pI,pR,cst_pI);
  extern void dsbasis_	(cst_pI ,cst_pI ,cst_pR ,cst_pI,pR,cst_pI);
  extern void obasis_	(cst_pI ,cst_pI ,cst_pR ,cst_pI,pR,cst_pI);
  extern void drobasis_	(cst_pI ,cst_pI ,cst_pR ,cst_pI,pR,cst_pI);
  extern void dsobasis_	(cst_pI ,cst_pI ,cst_pR ,cst_pI,pR,cst_pI);
#endif
  cst_eElement element = fem_->elm; 
  const I negal3=3;  
  
#define define_basis(_a,_v) ensBASIS_b(_a,element,"T","T",qn_,qp_,qoff_,fem_->shape##_v,&n##_v); \
  ensBASIS_dr(_a,element,"T","T",qn_,qp_,qoff_,fem_->drshape##_v,&n##_v); \
  ensBASIS_ds(_a,element,"T","T",qn_,qp_,qoff_,fem_->dsshape##_v,&n##_v)  

  ensBASIS_b(__ensBASIS_LAGRANGE_1,element,"T","T",qn_,qp_,qoff_,fem_->shaper,&negal3);
  define_basis(_shape_u,u);
  define_basis(_shape_p,p);
  define_basis(_shape_dg,dg);
  
#if 0
  basis_(&negal3,qn_,qp_,qoff_,fem_->shaper,&negal3);

#define define_basis(_a,_v) basis_(_a,qn_,qp_,qoff_,fem_->shape##_v,&n##_v); \
  drbasis_(_a,qn_,qp_,qoff_,fem_->drshape##_v,&n##_v);			\
  dsbasis_(_a,qn_,qp_,qoff_,fem_->dsshape##_v,&n##_v)

#define define_obasis(_a,_v) obasis_(_a,qn_,qp_,qoff_,fem_->shape##_v,&n##_v); \
  drobasis_(_a,qn_,qp_,qoff_,fem_->drshape##_v,&n##_v);			\
  dsobasis_(_a,qn_,qp_,qoff_,fem_->dsshape##_v,&n##_v)

  define_basis	(&nu,u);
  define_basis	(&np,p);
  define_obasis	(&ndg,dg);
#endif
  { I i,j,k;
    for (k=0;k<3;++k)
      {
	for (j=0;j<_ndg;++j)
	  {
	    for (i=0;i<qn_[0];++i)
	      {
		const R a = qw_[i] * fem_->shapedg[i*_ndg+j] * fem_->shaper[i*3+k];
		dger(&nu,&nu,&a,&fem_->drshapeu[i*_nu],&negal1,&fem_->drshapeu[i*_nu],&negal1,&fem_->qref_radius_visc_dru_dru[_ndg*nuxnu*k+j*_nuxnu],&nu);
		dger(&nu,&nu,&a,&fem_->dsshapeu[i*_nu],&negal1,&fem_->dsshapeu[i*_nu],&negal1,&fem_->qref_radius_visc_dsu_dsu[_ndg*nuxnu*k+j*_nuxnu],&nu); 	
		dger(&nu,&nu,&a,&fem_->dsshapeu[i*_nu],&negal1,&fem_->drshapeu[i*_nu],&negal1,&fem_->qref_radius_visc_dru_dsu[_ndg*nuxnu*k+j*_nuxnu],&nu); 	
		dger(&nu,&nu,&a,&fem_->drshapeu[i*_nu],&negal1,&fem_->dsshapeu[i*_nu],&negal1,&fem_->qref_radius_visc_dsu_dru[_ndg*nuxnu*k+j*_nuxnu],&nu); 	  
	      }
	  }        	
	for (j=0;j<_ndg;++j)
	  {
	    for (i=0;i<qn_[0];++i)
	      {
		const R a = qw_[i] * fem_->shapedg[i*_ndg+j] *fem_->shaper[i*3+k];
		dger(&nu,&nu,&a,&fem_->shapeu[i*_nu],&negal1,&fem_->shapeu[i*_nu],&negal1,&fem_->qref_radius_dens_u_u[_ndg*nuxnu*k+j*_nuxnu],&nu);
	      }
	  }	
	for (i=0;i<qn_[0];++i)
	  {
	    const R a = qw_[i] * fem_->shaper[i*3+k];
	    dger(&np,&nu,&a,&fem_->shapep[i*_np],&negal1,&fem_->drshapeu[i*_nu],&negal1,&fem_->qref_radius_dru_p[k*_nu*_np],&np);
	    dger(&np,&nu,&a,&fem_->shapep[i*_np],&negal1,&fem_->dsshapeu[i*_nu],&negal1,&fem_->qref_radius_dsu_p[k*_nu*_np],&np); 	
	    dger(&nu,&np,&a,&fem_->drshapeu[i*_nu],&negal1,&fem_->shapep[i*_np],&negal1,&fem_->qref_radius_p_dru[k*_nu*_np],&nu);
	    dger(&nu,&np,&a,&fem_->dsshapeu[i*_nu],&negal1,&fem_->shapep[i*_np],&negal1,&fem_->qref_radius_p_dsu[k*_nu*_np],&nu); 	
	  }  
      } }
  return err;
}


Err nsFEM_DATA_def(nsFEM_DATA_ST*const fem_,cst_eElement element_  )
{
  Err err = __eErr_no;
  I dec = 0;
  I n;
  fem_->elm = element_;
  /* L ORDRE EST IMPORTANT */
  fem_->exact_dyu_p			= &fem_->build_partialr[_build_partialoff+0];
  fem_->exact_dens_dyu_u_u 		= &fem_->build_partialr[_build_partialoff+_nu*_np];
  fem_->exact_dens_u_dyu_u 		= &fem_->build_partialr[_build_partialoff+_nu*_np+_ndg*_nu*_nu*_nu];
  fem_->exact_u_dyu_dens_u 		= &fem_->build_partialr[_build_partialoff+_nu*_np+_ndg*_nu*_nu*_nu*2];
  fem_->exact_u_dyf_f 			= &fem_->build_partialr[_build_partialoff+_nu*_np+_ndg*_nu*_nu*_nu*3];
  fem_->exact_dyf_u_f 			= &fem_->build_partialr[_build_partialoff+_nu*_np+_ndg*_nu*_nu*_nu*3+_nu*_ndg*_ndg];
  fem_->exact_u_f_dyf 			= &fem_->build_partialr[_build_partialoff+_nu*_np+_ndg*_nu*_nu*_nu*3+_nu*_ndg*_ndg*2];
  fem_->exact_f_u_dyf 			= &fem_->build_partialr[_build_partialoff+_nu*_np+_ndg*_nu*_nu*_nu*3+_nu*_ndg*_ndg*3];
  fem_->exact_dxu_p			= &fem_->build_partialr[0];
  fem_->exact_dens_dxu_u_u 		= &fem_->build_partialr[_nu*_np];
  fem_->exact_dens_u_dxu_u 		= &fem_->build_partialr[_nu*_np+_ndg*_nu*_nu*_nu];
  fem_->exact_u_dxu_dens_u 		= &fem_->build_partialr[_nu*_np+_ndg*_nu*_nu*_nu*2];
  fem_->exact_u_dxf_f 			= &fem_->build_partialr[_nu*_np+_ndg*_nu*_nu*_nu*3];
  fem_->exact_dxf_u_f 			= &fem_->build_partialr[_nu*_np+_ndg*_nu*_nu*_nu*3+_nu*_ndg*_ndg];
  fem_->exact_u_f_dxf 			= &fem_->build_partialr[_nu*_np+_ndg*_nu*_nu*_nu*3+_nu*_ndg*_ndg*2];
  fem_->exact_f_u_dxf 			= &fem_->build_partialr[_nu*_np+_ndg*_nu*_nu*_nu*3+_nu*_ndg*_ndg*3];
  fem_->exact_visc_dxu_dxu 		= &fem_->build_twice_partialr[0];
  fem_->exact_dxu_visc_dxu 		= &fem_->build_twice_partialr[_nu*_nu*_ndg];
  fem_->exact_u_dxf_dxf 		= &fem_->build_twice_partialr[_nu*_nu*_ndg+_nu*_nu*_ndg];
  fem_->exact_dxf_u_dxf 		= &fem_->build_twice_partialr[_nu*_nu*_ndg*2 + _nu*_ndg*_ndg];
  fem_->exact_visc_dxu_dyu 		= &fem_->build_twice_partialr[_build_twice_partialoff];
  fem_->exact_dxu_visc_dyu 		= &fem_->build_twice_partialr[_build_twice_partialoff*1 + _nu*_nu*_ndg];
  fem_->exact_u_dxf_dyf 		= &fem_->build_twice_partialr[_build_twice_partialoff+_nu*_nu*_ndg*2];
  fem_->exact_dxf_u_dyf 		= &fem_->build_twice_partialr[_build_twice_partialoff*1+_nu*_nu*_ndg*2 + _nu*_ndg*_ndg];
  fem_->exact_visc_dyu_dxu 		= &fem_->build_twice_partialr[_build_twice_partialoff*2];
  fem_->exact_dyu_visc_dxu 		= &fem_->build_twice_partialr[_build_twice_partialoff*2 + _nu*_nu*_ndg];
  fem_->exact_u_dyf_dxf 		= &fem_->build_twice_partialr[_build_twice_partialoff*2+_nu*_nu*_ndg*2];
  fem_->exact_dyf_u_dxf 		= &fem_->build_twice_partialr[_build_twice_partialoff*2+_nu*_nu*_ndg*2 + _nu*_ndg*_ndg];
  fem_->exact_visc_dyu_dyu 		= &fem_->build_twice_partialr[_build_twice_partialoff*3];
  fem_->exact_dyu_visc_dyu 		= &fem_->build_twice_partialr[_build_twice_partialoff*3 + _nu*_nu*_ndg];
  fem_->exact_u_dyf_dyf 		= &fem_->build_twice_partialr[_build_twice_partialoff*3+_nu*_nu*_ndg*2];
  fem_->exact_dyf_u_dyf 		= &fem_->build_twice_partialr[_build_twice_partialoff*3+_nu*_nu*_ndg*2 + _nu*_ndg*_ndg];
  fem_->exact_u_dens_u 			= exactref_u_dens_u;
  fem_->exact_dens_u_u 			= exactref_dens_u_u;
  dec = 0;
  n = _nu*_np;memcpy(&fem_->build_partialx[dec],exactref_dru_p,sizeof(R)*n);dec+=n;
  n = _nu*_nu*_nu*_ndg;memcpy(&fem_->build_partialx[dec],exactref_dens_dru_u_u,sizeof(R)*n);dec+=n;
  n = _nu*_nu*_nu*_ndg;memcpy(&fem_->build_partialx[dec],exactref_dens_u_dru_u,sizeof(R)*n);dec+=n;
  n = _nu*_nu*_nu*_ndg;memcpy(&fem_->build_partialx[dec],exactref_u_dru_dens_u,sizeof(R)*n);dec+=n;
  n = _nu*_ndg*_ndg;memcpy(&fem_->build_partialx[dec],exactref_u_drf_f,sizeof(R)*n);dec+=n;
  n = _nu*_ndg*_ndg;memcpy(&fem_->build_partialx[dec],exactref_drf_u_f,sizeof(R)*n);dec+=n;
  n = _nu*_ndg*_ndg;memcpy(&fem_->build_partialx[dec],exactref_u_f_drf,sizeof(R)*n);dec+=n;
  n = _nu*_ndg*_ndg;memcpy(&fem_->build_partialx[dec],exactref_f_u_drf,sizeof(R)*n);dec+=n;
  dec = 0;
  n = _nu*_np;memcpy(&fem_->build_partialx[build_partialoff+dec],exactref_dsu_p,sizeof(R)*n);dec+=n;
  n = _nu*_nu*_nu*_ndg;memcpy(&fem_->build_partialx[build_partialoff+dec],exactref_dens_dsu_u_u,sizeof(R)*n);dec+=n;
  n = _nu*_nu*_nu*_ndg;memcpy(&fem_->build_partialx[build_partialoff+dec],exactref_dens_u_dsu_u,sizeof(R)*n);dec+=n;
  n = _nu*_nu*_nu*_ndg;memcpy(&fem_->build_partialx[build_partialoff+dec],exactref_u_dsu_dens_u,sizeof(R)*n);dec+=n;
  n = _nu*_ndg*_ndg;memcpy(&fem_->build_partialx[build_partialoff+dec],exactref_u_dsf_f,sizeof(R)*n);dec+=n;
  n = _nu*_ndg*_ndg;memcpy(&fem_->build_partialx[build_partialoff+dec],exactref_dsf_u_f,sizeof(R)*n);dec+=n;
  n = _nu*_ndg*_ndg;memcpy(&fem_->build_partialx[build_partialoff+dec],exactref_u_f_dsf,sizeof(R)*n);dec+=n;
  n = _nu*_ndg*_ndg;memcpy(&fem_->build_partialx[build_partialoff+dec],exactref_f_u_dsf,sizeof(R)*n);dec+=n;
  dec = 0;
  n = _nu*_nu*_ndg;memcpy(&fem_->build_twice_partialx[_build_twice_partialoff*0+dec],exactref_visc_dru_dru,sizeof(R)*n);dec+=n;
  n = _nu*_nu*_ndg;memcpy(&fem_->build_twice_partialx[_build_twice_partialoff*0+dec],exactref_dru_visc_dru,sizeof(R)*n);dec+=n;
  n = _nu*_ndg*_ndg;memcpy(&fem_->build_twice_partialx[_build_twice_partialoff*0+dec],exactref_u_drf_drf,sizeof(R)*n);dec+=n;
  n = _nu*_ndg*_ndg;memcpy(&fem_->build_twice_partialx[_build_twice_partialoff*0+dec],exactref_drf_u_drf,sizeof(R)*n);dec+=n;
  dec = 0;
  n = _nu*_nu*_ndg;memcpy(&fem_->build_twice_partialx[_build_twice_partialoff*1+dec],exactref_visc_dru_dsu,sizeof(R)*n);dec+=n;
  n = _nu*_nu*_ndg;memcpy(&fem_->build_twice_partialx[_build_twice_partialoff*1+dec],exactref_dru_visc_dru,sizeof(R)*n);dec+=n;
  n = _nu*_ndg*_ndg;memcpy(&fem_->build_twice_partialx[_build_twice_partialoff*1+dec],exactref_u_drf_dsf,sizeof(R)*n);dec+=n;
  n = _nu*_ndg*_ndg;memcpy(&fem_->build_twice_partialx[_build_twice_partialoff*1+dec],exactref_drf_u_dsf,sizeof(R)*n);dec+=n;
  dec = 0;
  n = _nu*_nu*_ndg;memcpy(&fem_->build_twice_partialx[_build_twice_partialoff*2+dec],exactref_visc_dsu_dru,sizeof(R)*n);dec+=n;
  n = _nu*_nu*_ndg;memcpy(&fem_->build_twice_partialx[_build_twice_partialoff*2+dec],exactref_dsu_visc_dru,sizeof(R)*n);dec+=n;
  n = _nu*_ndg*_ndg;memcpy(&fem_->build_twice_partialx[_build_twice_partialoff*2+dec],exactref_u_dsf_drf,sizeof(R)*n);dec+=n;
  n = _nu*_ndg*_ndg;memcpy(&fem_->build_twice_partialx[_build_twice_partialoff*2+dec],exactref_dsf_u_drf,sizeof(R)*n);dec+=n;
  dec = 0;
  n = _nu*_nu*_ndg;memcpy(&fem_->build_twice_partialx[_build_twice_partialoff*3+dec],exactref_visc_dsu_dsu,sizeof(R)*n);dec+=n;
  n = _nu*_nu*_ndg;memcpy(&fem_->build_twice_partialx[_build_twice_partialoff*3+dec],exactref_dsu_visc_dsu,sizeof(R)*n);dec+=n;
  n = _nu*_ndg*_ndg;memcpy(&fem_->build_twice_partialx[_build_twice_partialoff*3+dec],exactref_u_dsf_dsf,sizeof(R)*n);dec+=n;
  n = _nu*_ndg*_ndg;memcpy(&fem_->build_twice_partialx[_build_twice_partialoff*3+dec],exactref_dsf_u_dsf,sizeof(R)*n);dec+=n;
  return err;
}


void nsFEM_DATA_setelm(nsFEM_DATA_ST*const 	fem_,
		       nsWORKELM_ST*const 	elm_)
{
  R sbelm2[_dim*_dim*_dim*_dim];  

  /*
    dx = (a0 dr + a2 ds)
    dy = (a1 dr + a3 ds)
    (a0 dr + a2 ds)*(a1 dr + a3 ds) = a0a1 drdr +  a2a3 dsds + a2a1 dsdr  + a0a3 drds 
    (a1 dr + a3 ds)*(a0 dr + a2 ds) = a0a1 drdr +  a2a3 dsds + a0a3 dsdr  + a2a1 drds
  */  
  sbelm2[0] 	= elm_->sbelm[0]*elm_->belm[0];
  sbelm2[1] 	= elm_->sbelm[0]*elm_->belm[2];
  sbelm2[2] 	= elm_->sbelm[0]*elm_->belm[2];
  sbelm2[3] 	= elm_->sbelm[2]*elm_->belm[2];
  
  sbelm2[4+0] 	= elm_->sbelm[0]*elm_->belm[1];
  sbelm2[4+1] 	= elm_->sbelm[0]*elm_->belm[3];
  sbelm2[4+2] 	= elm_->sbelm[1]*elm_->belm[2];
  sbelm2[4+3] 	= elm_->sbelm[2]*elm_->belm[3];
  sbelm2[8+0] 	= elm_->sbelm[0]*elm_->belm[1];
  sbelm2[8+2] 	= elm_->sbelm[0]*elm_->belm[3];/*ATTENTION C EST INVERSE PAR RAPPORT A LA COLONNE PRECEDENTE */
  sbelm2[8+1] 	= elm_->sbelm[1]*elm_->belm[2];
  sbelm2[8+3] 	= elm_->sbelm[2]*elm_->belm[3];
  sbelm2[12+0] 	= elm_->sbelm[1]*elm_->belm[1];
  sbelm2[12+1] 	= elm_->sbelm[1]*elm_->belm[3];
  sbelm2[12+2] 	= elm_->sbelm[1]*elm_->belm[3];
  sbelm2[12+3] 	= elm_->sbelm[3]*elm_->belm[3];

  /* CALCUL DES MATRICES AVEC UNE DERIVEE PREMIERE */
  nsblas_dgemm	(transN,
		 transT,
		 &build_partialoff,
		 &dim,
		 &dim,
		 &regal1,
		 fem_->build_partialx,
		 &build_partialoff,
		 elm_->sbelm,
		 &dim,
		 &regal0,
		 fem_->build_partialr,
		 &build_partialoff);  

  /* CALCUL DES MATRICES AVEC DEUX DERIVEE PREMIERE */
  nsblas_dgemm	(transN,
		 transN,
		 &build_twice_partialoff,
		 &dimxdim,
		 &dimxdim,
		 &regal1,
		 fem_->build_twice_partialx,
		 &build_twice_partialoff,
		 sbelm2,
		 &dimxdim,
		 &regal0,
		 fem_->build_twice_partialr,
		 &build_twice_partialoff);    


}
