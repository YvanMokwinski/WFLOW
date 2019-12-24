#include  "ns_sys.h"
#include "ensBASIS.h"
#include "eHeaviside.h"
#include "eDirac.h"
// #include "Cmdline.h"
#include "ns_config_lapack.h"
#include "ns_constantes.h"
#include "ens_method_capturing.h"
#include "ns_var.h"
#include "ns_config.h"
#include "nsGLOBAL.h"
#include "ns_automatic_constantes.h"
#include "Monitor.h"



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



void 	nsGLOBAL_get_ddlcnc_u		(nsGLOBAL_ST*const 	self_,
					 const I 		ielm_,
					 pI 			ddl_)
{
  nsSPACE_cncelm(self_->space_u,
		 &ielm_,
		 ddl_);
}

void 	nsGLOBAL_get_ddlcnc_v		(nsGLOBAL_ST*const 	self_,
					 const I 		ielm_,
					 pI 			ddl_)
{
  nsSPACE_cncelm(self_->space_u,
		 &ielm_,
		 ddl_);
  { I i;
    for (i=0;i<_nu;++i)
      {
	ddl_[i] += self_->space_u->nddl;
      } }
}

void 	nsGLOBAL_get_ddlcnc_p		(nsGLOBAL_ST*const 	self_,
					 const I 		ielm_,
					 pI 			ddl_)
{
  nsSPACE_cncelm	(self_->space_p,
			 &ielm_,
			 ddl_);  
  { I i;
    for (i=0;i<_np;++i)
      {
	ddl_[i] += _dim*self_->space_u->nddl;
      } }
}


void nsGLOBAL_get_ddlcnc(nsGLOBAL_ST*const 	self_,
			 const I 		ielm_,
			 pI 			ddl_)
{
  nsGLOBAL_get_ddlcnc_u	(self_,ielm_,ddl_ + _ju);
  nsGLOBAL_get_ddlcnc_v	(self_,ielm_,ddl_ + _jv);
  nsGLOBAL_get_ddlcnc_p	(self_,ielm_,ddl_ + _jp);
}



#define	clr_dvect(_n,_x) { I _k; for (_k=0;_k<(_n);++_k) { (_x)[_k] = regal0;  }}
#define	cpy_dvect(_n,_x,_y) { I _k; for (_k=0;_k<(_n);++_k) { (_y)[_k] = (_x)[_k];  }}
#define	cpy_dvect_indirect(_n,_x,_b,_y) { I _k; for (_k=0;_k<(_n);++_k) { (_y)[_k] = (_x)[(_b)[_k]];  }}

void integrator_compute_viscosity_density_nodalelm(pGlobal const self_,const I iproc_,cst_pR felm_);


#define mvass(_n,_m,_ax,_a,_x,_y)  dgemv("N",_n,_m,(_ax),(_a),_n,(_x),&negal1,&regal1,(_y),&negal1)
#if 0
/*\brief data structure for mesh */
#define vheaviside(_self,_rn,_r) 		ns_heaviside	((_self)->params.iinfo[__ens_iinfo_heaviside],_rn,_r,&(_self)->rv[__ens_rv_eps]);
#endif

void ns_setelm			(nsGLOBAL_ST*const 	Z,
				 const I 		jelm_,
				 cst_pR 		x_,
				 cst_pR 		xi_,
				 cst_pR 		xii_)
{
#if 0
  I cncelm[256];
  I cellToNodes[16];
#endif
  I cellToEdges[16];
  //  nsGLOBAL_ST*const Z 		= &global[iproc_];
  pWorkelm const elm 		= Global_get_Workelm(Z);
  pParameters const gParameters	= Global_get_Parameters(Z);
  const I itime 		= TimeReadOnly_get_itime(GlobalReadOnly_get_Time(Z));
#if 0
  const I nvertex 		= Z->nvertex;
#endif
  
  ns_mesh_cooelm		((ns_mesh*)Z->mesh_usrptr,
				 &jelm_,
				 elm->cooelm);
  
#if 0
  ns_mesh_cncelm		((ns_mesh*)Z->mesh_usrptr,
				 &jelm_,
				 cncelm);
  ns_mesh_get_cellToNodes	((ns_mesh*)Z->mesh_usrptr,
				 &jelm_,
				 cellToNodes);
#endif

  ns_mesh_get_cellToEdges	((ns_mesh*)Z->mesh_usrptr,
				 &jelm_,
				 cellToEdges);

  ns_mesh_get_trelm		((ns_mesh*)Z->mesh_usrptr,
				 &jelm_,
				 elm->belm);
  
  elm->alpha2 		= Z->jacelm[jelm_];
  elm->ialpha2 		= regal1/elm->alpha2;
  elm->edgelm[0]     	= elm->cooelm[1]-elm->cooelm[0];
  elm->edgelm[1]     	= elm->cooelm[4]-elm->cooelm[3];
  elm->edgelm[2]     	= elm->cooelm[2]-elm->cooelm[1];
  elm->edgelm[3]     	= elm->cooelm[5]-elm->cooelm[4];
  elm->edgelm[4]     	= elm->cooelm[0]-elm->cooelm[2];
  elm->edgelm[5]     	= elm->cooelm[3]-elm->cooelm[5];

  elm->sbelm[3]   	= elm->belm[3]*elm->alpha2;
  elm->sbelm[1]   	= elm->belm[1]*elm->alpha2;
  elm->sbelm[2]   	= elm->belm[2]*elm->alpha2;
  elm->sbelm[0]   	= elm->belm[0]*elm->alpha2;

  
  {
#if 0
    I jedges[3];    
    jedges[0] 		= cncelm[3] - nvertex;
    jedges[1] 		= cncelm[4] - nvertex;
    jedges[2] 		= cncelm[5] - nvertex;
#endif
    elm->normalelmx[0] 	= Z->normaledge[_dim*cellToEdges[0]+0];
    elm->normalelmy[0] 	= Z->normaledge[_dim*cellToEdges[0]+1];
    elm->normalelmx[1] 	= Z->normaledge[_dim*cellToEdges[1]+0];
    elm->normalelmy[1] 	= Z->normaledge[_dim*cellToEdges[1]+1];
    elm->normalelmx[2] 	= Z->normaledge[_dim*cellToEdges[2]+0];
    elm->normalelmy[2] 	= Z->normaledge[_dim*cellToEdges[2]+1];

#if 0
#ifndef NDEBUG
    printf("-> "ifmt","ifmt","ifmt","ifmt"\n",Z->nedge,cellToEdges[0],cellToEdges[1],cellToEdges[2]);
    { I i;
      for (i=0;i<3;++i)
	{
	  if (cellToEdges[i] >= Z->nedge)
	    {
	      printf("problem cellToEdges["ifmt"] >= "ifmt"\n",
		     cellToEdges[i],
		     Z->nedge);
	    }
	} }
#endif
#endif
    elm->alpha3[0]     	= Z->jacedge[cellToEdges[0]];
    elm->alpha3[1]     	= Z->jacedge[cellToEdges[1]];
    elm->alpha3[2]     	= Z->jacedge[cellToEdges[2]];
  }

  nsGLOBAL_get_ddlcnc(Z,jelm_,elm->locnumer);
  cpy_dvect_indirect(_total_nddlelm,x_,elm->locnumer,elm->ddlelm);
  if (Parameters_geti(gParameters,__ens_iinfo_ntime)>1)	
    cpy_dvect_indirect	(_total_nddlelm,xi_,elm->locnumer,elm->ddlelmi);  
  if ( (Parameters_geti(gParameters,__ens_iinfo_ntime)>1) AND (itime>0) )					
    cpy_dvect_indirect	(_total_nddlelm,xii_,elm->locnumer,elm->ddlelmii);

}	

#ifndef NDEBUG

static void ns_rhselm2rhs(pGlobal const self_,
			  const I jelm_,
			  pR global_rhs_)
{ 
  nsGLOBAL_get_ddlcnc(self_,jelm_,self_->m_workelm.locnumer);
  { I i;
    for (i=0;i<_total_nddlelm;++i)
      {
	global_rhs_[self_->m_workelm.locnumer[i]] += self_->m_workelm.locresidu[i];
      } }    

}
#else

#define ns_rhselm2rhs(_self,_jelm,_global_rhs)  nsGLOBAL_get_ddlcnc((_self),_jelm,(_self)->m_workelm.locnumer); {I i;for (i=0;i<_total_nddlelm;++i){_global_rhs[(_self)->m_workelm.locnumer[i]] += (_self)->m_workelm.locresidu[i];}}

#endif

//#include "libns_exact.h"
#include "libns_exact_f3_u6_p3.c"
//#include "PredefinedData/libns_exact.c"

void ns_residu_exact_dens_time_derivative(pGlobal 	self_,
					  const I 	jelm_,
					  cst_pR 	xa_,
					  cst_pR 	denselm_,
					  cst_pR 	uelm_,
					  cst_pR 	uelmi_,
					  cst_pR 	uelmii_,
					  pR 		residu_,
					  const eTransientMethod scheme)
{

  pTimeReadOnly 		const gTimeInfo 	= GlobalReadOnly_get_Time(self_);

  pFem const fem					= Global_get_Fem(self_);
  
#if 0
#ifndef NDEBUG
  wrong_param(NOT xa_,1,ns_residu_exact_dens_time_derivative);
  wrong_param(NOT denselm_,2,ns_residu_exact_dens_time_derivative);
  wrong_param(NOT uelm_,3,ns_residu_exact_dens_time_derivative);
  wrong_param(NOT uelmi_,4,ns_residu_exact_dens_time_derivative);
  wrong_param(NOT uelmii_,5,ns_residu_exact_dens_time_derivative);
  wrong_param(NOT residu_,6,ns_residu_exact_dens_time_derivative);
#endif
#endif
  
  const R aa = xa_[0] * self_->jacelm[jelm_];
  switch(scheme)
    {
      
    case __eTransientMethod_NO:
      {
	break;
      }
      
    case __eTransientMethod_EULER:
      {
	
	const R a = self_->rv[__ens_rv_idt] * aa;

	R tmp[_nu];
	{ I k;
	  for (k=0;k<_nu;++k)
	    {		
	      tmp[k] = (uelm_[k]-uelmi_[k]);
	    } }
	
	nsblas_dgemv	(transN,
			 &nuxnu,
			 &ndg,
			 &regal1,
			 exactref_dens_u_u,
			 &nuxnu,
			 denselm_,
			 &negal1,
			 &regal0,
			 Fem_get_dens_u_u(fem),
			 &negal1);
	
	mvass(&nu,
	      &nu,
	      &a,
	      FemReadOnly_get_dens_u_u(fem),
	      tmp,
	      residu_);
	
	break;
      }
      
    case __eTransientMethod_IMR:
      {
	const R a = self_->rv[__ens_rv_idt]*((R)2.0);

	R tmp[_nu];
	{ I k;
	  for (k=0;k<_nu;++k)
	    {		
	      tmp[k] = (uelm_[k]-uelmi_[k]);
	    } }
	
	nsblas_dgemv	(transN,
			 &nuxnu,
			 &ndg,
			 &regal1,
			 exactref_dens_u_u,
			 &nuxnu,
			 denselm_,
			 &negal1,
			 &regal0,
			 Fem_get_dens_u_u(fem),
			 &negal1);

	mvass		(&nu,
			 &nu,
			 &a,
			 FemReadOnly_get_dens_u_u(fem),
			 tmp,
			 residu_);
	
	break;
      }
      
    case __eTransientMethod_GEAREULER:
      {
	if (gTimeInfo->itime>0)
	  {
	    const R a = aa;
	    
	    R tmp[_nu];
	    { I k;
	      for (k=0;k<_nu;++k)
		{		
		  tmp[k] = ( (regal2+self_->rv[__ens_rv_ratio_dti])*uelm_[k]-(regal2 + self_->rv[__ens_rv_ratio_dti] + self_->rv[__ens_rv_iratio_dti])*uelmi_[k] + self_->rv[__ens_rv_iratio_dti]*uelmii_[k] )/(self_->rv[__ens_rv_dt]+self_->rv[__ens_rv_dti]);
		} }

	    nsblas_dgemv	(transN,
				 &nuxnu,
				 &ndg,
				 &regal1,
				 exactref_dens_u_u,
				 &nuxnu,
				 denselm_,
				 &negal1,
				 &regal0,
				 Fem_get_dens_u_u(fem),
				 &negal1);
	    
	    mvass		(&nu,
				 &nu,
				 &a,
				 FemReadOnly_get_dens_u_u(fem),
				 tmp,
				 residu_);
	    
	  }
	else
	  {
	    const R a = self_->rv[__ens_rv_idt]*aa;
	    R tmp[_nu];
	    { I k;
	      for (k=0;k<_nu;++k)
		{		
		  tmp[k] = (uelm_[k]-uelmi_[k]);
		} }
	    nsblas_dgemv(transN,&nuxnu,&ndg,&regal1,exactref_dens_u_u,&nuxnu,denselm_,&negal1,&regal0,Fem_get_dens_u_u(fem),&negal1);
	    mvass(&nu,&nu,&a,FemReadOnly_get_dens_u_u(fem),tmp,residu_);
	  }
	break;
      }

    case __eTransientMethod_GEARIMR:
      {
	if (gTimeInfo->itime>0)
	  {
	    const R a = aa;
	    R tmp[_nu];
	    { I k;
	      for (k=0;k<_nu;++k)
		{		
		  tmp[k] = ( (regal2+self_->rv[__ens_rv_ratio_dti])*uelm_[k]-(regal2 + self_->rv[__ens_rv_ratio_dti] + self_->rv[__ens_rv_iratio_dti])*uelmi_[k] + self_->rv[__ens_rv_iratio_dti]*uelmii_[k] )/(self_->rv[__ens_rv_dt]+self_->rv[__ens_rv_dti]);
		} }
	    
	    nsblas_dgemv(transN,
			 &nuxnu,
			 &ndg,
			 &regal1,
			 exactref_dens_u_u,
			 &nuxnu,
			 denselm_,
			 &negal1,
			 &regal0,
			 Fem_get_dens_u_u(fem),
			 &negal1);

	    mvass(&nu,
		  &nu,
		  &a,
		  FemReadOnly_get_dens_u_u(fem),
		  tmp,
		  residu_);	    
	  }
	else
	  {
	    const R a = aa*self_->rv[__ens_rv_idt]*((R)2.0);
	    R tmp[_nu];
	    { I k;
	      for (k=0;k<_nu;++k)
		{		
		  tmp[k] = (uelm_[k]-uelmi_[k]);
		} }
	    
	    nsblas_dgemv(transN,
			 &nuxnu,
			 &ndg,
			 &regal1,
			 exactref_dens_u_u,
			 &nuxnu,
			 denselm_,
			 &negal1,
			 &regal0,
			 Fem_get_dens_u_u(fem),
			 &negal1);
	    
	    mvass(&nu,
		  &nu,
		  &a,
		  FemReadOnly_get_dens_u_u(fem),
		  tmp,
		  residu_);
	    
	  }
	break;
      }
      
    case __eTransientMethod_TRAPEZE:
      {
	Monitor_errmsg(self_->iproc,
		       "__eTransientMethod_TRAPEZE not yet");
	break;
      }
      
    case __eTransientMethod_ERROR:
    case __eTransientMethod_ALL:
      {
	Monitor_errmsg(self_->iproc,
		       "ns_residu_exact_dens_time_derivative:wrong case");
	break;
      }
      
    }
}






void ns_build_residuelm_u(nsGLOBAL_ST*const Z,
			  const I 	jelm_,
			  cst_pR 	denselm_,
			  cst_pR 	viscelm_,
			  cst_pR 	uelm_,
			  cst_pR 	velm_,
			  cst_pR 	pelm_,
			  pR 		residu_u_,
			  pR 		residu_v_)
{
  
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
 
 //  nsGLOBAL_ST*const Z 	= &global[iproc_];
  pFem const fem	= Global_get_Fem(Z);

  /* terme de gravite (0,-1) */
  { const R a = Z->jacelm[jelm_]*Z->rv[__ens_rv_ifroude];
    nsblas_dgemv(transN,
		 &nu,
		 &ndg,
		 &a,
		 gravity_term,
		 &nu,
		 denselm_,
		 &negal1,
		 &regal1,
		 residu_v_,
		 &negal1); }
  
  /* residu pressure -/- u */
  { const R a = ((R)-1.0);
    nsblas_dgemv(transT,
		 &np,
		 &nu,
		 &a,
		 fem->exact_dxu_p,
		 &np,pelm_,
		 &negal1,
		 &regal1,
		 residu_u_,
		 &negal1);
    
    nsblas_dgemv(transT,
		 &np,
		 &nu,
		 &a,
		 fem->exact_dyu_p,
		 &np,
		 pelm_,
		 &negal1,
		 &regal1,
		 residu_v_,
		 &negal1); }
  
  /* stress  */
  nsblas_dgemv(transN,
	       &nuxnu,
	       &ndg,
	       &regal1,
	       fem->exact_visc_dxu_dxu,
	       &nuxnu,viscelm_,
	       &negal1,
	       &regal0,
	       fem->visc_dxu_dxu,
	       &negal1);
  
  { const R a = Z->rv[__ens_rv_ireynold] * regal2;
    nsblas_dgemv(transN,&nu,&nu,&a,fem->visc_dxu_dxu,&nu,uelm_,&negal1,&regal1,residu_u_,&negal1); }
  { const R a = Z->rv[__ens_rv_ireynold];
    nsblas_dgemv(transN,&nu,&nu,&a,fem->visc_dxu_dxu,&nu,velm_,&negal1,&regal1,residu_v_,&negal1); }
  
  nsblas_dgemv(transN,&nuxnu,&ndg,&regal1,fem->exact_visc_dyu_dyu,&nuxnu,viscelm_,&negal1,&regal0,fem->visc_dyu_dyu,&negal1);
  { const R a = Z->rv[__ens_rv_ireynold] * regal2;
    nsblas_dgemv(transN,&nu,&nu,&a,fem->visc_dyu_dyu,&nu,velm_,&negal1,&regal1,residu_v_,&negal1); }
  { const R a = Z->rv[__ens_rv_ireynold];
    nsblas_dgemv(transN,&nu,&nu,&a,fem->visc_dyu_dyu,&nu,uelm_,&negal1,&regal1,residu_u_,&negal1); }
  
  { const R a = Z->rv[__ens_rv_ireynold];
    nsblas_dgemv(transN,&nuxnu,&ndg,&regal1,fem->exact_visc_dyu_dxu,&nuxnu,viscelm_,&negal1,&regal0,fem->visc_dyu_dxu,&negal1);
    nsblas_dgemv(transN,&nu,&nu,&a,fem->visc_dyu_dxu,&nu,uelm_,&negal1,&regal1,residu_v_,&negal1); 

    nsblas_dgemv(transN,&nuxnu,&ndg,&regal1,fem->exact_visc_dxu_dyu,&nuxnu,viscelm_,&negal1,&regal0,fem->visc_dxu_dyu,&negal1);
    nsblas_dgemv(transN,&nu,&nu,&a,fem->visc_dxu_dyu,&nu,velm_,&negal1,&regal1,residu_u_,&negal1); }
  
  /* residu d'acceleration */

#if 0
  nsblas_dgemv(transN,&nuxnuxnu,&ndg,&regal1,fem->exact_dens_u_dxu_u,&nuxnuxnu,denselm_,&negal1,&regal0,fem->dens_u_dxu_u,&negal1);
  nsblas_dgemv(transN,&nuxnuxnu,&ndg,&regal1,fem->exact_dens_u_dyu_u,&nuxnuxnu,denselm_,&negal1,&regal0,fem->dens_u_dyu_u,&negal1);
  /**/
  nsblas_dgemv(transN,&nuxnu,&nu,&regal1,fem->dens_u_dxu_u,&nuxnu,uelm_,&negal1,&regal0,fem->dens_dxu_u,&negal1);
  nsblas_dgemv(transN,&nuxnu,&nu,&regal1,fem->dens_u_dyu_u,&nuxnu,velm_,&negal1,&regal0,fem->dens_dyu_u,&negal1);
  /**/	  
  nsblas_dgemv(transN,&nu,&nu,&regal1,fem->dens_dxu_u,&nu,uelm_,&negal1,&regal1,residu_u_,&negal1);
  nsblas_dgemv(transN,&nu,&nu,&regal1,fem->dens_dyu_u,&nu,uelm_,&negal1,&regal1,residu_u_,&negal1);
  nsblas_dgemv(transN,&nu,&nu,&regal1,fem->dens_dxu_u,&nu,velm_,&negal1,&regal1,residu_v_,&negal1);
  nsblas_dgemv(transN,&nu,&nu,&regal1,fem->dens_dyu_u,&nu,velm_,&negal1,&regal1,residu_v_,&negal1);
#endif
  
}


Err ns_build_residu	(pGlobal 	self_,
			 pR 		global_rhs_,
			 cst_pR 	x_,
			 cst_pR 	xi_,
			 cst_pR 	xii_)
{
  pFem const 		fem		= Global_get_Fem(self_);
  pParameters const 	gParameters 	= Global_get_Parameters(self_);
  pWorkelm const 	elm		= Global_get_Workelm(self_);
  pSparseStokes const  	matrix 		= self_->sparseStokes;

  const eTransientMethod scheme		= Parameters_get_eTransientScheme(gParameters,__eKindEquation_VELOCITY);

  
#ifndef NDEBUG
  const L  verbose 			= Parameters_getl(gParameters,__ens_linfo_verbose);
#endif

  
  Err err  = __eErr_no;

#ifndef NDEBUG
  if (verbose)
    {
      Monitor_msg(self_->iproc,"ns_build_residu:compute linear pressure residual ...");
    }
#endif
  printf("call pressure rhs\n");
  SparseStokes_pressure_rhs(matrix,
			    x_,
			    global_rhs_);
  
  /*
    on calcule residu en fuv
  */
  { const I NELM = self_->nelm;
    I jelm;
    for (jelm=0;jelm<NELM;++jelm)
      {      

	ns_setelm(self_,
		  jelm,
		  x_,
		  xi_,
		  xii_);

	Fem_setelm(fem,
		   elm);	
	
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
	pWorkelm const workelm 		= Global_get_Workelm(self_);

	{ I k;
	  for (k=0;k<_ndg;++k)						
	    {
	      R hh 		= workelm->viscelm[k];
		if (hh<0.0)
		  hh=0.0;
		if (hh>1.0)
		  hh=1.0;
	      workelm->viscelm[k] 	= 1.0+hh*(self_->rv[__ens_rv_coeff_ratio_viscosity]-1.0);
	      workelm->denselm[k] 	= 1.0+hh*(self_->rv[__ens_rv_coeff_ratio_density]-1.0);
	    } }
	
	/*
	  v1 + F(v2-v1) 
	  r + F(1-r)
	  printf("%e\n",Z->rv[__ens_rinfo_ratio_viscosity]); 
	*/
	
	clr_dvect(_total_nddlelm,
		  elm->locresidu);
	
	ns_build_residuelm_u			(self_,	
						 jelm,
						 elm->denselm,
						 elm->viscelm,
						 elm->uelm,
						 elm->velm,
						 elm->pelm,
						 &elm->locresidu[offu],
						 &elm->locresidu[offv]); 

	ns_residu_exact_dens_time_derivative	(self_,
						 jelm,
						 &regal1,
						 elm->denselm,
						 elm->uelm,
						 elm->uelmi,
						 elm->uelmii,
						 &elm->locresidu[offu],
						 scheme);

	ns_residu_exact_dens_time_derivative	(self_,
						 jelm,
						 &regal1,
						 elm->denselm,
						 elm->velm,
						 elm->velmi,
						 elm->velmii,
						 &elm->locresidu[offv],
						 scheme);
	
	ns_rhselm2rhs				(self_,
						 jelm,
						 global_rhs_);
	
      } }
  
  //  printf("putain de ta race\n");
  //  exit(1);

  /* RESIDU SLIP CONDITION */
  if (Parameters_getl(gParameters,__ens_linfo_slip))
    {
#ifndef NDEBUG
      if (verbose)
	{
	  Monitor_msg(self_->iproc,"ns_build_residu:compute slip boundary condition ...");
	}
#endif
      SparseStokes_slip_rhs(matrix,
			    x_,
			    global_rhs_);
    }   
  else
    {
#ifndef NDEBUG
      if (verbose)
	{
	  Monitor_msg(self_->iproc,"ns_build_residu:no slip boundary condition to compute ");
	}
#endif
    }
  
  /* IMPOSITION DIRICHLET */
#ifndef NDEBUG
  if (verbose)
    {
      Monitor_msg(self_->iproc,"ns_build_residu:compute boundary dirichlet condition");
    }
#endif
  
  /* DIRICHLET CONDITION */
  SparseStokes_dirichlet_rhs(matrix,
			     global_rhs_);
  
  return err;
}

#if 0
void  SparseBlock_ass(const I 		iproc_,
		      pSparseBlock 	block_,
		      cst_pI 		ni_,
		      cst_pI 		ddli_,
		      cst_pI 		nj_,
		      cst_pI 		ddlj_,
		      cst_pR 		locmat_,
		      cst_pI 		locmatoff_)
{
  SparseBlock_ass_blank(block_,
			ni_,
			ddli_,
			nj_,
			ddlj_,
			locmat_,
			locmatoff_,
			global[iproc_].blank);
}
#endif

void ns_build_system_divergence(pGlobal	const self_)
{ 
  const I N 		= self_->nelm;
  const I nuxdim	= _nuxdim;
  R 	locmat_BU[_nuxdim*_np];
  R 	locmat_UB[_nuxdim*_np];
  I 	locnumeru[_nuxdim];
  I 	locnumerp[_np];
  
  pWorkelm const elm = Global_get_Workelm(self_);
  { I jelm;
    for (jelm=0;jelm<N;++jelm)
      {            
	clr_dvect(total_nddlelmxtotal_nddlelm,
		  elm->locmat);
	
	clr_dvect(_total_nddlelm,
		  elm->locresidu);      

	nsGLOBAL_get_ddlcnc_u(self_,
			      jelm,
			      locnumeru);
	
	nsGLOBAL_get_ddlcnc_v(self_,
			      jelm,
			      &locnumeru[_nu]);

	nsGLOBAL_get_ddlcnc_p(self_,
			      jelm,
			      locnumerp);    

	elm->alpha2	= self_->jacelm[jelm];
	elm->belm[3]   	= self_->trelm[4*jelm+3];
	elm->belm[1]   	= self_->trelm[4*jelm+1];
	elm->belm[2]   	= self_->trelm[4*jelm+2];
	elm->belm[0]   	= self_->trelm[4*jelm+0];
	elm->sbelm[3]   = elm->belm[3]*elm->alpha2;
	elm->sbelm[1]   = elm->belm[1]*elm->alpha2;
	elm->sbelm[2]   = elm->belm[2]*elm->alpha2;
	elm->sbelm[0]   = elm->belm[0]*elm->alpha2;
	
	{ I i;
	  for (i=0;i<_np;++i)
	    {
	      { I j;
		for (j=0;j<_nu;++j)
		  {
		    locmat_BU[j*_np+i]    = elm->sbelm[0] * exactref_dru_p[j*_np+i]+elm->sbelm[2] * exactref_dsu_p[j*_np+i];
		    locmat_UB[i*nuxdim+j] = -locmat_BU[j*_np+i];		
		  } }
	      { I j;
		for (j=_nu;j<_nuxdim;++j)
		  {
		    locmat_BU[j*_np+i] = elm->sbelm[1] * exactref_dru_p[(j-_nu)*_np+i]+elm->sbelm[3] * exactref_dsu_p[(j-_nu)*_np+i];		
		    locmat_UB[i*nuxdim+j] = -locmat_BU[j*_np+i];		
		  } }
	    } }         
  
	SparseBlock_ass_blank(&self_->sparseStokes->block_UB,
			      &nuxdim,
			      locnumeru,
			      &np,
			      locnumerp,
			      locmat_UB,
			      &nuxdim,
			      self_->blank);
  
	SparseBlock_ass_blank(&self_->sparseStokes->block_BU,
			      &np,
			      locnumerp,
			      &nuxdim,
			      locnumeru,
			      locmat_BU,
			      &np,
			      self_->blank);
  
      } }
}


#if 0
void nsdg_csf(pR trelm_,
	      pR itrelm_)
{


  pR wei 		= &global[0].var_dg.x[ielm*_ndg+0];
  const R vol 	= wei[0]*nsSQRT(2.0);
  const R dr  	= wei[1]*6.0+wei[2]*3.464101615137755e0;
  const R ds  	= 6.928203230275509e0*wei[2];      
  const R dx 	= itrelm_[0]*dr + itrelm_[2]*ds;
  const R dy 	= itrelm_[1]*dr + itrelm_[3]*ds;
  const R len	= dx*dx+dy*dy;
  if (len>0.5)
    {
      const R nx 	= dx/nsSQRT(len);
      const R ny 	= dy/nsSQRT(len);
      wei[1] 		= trelm_[0]*nx + trelm_[0]*ny;
      wei[2] 		= trelm_[1]*nx + trelm_[2]*ny;
    }
}
#endif


void integrator_assmatelm_full(pGlobal 	self_,
			       const I 	jelm_)
{ 
  //  nsGLOBAL_ST*const Z = &global[iproc_];
  pParametersReadOnly const gParameters 	= GlobalReadOnly_get_Parameters(self_);
  pWorkelm const elm 		= Global_get_Workelm(self_);

  nsGLOBAL_get_ddlcnc(self_,
		      jelm_,
		      elm->locnumer);
  
  I N = 0;  

  { I i;
    for (i=0;i<_total_nddlelm;++i)
      {
	if (NOT self_->blank[elm->locnumer[i]])
	  {
	    self_->blank[elm->locnumer[i]] 	= ++N;
	    self_->iwork[N-1] 	  	= elm->locnumer[i];
	  }
      } }

  cst_pI Ab = Sparse_get_b(self_->sparseStokes->A);
  cst_pI Ai = Sparse_get_i(self_->sparseStokes->A);
  pR Ax = Sparse_get_ownx(self_->sparseStokes->A);

  { I j,k;
    for (j=0;j<_total_nddlelm-_np;++j)
      {
	I jddl 	= elm->locnumer[j];
	const I n		= Ab[jddl+1]-Ab[jddl];
	cst_pI it 		= &Ai[Ab[jddl]-1];
	pR rt 			= &Ax[Ab[jddl]-1];
	{ I i;
	  for (i=0;i<n;++i)
	    {
	      if ( (k=self_->blank[it[i]-1]) )
		{
		  if (Parameters_getl(gParameters,__ens_linfo_pspg))
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
	    } }
      } }  

  { I j,k;
    for (j=_nuxdim;j<_total_nddlelm;++j)
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
	{ I i;
	  for (i=0;i<n;++i)
	    {
	      if ( (k=self_->blank[it[i]-1]) )
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
		      
		      if (Parameters_getl(gParameters,__ens_linfo_pspg))
			{
			  rt[i] += elm->locmat[j + (k-1)*off];
			}
		    } 
		  else
		    {
		      if (Parameters_getl(gParameters,__ens_linfo_pspg))
			rt[i] += elm->locmat[j + (k-1)*off];
		    }
		}
	    } }
      } }   
  { I i;
    for (i=0;i<N;++i) 
      {
	self_->blank[self_->iwork[i]] = 0;
      } }
}




void ns_jacobian_exact_dens_time_derivative(pGlobal 	self_,
					    const I 	jelm_)
{
#ifndef NDEBUG
  DebugVerif(self_);
  DebugVerif(jelm_ >= 0);
#endif
  
  pParametersReadOnly 		const gParameters 	= GlobalReadOnly_get_Parameters(self_);
  pTimeReadOnly 		const gTimeInfo		= GlobalReadOnly_get_Time(self_);
  
  pWorkelmReadOnly const elm 	= GlobalReadOnly_get_Workelm(self_);
  pFem const fem		= Global_get_Fem(self_);

  pR J_UU = Fem_getJ_UU(fem);
  pR J_VV = Fem_getJ_VV(fem);

  const eTransientMethod scheme 			= Parameters_get_eTransientScheme(gParameters,__eKindEquation_VELOCITY);
  switch(scheme)
    {
    case __eTransientMethod_NO:
      {
	break;
      }

    case __eTransientMethod_EULER:
      {
	{ const R a = self_->rv[__ens_rv_idt]*self_->jacelm[jelm_];
	  mvass(&nuxnu,&ndg,&a,exactref_dens_u_u,elm->denselm,J_UU);
	  mvass(&nuxnu,&ndg,&a,exactref_dens_u_u,elm->denselm,J_VV); }
#if 0
	{ const R a = self_->rv[__ens_rv_idt]*self_->jacelm[jelm_]* (-self_->rv[__ens_rv_coeff_ratio_density]);
	  mvass(&nuxndg,&nu,&a,exactref_u_dens_u,&elm->ddlelm[_ju],fem->J_UF);
	  mvass(&nuxndg,&nu,&a,exactref_u_dens_u,&elm->ddlelm[_jv],fem->J_VF); }
#endif
	break;
      }

    case __eTransientMethod_IMR:
      {
	const R a = self_->rv[__ens_rv_idt]*((R)2.0)*self_->jacelm[jelm_];
	mvass(&nuxnu,&ndg,&a,exactref_dens_u_u,elm->denselm,J_UU);
	mvass(&nuxnu,&ndg,&a,exactref_dens_u_u,elm->denselm,J_VV);
	break;
      }

    case __eTransientMethod_GEAREULER:
      {
	if (gTimeInfo->itime>0)
	  {
	    const R a = (regal2+self_->rv[__ens_rv_ratio_dti])*self_->jacelm[jelm_];
	    mvass(&nuxnu,&ndg,&a,exactref_dens_u_u,elm->denselm,J_UU);
	    mvass(&nuxnu,&ndg,&a,exactref_dens_u_u,elm->denselm,J_VV);
	  }
	else
	  {
	    const R a = self_->rv[__ens_rv_idt]*self_->jacelm[jelm_];
	    mvass(&nuxnu,&ndg,&a,exactref_dens_u_u,elm->denselm,J_UU);
	    mvass(&nuxnu,&ndg,&a,exactref_dens_u_u,elm->denselm,J_VV);
	  }
	break;
      }

    case __eTransientMethod_GEARIMR:
      {
	if (gTimeInfo->itime>0)
	  {
	    const R a = (regal2+self_->rv[__ens_rv_ratio_dti])*self_->jacelm[jelm_];
	    mvass(&nuxnu,&ndg,&a,exactref_dens_u_u,elm->denselm,J_UU);
	    mvass(&nuxnu,&ndg,&a,exactref_dens_u_u,elm->denselm,J_VV);
	  }
	else
	  {
	    const R a = self_->rv[__ens_rv_idt]*((R)2.0)*self_->jacelm[jelm_];
	    mvass(&nuxnu,&ndg,&a,exactref_dens_u_u,elm->denselm,J_UU);
	    mvass(&nuxnu,&ndg,&a,exactref_dens_u_u,elm->denselm,J_VV);
	  }
	break;
      }

    case __eTransientMethod_TRAPEZE:
      {
	Monitor_errmsg(self_->iproc,"__eTransientMethod_TRAPEZE not yet");
	break;
      }

    case __eTransientMethod_ERROR:
    case __eTransientMethod_ALL:
      {
	break;
      }
    }
}

#define changement_variable(_ratio,_f)   ((_ratio)+(((R)1.0)-(_ratio))*(_f))
#define changement_variable2(_ratio,_f)  ((_ratio)+(((R)1.0)-(_ratio))*(_f))

void integrator_compute_viscosity_density_nodalelm(pGlobal 	self_,
						   const I 	iproc_,
						   cst_pR 	felm_)
{
  //  nsGLOBAL_ST*const Z 		= &global[iproc_];
  pWorkelm const workelm 	= Global_get_Workelm(self_);
  { I k;							
    for (k=0;k<_ndg;++k)						
      {
	workelm->denselm[k]=1.0;
	workelm->viscelm[k]=1.0;
      } }

#if 0
  if (Parameters_getl(gParameters,__ens_linfo_transport))
    {
      if (NOT Parameters_getl(gParameters,__ens_linfo_transport_uncoupled))
	{
	  if (Parameters_getl(gParameters,__ens_linfo_color))							
	    {									
	      { I k;							
		for (k=0;k<_nf;++k)						
		  {
		    workelm->denselm[k] = changement_variable(self_->rv[__ens_rv_coeff_ratio_density],felm_[k]);	
		    workelm->viscelm[k] = changement_variable(self_->rv[__ens_rv_coeff_ratio_viscosity],felm_[k]);	
		  } }
	    }								
	  else									
	    {			
	      cpy_dvect(_nf,felm_,workelm->viscelm);
	      vheaviside(self_,_nf,workelm->viscelm);
	      { I k;
		for (k=0;k<_nf;++k)
		  {
		    const R aa = workelm->viscelm[k];

		    workelm->viscelm[k] = changement_variable(self_->rv[__ens_rv_coeff_ratio_viscosity],
							      aa);
		    
		    workelm->denselm[k] = changement_variable(self_->rv[__ens_rv_coeff_ratio_density],
							      aa);
		    
		  } }								
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
}


void integrator_compute_viscosity_density_nodalelm2(pGlobal self_,
						    cst_pR felm_)
{

  pWorkelm const workelm = Global_get_Workelm(self_);

  { I k;							
    for (k=0;k<_ndg;++k)						
      {
	workelm->denselm[k] = 1.0;
	workelm->viscelm[k] = 1.0;
      } }
  
#if 0

  /*
    pParameters const gParameters = Global_get_Parameters(Z);
  */
  if (Parameters_getl(gParameters,__ens_linfo_transport))
    {
      if (NOT Parameters_getl(gParameters,__ens_linfo_transport_uncoupled))
	{
	  if (Parameters_getl(gParameters,__ens_linfo_color))							
	    {									
	      { I k;							
		for (k=0;k<_nf;++k)						
		  {
		    workelm->denselm[k] = changement_variable2(self_->rv[__ens_rv_coeff_ratio_density],
							       felm_[k]);
		    workelm->viscelm[k] = changement_variable2(self_->rv[__ens_rv_coeff_ratio_viscosity],
							       felm_[k]);
		  } }
	    }								
	  else									
	    {			
	      cpy_dvect(_nf,felm_,workelm->viscelm);
	      vheaviside(self_,_nf,workelm->viscelm);
	      { I k;
		for (k=0;k<_nf;++k)
		  {
		    const R aa = workelm->viscelm[k];
		    workelm->viscelm[k] = changement_variable2(self_->rv[__ens_rv_coeff_ratio_viscosity],
							       aa);
		    workelm->denselm[k] = changement_variable2(self_->rv[__ens_rv_coeff_ratio_density],
							       aa);
		  } }								
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

}




void ns_build_system_stokes_quadrature_free(pGlobal 	self_,
					    cst_pR 	x_,
					    cst_pR 	xi_,
					    cst_pR 	xii_)
{ 
  pParameters const gParameters = Global_get_Parameters(self_);
  pFem const fem		= Global_get_Fem(self_);

  pR J_UU = Fem_getJ_UU(fem);
  pR J_VV = Fem_getJ_VV(fem);
  pR J_UV = Fem_getJ_UV(fem);
  pR J_VU = Fem_getJ_VU(fem);
  
#if 0
#ifndef NDEBUG
  wrong_param(NOT x_,1,ns_build_system_stokes_quadrature_free);
  wrong_param(NOT xi_,2,ns_build_system_stokes_quadrature_free);
  wrong_param(NOT xii_,3,ns_build_system_stokes_quadrature_free);
#endif
#endif

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

  /**OR (Z->newton_iter>0)*/
  if ( (NOT Parameters_getl(gParameters,__ens_linfo_time_adaptivity))  )
    {   
      pWorkelm const elm 	= Global_get_Workelm(self_);
      
      I jelm=0;
      const I N = self_->nelm;
      for (jelm=0;jelm<N;++jelm)
	{      	  
	  clr_dvect(total_nddlelmxtotal_nddlelm,elm->locmat);
	  ns_setelm(self_,
		    jelm,
		    x_,
		    xi_,
		    xii_);
	  
	  Fem_setelm(fem,
		     elm);
	  
#if 0
	  if (this->transport_usrptr)
	    {
	      transport_get_viscelm(this->transport_usrptr,jelm,elm->viscelm);
	      transport_get_denselm(this->transport_usrptr,jelm,elm->denselm);
	    }
	  else
	    {
	      elm->viscelm[0]=((R)1.0)/nsSQRT(2.0);
	      for (i=1;i<_ndg;++i)
		elm->viscelm[i]=((R)0.0);
	      elm->denselm[0]=((R)1.0)/nsSQRT(2.0);
	      for (i=1;i<_ndg;++i)
		elm->denselm[i]=((R)0.0);
	    }
	  if (Parameters_getl(gParameters,__ens_linfo_color))							
	    {									
	      { I k;							
		for (k=0;k<_nf;++k)						
		  {
		    workelm->denselm[k] = changement_variable2(self_->rv[__ens_rv_coeff_ratio_density],felm_[k]);	
		    workelm->viscelm[k] = changement_variable2(self_->rv[__ens_rv_coeff_ratio_viscosity],felm_[k]);	
		  } }
	    }								
	  else									
	    {			
	      cpy_dvect(_nf,felm_,workelm->viscelm);
	      vheaviside(iproc_,_nf,workelm->viscelm);
	      { I k;							
		for (k=0;k<_nf;++k)						
		  {								
		    const R aa = workelm->viscelm[k];
		    workelm->viscelm[k] = changement_variable2(self_->rv[__ens_rv_coeff_ratio_viscosity],aa);	
		    workelm->denselm[k] = changement_variable2(self_->rv[__ens_rv_coeff_ratio_density],aa);	
		  } }								
	    }
#endif

	  { I k;							
	    for (k=0;k<_ndg;++k)						
	      {
		elm->denselm[k] = 1.0;
	      } }
	  { I k;							
	    for (k=0;k<_ndg;++k)						
	      {
		elm->viscelm[k] = 1.0;
	      } }

#if 0
	  nsDG_transport_get_nodalvalues(self_->transport_usrptr,&jelm,elm->denselm);
	  nsDG_transport_get_nodalvalues(self_->transport_usrptr,&jelm,elm->viscelm);
#endif
	  
	  { I k;
	    for (k=0;k<_ndg;++k)						
	      {
		R hh 		= elm->viscelm[k];
		if (hh<0.0)
		  hh=0.0;
		if (hh>1.0)
		  hh=1.0;
		elm->viscelm[k] 	= 1.0+hh*(self_->rv[__ens_rv_coeff_ratio_viscosity]-1.0);
		elm->denselm[k] 	= 1.0+hh*(self_->rv[__ens_rv_coeff_ratio_density]-1.0);
#if 0
		elm->viscelm[k] 	= self_->rv[__ens_rv_coeff_ratio_viscosity]+hh*(1.0-self_->rv[__ens_rv_coeff_ratio_viscosity]);
		elm->denselm[k] 	= self_->rv[__ens_rv_coeff_ratio_viscosity]+hh*(1.0-self_->rv[__ens_rv_coeff_ratio_density]);
#endif
	      } }
	  
	  
	  clr_dvect(_nu*_nu,J_VU);
	  clr_dvect(_nu*_nu,J_UV);
	  clr_dvect(_nu*_nu,J_VV);
	  clr_dvect(_nu*_nu,J_UU);


	  /* CONSTRUCTION DES MATRICES */
	  /* jacobien exact terme de stress -/- u */
	  nsblas_dgemv(transN,&nuxnu,&ndg,&regal1,fem->exact_visc_dxu_dxu,&nuxnu,elm->viscelm,&negal1,&regal0,fem->visc_dxu_dxu,&negal1);
	  nsblas_dgemv(transN,&nuxnu,&ndg,&regal1,fem->exact_visc_dxu_dyu,&nuxnu,elm->viscelm,&negal1,&regal0,fem->visc_dxu_dyu,&negal1);
	  nsblas_dgemv(transN,&nuxnu,&ndg,&regal1,fem->exact_visc_dyu_dxu,&nuxnu,elm->viscelm,&negal1,&regal0,fem->visc_dyu_dxu,&negal1);
	  nsblas_dgemv(transN,&nuxnu,&ndg,&regal1,fem->exact_visc_dyu_dyu,&nuxnu,elm->viscelm,&negal1,&regal0,fem->visc_dyu_dyu,&negal1);	  

	  { const R a = self_->rv[__ens_rv_ireynold] * regal2;
	    nsblas_daxpy(&nuxnu,&a,fem->visc_dxu_dxu,&negal1,J_UU,&negal1); }
	  { const R a = self_->rv[__ens_rv_ireynold];
	    nsblas_daxpy(&nuxnu,&a,fem->visc_dyu_dyu,&negal1,J_UU,&negal1); }
	  { const R a = self_->rv[__ens_rv_ireynold];
	    nsblas_daxpy(&nuxnu,&a,fem->visc_dxu_dyu,&negal1,J_UV,&negal1); }
	  { const R a = self_->rv[__ens_rv_ireynold] * regal2;
	    nsblas_daxpy(&nuxnu,&a,fem->visc_dyu_dyu,&negal1,J_VV,&negal1); }
	  { const R a = self_->rv[__ens_rv_ireynold];
	    nsblas_daxpy(&nuxnu,&a,fem->visc_dxu_dxu,&negal1,J_VV,&negal1); }
	  { const R a = self_->rv[__ens_rv_ireynold];
	    nsblas_daxpy(&nuxnu,&a,fem->visc_dyu_dxu,&negal1,J_VU,&negal1); }  

#if 0
	  ns_jacobian_exact_dens_time_derivative(iproc_,jelm);
	  nsblas_dgemv(transN,&nuxnuxnu,&ndg,&regal1,fem->exact_dens_u_dxu_u,&nuxnuxnu,elm->denselm,&negal1,&regal0,fem->dens_u_dxu_u,&negal1);
	  nsblas_dgemv(transN,&nuxnuxnu,&ndg,&regal1,fem->exact_dens_u_dyu_u,&nuxnuxnu,elm->denselm,&negal1,&regal0,fem->dens_u_dyu_u,&negal1);
	  nsblas_dgemv(transN,&nuxnuxnu,&ndg,&regal1,fem->exact_dens_dxu_u_u,&nuxnuxnu,elm->denselm,&negal1,&regal0,fem->dens_dxu_u_u,&negal1);
	  nsblas_dgemv(transN,&nuxnuxnu,&ndg,&regal1,fem->exact_dens_dyu_u_u,&nuxnuxnu,elm->denselm,&negal1,&regal0,fem->dens_dyu_u_u,&negal1);
	  { const R a = self_->rv[__ens_rv_ireynold] * regal2;
	    nsblas_daxpy(&nuxnu,&a,fem->visc_dxu_dxu,&negal1,J_UU,&negal1); }
	  { const R a = self_->rv[__ens_rv_ireynold];
	    nsblas_daxpy(&nuxnu,&a,fem->visc_dyu_dyu,&negal1,J_UU,&negal1); }
	  { const R a = self_->rv[__ens_rv_ireynold];
	    nsblas_daxpy(&nuxnu,&a,fem->visc_dxu_dyu,&negal1,J_UV,&negal1); }
	  { const R a = self_->rv[__ens_rv_ireynold] * regal2;
	    nsblas_daxpy(&nuxnu,&a,fem->visc_dyu_dyu,&negal1,J_VV,&negal1); }
	  { const R a = self_->rv[__ens_rv_ireynold];
	    nsblas_daxpy(&nuxnu,&a,fem->visc_dxu_dxu,&negal1,J_VV,&negal1); }
	  { const R a = self_->rv[__ens_rv_ireynold];
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
	  nsblas_dgemv(transN,&nuxndg,&nu,&regal1,fem->exact_dxu_visc_dxu,&nuxndg,uelm,&negal1,&regal0,dxu_visc_dxu,&negal1);
	  nsblas_dgemv(transN,&nuxndg,&nu,&regal1,fem->exact_dyu_visc_dyu,&nuxndg,uelm,&negal1,&regal0,dyu_visc_dyu,&negal1);
	  nsblas_dgemv(transN,&nuxndg,&nu,&regal1,fem->exact_dxu_visc_dyu,&nuxndg,velm,&negal1,&regal0,dxu_visc_dyu,&negal1);
	  { const R a = self_->rv[__ens_rv_ireynold] * regal2* (-self_->rv[__ens_rv_coeff_ratio_viscosity]);
	    nsblas_daxpy(&nuxndg,&a,dxu_visc_dxu,&negal1,J_UF,&negal1); }
	  { const R a = self_->rv[__ens_rv_ireynold]* (-self_->rv[__ens_rv_coeff_ratio_viscosity]);
	    nsblas_daxpy(&nuxndg,&a,dyu_visc_dyu,&negal1,J_UF,&negal1); }
	  { const R a = self_->rv[__ens_rv_ireynold]* (-self_->rv[__ens_rv_coeff_ratio_viscosity]);
	    nsblas_daxpy(&nuxndg,&a,dxu_visc_dyu,&negal1,J_UF,&negal1); }      
	  nsblas_dgemv(transN,&nuxndg,&nu,&regal1,fem->exact_dxu_visc_dxu,&nuxndg,velm,&negal1,&regal0,dxu_visc_dxu,&negal1);
	  nsblas_dgemv(transN,&nuxndg,&nu,&regal1,fem->exact_dyu_visc_dyu,&nuxndg,velm,&negal1,&regal0,dyu_visc_dyu,&negal1);
	  nsblas_dgemv(transN,&nuxndg,&nu,&regal1,fem->exact_dyu_visc_dxu,&nuxndg,uelm,&negal1,&regal0,dyu_visc_dxu,&negal1);
	  { const R a = self_->rv[__ens_rv_ireynold] * regal2* (-self_->rv[__ens_rv_coeff_ratio_viscosity]);
	    nsblas_daxpy(&nuxndg,&a,dyu_visc_dyu,&negal1,J_VF,&negal1); }
	  { const R a = self_->rv[__ens_rv_ireynold]* (-self_->rv[__ens_rv_coeff_ratio_viscosity]);
	    nsblas_daxpy(&nuxndg,&a,dxu_visc_dxu,&negal1,J_VF,&negal1); }
	  { const R a = self_->rv[__ens_rv_ireynold]* (-self_->rv[__ens_rv_coeff_ratio_viscosity]);
	    nsblas_daxpy(&nuxndg,&a,dyu_visc_dxu,&negal1,J_VF,&negal1); }
#endif
	  /* 
	     jacobien gravite -/- denselm
	  */
#if 0
	  { const R a = self_->rv[__ens_rv_ifroude]*self_->jacelm[jelm] * (-self_->rv[__ens_rv_coeff_ratio_density]);
	    nsblas_daxpy(&nuxndg,&a,gravity_term,&negal1,J_VF,&negal1); }
#endif

	  { I i,j;
	    for (j=0;j<_nu;++j)
	      for (i=0;i<_nu;++i)
		{
		  elm->locmat[(offu+j)*off+(offu+i)] += J_UU[j*_nu+i];
		  elm->locmat[(offv+j)*off+(offu+i)] += J_UV[j*_nu+i];
		  elm->locmat[(offu+j)*off+(offv+i)] += J_VU[j*_nu+i];
		  elm->locmat[(offv+j)*off+(offv+i)] += J_VV[j*_nu+i];
		} } 

	  integrator_assmatelm_full(self_,jelm);	  

	} 
    }
  else
    {
      pWorkelm const elm 	= Global_get_Workelm(self_);
      
      I jelm;
      const I N = self_->nelm;
      for (jelm=0;jelm<N;++jelm)
	{      
	  clr_dvect(total_nddlelmxtotal_nddlelm,elm->locmat);
	  
	  ns_setelm(self_,
		    jelm,
		    x_,
		    xi_,
		    xii_);
	  
	  Fem_setelm(&self_->m_fem,elm);
	  integrator_compute_viscosity_density_nodalelm2(self_,elm->ddlelm);
	  clr_dvect(_nu*_nu,J_VU);
	  clr_dvect(_nu*_nu,J_UV);
	  clr_dvect(_nu*_nu,J_VV);
	  clr_dvect(_nu*_nu,J_UU);
	  { I i,j;
	    for (j=0;j<_nu;++j)
	      for (i=0;i<_nu;++i)
		{
		  elm->locmat[(offu+j)*off+(offu+i)] += J_UU[j*_nu+i];
		  elm->locmat[(offv+j)*off+(offu+i)] += J_UV[j*_nu+i];
		  elm->locmat[(offu+j)*off+(offv+i)] += J_VU[j*_nu+i];
		  elm->locmat[(offv+j)*off+(offv+i)] += J_VV[j*_nu+i];
		} } 
	  integrator_assmatelm_full(self_,jelm);
	} 
    }
}






Err ns_build_system(nsGLOBAL_ST*const Z,
		    cst_pI 	nx_,
		    cst_pR 	x_,
		    cst_pI 	xoff_)
{ 
  //  nsGLOBAL_ST*const Z 			= &global[iproc_];
  pParameters const gParameters		= Global_get_Parameters(Z);
  const L pressure_uncoupled 	= Parameters_getl(gParameters,__ens_linfo_pressure_uncoupled);  
  const L pressure_freematrix 	= Parameters_getl(gParameters,__ens_linfo_pressure_freematrix);
  const L slip 			= Parameters_getl(gParameters,__ens_linfo_slip);
  Err err 				= __eErr_no;

  /* clean jacobian */
  nsGLOBAL_linsys_clear(Z);

  /* build jacobian */
  ns_build_system_stokes_quadrature_free(Z,
					 x_,
					 &x_[xoff_[0]],
					 &x_[2*xoff_[0]]);


  /* dirichlet condition velocity */
  Sparse_dirichlet(Z->sparseStokes->A,
		   Z->sparseStokes->nddlu_dirichlet,
		   Z->sparseStokes->ddlu_dirichlet,
		   Z->dec_ddlu);
  
  Sparse_dirichlet(Z->sparseStokes->A,
		   Z->sparseStokes->nddlv_dirichlet,
		   Z->sparseStokes->ddlv_dirichlet,
		   Z->dec_ddlu+Z->nddlu);


  /* 
     TANT QUE LES CONDITIONS DE DIICHLETS NE CHANGENT PAS SYMBOLIQUEMENT EN TEMPS
     CE QUI EST CI DESSOUS EST FAIT A CHAQUE PAS DE TEMPS POUR RIEN
     SI LES CONDITIONS CHANGENT IL FAUDRA RECALCULER LA MATRICE ECRASEE PRECEDENTE
     SUR LES MATRICES LINEAIRES
  */
  

  /* CONDITION DIRICHLET EN P */
  if (NOT pressure_uncoupled)
    {
      Sparse_dirichlet(Z->sparseStokes->A,
		       Z->sparseStokes->nddlp_dirichlet,
		       Z->sparseStokes->ddlp_dirichlet,Z->dec_ddlp);	
    }
  else if ( (pressure_uncoupled)
	    AND
	    (NOT pressure_freematrix) )
    {
      Sparse_dirichlet(Z->sparseB,
		       Z->sparseStokes->nddlp_dirichlet,
		       Z->sparseStokes->ddlp_dirichlet,
		       Z->dec_ddlp);
    }

  /* conditions de slip*/
  if (slip)
    {	    

      cst_pI 	Ab = Sparse_get_b(Z->sparseStokes->A);
      cst_pI 	Ai = Sparse_get_i(Z->sparseStokes->A);
      pR 	Ax = Sparse_get_ownx(Z->sparseStokes->A);
      
      { I i;
	for (i=0;i<Z->sparseStokes->nddls_dirichlet;++i)
	  {		      
	    const I k = Z->slip_iperm[Z->sparseStokes->ddls_dirichlet[i]];
	    if (NOT k)
	      {
		Monitor_errmsg(Z->iproc,"erreur slip iperm, iddl n est pas frontiere ...");
	      }
	    else
	      {
		const I k2=Z->total_nddl+k-1;
		{ I j;
		  for (j=Ab[k2];j<Ab[k2+1];++j)
		    Ax[j-1] = ((R)0.0);
		  for (j=Ab[k2];j<Ab[k2+1];++j)
		    if (Ai[j-1]!=k2+1)  { Ax[j-1] =  ((R)1.0); break; }
		  if (j>=Ab[k2+1])
		    {
		      printf("wrong algo build system missing diagonal\n");
		    }
		}
	      } 
	  } }
    }


  
  return err;
}







void ns_build_system_slip_ass(nsGLOBAL_ST*const Z,
			      cst_pR 	s_,
			      cst_pI 	soff_,
			      cst_pI 	n_,
			      cst_pI 	numeri_,
			      cst_pI 	m_,
			      cst_pI 	numerj_,
			      cst_pI 	dec_s_,
			      pI     	blank_,
			      pI 	iwork_)
{ 

  I N  		= 0;     

  { I i;
    for (i=0;i<m_[0];++i)
      {
	if (NOT blank_[numerj_[i]])
	  {
	    blank_[numerj_[i]] 		= ++N;
	    iwork_[N-1] 		= numerj_[i];
	  }
      } }  
  
  cst_pI Ab 	= Sparse_get_b(Z->sparseStokes->A);
  cst_pI Ai 	= Sparse_get_i(Z->sparseStokes->A);
  pR Ax 	= Sparse_get_ownx(Z->sparseStokes->A);
      
  { I j,k;
    for (j=0;j<n_[0];++j)
      {
	const I jddl 	= dec_s_[0]+numeri_[j];
	const I n		= Ab[jddl+1]-Ab[jddl];
	cst_pI it 		= &Ai[ Ab[jddl]-1 ];
	pR rt 			= &Ax[ Ab[jddl]-1 ];
	{ I i;
	  for (i=0;i<n;++i)
	    {
	      if ( (k=blank_[it[i]-1]) )
		{
		  Sparse_ass(Z->sparseStokes->A,it[i]-1,jddl,s_[j + (k-1)*soff_[0]]);
		  rt[i] -= s_[j + (k-1)*soff_[0]];			  
		}		      
	    } }	    
      } }
  
  { I i;
    for (i=0;i<N;++i) 
      {
	blank_[iwork_[i]] = 0;
      } }
}






void ns_build_system_slip(nsGLOBAL_ST*const Z,
			  cst_pI slip_iperm_,
			  cst_pI slip_nddlelm_)
{  

  const R f1[6] = {((R)1.0)/((R)3.0),
			((R)0.0),
			
			((R)2.0)/((R)3.0),
			((R)2.0)/((R)3.0),
			
			((R)0.0),
			((R)1.0)/((R)3.0)};
  
  const R f2[9] = {((R)4.0)/((R)15.0),((R)2.0)/((R)15.0),((R)-1.0)/((R)15.0),
			((R)2.0)/((R)15.0),((R)16.0)/((R)15.0),((R)2.0)/((R)15.0),
			((R)-1.0)/((R)15.0),((R)2.0)/((R)15.0),((R)4.0)/((R)15.0)};
  
  
  I 	loci[3];
  I 	locu[6];
  I 	locs[3];
  I 	locc[3];  
  R 	s[18];
  const I N = Z->nelm;
  const I Nvertex= Z->nvertex;
  { I i;
    for (i=0;i<N;++i)
      {
	{ I j;
	  for (j=0;j<3;++j)
	    {
	      if (NOT Z->adj[3*i+j])
		{
		  const I jedge = Z->ddlcnc[_nu*i+3+j]-1-Nvertex;
		  		  
		  loci[0] 	= j;
		  loci[1] 	= 3+j;
		  loci[2] 	= (1+j)%3;
		  
		  locc[0] 	= Z->ddlcnc[_nu*i+loci[0]]-1;
		  locc[1] 	= Z->ddlcnc[_nu*i+loci[1]]-1;
		  locc[2] 	= Z->ddlcnc[_nu*i+loci[2]]-1;
		  
		  locu[0] 	= Z->dec_ddlu + locc[0];
		  locu[1] 	= Z->dec_ddlu + locc[1];
		  locu[2] 	= Z->dec_ddlu + locc[2];
		  locu[3] 	= locu[0]  + Z->nddlu;
		  locu[4] 	= locu[1]  + Z->nddlu;
		  locu[5] 	= locu[2]  + Z->nddlu;

		  if (slip_nddlelm_[0]==6)
		    {
		      locs[0] 	= slip_iperm_[locc[0]]-1;
		      locs[1] 	= slip_iperm_[locc[1]]-1;
		      locs[2] 	= slip_iperm_[locc[2]]-1;	
		      if ( (locs[0]<0) OR (locs[1]<0) OR (locs[2]<0) )
			{
			  Monitor_errmsg(Z->iproc,"integrator_slip_constraint failed\n");
			  
			}
		    }
		  else if (slip_nddlelm_[0]==3)
		    {		  
		      locs[0] = slip_iperm_[locc[0]]-1;
		      locs[1] = slip_iperm_[locc[2]]-1;		     
		      if ( (locs[0]<0) OR (locs[1]<0)  )
			{
			  Monitor_errmsg(Z->iproc,"integrator_slip_constraint failed\n");
			}
		    }
		  if (slip_nddlelm_[0]==6)
		    {
		      clr_dvect(18,s);		    
		      const I n = 9;
		      R a = Z->jacedge[jedge]*Z->normaledge[2*jedge];
		      nsblas_daxpy(&n,&a,f2,&negal1,s,&negal1);
		      a = Z->jacedge[jedge]*Z->normaledge[2*jedge+1];
		      nsblas_daxpy(&n,&a,f2,&negal1,&s[n],&negal1);
		      ns_build_system_slip_ass		(Z,
							 s,
							 &negal3,
							 &negal3,
							 locs,
							 &nu,
							 locu,
							 &Z->total_nddl,
							 Z->blank,Z->iwork);
		    }
		  else
		    {
		      clr_dvect(12,s);		    
		      const I n 	= 6;

#if 0
		      R a 		= Z->jacedge[jedge]*Z->normaledge[2*jedge];
		      a = Z->jacedge[jedge]*Z->normaledge[2*jedge+1];
#endif
		      nsblas_daxpy	(&n,&Z->normaledge[2*jedge],f1,&negal1,s,&negal1);
		      nsblas_daxpy	(&n,&Z->normaledge[2*jedge+1],f1,&negal1,&s[n],&negal1);
		      ns_build_system_slip_ass		(Z,
							 s,
							 &negal2,
							 &negal2,
							 locs,
							 &nu,
							 locu,
							 &Z->total_nddl,
							 Z->blank,Z->iwork);
		    }
		}
	    } }
      } }
}

I line_symbolicu_u( nsGLOBAL_ST*const Z,
			const I 	iddl_,
			pI 		ddl_,
			cst_pI		bpatch_,
			cst_pI		patch_,
			pI		blank_,
			pI 		symbolic_cnc_)
{
  I i,j,N=0;
  for (i=bpatch_[iddl_];i<bpatch_[iddl_+1];++i)
    {
      const I jelm = patch_[i];
      nsGLOBAL_get_ddlcnc_u(Z,jelm,symbolic_cnc_);
      nsGLOBAL_get_ddlcnc_v(Z,jelm,&symbolic_cnc_[_nu]);
      for (j=0;j<_nu*2;++j)
	{
	  if (NOT blank_[symbolic_cnc_[j]])
	    {
	      blank_[symbolic_cnc_[j]] = ++N;
	      ddl_[N-1]=symbolic_cnc_[j];
	    }
	}
    }
  for (i=0;i<N;++i)
    blank_[ddl_[i]]=(I)0;
  return N;
}

I line_symbolicu_up( nsGLOBAL_ST*const Z,
			const I 	iddl_,
			 pI 		ddl_,
			 cst_pI	bpatch_,
			 cst_pI	patch_,
			 pI		blank_,
			 pI 		symbolic_cnc_)
{

  pParameters const gParameters = Global_get_Parameters(Z);
  I N2=0,N=0;
  { I i;
    for (i=bpatch_[iddl_];i<bpatch_[iddl_+1];++i)
      {
	const I jelm = patch_[i];
	nsGLOBAL_get_ddlcnc(Z,jelm,symbolic_cnc_);
	{ I j;
	  for (j=0;j<_total_nddlelm;++j)
	    {
	      if (NOT blank_[symbolic_cnc_[j]])
		{
		  blank_[symbolic_cnc_[j]] = ++N;
		  ddl_[N2+N-1]=symbolic_cnc_[j];		
		}
	    } }      
	if (Parameters_getl(gParameters,__ens_linfo_slip))
	  {
	    { I j;
	      for (j=0;j<_nu;++j)
		{
		  if (Z->ddlcod[symbolic_cnc_[j]-Z->dec_ddlu]!=Parameters_getl(gParameters,__ens_iinfo_noboundary_vcod))
		    {
		      if (NOT blank_[Z->slip_iperm[symbolic_cnc_[j]-Z->dec_ddlu]-1+Z->total_nddl])
			{
			  blank_[Z->slip_iperm[symbolic_cnc_[j]-Z->dec_ddlu]-1+Z->total_nddl]	= ++N;
			  ddl_[N-1] = (Z->slip_iperm[symbolic_cnc_[j]-Z->dec_ddlu]-1)+Z->total_nddl;		      
			}
		    }
		} }
	  }
      } }
  { I i;
    for (i=0;i<N;++i)
      {
	blank_[ddl_[N2+i]]=(I)0;
      } }
  return N+N2;
}



I line_symbolics( nsGLOBAL_ST*const Z,
		  const I 	iddl_,
		  pI 		ddl_,
		  cst_pI	bpatch_,
		  cst_pI	patch_,
		  pI		blank_,
		  pI 		symbolic_cnc_)
{

  pParameters const gParameters 	= Global_get_Parameters(Z);
  I i,j,N=0;
  for (i=bpatch_[iddl_];i<bpatch_[iddl_+1];++i)
    {
      const I jelm = patch_[i];
      nsGLOBAL_get_ddlcnc_u(Z,jelm,&symbolic_cnc_[0]);
      nsGLOBAL_get_ddlcnc_v(Z,jelm,&symbolic_cnc_[_nu]);
      for (j=0;j<_nu;++j)
	{
	  if (Z->ddlcod[symbolic_cnc_[j]-Z->dec_ddlu]!=Parameters_geti(gParameters,__ens_iinfo_noboundary_vcod))
	    {
	      if (NOT blank_[symbolic_cnc_[j]])
		{
		  blank_[symbolic_cnc_[j]] 	= ++N;
		  ddl_[N-1]   		= symbolic_cnc_[j];
		}
	    }
	}
      for (j=0;j<_nu;++j)
	{
	  if (Z->ddlcod[symbolic_cnc_[_nu+j]-Z->dec_ddlu-Z->nddlu]!=Parameters_geti(gParameters,__ens_iinfo_noboundary_vcod))
	    {
	      if (NOT blank_[symbolic_cnc_[_nu+j]])
		{

		  blank_[symbolic_cnc_[_nu+j]] = ++N;
		  ddl_[N-1] = symbolic_cnc_[_nu+j];
		}
	    }
	}
    }
  for (i=0;i<N;++i)
    {
      blank_[ddl_[i]]=(I)0;
    }
  /*
    on rajoute un element diagonal pour pouvoir annuler la contrainte
  */
#if 0
  ddl_[N] = iddl_;
  ++N;
#endif
  return N;
}

I line_symbolicp( nsGLOBAL_ST*const Z,
		      const I 	iddl_,
		      pI 		ddl_,
		      cst_pI		bpatch_,
		      cst_pI		patch_,
		      pI		blank_,
		      pI 		symbolic_cnc_)
{
  I N=0;
  { I i;
    for (i=bpatch_[iddl_];i<bpatch_[iddl_+1];++i)
      {
	const I jelm = patch_[i];
	nsGLOBAL_get_ddlcnc_u	(Z,jelm,&symbolic_cnc_[0]);
	nsGLOBAL_get_ddlcnc_v	(Z,jelm,&symbolic_cnc_[_nu]);
	nsGLOBAL_get_ddlcnc_p	(Z,jelm,&symbolic_cnc_[_nu*2]);
	{ I j;
	  for (j=0;j<_nuxdim+_np;++j)
	    {
	      if (NOT blank_[symbolic_cnc_[j]])
		{
		  blank_[symbolic_cnc_[j]] = ++N;
		  ddl_[N-1]=symbolic_cnc_[j];
		}
	    } }
      } }
  { I i;
    for (i=0;i<N;++i)
      {
	blank_[ddl_[i]]=(I)0;
      } }
  return N;
}

Err sparse_A_init(nsGLOBAL_ST*const Z,
		  cst_eKindSystem 	kind_)
{

  pParameters const gParameters 	= Global_get_Parameters(Z);
  const L  slip 		= Parameters_getl(gParameters,__ens_linfo_slip);
  const I nedge_boundary 	= Z->nedge_boundary;
  const I nelm 		= Z->nelm;
  /*  const I nvertex 		= Z->nvertex;*/
  Err 	err 		= __eErr_no;
  I 	An 	= ((I)0);
  I 	Am 	= ((I)0);
  I 	Ancoeff = ((I)0);
  pI 		Ab 	= NULL;
  pI 		Ai 	= NULL;

  if (slip)
    {
      Z->slip_iperm = (pI)calloc(Z->nddlu,sizeof(I));
      I iddl,jddl=0;
      for (iddl=0;iddl<Z->nddlu;++iddl)
	{      
	  if (Z->ddlcod[iddl]!=Parameters_geti(gParameters,__ens_iinfo_noboundary_vcod))
	    {
	      Z->slip_iperm[iddl] = jddl+1;
	      ++jddl;
	    }
	}
    }
  else
    {
      Z->slip_iperm= NULL;
    }
  I	decu=0,decv=0,decp=0;
  switch(kind_)
    {
    case __eKindSystem_up:
      {
	An 		= dim*Z->nddlu+Z->nddlp;
	decu      	= 0;
	decv      	= Z->nddlu;
	decp      	= Z->nddlu*2;
	break;
      }
    case __eKindSystem_u:
      {
	An 		= dim*Z->nddlu;
	decu    	= 0;
	decv    	= Z->nddlu;
	break;
      }

    case __eKindSystem_ALL:
    case __eKindSystem_ERROR:
      {
	err = __eErr_switch;
	Monitor_errmsg(Z->iproc,"sparse_A_init:switch failed on __eKindSystem");
	break;
      }
    }

  if (err)
    return err;

  if (slip)
    {
      An += ( ensBASIS_degree(_shape_u)-1 );
    }
  
  Am		= An;
  Ab 		= (pI)malloc((An+1)*sizeof(I));
  if (slip)
    {
      Ai 	= (pI)malloc((_total_nddlelm*_total_nddlelm*nelm+2*nedge_boundary*(12+1))*sizeof(I));     /*+1 pour l'element diagonal du slip */
    }
  else
    {
      Ai 	= (pI)malloc(_total_nddlelm*_total_nddlelm*nelm*sizeof(I));
    }
 
  Ab[0] 	= 0;
  Ancoeff 	= 0;
  I symbolic_cnc[256];
  switch(kind_)
    {
    case __eKindSystem_up:
      {
	{ I iddl;
	  for (iddl=0;iddl<Z->nddlu;++iddl)
	    {
	      const I nline = line_symbolicu_up(Z,iddl,&Ai[Ancoeff],Z->bpatch,Z->patch,Z->blank,symbolic_cnc);
	      Ancoeff +=nline;
	      Ab[decu+iddl+1] = Ab[decu+iddl] + nline;
	    }  }
	{ I iddl;
	  for (iddl=0;iddl<Z->nddlu;++iddl)
	    {
	      const I nline = line_symbolicu_up(Z,iddl,&Ai[Ancoeff],Z->bpatch,Z->patch,Z->blank,symbolic_cnc);
	      Ancoeff +=nline;
	      Ab[decv+iddl+1] = Ab[decv+iddl] + nline;
	    }    }
	{ I iddl;
	  for (iddl=0;iddl<Z->nddlp;++iddl)
	    {
	      const I nline = line_symbolicp(Z,iddl,&Ai[Ancoeff],Z->bpatch,Z->patch,Z->blank,symbolic_cnc);
	      Ancoeff +=nline;
	      Ab[decp+iddl+1] = Ab[decp+iddl] + nline;
	    } }
	break;
      }
    case __eKindSystem_u:
      {
	{ I iddl;
	  for (iddl=0;iddl<Z->nddlu;++iddl)
	    {
	      const I nline = line_symbolicu_u(Z,iddl,&Ai[Ancoeff],Z->bpatch,Z->patch,Z->blank,symbolic_cnc);
	      Ancoeff +=nline;
	      Ab[decu+iddl+1] = Ab[decu+iddl] + nline;
	    } }
	{ I iddl;
	  for (iddl=0;iddl<Z->nddlu;++iddl)
	    {
	      const I nline = line_symbolicu_u(Z,iddl,&Ai[Ancoeff],Z->bpatch,Z->patch,Z->blank,symbolic_cnc);
	      Ancoeff +=nline;
	      Ab[decv+iddl+1] = Ab[decv+iddl] + nline;
	    }  }
	break;
      }
    case __eKindSystem_ALL:
    case __eKindSystem_ERROR:
      {
	err = __eErr_switch;
	Monitor_errmsg(Z->iproc,"sparse_A_init:switch failed on __eKindSystem");
	break;
      }
    }
  if (slip)
    {
      I jddl = 0;
      { I iddl;
	for (iddl=0;iddl<Z->nddlu;++iddl)
	  {
	    if (Z->ddlcod[iddl]!=Parameters_geti(gParameters,__ens_iinfo_noboundary_vcod))
	      {
		const I nline = line_symbolics(Z,iddl,&Ai[Ancoeff],Z->bpatch,Z->patch,Z->blank,symbolic_cnc);
		Ancoeff +=nline;
		/* on ajoute un element diagonal pour pouvoir annuler la condition */
		/* s'il n'y est pas normalement ca doit guEULER quand on vuet appliquer du dirichlet nul sur le slip
		   puisqu'on veut rajouter un element sur la diagonale */
		Ai[Ancoeff]=Z->total_nddl+jddl;/* la notation est C pas fortran, la version fortran vient apres*/
		++Ancoeff;
		Ab[Z->total_nddl+jddl+1] = Ab[Z->total_nddl+jddl] + (nline+1);
		++jddl;
	      }
	  } }
      
      Z->nddl_slip	= An-Z->total_nddl;
    }
  else
    {
      Z->nddl_slip = 0;
    }
  Ai 				= (pI)realloc(Ai,Ancoeff*sizeof(I));  


  
  Z->sparseStokes->A = Sparse_build(An,Am,Ancoeff,Ab,Ai,NULL);
#if 0
  Z->sparse_stokes.A.own_i   	= Ai;
  Z->sparse_stokes.A.own_x   	= (pR)malloc(Ancoeff*sizeof(R));
  Z->sparse_stokes.A.own_b   	= Ab;
  Z->sparse_stokes.A.n    	= An;
  Z->sparse_stokes.A.m    	= Am;  
  Z->sparse_stokes.A.nc   	= Ancoeff;
  Z->sparse_stokes.A.b 		= Z->sparse_stokes.A.own_b;
  Z->sparse_stokes.A.x 		= Z->sparse_stokes.A.own_x;
  Z->sparse_stokes.A.i 		= Z->sparse_stokes.A.own_i;
#endif
  Sparse_sort(Z->sparseStokes->A);
  return __eErr_no;
}









Err sparse_B_init(pGlobal self_)
{  
  const I nelm 		= self_->nelm;
  const I Bn 		= self_->nddlp;
  const I Bm 		= 2*self_->nddlu;
  pI own_b		= (pI)malloc((Bn+1)*sizeof(I));
  pI own_i		= (pI)malloc(_nu*_nu*nelm*sizeof(I));
  own_b[0] 		= 0;
  I Bnc 		= 0;
  
#if 0
  self_->sparse_B.n		= self_->nddlp;
  self_->sparse_B.m    		= 2*self_->nddlu;
  self_->sparse_B.own_b		= (pI)malloc((self_->sparse_B.n+1)*sizeof(I));
  self_->sparse_B.own_i		= (pI)malloc(_nu*_nu*nelm*sizeof(I));
  self_->sparse_B.own_b[0] 	= 0;
  self_->sparse_B.nc 		= 0;

  I symbolic_cnc[128];
  I iddl;
  for (iddl=0;iddl<self_->nddlp;++iddl)
    {
      const I nline = line_symbolicp(Z,iddl,&self_->sparse_B.own_i[self_->sparse_B.nc],self_->bpatch,self_->patch,self_->blank,symbolic_cnc);
      self_->sparse_B.nc +=nline;
      self_->sparse_B.own_b[iddl+1] = self_->sparse_B.own_b[iddl] + nline;
    }
#endif

  I symbolic_cnc[128];
  I iddl;
  for (iddl=0;iddl<self_->nddlp;++iddl)
    {
      const I nline = line_symbolicp(self_,iddl,&own_i[Bnc],self_->bpatch,self_->patch,self_->blank,symbolic_cnc);
      Bnc +=nline;
      own_b[iddl+1] = own_b[iddl] + nline;
    }

  printf("realll "ifmt"\n",_nu*_nu*nelm-Bnc);
  own_i		= (pI)realloc(own_i,Bnc*sizeof(I));  


#if 0
  own_x		= (pR)malloc(Bnc*sizeof(R));  
#endif

  self_->sparseB = Sparse_build(Bn,Bm,Bnc,own_b,own_i,NULL);
#if 0
  self_->sparse_B.b 		= self_->sparse_B.own_b;
  self_->sparse_B.x 		= self_->sparse_B.own_x;
  self_->sparse_B.i 		= self_->sparse_B.own_i;
#endif
  return __eErr_no;
}


#ifdef WAS_VALID_ROUTINE
void LinsysVectors_update_previous_x(pGlobal 		self_,
				     pLinsysVectors 	v_)
{
#ifndef NDEBUG
  DebugVerif(self_);
  DebugVerif(v_);
#endif
  // nsGLOBAL_ST*const Z 		= &global[iproc_];
  pParametersReadOnly 	const 	gParameters	= GlobalReadOnly_get_Parameters(self_);
  pTimeReadOnly 	const 	gTimeInfo	= GlobalReadOnly_get_Time(self_);
  
  const I vN = LinsysVectors_get_n(v_);
  cst_pR x = (cst_pR)LinsysVectors_get_x(v_);
  
  switch(Parameters_get_eTransientScheme(gParameters,__eKindEquation_VELOCITY])
    {
    case __eTransientMethod_NO:
      {
	break;
      }
    case __eTransientMethod_EULER:
      {
	static const I istep = 1;
	cst_pR xi = (cst_pR)LinsysVectors_get_friend_xi(v_,&istep);
	nsblas_dcopy(&vN,x,&negal1,xi,&negal1);	
	break;
      }
    case __eTransientMethod_IMR:
      {
	static const I istep = 1;
	cst_pR xi = (cst_pR)LinsysVectors_get_friend_xi(v_,&istep);
	
	const I N  = _dim*self_->nddlu;
	const I N2 = vN-N;
	/* 
	   mise a jour imr
	*/
	nsblas_dscal	(&N,&mregal1,xi,&negal1);
	nsblas_daxpy	(&N,&regal2,x,&negal1,xi,&negal1);	
	nsblas_dcopy	(&N2,&x[N],&negal1,&xi[N],&negal1);
	break;
      }
    case __eTransientMethod_GEAREULER:
      {
	static const I istep1 = 1;
	cst_pR xi = (cst_pR)LinsysVectors_get_friend_xi(v_,&istep1);

	static const I istep = 2;
	cst_pR xii = (cst_pR)LinsysVectors_get_friend_xi(v_,&istep);

	nsblas_dcopy(&vN,xi,&negal1,xii,&negal1);
	nsblas_dcopy(&vN,x,&negal1,xi,&negal1);		       	
	break;
      }
    case __eTransientMethod_GEARIMR:
      {
	static const I istep1 = 1;
	cst_pR xi = (cst_pR)LinsysVectors_get_friend_xi(v_,&istep1);

	static const I istep = 2;
	cst_pR xii = (cst_pR)LinsysVectors_get_friend_xi(v_,&istep);

	if (gTimeInfo->itime>0)
	  {
	    nsblas_dcopy(&vN,xi,&negal1,xii,&negal1);
	    nsblas_dcopy(&vN,x,&negal1,xi,&negal1);		       	
	  }
	else
	  {
	    nsblas_dcopy(&vN,xi,&negal1,xii,&negal1);
	    const I N  =  _dim*self_->nddlu;
	    const I N2 = vN-N;
	    /* 
	       mise a jour imr
	    */
	    nsblas_dscal	(&N,&mregal1,xi,&negal1);
	    nsblas_daxpy	(&N,&regal2,x,&negal1,xi,&negal1);	
	    nsblas_dcopy	(&N2,&x[N],&negal1,&xi[N],&negal1);
	  }
	break;
      }
    case __eTransientMethod_TRAPEZE:
      {
	Monitor_errmsg(self_->iproc,"__eTransientMethod_TRAPEZE not yet");
	break;
      }
    case __eTransientMethod_ERROR:
    case __eTransientMethod_ALL:
      {
	Monitor_errmsg(self_->iproc,"nsLINSYS_VECTORS_update_previous_x:switch failed on __eTransientMethod \n");
	break;
      }
    }	      
}
#endif


