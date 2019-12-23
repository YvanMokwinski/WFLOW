#ifdef __cplusplus
extern "C"
{
#endif

#define GmfStrSiz 1024
#define GmfMaxTyp 20
#define GmfMaxKwd 79
#define GmfMshVer 1
#define GmfRead 1
#define GmfWrite 2
#define GmfSca 1
#define GmfVec 2
#define GmfSymMat 3
#define GmfMat 4
#define GmfFloat 1
#define GmfDouble 2

enum GmfKwdCod
{
  GmfReserved1, \
  GmfVersionFormatted, \
  GmfReserved2, \
  GmfDimension, \
  GmfVertices, \
  GmfEdges, \
  GmfTriangles, \
  GmfQuadrilaterals, \
  GmfTetrahedra, \
  GmfPentahedra, \
  GmfHexahedra, \
  GmfWedges, \
  GmfPyramids, \
  GmfCorners, \
  GmfRidges, \
  GmfRequiredVertices, \
  GmfRequiredEdges, \
  GmfRequiredTriangles, \
  GmfRequiredQuadrilaterals, \
  GmfTangentAtEdgeVertices, \
  GmfNormalAtVertices, \
  GmfNormalAtTriangleVertices, \
  GmfNormalAtQuadrilateralVertices, \
  GmfAngleOfCornerBound, \
  GmfTrianglesP2, \
  GmfTrianglesP3, \
  GmfTrianglesP4, \
  GmfQuadrilateralsP2, \
  GmfQuadrilateralsP3, \
  GmfQuadrilateralsP4, \
  GmfTetrahedraP2, \
  GmfTetrahedraP3, \
  GmfTetrahedraP4, \
  GmfHexahedraP2, \
  GmfHexahedraP3, \
  GmfHexahedraP4, \
  GmfReserved17, \
  GmfReserved18, \
  GmfReserved19, \
  GmfReserved20, \
  GmfReserved21, \
  GmfReserved22, \
  GmfReserved23, \
  GmfReserved24, \
  GmfReserved25, \
  GmfReserved26, \
  GmfReserved27, \
  GmfReserved28, \
  GmfReserved29, \
  GmfReserved30, \
  GmfBoundingBox, \
  GmfReserved31, \
  GmfReserved32, \
  GmfReserved33, \
  GmfEnd, \
  GmfReserved34, \
  GmfReserved35, \
  GmfReserved36, \
  GmfReserved37, \
  GmfTangents, \
  GmfNormals, \
  GmfTangentAtVertices, \
  GmfSolAtVertices, \
  GmfSolAtEdges, \
  GmfSolAtTriangles, \
  GmfSolAtQuadrilaterals, \
  GmfSolAtTetrahedra, \
  GmfSolAtPentahedra, \
  GmfSolAtHexahedra, \
  GmfDSolAtVertices, \
  GmfISolAtVertices, \
  GmfISolAtEdges, \
  GmfISolAtTriangles, \
  GmfISolAtQuadrilaterals, \
  GmfISolAtTetrahedra, \
  GmfISolAtPentahedra, \
  GmfISolAtHexahedra, \
  GmfIterations, \
  GmfTime, \
  GmfReserved38
};


extern int 	GmfOpenMesh(char *, int, ...);
extern int 	GmfCloseMesh(int);
extern int 	GmfStatKwd(int, int, ...);
extern int 	GmfGotoKwd(int, int);
extern int 	GmfSetKwd(int, int, ...);
extern void 	GmfGetLin(int, int, ...);
extern void 	GmfSetLin(int, int, ...);

#if defined(F77_NO_UNDER_SCORE)
#define call(x) x
#else
#define call(x) x ## _
#endif

#ifdef TRANSMESH

extern char *GmfKwdFmt[ GmfMaxKwd + 1 ][4];
extern int GmfCpyLin(int, int, int);

#endif

#ifdef __cplusplus
};
#endif
