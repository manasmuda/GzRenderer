#include    "gz.h"
#ifndef GZRENDER_
#define GZRENDER_


/* Camera defaults */
#define DEFAULT_FOV     35.0
#define DEFAULT_IM_Z    (-10.0)  /* world coords for image plane origin */
#define DEFAULT_IM_Y    (5.0)    /* default look-at point = 0,0,0 */
#define DEFAULT_IM_X    (-10.0)

#define DEFAULT_AMBIENT {0.1, 0.1, 0.1}
#define DEFAULT_DIFFUSE {0.7, 0.6, 0.5}
#define DEFAULT_SPECULAR    {0.2, 0.3, 0.4}
#define DEFAULT_SPEC        32

#define MATLEVELS   100     /* how many matrix pushes allowed */
#define MAX_LIGHTS  10      /* how many lights allowed */

class Vertex;
class Tris;

class GzRender {         /* define a renderer */


public:
    unsigned short  xres;
    unsigned short  yres;
    GzPixel* pixelbuffer;       /* frame buffer array */
    char* framebuffer;

    GzCamera        m_camera;
    short           matlevel;           /* top of stack - current xform */
    GzMatrix        Ximage[MATLEVELS];  /* stack of xforms (Xsm) */
    GzMatrix        Xnorm[MATLEVELS];   /* xforms for norms (Xim) */
    GzMatrix        Xsp;                /* NDC to screen (pers-to-screen) */
    GzColor     flatcolor;          /* color state for flat shaded triangles */
    int         interp_mode;
    int         numlights;
    GzLight     lights[MAX_LIGHTS];
    GzLight     ambientlight;
    GzColor     Ka, Kd, Ks;
    float           spec;       /* specular power */
    GzTexture       tex_fun;    /* tex_fun(float u, float v, GzColor color) */

    float xOffset, yOffset;

    // Constructors
    GzRender(int xRes, int yRes);
    GzRender();
    ~GzRender();

public:
    // Function declaration
    void CopyRenderer(GzRender* rend);

    // HW1: Display methods
    int GzDefault();
    int GzBeginRender();
    int GzPut(int i, int j, GzIntensity r, GzIntensity g, GzIntensity b, GzIntensity a, GzDepth z);
    int GzGet(int i, int j, GzIntensity* r, GzIntensity* g, GzIntensity* b, GzIntensity* a, GzDepth* z);

    int GzFlushDisplay2File(FILE* outfile);
    int GzFlushDisplay2FrameBuffer();

    // HW2: Render methods
    int GzPutAttribute(int numAttributes, GzToken* nameList, GzPointer* valueList);
    int GzPutTriangle(int numParts, GzToken* nameList, GzPointer* valueList);

    // HW3
    int GzDefaultCamera();
    int GzPutCamera(GzCamera camera);
    int GzPushMatrix(GzMatrix matrix);
    int GzPopMatrix();

    // Extra methods: NOT part of API - just for general assistance */
    inline int ARRAY(int x, int y) { return (x + y * xres); }   /* simplify fbuf indexing */
    inline short ctoi(float color) { return(short)((int)(color * ((1 << 12) - 1))); }        /* convert float color to GzIntensity short */
    inline bool IsUVinRange(GzTextureIndex uv) { return uv[0] >= 0 && uv[1] >= 0 && uv[0] <= 1 && uv[1] <= 1; }

    // Object Translation
    int GzRotXMat(float degree, GzMatrix mat);
    int GzRotYMat(float degree, GzMatrix mat);
    int GzRotZMat(float degree, GzMatrix mat);
    int GzTrxMat(GzCoord translate, GzMatrix mat);
    int GzScaleMat(GzCoord scale, GzMatrix mat);

private:
    GzIntensity ClampGzIntensity(GzIntensity x);
    char DownsizeGzIntesity(GzIntensity x);

    void SetPixelColor(GzCoord pixelPos, GzColor pixelColor);

    void DoScanLineRaterization(Tris& tris);
    void FormLinesAndDoRaterization(float y1, float y2, GzCoord v, float slopex1, float slopez1, float slopex2, float slopez2);
    void DoLineRasterization(GzCoord v1, GzCoord v2);

    void DoScanLineRaterizationGourardShading(Tris& tris);
    void FormLinesAndDoRaterizationGourardShading(float y1, float y2, Vertex vc, Vertex voc, Vertex vm, float slopex1, float slopez1, float slopex2, float slopez2, bool useTexture);
    void DoLineRasterizationGourardShading(Vertex v1, Vertex v2, bool useTexture);
    void InterpolateColor(GzCoord start, GzCoord end, GzColor startColor, GzColor endColor, GzCoord interpolatePoint, GzColor interpolateColor);

    void DoScanLineRaterizationPhongShading(Tris& tris);
    void FormLinesAndDoRaterizationPhongShading(float y1, float y2, Vertex vc, Vertex voc, Vertex vm, float slopex1, float slopez1, float slopex2, float slopez2);
    void DoLineRasterizationPhongShading(Vertex v1, Vertex v2);
    void InterpolateNormal(GzCoord start, GzCoord end, GzCoord startNormal, GzCoord endNormal, GzCoord interpolatePoint, GzCoord interpolateNormal);

    void InterpolateUV(GzCoord start, GzCoord end, GzTextureIndex startUV, GzTextureIndex endUV, GzCoord interpolatePoint, GzTextureIndex interpolateUV);

    void DoCrossProduct(GzCoord v1, GzCoord v2, GzCoord rv);
    void MultiplyVector(GzCoord v1, float scale, GzCoord rv);
    void MultiplyVectors(GzCoord v1, GzCoord v2, GzCoord rv);
    void AddVectors(GzCoord v1, GzCoord v2, GzCoord rv);
    void SubtractVectors(GzCoord v1, GzCoord v2, GzCoord rv);
    float GetVectorsDistance(GzCoord v1, GzCoord v2);
    void GetUnitVector(GzCoord vector, GzCoord unitvector);
    float GetVectorMagnitude(GzCoord v1);
    float DoDotProduct(GzCoord v1, GzCoord v2);
    void ConvertMatrixToUnitaryRotation(GzMatrix matrix);
    void Convert4dCoordTo3dCoord(GzMatrix matrix, int matrixrow, GzCoord coord);

    void MultiplyGzMatrix(GzMatrix m1, GzMatrix m2, GzMatrix rm);
    void CreateXspMatrix(GzMatrix xsp);
    void CreateXpiMatrix(GzCamera camera, GzMatrix xpi);
    void CreateXiwMatrix(GzCamera camera, GzMatrix xiw);

    void GetTransformedVertices(GzMatrix verts, GzMatrix transformedVerts);
    void GetTransformedNormals(GzMatrix norms, GzMatrix transformedNorms);

    void ShadePixel(GzCoord norm, GzColor Ks, GzColor Kd, GzColor Ka, GzColor color);
    void ShadePixel(GzCoord norm, GzColor color);

    float Clamp(float org, float min, float max) {
        if (org < min) {
            org = min;
        }
        else if (org > max) {
            org = max;
        }
        return org;
    }

};

class Vertex {
public:
    GzCoord pos, norm;
    GzTextureIndex uv;
    GzColor shadedColor;
    GzColor textureColor;

public:
    Vertex(GzCoord pos, GzCoord norm, GzTextureIndex uv);
    Vertex(GzCoord pos, GzTextureIndex uv);
    ~Vertex();

    void SetShadedColor(GzColor shadedColor);
};

class Tris
{

public:
    Vertex *v1, *v2, *v3;
    Vertex *fn; //flat shading vertex
    float slope12X, slope12Z, slope23X, slope23Z, slope13X, slope13Z;
    float y1, y21, y23, y3;

public:
    Tris::Tris(Vertex* v1, Vertex* v2, Vertex* v3);
    ~Tris();

private:
    float CalculateSlopeX(GzCoord a, GzCoord b) {
        return (b[0] - a[0]) / (b[1] - a[1]);
    }
    float CalculateSlopeZ(GzCoord a, GzCoord b) {
        return (b[2] - a[2]) / (b[1] - a[1]);
    }

};



#endif