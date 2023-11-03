#include "stdafx.h"
#include "stdio.h"
#include "math.h"
#include "Gz.h"
#include "rend.h"
#include <algorithm>
#include <cstring>
#include <cmath>

#define PI (float) 3.14159265358979323846


int GzRender::GzRotXMat(float degree, GzMatrix mat)
{
	float rad = degree * PI / 180;
	GzMatrix rotationMatrix =
	{
	1,0,0,0,
	0,cos(rad),-1 * sin(rad),0,
	0,sin(rad),cos(rad),0,
	0,0,0,1
	};
	MultiplyGzMatrix(rotationMatrix, mat, mat);
	return GZ_SUCCESS;
}

int GzRender::GzRotYMat(float degree, GzMatrix mat)
{
	float rad = degree * PI / 180;
	GzMatrix rotationMatrix =
	{
	cos(rad),0,sin(rad),0,
	0,1,0,0,
	-1 * sin(rad),0,cos(rad),0,
	0,0,0,1
	};
	MultiplyGzMatrix(rotationMatrix, mat, mat);
	return GZ_SUCCESS;
}

int GzRender::GzRotZMat(float degree, GzMatrix mat)
{
	float rad = degree * PI / 180;
	GzMatrix rotationMatrix =
	{
	cos(rad),-1 * sin(rad),0,0,
	sin(rad),cos(rad),0,0,
	0,0,1,0,
	0,0,0,1
	};
	MultiplyGzMatrix(rotationMatrix, mat, mat);
	return GZ_SUCCESS;
}

int GzRender::GzTrxMat(GzCoord translate, GzMatrix mat)
{
	GzMatrix translateMatrix =
	{
	1,0,0,translate[0],
	0,1,0,translate[1],
	0,0,1,translate[2],
	0,0,0,1
	};
	MultiplyGzMatrix(translateMatrix, mat, mat);
	return GZ_SUCCESS;
}


int GzRender::GzScaleMat(GzCoord scale, GzMatrix mat)
{
	GzMatrix scaleMatrix =
	{
	scale[0],0,0,0,
	0,scale[1],0,0,
	0,0,scale[2],0,
	0,0,0,1
	};
	MultiplyGzMatrix(scaleMatrix, mat, mat);
	return GZ_SUCCESS;
}

void GzRender::MultiplyVector(GzCoord v1, float scale, GzCoord rv) {
	rv[0] = v1[0] * scale;
	rv[1] = v1[1] * scale;
	rv[2] = v1[2] * scale;
}

void GzRender::MultiplyVectors(GzCoord v1, GzCoord v2, GzCoord rv) {
	rv[0] = v1[0] * v2[0];
	rv[1] = v1[1] * v2[1];
	rv[2] = v1[2] * v2[2];
}

void GzRender::AddVectors(GzCoord v1, GzCoord v2, GzCoord rv) {
	rv[0] = v1[0] + v2[0];
	rv[1] = v1[1] + v2[1];
	rv[2] = v1[2] + v2[2];
}

void GzRender::SubtractVectors(GzCoord v1, GzCoord v2, GzCoord rv) {
	rv[0] = v1[0] - v2[0];
	rv[1] = v1[1] - v2[1];
	rv[2] = v1[2] - v2[2];
}

float GzRender::GetVectorsDistance(GzCoord v1, GzCoord v2) {
	return sqrt(pow(v1[0] - v2[0], 2) + pow(v1[1] - v2[1], 2) + pow(v1[2] - v2[2], 2));
}

void GzRender::GetUnitVector(GzCoord vector, GzCoord unitvector) {
	float magnitude = GetVectorMagnitude(vector);
	unitvector[0] = vector[0] / magnitude;
	unitvector[1] = vector[1] / magnitude;
	unitvector[2] = vector[2] / magnitude;
}

float GzRender::GetVectorMagnitude(GzCoord v1) {
	return sqrt(pow(v1[0], 2) + pow(v1[1], 2) + pow(v1[2], 2));
}

void GzRender::DoCrossProduct(GzCoord v1, GzCoord v2, GzCoord rv) {
	rv[0] = v1[1] * v2[2] - v1[2] * v2[1];
	rv[1] = v1[2] * v2[0] - v1[0] * v2[2];
	rv[2] = v1[0] * v2[1] - v1[1] * v2[0];
}

float GzRender::DoDotProduct(GzCoord v1, GzCoord v2) {
	return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

void GzRender::ConvertMatrixToUnitaryRotation(GzMatrix matrix) {
	matrix[0][3] = 0;
	matrix[1][3] = 0;
	matrix[2][3] = 0;

	float magnitude = sqrt(pow(matrix[0][0], 2) + pow(matrix[0][1], 2) + pow(matrix[0][2], 2));

	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			matrix[i][j] = matrix[i][j] / magnitude;
		}
	}

	matrix[3][0] = 0;
	matrix[3][1] = 0;
	matrix[3][2] = 0;
	matrix[3][3] = 1;
}

void GzRender::Convert4dCoordTo3dCoord(GzMatrix matrix, int matrixrow, GzCoord coord) {
	coord[0] = matrix[matrixrow][0] / matrix[matrixrow][3];
	coord[1] = matrix[matrixrow][1] / matrix[matrixrow][3];
	coord[2] = matrix[matrixrow][2] / matrix[matrixrow][3];
}

GzRender::GzRender() {

}

GzRender::GzRender(int xRes, int yRes)
{
	xres = xRes;
	yres = yRes;
	matlevel = -1;
	numlights = 0;
	pixelbuffer = new GzPixel[xres * yres];
	framebuffer = new char[3 * xres * yres];

	CreateXspMatrix(Xsp);
	GzDefaultCamera();
}

void GzRender::CopyRenderer(GzRender* rend) {
	xres = rend->xres;
	yres = rend->yres;

	pixelbuffer = new GzPixel[xres * yres];
	framebuffer = new char[3 * xres * yres];

	m_camera = rend->m_camera;

	matlevel = rend->matlevel;
	std::memcpy(Ximage, rend->Ximage, (matlevel+1) * sizeof(GzMatrix));
	std::memcpy(Xnorm, rend->Xnorm, (matlevel+1) * sizeof(GzMatrix));

	std::memcpy(Xsp, rend->Xsp, sizeof(GzMatrix));
	std::memcpy(flatcolor, rend->flatcolor, sizeof(GzColor));

	interp_mode = rend->interp_mode;

	numlights = rend->numlights;
	std::memcpy(lights, rend->lights, numlights * sizeof(GzLight));
	ambientlight = rend->ambientlight;

	std::memcpy(Ka, rend->Ka, sizeof(GzColor));
	std::memcpy(Kd, rend->Kd, sizeof(GzColor));
	std::memcpy(Ks, rend->Ks, sizeof(GzColor));
	spec = rend->spec;

	tex_fun = rend->tex_fun;
}

GzRender::~GzRender()
{
	delete[] pixelbuffer;
	delete[] framebuffer;
}

int GzRender::GzDefault()
{
	for (int j = 0; j < yres; j++)
	{
		for (int i = 0; i < xres; i++)
		{
			int index = ARRAY(i, j);
			GzPixel* curpixel = &pixelbuffer[index];
			curpixel->red = 3000;
			curpixel->green = 3000;
			curpixel->blue = 3000;
			pixelbuffer[index].alpha = 4000;
			pixelbuffer[index].z = MAXINT;
		}
	}
	return GZ_SUCCESS;
}

int GzRender::GzDefaultCamera() {
	GzCamera defaultCamera;
	defaultCamera.FOV = DEFAULT_FOV;
	defaultCamera.lookat[0] = 0;
	defaultCamera.lookat[1] = 0;
	defaultCamera.lookat[2] = 0;
	defaultCamera.position[0] = DEFAULT_IM_X;
	defaultCamera.position[1] = DEFAULT_IM_Y;
	defaultCamera.position[2] = DEFAULT_IM_Z;
	defaultCamera.worldup[0] = 0;
	defaultCamera.worldup[1] = 1;
	defaultCamera.worldup[2] = 0;
	GzPutCamera(defaultCamera);
	return GZ_SUCCESS;
}

int GzRender::GzBeginRender()
{
	matlevel = -1;
	GzPushMatrix(Xsp);
	GzPushMatrix(m_camera.Xpi);
	GzPushMatrix(m_camera.Xiw);

	return GZ_SUCCESS;
}

int GzRender::GzPutCamera(GzCamera camera)
{
	GzMatrix xpi, xiw;
	CreateXpiMatrix(camera, xpi);
	CreateXiwMatrix(camera, xiw);
	std::memcpy(camera.Xpi, xpi, sizeof(GzMatrix));
	std::memcpy(camera.Xiw, xiw, sizeof(GzMatrix));
	this->m_camera = camera;
	return GZ_SUCCESS;
}

void GzRender::CreateXspMatrix(GzMatrix xsp) {
	float halfX = xres / 2;
	float halfY = yres / 2;
	GzMatrix result =
	{
	halfX, 0, 0, halfX,
	0, -1 * halfY, 0, halfY,
	0, 0, MAXINT, 0,
	0, 0, 0, 1
	};
	std::memcpy(xsp, result, sizeof(GzMatrix));
}

void GzRender::CreateXpiMatrix(GzCamera camera, GzMatrix xpi)
{
	float rad = (camera.FOV / 2) * PI / 180;
	GzMatrix result =
	{
	1,0,0,0,
	0,1,0,0,
	0,0,tan(rad),0,
	0,0,tan(rad),1
	};
	std::memcpy(xpi, result, sizeof(GzMatrix));
}

void GzRender::CreateXiwMatrix(GzCamera camera, GzMatrix xiw)
{
	GzCoord cl, x, y, z;
	SubtractVectors(camera.lookat, camera.position, cl);
	MultiplyVector(cl, 1 / GetVectorMagnitude(cl), z);
	float ts = GetVectorMagnitude(z);
	GzCoord zUnit;
	MultiplyVector(z, DoDotProduct(z, camera.worldup), zUnit);
	SubtractVectors(camera.worldup, zUnit, y);
	GetUnitVector(y, y);
	DoCrossProduct(y, z, x);
	GetUnitVector(x, x);
	GzMatrix result =
	{
	x[0],x[1],x[2],-1 * DoDotProduct(x,camera.position),
	y[0],y[1],y[2],-1 * DoDotProduct(y,camera.position),
	z[0],z[1],z[2],-1 * DoDotProduct(z,camera.position),
	0,0,0,1
	};
	std::memcpy(xiw, result, sizeof(GzMatrix));
}



void GzRender::MultiplyGzMatrix(GzMatrix m1, GzMatrix m2, GzMatrix rm) {
	GzMatrix result;
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			result[i][j] = 0;
			for (int k = 0; k < 4; ++k) {
				result[i][j] += m1[i][k] * m2[k][j];
			}
		}
	}
	std::memcpy(rm, result, sizeof(GzMatrix));
}


int GzRender::GzPushMatrix(GzMatrix matrix)
{
	GzMatrix orgMatrix, newMatrix, newMatrixNorm;
	std::memcpy(orgMatrix, matrix, sizeof(GzMatrix));

	if (matlevel > -1) {
		MultiplyGzMatrix(Ximage[matlevel], orgMatrix, newMatrix);
	}
	else
	{
		std::memcpy(newMatrix, orgMatrix, sizeof(GzMatrix));
	}

	ConvertMatrixToUnitaryRotation(orgMatrix);

	if (matlevel > 1) {
		MultiplyGzMatrix(Xnorm[matlevel], orgMatrix, newMatrixNorm);
	}
	else
	{
		std::memcpy(newMatrixNorm, orgMatrix, sizeof(GzMatrix));
	}

	if (matlevel < MATLEVELS) {
		matlevel++;
		std::memcpy(Ximage[matlevel], newMatrix, sizeof(GzMatrix));
		if (matlevel > 1) {
			std::memcpy(Xnorm[matlevel], newMatrixNorm, sizeof(GzMatrix));
		}
		return GZ_SUCCESS;
	}
	return GZ_FAILURE;
}



int GzRender::GzPopMatrix()
{
	if (matlevel > 0) {
		matlevel--;
		return GZ_SUCCESS;
	}
	return GZ_FAILURE;
}

GzIntensity GzRender::ClampGzIntensity(GzIntensity x) {
	if (x < 0) {
		return 0;
	}
	else if (x > 4095) {
		return 4095;
	}
	return x;
}

char GzRender::DownsizeGzIntesity(GzIntensity x) {
	x = x >> 4;
	return x;
}

int GzRender::GzPut(int i, int j, GzIntensity r, GzIntensity g, GzIntensity b, GzIntensity a, GzDepth z)
{
	if (i > 0 && i < xres && j >0 && j < yres) {
		int index = ARRAY(i, j);
		GzPixel* curpixel = &pixelbuffer[index];
		curpixel->red = ClampGzIntensity(r);
		curpixel->green = ClampGzIntensity(g);
		curpixel->blue = ClampGzIntensity(b);
		curpixel->alpha = ClampGzIntensity(a);
		curpixel->z = z;
		return GZ_SUCCESS;
	}

	return GZ_FAILURE;
}


int GzRender::GzGet(int i, int j, GzIntensity* r, GzIntensity* g, GzIntensity* b, GzIntensity* a, GzDepth* z)
{
	if (i > 0 && i < xres && j >0 && j < yres) {
		int index = ARRAY(i, j);
		GzPixel curpixel = pixelbuffer[index];
		*r = curpixel.red;
		*g = curpixel.green;
		*b = curpixel.blue;
		*a = curpixel.alpha;
		*z = curpixel.z;
		return GZ_SUCCESS;
	}
	return GZ_FAILURE;
}


int GzRender::GzFlushDisplay2File(FILE* outfile)
{
	fprintf(outfile, "P6 %d %d 255\r", xres, yres); // PPM header

	for (int j = 0; j < yres; j++)
	{
		for (int i = 0; i < xres; i++)
		{
			int index = ARRAY(i, j);
			GzPixel curpixel = pixelbuffer[index];
			framebuffer[3 * index + RED] = DownsizeGzIntesity(curpixel.red);
			framebuffer[3 * index + GREEN] = DownsizeGzIntesity(curpixel.green);
			framebuffer[3 * index + BLUE] = DownsizeGzIntesity(curpixel.blue);
		}
	}

	// Write the image data
	fwrite(framebuffer, sizeof(char), xres * yres * 3, outfile);
	return GZ_SUCCESS;
}

int GzRender::GzFlushDisplay2FrameBuffer()
{
	for (int j = 0; j < yres; j++)
	{
		for (int i = 0; i < xres; i++)
		{
			int index = ARRAY(i, j);
			GzPixel curpixel = pixelbuffer[index];
			framebuffer[3 * index + 2 - RED] = DownsizeGzIntesity(curpixel.red);
			framebuffer[3 * index + 2 - GREEN] = DownsizeGzIntesity(curpixel.green);
			framebuffer[3 * index + 2 - BLUE] = DownsizeGzIntesity(curpixel.blue);
			char a = framebuffer[3 * index + 2 - RED];
			char b = framebuffer[3 * index + 2 - RED];
			char c = framebuffer[3 * index + 2 - RED];
		}
	}

	return GZ_SUCCESS;
}

int GzRender::GzPutAttribute(int numAttributes, GzToken* nameList, GzPointer* valueList)
{
	for (int i = 0; i < numAttributes; i++) {
		if (nameList[i] == GZ_RGB_COLOR) {
			GzColor* colorPtr = (GzColor*)valueList[i];
			flatcolor[0] = ctoi((*colorPtr)[0]);
			flatcolor[1] = ctoi((*colorPtr)[1]);
			flatcolor[2] = ctoi((*colorPtr)[2]);
		}
		else if (nameList[i] == GZ_DIRECTIONAL_LIGHT) {
			GzLight* lightPtr = (GzLight*)valueList[i];
			GzLight light = *lightPtr;
			lights[numlights] = light;
			numlights++;
		}
		else if (nameList[i] == GZ_AMBIENT_LIGHT) {
			GzLight* lightPtr = (GzLight*)valueList[i];
			ambientlight = *lightPtr;
		}
		else if (nameList[i] == GZ_DIFFUSE_COEFFICIENT) {
			GzColor* colorPtr = (GzColor*)valueList[i];
			Kd[0] = (*colorPtr)[0];
			Kd[1] = (*colorPtr)[1];
			Kd[2] = (*colorPtr)[2];
		}
		else if (nameList[i] == GZ_AMBIENT_COEFFICIENT) {
			GzColor* colorPtr = (GzColor*)valueList[i];
			Ka[0] = (*colorPtr)[0];
			Ka[1] = (*colorPtr)[1];
			Ka[2] = (*colorPtr)[2];
		}
		else if (nameList[i] == GZ_SPECULAR_COEFFICIENT) {
			GzColor* colorPtr = (GzColor*)valueList[i];
			Ks[0] = (*colorPtr)[0];
			Ks[1] = (*colorPtr)[1];
			Ks[2] = (*colorPtr)[2];
		}
		else if (nameList[i] == GZ_DISTRIBUTION_COEFFICIENT) {
			float* specPtr = (float*)valueList[i];
			spec = *specPtr;
		}
		else if (nameList[i] == GZ_INTERPOLATE) {
			int* imPtr = (int*)valueList[i];
			interp_mode = *imPtr;
		}
		else if(nameList[i] == GZ_TEXTURE_MAP)
		{
			tex_fun = (GzTexture)valueList[i];
		}
	}
	return GZ_SUCCESS;
}

int GzRender::GzPutTriangle(int numParts, GzToken* nameList, GzPointer* valueList)
{
	GzCoord* vertsPtr = (GzCoord*)valueList[0];
	GzCoord* normsPtr = (GzCoord*)valueList[1];
	GzTextureIndex* uvPtr = (GzTextureIndex*)valueList[2];
	GzMatrix verts =
	{
	vertsPtr[0][0], vertsPtr[0][1], vertsPtr[0][2], 1,
	vertsPtr[1][0], vertsPtr[1][1], vertsPtr[1][2], 1,
	vertsPtr[2][0], vertsPtr[2][1], vertsPtr[2][2], 1,
	0,0,0,0
	};
	GzMatrix norms =
	{
	normsPtr[0][0], normsPtr[0][1], normsPtr[0][2], 1,
	normsPtr[1][0], normsPtr[1][1], normsPtr[1][2], 1,
	normsPtr[2][0], normsPtr[2][1], normsPtr[2][2], 1,
	0,0,0,0
	};

	GzMatrix transformedVerts, transformedNorms;
	GetTransformedVertices(verts, transformedVerts);
	GetTransformedNormals(norms, transformedNorms);

	GzCoord p1, p2, p3;
	Convert4dCoordTo3dCoord(transformedVerts, 0, p1);
	Convert4dCoordTo3dCoord(transformedVerts, 1, p2);
	Convert4dCoordTo3dCoord(transformedVerts, 2, p3);

	GzCoord n1, n2, n3;
	Convert4dCoordTo3dCoord(transformedNorms, 0, n1);
	Convert4dCoordTo3dCoord(transformedNorms, 1, n2);
	Convert4dCoordTo3dCoord(transformedNorms, 2, n3);

	GzCoord offset = { xOffset, yOffset, 0 };
	SubtractVectors(p1, offset, p1);
	SubtractVectors(p2, offset, p2);
	SubtractVectors(p3, offset, p3);

	Vertex* v1, * v2, * v3;
	v1 = new Vertex(p1, n1, uvPtr[0]);
	v2 = new Vertex(p2, n2, uvPtr[1]);
	v3 = new Vertex(p3, n3, uvPtr[2]);

	Tris tris(v1, v2, v3);

	if (interp_mode == GZ_COLOR) { //Gourard shading
		DoScanLineRaterizationGourardShading(tris);
	} else if (interp_mode == GZ_NORMALS) { //Phong shading
		DoScanLineRaterizationPhongShading(tris);
	} else { //Flat shading
		DoScanLineRaterization(tris);
	}

	return GZ_SUCCESS;
}

void GzRender::GetTransformedVertices(GzMatrix verts, GzMatrix transformedVerts) {
	GzMatrix transformMatrix;
	std::memcpy(transformMatrix, Ximage[matlevel], sizeof(GzMatrix));
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			transformedVerts[i][j] = 0;
			for (int k = 0; k < 4; k++)
			{
				transformedVerts[i][j] += transformMatrix[j][k] * verts[i][k];
			}
		}
	}
}


void GzRender::GetTransformedNormals(GzMatrix verts, GzMatrix transformedNorms) {
	GzMatrix transformMatrix;
	std::memcpy(transformMatrix, Xnorm[matlevel], sizeof(GzMatrix));
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			transformedNorms[i][j] = 0;
			for (int k = 0; k < 4; k++)
			{
				transformedNorms[i][j] += transformMatrix[j][k] * verts[i][k];
			}
		}
	}
}

void GzRender::SetPixelColor(GzCoord pixelPos, GzColor pixelColor) {
	int index = ARRAY(pixelPos[0], pixelPos[1]);
	GzPixel* curpixel = &pixelbuffer[index];
	if (curpixel->z > pixelPos[2] && pixelPos[2]<MAXINT) {
		curpixel->red = ClampGzIntensity(pixelColor[0]);
		curpixel->green = ClampGzIntensity(pixelColor[1]);
		curpixel->blue = ClampGzIntensity(pixelColor[2]);
		curpixel->z = pixelPos[2];
	}
}

#pragma region Flat Shading


void GzRender::DoScanLineRaterization(Tris& tris) {
	ShadePixel(tris.fn->norm, flatcolor);
	FormLinesAndDoRaterization(tris.y1, tris.y21, tris.v1->pos, tris.slope12X, tris.slope12Z, tris.slope13X, tris.slope13Z);
	FormLinesAndDoRaterization(tris.y23, tris.y3, tris.v3->pos, tris.slope23X, tris.slope23Z, tris.slope13X, tris.slope13Z);
}

void GzRender::FormLinesAndDoRaterization(float y1, float y2, GzCoord v, float slopex1, float slopez1, float slopex2, float slopez2) {
	for (int i = y1; i <= y2; i++)
	{
		if (i >= 0 && i < yres) {
			GzCoord v1, v2;
			float dy = i - v[1];
			v1[0] = v[0] + slopex1 * dy;
			v1[1] = i;
			v1[2] = v[2] + slopez1 * dy;

			v2[0] = v[0] + slopex2 * dy;
			v2[1] = i;
			v2[2] = v[2] + slopez2 * dy;

			DoLineRasterization(v1, v2);
		}
	}
}

void GzRender::DoLineRasterization(GzCoord v1, GzCoord v2) {
	float startX = ceil(v1[0]);
	float endX = ceil(v2[0]);
	if (startX > endX) {
		std::swap(startX, endX);
	}
	float y = v1[1];
	float slope = (v2[2] - v1[2]) / (v2[0] - v1[0]);
	for (int x = startX; x < endX; x++)
	{
		if (x >= 0 && x < xres) {
			float z = slope * (x - v1[0]) + v1[2];
			GzCoord point = { x, y, z };
			SetPixelColor(point, flatcolor);
		}
	}
}

#pragma endregion

#pragma region Gourard Shading


void GzRender::DoScanLineRaterizationGourardShading(Tris& tris) {
	GzColor c1, c2, c3;
	bool useTexture = tex_fun != nullptr && IsUVinRange(tris.v1->uv) && IsUVinRange(tris.v2->uv) && IsUVinRange(tris.v3->uv);
	if (useTexture) {
		GzCoord unitVector = { 1,1,1 };
		ShadePixel(tris.v1->norm, unitVector, unitVector, unitVector, tris.v1->shadedColor);
		ShadePixel(tris.v2->norm, unitVector, unitVector, unitVector, tris.v2->shadedColor);
		ShadePixel(tris.v3->norm, unitVector, unitVector, unitVector, tris.v3->shadedColor);
	} else	{
		ShadePixel(tris.v1->norm, tris.v1->shadedColor);
		ShadePixel(tris.v2->norm, tris.v2->shadedColor);
		ShadePixel(tris.v3->norm, tris.v3->shadedColor);
	}
	FormLinesAndDoRaterizationGourardShading(tris.y1, tris.y21, *tris.v1, *tris.v3, *tris.v2, tris.slope12X, tris.slope12Z, tris.slope13X, tris.slope13Z, useTexture);
	FormLinesAndDoRaterizationGourardShading(tris.y23, tris.y3, *tris.v3, *tris.v1, *tris.v2, tris.slope23X, tris.slope23Z, tris.slope13X, tris.slope13Z, useTexture);
}

void GzRender::FormLinesAndDoRaterizationGourardShading(float y1, float y2, Vertex vc, Vertex voc, Vertex vm, float slopex1, float slopez1, float slopex2, float slopez2, bool useTexture) {
	for (int i = y1; i <= y2; i++)
	{
		if (i >= 0 && i < yres) {
			GzCoord p1, p2;
			GzColor c1, c2;
			GzTextureIndex uv1, uv2;
			float dy = i - vc.pos[1];

			//Shorter Side
			p1[0] = vc.pos[0] + slopex1 * dy;
			p1[1] = i;
			p1[2] = vc.pos[2] + slopez1 * dy;
			InterpolateColor(vc.pos, vm.pos, vc.shadedColor, vm.shadedColor, p1, c1);
			InterpolateUV(vc.pos, vm.pos, vc.uv, vm.uv, p1, uv1);

			//Longer Side
			p2[0] = vc.pos[0] + slopex2 * dy;
			p2[1] = i;
			p2[2] = vc.pos[2] + slopez2 * dy;
			InterpolateColor(vc.pos, voc.pos, vc.shadedColor, voc.shadedColor, p2, c2);
			InterpolateUV(vc.pos, voc.pos, vc.uv, voc.uv, p2, uv2);

			Vertex v1(p1, uv1);
			Vertex v2(p2, uv2);
			v1.SetShadedColor(c1);
			v2.SetShadedColor(c2);

			DoLineRasterizationGourardShading(v1, v2, useTexture);
		}
	}
}

void GzRender::DoLineRasterizationGourardShading(Vertex v1, Vertex v2, bool useTexture) {
	GzColor pixelColor;
	float startX = ceil(v1.pos[0]);
	float endX = ceil(v2.pos[0]);
	if (startX > endX) {
		std::swap(startX, endX);
	}
	float y = v1.pos[1];
	float slope = (v2.pos[2] - v1.pos[2]) / (v2.pos[0] - v1.pos[0]);
	for (int x = startX; x < endX; x++)
	{
		if (x >= 0 && x < xres) {
			float z = slope * (x - v1.pos[0]) + v1.pos[2];
			GzCoord point = { x, y, z };
			InterpolateColor(v1.pos, v2.pos, v1.shadedColor, v2.shadedColor, point, pixelColor);
			GzTextureIndex interpolatedUV;
			InterpolateUV(v1.pos, v2.pos, v1.uv, v2.uv, point, interpolatedUV);
			GzColor textureColor;
			if (useTexture && tex_fun(interpolatedUV[0], interpolatedUV[1], textureColor) == GZ_SUCCESS) {
				MultiplyVectors(pixelColor, textureColor, pixelColor);
			}
			SetPixelColor(point, pixelColor);
		}
	}
}

void GzRender::InterpolateColor(GzCoord start, GzCoord end, GzColor startColor, GzColor endColor, GzCoord interpolatePoint, GzColor interpolateColor) {
	float totalDistance = GetVectorsDistance(start, end);
	float startToInterpolateDistance = GetVectorsDistance(start, interpolatePoint);
	for (int i = 0; i < 3; i++)
	{
		interpolateColor[i] = startColor[i] + ((endColor[i] - startColor[i]) * startToInterpolateDistance) / totalDistance;
	}
}

#pragma endregion

#pragma region Phong Shading

void GzRender::DoScanLineRaterizationPhongShading(Tris& tris) {
	FormLinesAndDoRaterizationPhongShading(tris.y1, tris.y21, *tris.v1, *tris.v3, *tris.v2, tris.slope12X, tris.slope12Z, tris.slope13X, tris.slope13Z);
	FormLinesAndDoRaterizationPhongShading(tris.y23, tris.y3, *tris.v3, *tris.v1, *tris.v2, tris.slope23X, tris.slope23Z, tris.slope13X, tris.slope13Z);
}

void GzRender::FormLinesAndDoRaterizationPhongShading(float y1, float y2, Vertex vc,  Vertex voc, Vertex vm,  float slopex1, float slopez1, float slopex2, float slopez2) {
	for (int i = y1; i <= y2; i++)
	{
		if (i >= 0 && i < yres) {
			GzCoord p1, p2;
			GzCoord n1, n2;
			GzTextureIndex uv1, uv2;
			float dy = i - vc.pos[1];

			//Shorter Side
			p1[0] = vc.pos[0] + slopex1 * dy;
			p1[1] = i;
			p1[2] = vc.pos[2] + slopez1 * dy;
			InterpolateNormal(vc.pos, vm.pos, vc.norm, vm.norm, p1, n1);
			InterpolateUV(vc.pos, vm.pos, vc.uv, vm.uv, p1, uv1);

			//Longer Side
			p2[0] = vc.pos[0] + slopex2 * dy;
			p2[1] = i;
			p2[2] = vc.pos[2] + slopez2 * dy;
			InterpolateNormal(vc.pos, voc.pos, vc.norm, voc.norm, p2, n2);
			InterpolateUV(vc.pos, voc.pos, vc.uv, voc.uv, p2, uv2);

			Vertex v1(p1, n1, uv1);
			Vertex v2(p2, n2, uv2);

			DoLineRasterizationPhongShading(v1, v2);
		}
	}
}

void GzRender::DoLineRasterizationPhongShading(Vertex v1, Vertex v2) {
	GzColor pixelColor;
	GzCoord interpolatedNormal;
	float startX = ceil(v1.pos[0]);
	float endX = ceil(v2.pos[0]);
	if (startX > endX) {
		std::swap(startX, endX);
	}
	float y = v1.pos[1];
	float slope = (v2.pos[2] - v1.pos[2]) / (v2.pos[0] - v1.pos[0]);
	for (int x = startX; x < endX; x++)
	{
		if (x >= 0 && x < xres) {
			float z = slope * (x - v1.pos[0]) + v1.pos[2];
			GzCoord point = { x, y, z };
			GzTextureIndex interpolatedUV;
			InterpolateNormal(v1.pos, v2.pos, v1.norm, v2.norm, point, interpolatedNormal);
			InterpolateUV(v1.pos, v2.pos, v1.uv, v2.uv, point, interpolatedUV);
			GzColor textureColor;
			if (interpolatedUV[0] > 0.9 && interpolatedUV[1] > 0.9) {
				float x = interpolatedUV[1];
			}
			if (tex_fun!=nullptr && tex_fun(interpolatedUV[0], interpolatedUV[1], textureColor) == GZ_SUCCESS) {
				ShadePixel(interpolatedNormal, Ks, textureColor, textureColor, pixelColor);
			} else {
				ShadePixel(interpolatedNormal, pixelColor);
			}
			SetPixelColor(point, pixelColor);
		}
	}
}

void GzRender::InterpolateNormal(GzCoord start, GzCoord end, GzCoord startNormal, GzCoord endNormal, GzCoord interpolatePoint, GzCoord interpolateNormal) {
	float totalDistance = GetVectorsDistance(start, end);
	float startToInterpolateDistance = GetVectorsDistance(start, interpolatePoint);
	for (int i = 0; i < 3; i++)
	{
		interpolateNormal[i] = startNormal[i] + ((endNormal[i] - startNormal[i]) * startToInterpolateDistance) / totalDistance;
	}
}

#pragma endregion 

void GzRender::InterpolateUV(GzCoord start, GzCoord end, GzTextureIndex startUV, GzTextureIndex endUV, GzCoord interpolatePoint, GzTextureIndex interpolateUV) {
	GzTextureIndex tempStartUV, tempEndUV;
	float modifier = (interpolatePoint[2] / ((float)MAXINT32 - interpolatePoint[2])) + 1;
	float totalDistance = GetVectorsDistance(start, end);
	float startToInterpolateDistance = GetVectorsDistance(start, interpolatePoint);
	tempStartUV[0] = startUV[0] / ((start[2] / ((float)MAXINT - start[2])) + 1);
	tempStartUV[1] = startUV[1] / ((start[2] / ((float)MAXINT - start[2])) + 1);
	tempEndUV[0] = endUV[0] / ((end[2] / ((float)MAXINT - end[2])) + 1);
	tempEndUV[1] = endUV[1] / ((end[2] / ((float)MAXINT - end[2])) + 1);
	for (int i = 0; i < 2; i++)
	{
		interpolateUV[i] = (tempStartUV[i] + (((tempEndUV[i] - tempStartUV[i])) * startToInterpolateDistance) / totalDistance)*((interpolatePoint[2] / ((float)MAXINT - interpolatePoint[2])) + 1);
	}	 
}

void GzRender::ShadePixel(GzCoord norm, GzColor color) {
	ShadePixel(norm, Ks, Kd, Ka, color);
}

void GzRender::ShadePixel(GzCoord norm, GzColor Ks, GzColor Kd, GzColor Ka, GzColor color)
{
	GzCoord reflectedLight;
	GzCoord light;
	GzCoord colorVector = { 0, 0, 0 };
	GzCoord eye = { 0,0,-1 };
	GetUnitVector(norm, norm);
	for (int i = 0; i < numlights; i++)
	{
		GzCoord specReflectionColor, diffuseColor;

		GzCoord reflectedLight, alteredNorm, temp;
		MultiplyVector(norm, 2 * DoDotProduct(norm, lights[i].direction), alteredNorm);
		SubtractVectors(alteredNorm, lights[i].direction, reflectedLight);

		float x1, x2;
		x1 = DoDotProduct(norm, lights[i].direction);
		x2 = DoDotProduct(norm, eye);

		if ((x1 > 0 && x2 > 0) || (x1 < 0 && x2 < 0)) {
			GzCoord curNormal;
			std::memcpy(curNormal, norm, sizeof(GzCoord));
			if (x1 < 0) {
				MultiplyVector(curNormal, -1, curNormal);
			}
			MultiplyVector(lights[i].color, pow(DoDotProduct(reflectedLight, eye), spec), temp);
			MultiplyVectors(temp, Ks, specReflectionColor);
			AddVectors(specReflectionColor, colorVector, colorVector);

			MultiplyVector(lights[i].color, Clamp(DoDotProduct(curNormal, lights[i].direction), 0, 1), temp);
			MultiplyVectors(temp, Kd, diffuseColor);
			AddVectors(diffuseColor, colorVector, colorVector);
		}
	}

	GzCoord ambientColor;
	MultiplyVectors(ambientlight.color, Ka, ambientColor);

	AddVectors(ambientColor, colorVector, colorVector);

	for (int i = 0; i < 3; i++)
	{
		colorVector[i] = Clamp(colorVector[i], 0, 1);
	}

	color[0] = ctoi(colorVector[0]);
	color[1] = ctoi(colorVector[1]);
	color[2] = ctoi(colorVector[2]);
}

Vertex::Vertex(GzCoord pos, GzCoord norm, GzTextureIndex uv) {
	std::memcpy(this->pos, pos, sizeof(GzCoord));
	std::memcpy(this->norm, norm, sizeof(GzCoord));
	std::memcpy(this->uv, uv, sizeof(GzTextureIndex));
}

Vertex::Vertex(GzCoord pos, GzTextureIndex uv) {
	std::memcpy(this->pos, pos, sizeof(GzCoord));
	std::memcpy(this->uv, uv, sizeof(GzTextureIndex));
}

void Vertex::SetShadedColor(GzColor shadedColor) {
	std::memcpy(this->shadedColor, shadedColor, sizeof(GzColor));
}

Vertex::~Vertex() {
	
}

Tris::Tris(Vertex* v1, Vertex* v2, Vertex* v3)
{
	fn = v1;
	if (v1->pos[1] > v2->pos[1]) {
		std::swap(v1, v2);
	}
	if (v1->pos[1] > v3->pos[1]) {
		std::swap(v1, v3);
	}
	if (v2->pos[1] > v3->pos[1]) {
		std::swap(v2, v3);
	}

	this->v1 = v1;
	this->v2 = v2;
	this->v3 = v3;


	this->slope12X = CalculateSlopeX(this->v1->pos, this->v2->pos);
	this->slope12Z = CalculateSlopeZ(this->v1->pos, this->v2->pos);

	this->slope13X = CalculateSlopeX(this->v1->pos, this->v3->pos);
	this->slope13Z = CalculateSlopeZ(this->v1->pos, this->v3->pos);

	this->slope23X = CalculateSlopeX(this->v2->pos, this->v3->pos);
	this->slope23Z = CalculateSlopeZ(this->v2->pos, this->v3->pos);

	this->y1 = ceil(this->v1->pos[1]);
	this->y21 = floor(this->v2->pos[1]);
	this->y23 = ceil(this->v2->pos[1]);
	this->y3 = floor(this->v3->pos[1]);
}

Tris::~Tris() {
	delete v1, v2, v3;
}