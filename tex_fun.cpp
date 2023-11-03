/* Texture functions for cs580 GzLib  */
#include    "stdafx.h" 
#include  "stdio.h"
#include  "Gz.h"
#include <algorithm>
#include <cmath>

GzColor* image = NULL;
int xs, ys;
int reset = 1;

int GetIndex(int u, int v) {
    return v * (xs) + u;
}

void MultiplyColorVector(GzColor c1, float scale, GzColor rc) {
    rc[0] = c1[0] * scale;
    rc[1] = c1[1] * scale;
    rc[2] = c1[2] * scale;
}

void DoColorBilinearInterpolation(float u, float v, GzColor finalColor) {
    float scaledU = u * (xs-1);
    float scaledV = v * (ys - 1);
    int minU = floor(scaledU);
    int maxU = ceil(scaledU);
    int minV = floor(scaledV);
    int maxV = ceil(scaledV);

    float s = scaledU - minU;
    float t = scaledV - minV;

    GzColor A, B, C, D;

    std::memcpy(A, image[GetIndex(minU, minV)], sizeof(GzColor));
    std::memcpy(B, image[GetIndex(maxU, minV)], sizeof(GzColor));
    std::memcpy(C, image[GetIndex(maxU, maxV)], sizeof(GzColor));
    std::memcpy(D, image[GetIndex(minU, maxV)], sizeof(GzColor));

    MultiplyColorVector(A, (1 - s) * (1 - t), A);
    MultiplyColorVector(B, s * (1 - t), B);
    MultiplyColorVector(C, s * t, C);
    MultiplyColorVector(D, (1 - s) * t, D);

    finalColor[RED] = A[RED] + B[RED] + C[RED] + D[RED];
    finalColor[GREEN] = A[GREEN] + B[GREEN] + C[GREEN] + D[GREEN];
    finalColor[BLUE] = A[BLUE] + B[BLUE] + C[BLUE] + D[BLUE];

}

/* Image texture function */
int tex_fun(float u, float v, GzColor color)
{
    unsigned char   pixel[3];
    unsigned char     dummy;
    char      foo[8];
    int       i, j;
    FILE* fd;

    if (reset) {          /* open and load texture file */
        fd = fopen("texture", "rb");
        if (fd == NULL) {
            fprintf(stderr, "texture file not found\n");
            exit(-1);
        }
        fscanf(fd, "%s %d %d %c", foo, &xs, &ys, &dummy);
        image = (GzColor*)malloc(sizeof(GzColor) * (xs + 1) * (ys + 1));
        if (image == NULL) {
            fprintf(stderr, "malloc for texture image failed\n");
            exit(-1);
        }

        for (i = 0; i < xs * ys; i++) { /* create array of GzColor values */
            fread(pixel, sizeof(pixel), 1, fd);
            image[i][RED] = (float)((int)pixel[RED]) * (1.0 / 255.0);
            image[i][GREEN] = (float)((int)pixel[GREEN]) * (1.0 / 255.0);
            image[i][BLUE] = (float)((int)pixel[BLUE]) * (1.0 / 255.0);
        }

        reset = 0;          /* init is done */
        fclose(fd);
    }

    if (u >= 0 && u <= 1 && v >= 0 && v <= 1) {
        DoColorBilinearInterpolation(u, v, color);
        return (color[0]>=0 && color[0]<=1 && color[1]>=0 && color[1]<=1 && color[2]>=0 && color[2]<=1) ? GZ_SUCCESS : GZ_FAILURE;
    }
    return GZ_FAILURE;
}

/* Procedural texture function */
int ptex_fun(float u, float v, GzColor color)
{
    unsigned char   pixel[3];
    unsigned char     dummy;
    char      foo[8];
    int       i, j;

    if (reset) {          /* open and load texture file */
        xs = 100;
        ys = 100;
        image = (GzColor*)malloc(sizeof(GzColor) * (xs + 1) * (ys + 1));
        if (image == NULL) {
            fprintf(stderr, "malloc for texture image failed\n");
            exit(-1);
        }
        /*int max = xs>ys ? pow(xs, 2) : pow(ys, 2);
        int maxThird = max / 3;*/
        
        float maxNum = xs > ys ? xs : ys;
        float max = pow(maxNum/2, 2);
        float transformValue = maxNum / 2;
        for (i = 0; i < xs; i++) { /* create array of GzColor values */
            for (j = 0; j < ys; j++) {
                int index = GetIndex(i, j);

                float value = pow(i - transformValue, 2) + pow(j - transformValue, 2);
                if (value < max/9) {
                    image[index][RED] = 1;
                    image[index][GREEN] = 1;
                    image[index][BLUE] = (((sin(value) + 1.0) / 2.0) + (value/(max/9)))/2;
                }
                else if (value < 4 * max / 9) {
                    image[index][RED] = 1;
                    image[index][GREEN] = (((sin(value) + 1.0) / 2.0) + (value / (4 * max / 9))) / 2;
                    image[index][BLUE] = 1;
                }
                else {
                    image[index][RED] = (((sin(value) + 1.0) / 2.0) + (value / (max))) / 2;
                    image[index][GREEN] = 1;
                    image[index][BLUE] = 1;
                }

            }
        }

        reset = 0;
    }

    if (u >= 0 && u <= 1 && v >= 0 && v <= 1) {
        DoColorBilinearInterpolation(u, v, color);
        return (color[0] >= 0 && color[0] <= 1 && color[1] >= 0 && color[1] <= 1 && color[2] >= 0 && color[2] <= 1) ? GZ_SUCCESS : GZ_FAILURE;
    }
    return GZ_FAILURE;
}

/* Free texture memory */
int GzFreeTexture()
{
    if (image != NULL)
        free(image);
    return GZ_SUCCESS;
}