











#ifndef TEXTURE_H
#define TEXTURE_H

#include "includes.h"

struct texdata
{
	int channels;			// The channels in the image (3 = RGB : 4 = RGBA)
	int sizex;				// The width of the image in pixels
	int sizey;				// The height of the image in pixels
	unsigned char *data;	// The image pixel data
};

typedef struct texdata texdata;

void ltexinit(texdata *lt);
void ltexfree(texdata *lt);
int savepng(const char* full, unsigned char* data, int sizex, int sizey, int channels);
texdata *loadpng(const char *file);

#endif
