/*
 * binmorph.c - functions for morphological operations on binary image
 *
 * Copyright 2015 Edward V. Emelianoff <eddy@sao.ru>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <err.h>
#include <sys/time.h>
#include "binmorph.h"
#include "types.h"
#include "usefull_macros.h"
#include "linfilter.h"

// global arrays for erosion/dilation masks
static uint8_t *ER = NULL, *DIL = NULL;
static bool __Init_done = false; // == true if arrays are inited

/*
 * =================== AUXILIARY FUNCTIONS ===================>
 */

/**
 * This function inits masks arrays for erosion and dilation
 * You may call it yourself or it will be called when one of
 * `erosion` or `dilation` functions will be ran first time
 */
static void morph_init(){
	if(__Init_done) return;
	int i;
	ER = MALLOC(uint8_t, 256);
	DIL = MALLOC(uint8_t, 256);
	for(i = 0; i < 256; i++){
		ER[i]  = i & ((i << 1) | 1) & ((i >> 1) | (0x80)); // don't forget that << and >> set borders to zero
		DIL[i] = i | (i << 1) | (i >> 1);
	}
	__Init_done = true;
}
/*
 * <=================== AUXILIARY FUNCTIONS ===================
 */

/*
 * =================== CONVERT IMAGE TYPES ===================>
 */

/**
 * Convert boolean image into pseudo-packed (1 char == 8 pixels)
 * @param im (i)  - image to convert
 * @param W, H    - size of image im (must be larger than 1)
 * @param W_0 (o) - (stride) new width of image
 * @return allocated memory area with "packed" image
 */
uint8_t *u16tochar(uint16_t *im, int W, int H, int *W_0){
	if(W < 2 || H < 2) ERRX("image size too small");
	int y, W0 = (W + 7) / 8;
	uint8_t *ret = MALLOC(uint8_t, W0 * H);
	OMP_FOR()
	for(y = 0; y < H; y++){
		int x, i, X;
		uint16_t *ptr = &im[y*W];
		uint8_t *rptr = &ret[y*W0];
		for(x = 0, X = 0; x < W0; x++, rptr++){
			for(i = 7; i > -1 && X < W; i--, X++, ptr++){
				*rptr |= *ptr << i;
			}
		}
	}
	if(W_0) *W_0 = W0;
	return ret;
}

/**
 * Convert "packed" image into boolean
 * @param image (i) - input image
 * @param W, H, W_0 - size of image and width of "packed" image
 * @return allocated memory area with "unpacked" image
 *
bool *char2bool(uint8_t *image, int W, int H, int W_0){
	int y;
	bool *ret = MALLOC(bool, W * H);
	OMP_FOR()
	for(y = 0; y < H; y++){
		int x, X, i;
		bool *optr = &ret[y*W];
		uint8_t *iptr = &image[y*W_0];
		for(x = 0, X = 0; x < W_0; x++, iptr++)
			for(i = 7; i > -1 && X < W; i--, X++, optr++){
				*optr = (*iptr >> i) & 1;
			}
	}
	return ret;
}*/

/**
 * Convert "packed" image into size_t array for conncomp procedure
 * @param image (i) - input image
 * @param W, H, W_0 - size of image and width of "packed" image
 * @return allocated memory area with copy of an image
 */
uint16_t *chartou16(uint8_t *image, int W, int H, int W_0){
	int y;
	uint16_t *ret = MALLOC(uint16_t, W * H);
	OMP_FOR()
	for(y = 0; y < H; y++){
		int x, X, i;
		uint16_t *optr = &ret[y*W];
		uint8_t *iptr = &image[y*W_0];
		for(x = 0, X = 0; x < W_0; x++, iptr++)
			for(i = 7; i > -1 && X < W; i--, X++, optr++){
				*optr = (*iptr >> i) & 1;
			}
	}
	return ret;
}

/*
 * <=================== CONVERT IMAGE TYPES ===================
 */

/*
 * =================== MORPHOLOGICAL OPERATIONS ===================>
 */

/**
 * Remove all non-4-connected pixels
 * @param image (i) - input image
 * @param W, H      - size of image
 * @return allocated memory area with converted input image
 */
uint8_t *FC_filter(uint8_t *image, int W, int H){
	if(W < 1 || H < 2) errx(1, "4-connect: image size too small");
	uint8_t *ret = MALLOC(uint8_t, W*H);
	int y = 0, w = W-1, h = H-1;
	// top of image, y = 0
	#define IM_UP
	#include "fc_filter.h"
	#undef IM_UP
	// mid of image, y = 1..h-1
	#include "fc_filter.h"
	// image bottom, y = h
	y = h;
	#define IM_DOWN
	#include "fc_filter.h"
	#undef IM_DOWN
	return ret;
}

/**
 * Make morphological operation of dilation
 * @param image (i) - input image
 * @param W, H      - size of image
 * @return allocated memory area with dilation of input image
 */
uint8_t *dilation(uint8_t *image, int W, int H){
	if(W < 2 || H < 2) errx(1, "Dilation: image size too small");
	if(!__Init_done) morph_init();
	uint8_t *ret = MALLOC(uint8_t, W*H);
	int y = 0, w = W-1, h = H-1;
	// top of image, y = 0
	#define IM_UP
	#include "dilation.h"
	#undef IM_UP
	// mid of image, y = 1..h-1
	#include "dilation.h"
	// image bottom, y = h
	y = h;
	#define IM_DOWN
	#include "dilation.h"
	#undef IM_DOWN
	return ret;
}

/**
 * Make morphological operation of erosion
 * @param image (i) - input image
 * @param W, H      - size of image
 * @return allocated memory area with erosion of input image
 */
uint8_t *erosion(uint8_t *image, int W, int H){
	if(W < 2 || H < 2) errx(1, "Erosion: image size too small");
	if(!__Init_done) morph_init();
	uint8_t *ret = MALLOC(uint8_t, W*H);
	int y, w = W-1, h = H-1;
	OMP_FOR()
	for(y = 1; y < h; y++){ // reset first & last rows of image
		uint8_t *iptr = &image[W*y];
		uint8_t *optr = &ret[W*y];
		uint8_t p = ER[*iptr] & 0x7f & iptr[-W] & iptr[W];
		int x;
		if(!(iptr[1] & 0x80)) p &= 0xfe;
		*optr++ = p;
		iptr++;
		for(x = 1; x < w; x++, iptr++, optr++){
			p = ER[*iptr] & iptr[-W] & iptr[W];
			if(!(iptr[-1] & 1)) p &= 0x7f;
			if(!(iptr[1] & 0x80)) p &= 0xfe;
			*optr = p;
		}
		p = ER[*iptr] & 0xfe & iptr[-W] & iptr[W];
		if(!(iptr[-1] & 1)) p &= 0x7f;
		*optr++ = p;
		iptr++;
	}
	return ret;
}

/*
 * <=================== MORPHOLOGICAL OPERATIONS ===================
 */

/*
 * =================== LOGICAL OPERATIONS ===================>
 */

/**
 * Logical AND of two images
 * @param im1, im2 (i) - two images
 * @param W, H         - their size (of course, equal for both images)
 * @return allocated memory area with   image = (im1 AND im2)
 */
uint8_t *imand(uint8_t *im1, uint8_t *im2, int W, int H){
	uint8_t *ret = MALLOC(uint8_t, W*H);
	int y;
	OMP_FOR()
	for(y = 0; y < H; y++){
		int x, S = y*W;
		uint8_t *rptr = &ret[S], *p1 = &im1[S], *p2 = &im2[S];
		for(x = 0; x < W; x++)
			*rptr++ = *p1++ & *p2++;
	}
	return ret;
}

/**
 * Substitute image 2 from image 1: reset to zero all bits of image 1 which set to 1 on image 2
 * @param im1, im2 (i) - two images
 * @param W, H         - their size (of course, equal for both images)
 * @return allocated memory area with    image = (im1 AND (!im2))
 */
uint8_t *substim(uint8_t *im1, uint8_t *im2, int W, int H){
	uint8_t *ret = MALLOC(uint8_t, W*H);
	int y;
	OMP_FOR()
	for(y = 0; y < H; y++){
		int x, S = y*W;
		uint8_t *rptr = &ret[S], *p1 = &im1[S], *p2 = &im2[S];
		for(x = 0; x < W; x++)
			*rptr++ = *p1++ & (~*p2++);
	}
	return ret;
}

/*
 * <=================== LOGICAL OPERATIONS ===================
 */

/*
 * =================== CONNECTED COMPONENTS LABELING ===================>
 */



/**
 * label 4-connected components on image
 * (slow algorythm, but easy to parallel)
 *
 * @param I (i)    - image ("packed")
 * @param W,H,W_0  - size of the image (W - width in pixels, W_0 - width in octets)
 * @param Nobj (o) - number of objects found
 * @return an array of labeled components
 */
uint16_t *_cclabel4(uint8_t *Img, int W, int H, int W_0, size_t *Nobj){
	uint8_t *I = FC_filter(Img, W_0, H);
	uint16_t *labels = chartou16(I, W, H, W_0);
	FREE(I);
	#include "cclabling.h"
	return labels;
}

// label 8-connected components, look cclabel4
uint16_t *_cclabel8(uint16_t *labels, int W, int H, size_t *Nobj){
	#define LABEL_8
	#include "cclabling.h"
	#undef LABEL_8
	return labels;
}

/**
 * Make connection-component labeling
 * output image vould have uint16_t data
 * @param img (i)       - input image
 * @param threshold (i) - threshold level in value of dynamic range (0,1)
 * @param Nobj (o)      - amount of object found (or NULL if not needed)
 */
IMAGE *cclabel4(IMAGE *img, double threshold, size_t *Nobj){
	/*if(N != 4 || N != 8){
		/// Могу работать лишь для четырех- и восьмисвязных областей
		ERRX(_("Can work only for 4- or 8-connected components"));
	}*/
	double thrval;
	uint16_t *binary = binarize(img, threshold, &thrval);
	if(!binary) return NULL;
	int W_0;
	uint8_t *Ima = u16tochar(binary, img->width, img->height, &W_0);
	FREE(binary);
	size_t N;
	uint16_t *dat = _cclabel4(Ima, img->width, img->height, W_0, &N);
	if(Nobj) *Nobj = N;
	FREE(Ima);
	IMAGE *ret = buildFITSfromdat(img->height, img->width, SHORT_IMG, (uint8_t*)dat);
	FREE(dat);
	char buf[80];
	snprintf(buf, 80, "COMMENT found %zd 4-connected components, threshold value %g",
		N, (double)thrval);
	list_add_record(&ret->keylist, buf);
	snprintf(buf, 80, "COMMENT    (%g%% fromdata range%s)", fabs(threshold)*100.,
		(threshold < 0.) ? ", inverted" : "");
	list_add_record(&ret->keylist, buf);
	return ret;
}

IMAGE *cclabel8(IMAGE *img, double threshold, size_t *Nobj){
	double thrval;
	uint16_t *binary = binarize(img, threshold, &thrval);
	if(!binary) return NULL;
	size_t N;
	_cclabel8(binary, img->width, img->height, &N);
	if(Nobj) *Nobj = N;
	IMAGE *ret = buildFITSfromdat(img->height, img->width, SHORT_IMG, (uint8_t*)binary);
	FREE(binary);
	char buf[80];
	snprintf(buf, 80, "COMMENT found %zd 4-connected components, threshold value %g",
		N, (double)thrval);
	list_add_record(&ret->keylist, buf);
	snprintf(buf, 80, "COMMENT    (%g%% fromdata range%s)", fabs(threshold)*100.,
		(threshold < 0.) ? ", inverted" : "");
	list_add_record(&ret->keylist, buf);
	return ret;
}

/*
 * <=================== CONNECTED COMPONENTS LABELING ===================
 */


/*
 * <=================== template ===================>
 */
