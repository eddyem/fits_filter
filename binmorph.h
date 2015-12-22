/*
 * binmorph.h
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

#pragma once
#ifndef __EROSION_DILATION_H__
#define __EROSION_DILATION_H__

#include <stdbool.h>
#include <stdlib.h>
#include <stdint.h>
#include "fits.h"

// convert image types
uint8_t *u16tochar(uint16_t *im, int W, int H, int *stride);
//bool *char2bool(uint8_t *image, int W, int H, int W_0);
uint16_t *chartou16(uint8_t *image, int W, int H, int W_0);

// morphological operations
uint8_t *dilation(uint8_t *image, int W, int H);
uint8_t *erosion(uint8_t *image, int W, int H);
uint8_t *FC_filter(uint8_t *image, int W, int H);

// logical operations
uint8_t *imand(uint8_t *im1, uint8_t *im2, int W, int H);
uint8_t *substim(uint8_t *im1, uint8_t *im2, int W, int H);
/*
// conncomp
// this is a box structure containing one object; data is aligned by original image bytes!
typedef struct {
	uint8_t *data; // pattern data in "packed" format
	int x,   // x coordinate of LU-pixel of box in "unpacked" image (by pixels)
	y,       // y -//-
	x_0;     // x coordinate in "packed" image (morph operations should work with it)
	size_t N;// number of component, starting from 1
} CCbox;
*/

IMAGE *cclabel4(IMAGE *I, double threshold, size_t *Nobj);
IMAGE *cclabel8(IMAGE *I, double threshold, size_t *Nobj);
#endif // __EROSION_DILATION_H__


