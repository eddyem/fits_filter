/*
 * convfilter.h
 *
 * Copyright 2015 Edward V. Emelianov <eddy@sao.ru, edward.emelianoff@gmail.com>
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
#ifndef __GRADIENT_H__
#define __GRADIENT_H__

#include "types.h"
#include "fits.h"

// FilterType
typedef enum{
	 LAPGAUSS			// laplasian of gaussian
	,GAUSS				// gaussian
	,SOBELH				// Sobel horizontal
	,SOBELV				// -//- vertical
	,SIMPLEGRAD			// simple gradient (by Sobel)
	,PREWITTH			// Prewitt (horizontal) - simple derivative
	,PREWITTV			// -//- (vertical)
	,SCHARRH			// Scharr (modified Sobel)
	,SCHARRV
	,STEP				// "posterisation"
} FType;

typedef struct{
	FType FilterType;	// filter type
	int w;				// filter width
	int h;				// height
	double sx;			// x half-width
	double sy;			// y half-width (sx, sy - for Gaussian-type filters)
} Filter;

IMAGE *DiffFilter(IMAGE *img, Filter *f);
IMAGE *GradFilterSimple(IMAGE *img);

#endif // __GRADIENT_H__
