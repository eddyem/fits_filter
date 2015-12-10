/*
 * linfilter.h
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
#ifndef __LINFILTER_H__
#define __LINFILTER_H__

#include "types.h"
#include "fits.h"
#include "convfilter.h"

// f->h for STEP filter
typedef enum{
	 UNIFORM
	,LOG
	,EXP
	,SQRT
	,POW
} StepType;

void get_statictics(IMAGE *img, Item *min, Item *max,
					Item *mean, Item *std, Item *med);
IMAGE *StepFilter(IMAGE *img, Filter *f, Itmarray *scale);

#endif // __LINFILTER_H__
