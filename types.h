/*
 * types.h
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
#ifndef __TYPES_H__
#define __TYPES_H__

#ifndef THREAD_NUMBER
	#define THREAD_NUMBER 4		// default - 4 threads
#endif

#ifndef DBL_EPSILON
#define DBL_EPSILON        2.2204460492503131e-16
#endif

typedef double Item;
// ITM_EPSILON is for data comparing, set it to zero for integer types
#define ITM_EPSILON  DBL_EPSILON

#define OMP_NUM_THREADS THREAD_NUMBER
#define Stringify(x) #x
#define OMP_FOR(x) _Pragma(Stringify(omp parallel for x))
#ifndef MAX
#define MAX(x,y) ((x) > (y) ? (x) : (y))
#endif
#ifndef MIN
#define MIN(x,y) ((x) < (y) ? (x) : (y))
#endif

#endif // __TYPES_H__

