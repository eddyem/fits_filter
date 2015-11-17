/*
 * linfilter.c - linear filters
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

#include <math.h>
#include "usefull_macros.h"
#include "linfilter.h"
#include "median.h"

/**
 * calculate simple statistics by image
 * @param img (i) - input image
 * @param min, max, mean, std, med (o) - statistical values
 */
void get_statictics(IMAGE *img, Item *min, Item *max,
					Item *mean, Item *std, Item *med){
	size_t sz = img->width * img->height;
	if(min || max || mean || std){
		Item *idata = img->data;
		Item minval = *idata, maxval = minval, sum = minval, sum2 = minval*minval;
		size_t i;
		++idata;
		for(i = 1; i < sz; ++i, ++idata){
			Item val = *idata;
			if(val > maxval) maxval = val;
			else if(val < minval) minval = val;
			sum += val;
			sum2 += val*val;
		}
		if(min){
			*min = minval;
			DBG("minimum: %g", minval);
		}
		if(max){
			*max = maxval;
			DBG("maximum: %g", maxval);
		}
		if(mean){
			*mean = sum / sz;
			DBG("mean: %g", *mean);
		}
		if(std){
			sum /= sz;
			*std = sqrt(sum2/sz - sum*sum);
			DBG("std: %g", *std);
		}
	}
	if(med){
		*med = quick_select(img->data, sz);
		DBG("median: %g", *med);
	}
}


/*
 * Fill isolines' scale (an array)
 * Input:
 * 		f - filter for given method
 * 		min - minimum value of intensity
 * 		wd - max-min (dinamic range)
 * Output:
 * 		scale - a pointer to array (allocated in this function)
 */
void fillIsoScale(Filter *f, Item **scale, Item min, Item wd){
	int M = f->w, y;
	Item (*scalefn)(Item in);
	Item step, Nsteps = (Item)f->w; // amount of intervals
	Item Suniform(Item in){
		return step*in + min;
	}
	Item Slog(Item in){
		return min + step * log(in + 1.);
	}
	Item Sexp(Item in){
		return min - 1. + exp(step * in);
	}
	Item Ssqrt(Item in){
		return min + step * sqrt(in);
	}
	Item Spow(Item in){
		return min + step*in*in;
	}
	*scale = MALLOC(Item, M);
	switch(f->h){
		case LOG:
			scalefn = Slog; step = wd / log(Nsteps + 1.);
		break;
		case EXP:
			scalefn = Sexp; step = log(wd + 1.) / Nsteps;
		break;
		case SQRT:
			scalefn = Ssqrt; step = wd/sqrt(Nsteps);
		break;
		case POW:
			scalefn = Spow; step = wd/Nsteps/Nsteps;
		break;
		default:
			scalefn = Suniform; step = wd/Nsteps;
	}
	for(y = 0; y < M; y++){
		(*scale)[y] = scalefn(y+1);
		DBG("level %d: I=%g", y+1, (*scale)[y]);
	}
}
/*
 * Threshold filtering ("posterization")
 * Input:
 *		img - input image
 *		f - filter:
 *			f-> w - number of levels of posterization, [2.255]
 *			f-> h - type of posterization (StepType)
 * Output:
 *		result - filtered image, the memory is allocated in this procedure
 *		scale - the scale of intensities, the memory is allocated here (if the scale!=NULL)
 *
 * TODO: save the result in the char, not float; learn display function
 */
IMAGE *StepFilter(IMAGE *img, Filter *f, double **scale){
	if(f->w < 2 || f->w > 255) return FALSE;
	Item Nsteps = (Item)f->w; // amount of intervals
	Item step, max, min;
	get_statictics(img, &min, &max, NULL, NULL, NULL);
	Item wd = max - min;
	if(fabs(wd) < ITM_EPSILON) return FALSE;
	Item (*stepfn)(Item in);
	Item Funiform(Item in){
		return floor((in-min)/step);
	}
	Item Flog(Item in){
		return floor(exp((in-min)/step) - 1.);
	}
	Item Fexp(Item in){
		return floor(log(in - min + 1.) / step);
	}
	Item Fsqrt(Item in){
		return floor((in - min)*(in - min) / step);
	}
	Item Fpow(Item in){
		return floor(sqrt((in - min)/step));
	}
	#ifdef EBUG
	double t0 = dtime();
	#endif
	switch(f->h){
		case LOG: // I = Imin + step*ln(N+1)
			stepfn = Flog;
			step = wd / log(Nsteps + 1.);
		break;
		case EXP: // I = Imin - 1 + exp(step * N)
			stepfn = Fexp;
			step = log(wd + 1.) / Nsteps;
		break;
		case SQRT: // I = Imin + sqrt(step * N)
			stepfn = Fsqrt;
			step = wd*wd/Nsteps;
		break;
		case POW: // I = Imin + step * N^2
			stepfn = Fpow;
			step = wd/Nsteps/Nsteps;
		break;
		default: // I = Imin + step*N
			stepfn = Funiform;
			step = wd/Nsteps;
	}
	IMAGE *out = similarFITS(img, SHORT_IMG);
	Item *res = out->data, *inputima = img->data;
	size_t sizex = img->width, sizey = img->height;
	OMP_FOR(shared(res, inputima))
	for(size_t y = 0; y < sizey; ++y){
		Item *iout = &res[y*sizex], *iin = &inputima[y*sizex];
		for(size_t x = 0; x < sizex; ++x, ++iin, ++iout){
			*iout = stepfn(*iin);
		}
	}
	if(scale)
		fillIsoScale(f, scale, min, wd);
	DBG("SF: %d sublevels, step=%g, time=%f\n", f->w, step, dtime()-t0);
	return out;
}
