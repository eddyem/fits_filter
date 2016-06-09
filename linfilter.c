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

stepscalespairs scales[] = {
	{UNIFORM, "uniform"},
	{LOG,     "log"},
	{EXP,     "exp"},
	{SQRT,    "sqrt"},
	{POW,     "pow"},
	{0, NULL}
};


/**
 * calculate simple statistics by image
 * @param img (i) - input image
 * @param min, max, mean, std, med (o) - statistical values
 */
void get_statictics(IMAGE *img, Item *min, Item *max,
					Item *mean, Item *std, Item *med){
	if(!img) return;
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
void fillIsoScale(Filter *f, Itmarray *scale, Item min, Item wd){
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
	scale->data = MALLOC(Item, M);
	scale->size = M;
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
	Item *ar = scale->data;
	for(y = 0; y < M; y++){
		ar[y] = scalefn(y+1);
		DBG("level %d: I=%g", y+1, ar[y]);
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
IMAGE *StepFilter(IMAGE *img, Filter *f, Itmarray *scale){
	if(f->w < 2 || f->w > 255){
		// "Неправильное количество уровней: %d"
		WARNX(_("Wrong levels amount: %d"), f->w);
		return NULL;
	}
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
	if(img->keylist){ // remove BZERO & BSCALE for given image format
		char buf[80];
		list_modify_key(img->keylist, "BZERO", "0");
		list_modify_key(img->keylist, "BSCALE", "1");
		list_modify_key(img->keylist, "DATAMIN", "0");
		snprintf(buf, 21, "%d", f->w);
		list_modify_key(img->keylist, "DATAMAX", buf);
	//	snprintf(buf, 80, "HISTORY step filter with %d levels (%s distribution)",
	//		f->w, scales[f->h].name);
	//	list_add_record(&img->keylist, buf);
	}
	IMAGE *out = similarFITS(img, BYTE_IMG);
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

/**
 * set all values more than 'up' to 'up & less than 'low' to 'low'
 */
void cut_bounds(IMAGE *img, Item low, Item up){
	if(!(low < DBL_MAX - 1. || up < DBL_MAX - 1.)) return;
	Item min, max;
	bool lowct = FALSE, upct = FALSE;
	get_statictics(img, &min, &max, NULL, NULL, NULL);
	if(low < DBL_MAX - 1.)
		lowct = TRUE;
	if(up < DBL_MAX - 1.)
		upct = TRUE;
	int w = img->width, h = img->height, y;
	OMP_FOR(shared(img))
	for(y = 0; y < h; ++y){
		Item *data = &img->data[y * w];
		int x;
		for(x = 0; x < w; ++x, ++data){
			if(lowct && *data < low) *data = low;
			else if(upct && *data > up) *data = up;
		}
	}
	char buf[80];
	if(lowct && !upct)
		snprintf(buf, 80, "COMMENT cut lower bound to value %g", (double)low);
	else if(upct && !lowct)
		snprintf(buf, 80, "COMMENT cut upper bound to value %g", (double)up);
	else
		snprintf(buf, 80, "COMMENT cut lower bound to %g & upper to %g", (double)low, (double)up);
	list_add_record(&img->keylist, buf);
}

/**
 * Convert image to binary
 * image values [min, max] converted to [0, 1]
 * threshold \in (0, 1)
 * (I < threshold) = 0, (I >= threshold) = 1
 *
 * if threshold less than 0 image would be inverted!
 *
 * @param thrvalue (o) - threshold intensity level
 */
uint16_t *binarize(IMAGE *img, double threshold, Item *thrvalue){
	DBG("THRES: %g", threshold);
	if(threshold < -1. + DBL_EPSILON || threshold > 1. - DBL_EPSILON){
		/// Пороговое значение должно лежать в интервале (-1, 1)
		WARNX(_("The threshold value should be in interval (-1, 1)"));
		return NULL;
	}
	bool invert = FALSE;
	if(threshold < 0.){
		threshold = -threshold;
		invert = TRUE;
	}
	Item min, max;
	get_statictics(img, &min, &max, NULL, NULL, NULL);
	Item thrval = min + (max - min) * threshold;
	int w = img->width, h = img->height, y;
	uint16_t *ret = MALLOC(uint16_t, w*h);
	OMP_FOR(shared(img))
	for(y = 0; y < h; ++y){
		Item *idata = &img->data[y * w];
		uint16_t *odata = &ret[y * w];
		int x;
		for(x = 0; x < w; ++x, ++idata, ++odata){
			if(*idata < thrval) *odata = invert;
			else *odata = !invert;
		}
	}
	if(thrvalue) *thrvalue = thrval;
	return ret;
}

IMAGE *get_binary(IMAGE *img, double threshold){
	Item thrval;
	uint16_t *binary = binarize(img, threshold, &thrval);
	if(!binary) return NULL;
	IMAGE *ret = buildFITSfromdat(img->height, img->width, SHORT_IMG, (uint8_t*)binary);
	FREE(binary);
	char buf[80];
	snprintf(buf, 80, "COMMENT binarize image by threshold value %g", (double)thrval);
	list_add_record(&ret->keylist, buf);
	snprintf(buf, 80, "COMMENT    (%g%% fromdata range%s)", fabs(threshold)*100.,
		(threshold < 0.) ? ", inverted" : "");
	list_add_record(&ret->keylist, buf);
	return ret;
}
