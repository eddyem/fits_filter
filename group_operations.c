/*
 * group_operations.c - operations with a lot of files; if size differs, truncate to smallest!
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

#include "group_operations.h"
#include "fits.h"
#include "usefull_macros.h"
#include "median.h"
#include <omp.h>

// files is NULL-terminated list - array of images
typedef IMAGE * (*mathfuncptr)(IMAGE **files);

/**
 * calculate minimal sizes in list of images
 */
static void get_minsizes(int *minh, int *minw, IMAGE **infiles){
	FNAME();
	int hmin = (*infiles)->height, wmin = (*infiles)->width;
	DBG("initial min size: %dx%d", wmin, hmin);
	while(*(++infiles)){
		int h = (*infiles)->height, w = (*infiles)->width;
		DBG("process image with h=%d, w=%d", h, w);
		if(hmin > h) hmin = h;
		if(wmin > w) wmin = w;
	}
	if(minh) *minh = hmin;
	if(minw) *minw = wmin;
	DBG("minimal sizes: %dx%d", wmin, hmin);
}

/**
 * Calculate sum of images in list
 */
static IMAGE* math_sum(IMAGE **infiles){
	FNAME();
	if(!infiles || !*infiles) return NULL;
	int h, w;
	get_minsizes(&h, &w, infiles);
	IMAGE *outp = newFITS(h, w, DOUBLE_IMG);
	double *odata = outp->data;
	while(*infiles){
		IMAGE *in = *infiles;
		int y, oriW = in->width;
		DBG("process file with W=%d", oriW);
		OMP_FOR(shared(in, odata))
		for(y = 0; y < h; ++y){
			int x;
			double *iptr = &in->data[oriW*y], *optr = &odata[w*y];
			for(x = 0; x < w; ++x)
				*(optr++) += *(iptr++);
		}
		++infiles;
	}
	DBG("OK");
	return outp;
}

static IMAGE* math_mean(IMAGE **infiles){
	FNAME();
	IMAGE *out = math_sum(infiles);
	if(!out) return NULL;
	double images_amount = 1.;
	while(*(++infiles)) ++images_amount;
	double *odata = out->data;
	int y, h = out->height, w = out->width;
	OMP_FOR(shared(odata))
	for(y = 0; y < h; ++y){
		int x;
		double *optr = &odata[w*y];
		for(x = 0; x < w; ++x)
			*(optr++) /= images_amount;
	}
	return out;
}

static IMAGE* math_median(IMAGE **infiles){
	FNAME();
	if(!infiles || !*infiles) return NULL;
	int images_amount = 1;
	IMAGE **f = infiles;
	while(*(++f)) ++images_amount;
	if(images_amount < 2){
		/// Не могу вычислить медиану меньше чем для двух изображений
		ERRX(_("Can't calculate median for less than two images"));
	}
	int h, w, y;
	get_minsizes(&h, &w, infiles);
	IMAGE *out = newFITS(h, w, DOUBLE_IMG);
	double *odata = out->data;
	double *idata = MALLOC(double, OMP_NUM_THREADS * images_amount);
	OMP_FOR(shared(odata, idata))
	for(y = 0; y < h; ++y){
		int x, N;
		double *optr = &odata[w*y];
		double *inp = &idata[images_amount * omp_get_thread_num()];
		for(x = 0; x < w; ++x){
			for(N = 0; N < images_amount; ++N){
				inp[N] = infiles[N]->data[infiles[N]->width * y + x];
			}
			*(optr++) = calc_median(inp, images_amount);
		}
	}
	FREE(idata);
	return out;
}

IMAGE *make_group_operation(int names_amount, char **names, MathOper oper){
	FNAME();
	char buf[80];
	mathfuncptr anoper = NULL;
	KeyList *list = NULL;
	switch(oper){
		case MATH_SUM:
			anoper = math_sum;
			snprintf(buf, 80, "COMMENT math sum for next files:");
			list_add_record(&list, buf);
		break;
		case MATH_MEAN:
			anoper = math_mean;
			snprintf(buf, 80, "COMMENT math mean for next files:");
			list_add_record(&list, buf);
		break;
		case MATH_MEDIAN:
			anoper = math_median;
			snprintf(buf, 80, "COMMENT image-by-image median for next files:");
			list_add_record(&list, buf);
		break;
		case MATH_NONE:
		default:
			/// Неизвестная групповая операция
			ERRX(_("Unknown group operation"));
		break;
	}
	// now fill list of files
	IMAGE **filelist = MALLOC(IMAGE*, names_amount + 1); // +1 for terminated NULL
	int i, ctr = 0;
	for(i = 0; i < names_amount; ++i){
		if(!readFITS(names[i], &filelist[i])){
			/// Не могу прочесть файл %s из последовательности
			WARNX(_("Can't read file %s from sequence"), names[i]);
		}else{
			++ctr;
			snprintf(buf, 80, "COMMENT %d: %s", ctr, names[i]);
			list_add_record(&list, buf);
		}
	}
	if(ctr < 2){
		/// Количество доступных файлов меньше двух
		ERRX(_("The amount of available files less than two"));
	}
	IMAGE *ret = anoper(filelist);
	if(ret) ret->keylist = list;
	// free filelist
	for(i = 0; i < ctr; ++i) imfree(&filelist[i]);
	FREE(filelist);
	return ret;
}
