/*
 * convfilter.c - convolution-based filters
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
#include <stdio.h>
#include <fftw3.h>
#include <math.h>

#include "usefull_macros.h"
#include "convfilter.h"

// from http://graphics.stanford.edu/%7Eseander/bithacks.html#RoundUpPowerOf2
int nextpow2(int i, int j){
	inline int p2oi(int i){
		unsigned int v = (unsigned int)i - 1;
		v |= v >> 1;
		v |= v >> 2;
		v |= v >> 4;
		v |= v >> 8;
		v |= v >> 16;
		v++;
		return (int) v;
	}
	int p1 = p2oi(i), p2 = p2oi(j);
	return MAX(p1, p2);
}

void fftshift(int size, Item *m){
	int h = size/2, ss, ss1, i, j, k, l, p;
	ss = (2*h+1)*h;
	ss1 = (2*h-1)*h;
	for(j = 0; j < h; j++){
		k = j * size;
		l = k + h;
		for(i = 0; i < h; i++, k++, l++){
			register Item tmp;
			p = k + ss;
			tmp = m[k]; m[k] = m[p]; m[p] = tmp;
			p = l + ss1;
			tmp = m[l]; m[l] = m[p]; m[p] = tmp;
		}
	}
}

// Gaussian mask building/
Item *build_G_filter(int size, Filter *f){
	int y0 = 0, y1 = size, x0 = 0, x1 = size;
	if(f->sx < 1.){
		WARNX(_("sigma_x is too low, set to 1."));
		f->sx = 1.;
	}
	if(f->sy < 1.){
		WARNX(_("sigma_y is too low, set to 1."));
		f->sy = 1.;
	}
	double sx2 = f->sx * f->sx, sy2 = f->sy * f->sy;
	double hw, hh;
	Item *mask = MALLOC(Item, size * size);
	#ifdef EBUG
	double t0=dtime();
	#endif
	if(f->w < size && f->w > 0){
		x0 = (size - f->w + 1) / 2;
		x1 = x0 + f->w;
	}
	if(f->h < size && f->h > 0){
		y0 = (size - f->h + 1) / 2;
		y1 = y0 + f->h;
	}
	hh = -(double)size / 2.;
	hw = -(double)size / 2.+(double)x0;
	DBG("y0=%d, y1=%d, hw=%g, hh=%g",y0,y1,hw,hh);
	const double ss = 1./(2*M_PI*f->sx*f->sy) / size / size;
	//const double ss = 1./sqrt(2*M_PI*f->sx*f->sy) / size / size / 4.;
	OMP_FOR(shared(mask))
	for(int y = y0; y < y1; y++){
		double X, Y, y2, x2, R;
		int str, x;
		X = hw;
		str = y * size;
		Y = ((double)y) + hh;
		y2 = Y*Y/sy2;
		for(x = x0; x < x1; x++, X+=1.){
			x2 = X*X/sx2;
			R = x2 + y2;
			mask[str + x] = ss * exp(-R/2.);
		}
	}
	DBG("time=%f\n", dtime()-t0);
	return mask;
}


// Lapgauss mask building
Item *build_LG_filter(int size, Filter *f){
	int y0 = 0, y1 = size, x0 = 0, x1 = size;
	double sx2 = f->sx * f->sx, sy2 = f->sy * f->sy;
	double hw, hh;
	Item *mask = MALLOC(Item, size * size);
	#ifdef EBUG
	double t0=dtime();
	#endif
	if(f->w < size && f->w > 0){
		x0 = (size - f->w + 1) / 2;
		x1 = x0 + f->w;
	}
	if(f->h < size && f->h > 0){
		y0 = (size - f->h + 1) / 2;
		y1 = y0 + f->h;
	}
	hh = -(double)size / 2.;
	hw = -(double)size / 2.+(double)x0;
	DBG("y0=%d, y1=%d, x0=%d, x1=%d, hw=%g, hh=%g",y0,y1,x0,x1,hw,hh);
//	double ss = 3. / hh / hh / sqrt(-hh);
	const double ss = 1./sqrt(2*M_PI*f->sx*f->sy) / hh / hh / sqrt(-hh);
	OMP_FOR(shared(mask))
	for(int y = y0; y < y1; y++){
		double X, Y, y2, ys2, x2, R;
		int str, x;
		X = hw;
		str = y * size;
		Y = ((double)y) + hh;
		y2 = Y*Y/sy2;
		ys2 = (y2 - 1) / sy2;
		for(x = x0; x < x1; x++, X+=1.){
			x2 = X*X/sx2;
			R = x2 + y2;
			mask[str + x] = ss * ((x2-1.)/sx2 + ys2) * exp(-R/2.);
		}
	}
	DBG("time=%f\n", dtime()-t0);
	return mask;
}


// Elementary filter mask building
Item *build_S_filter(int size, Filter *f){
	int a0, a1;
	a0 = (size - 2) / 2;
	a1 = a0 + 3;
	//double hh = -(double)(size / 2);
	double ss = 1. / size / size / 8.;
	double Y = -1.;
	Item *mask = MALLOC(Item, size * size);
	#ifdef EBUG
	double t0 = dtime();
	#endif
	Item (*filtfun)(Item x, Item y);
	Item imcopy(Item x, Item y){return (x == 0 && y == 0) ? 1 : 0;}
	Item sobelh(Item x, Item y){return -x*(2.-fabs(y));}
	Item sobelv(Item x, Item y){return -y*(2.-fabs(x));}
	Item prewitth(Item x, _U_ Item y){return x;}
	Item prewittv(_U_ Item x, Item y){return y;}
	Item scharrh(Item x, Item y){return -x * (10. - 7.*fabs(y));}
	Item scharrv(Item x, Item y){return -y * (10. - 7.*fabs(x));}
	switch(f->FilterType){
		case SOBELH:
			filtfun = sobelh;
		break;
		case SOBELV:
			filtfun = sobelv;
		break;
		case PREWITTH:
			filtfun = prewitth;
		break;
		case PREWITTV:
			filtfun = prewittv;
		break;
		case SCHARRH:
			filtfun = scharrh;
		break;
		case SCHARRV:
			filtfun = scharrv;
		break;
		default:
			filtfun = imcopy;
	}
	for(int y = a0; y < a1; y++, Y+=1.){
		double X = -1.;
		int str, x;
		str = y * size;
		for(x = a0; x < a1; x++, X+=1.){
			mask[str + x] = ss*filtfun(X, Y);
		}
	}
	DBG("time=%f\n", dtime()-t0);
	return mask;
}

/*
 * Filtering by convolution with a filter
 * Input:
 *		ima - input image
 *		f - filter parameters
 * Output:
 * Returns NULL on error or converted image
 */
IMAGE *DiffFilter(IMAGE *img, Filter *f, _U_ Itmarray *u){
	int ssize;
	static int fftw_ini = 0;
	int sizex = img->width, sizey = img->height;
	size_t blklen = sizex * sizeof(Item);
	int size2 = nextpow2(sizex, sizey);
	#ifdef EBUG
	double t0 = dtime();
	#endif
	// build filter
	Item *mask;
	switch(f->FilterType){
		case LAPGAUSS:
			mask = build_LG_filter(size2, f);
		break;
		case GAUSS:
			mask = build_G_filter(size2, f);
		break;
		default:
			mask = build_S_filter(size2, f);
	}
	IMAGE *out = similarFITS(img, DOUBLE_IMG);
	Item *res = out->data, *inputima = img->data;
	ssize = size2 * size2; // FFT image size
	if(!fftw_ini)
		if(!(fftw_ini = fftw_init_threads())){
			WARN(_("FFTW error"));
			return NULL;
		}
	DBG("img (%d x %d) -> (%d x %d), time=%f\n", sizex,sizey, size2,size2, dtime()-t0);
	// allocate memory for objects
	Item *ima = MALLOC(Item, ssize);
	fftw_complex *Fimg = MALLOC(fftw_complex, ssize);
	// define direct FFT
	fftw_plan_with_nthreads(THREAD_NUMBER);
	fftw_plan fftimg = fftw_plan_dft_r2c_2d(size2, size2, ima, Fimg, FFTW_ESTIMATE);
	// copy ima -> img
	OMP_FOR(shared(ima, inputima))
	for(int j = 0; j < sizey; j++){
		int k = j * size2, l = j * sizex;
		memcpy(&ima[k], &inputima[l], blklen); // main image
	}
	int diff = size2 - sizey;
	if(diff){
		OMP_FOR(shared(ima, inputima))
		for(int j = 0; j < diff; ++j){
			int k = (j + sizey)* size2, l = (sizey - 1 - j) * sizex;
			memcpy(&ima[k], &inputima[l], blklen); // lower left subimage (mirrored vertically)
		}
	}
	// the rest: right subimages (mirrored by x)
	diff = size2 - sizex;
	if(diff){
		OMP_FOR(shared(ima))
		for(int j = 0; j < size2; ++j){
			int l = j * size2 + sizex, k = l - 1;
			for(int i = 0; i < diff; ++i, --k, ++l)
				ima[l] = ima[k];
		}
	}
	fftw_execute(fftimg); // now Fimg is fourier transform of input image
	FREE(ima); // we don't need in anymore
	// define filter FFT
	fftw_plan_with_nthreads(THREAD_NUMBER);
	fftw_complex *Fmask = MALLOC(fftw_complex, ssize);
	fftw_plan fftmask = fftw_plan_dft_r2c_2d(size2, size2, mask, Fmask, FFTW_ESTIMATE);
	fftw_execute(fftmask);
	FREE(mask);
	DBG("here");
	// filtered picture:
	DBG("filter image, time=%f\n", dtime()-t0);
	Item *resm = MALLOC(Item, ssize);
	OMP_FOR(shared(Fmask, Fimg))
	for(int i = 0; i < ssize; i++){ // convolution by multiplication in Fourier space
		Item a, b, c, d;
		a = Fimg[i][0]; c = Fmask[i][0];
		b = Fimg[i][1]; d = Fmask[i][1];
		Fimg[i][0] = a*c - b*d;
		Fimg[i][1] = b*c + a*d;
	}
	FREE(Fmask);
	// define inverse FFT
	fftw_plan_with_nthreads(THREAD_NUMBER);
	fftw_plan ifft = fftw_plan_dft_c2r_2d(size2, size2, Fimg, resm, FFTW_ESTIMATE);
	fftw_execute(ifft);
	FREE(Fimg);
	fftshift(size2, resm);
	OMP_FOR(shared(res, resm))
	for(int j = 0; j < sizey; j++){
		int k = j * size2, l = j * sizex;
		memcpy(&res[l], &resm[k], blklen);
	}
	FREE(resm);
	fftw_destroy_plan(fftimg);
	fftw_destroy_plan(ifft);
	fftw_destroy_plan(fftmask);
	fftw_cleanup_threads();
	DBG("time=%f\n", dtime()-t0);
	return out;
}

/*
 * Simple gradient filter based on two Sobel filters
 * output = sqrt(SobelH(input)^2+SobelV(input)^2)
 */
IMAGE *GradFilterSimple(IMAGE *img, _U_ Filter *fu, _U_ Itmarray *u){
	#ifdef EBUG
	double t0 = dtime();
	#endif
	Filter f;
	f.FilterType = SOBELH;
	IMAGE *horiz = DiffFilter(img, &f, NULL);
	f.FilterType = SOBELV;
	IMAGE *vert = DiffFilter(img, &f, NULL);
	int w = img->width, h = img->height;
	double *dst1 = horiz->data, *dst2 = vert->data;
	OMP_FOR(shared(dst1, dst2))
	for(int y = 0; y < h; y++){
		int x = y*w;
		Item *ptr2 = &dst2[x], *ptr1 = &dst1[x];
		for(x = 0; x < w; ++x, ++ptr1, ++ptr2)
			*ptr1 = sqrt((*ptr2)*(*ptr2) + (*ptr1)*(*ptr1));
	}
	imfree(&vert);
	DBG("time=%f\n", dtime()-t0);
	return horiz;
}
