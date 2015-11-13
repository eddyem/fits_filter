// FOR MEDIATOR:
// Copyright (c) 2011 ashelly.myopenid.com under <http://www.opensource.org/licenses/mit-license>
// FOR opt_medXX:
// Copyright (c) 1998 Nicolas Devillard. Public domain.
// FOR qickselect:
// "Numerical recipes in C", Second Edition,
//  Cambridge University Press, 1992, Section 8.5, ISBN 0-521-43108-5
//  Code by Nicolas Devillard - 1998. Public domain.

/*
 * median.c
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

// TODO: resolve problem with borders

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "median.h"
#include "usefull_macros.h"

// largest radius for adaptive median filter - 7x7 pix
#define LARGEST_ADPMED_RADIUS  (3)

#ifndef DBL_EPSILON
#define DBL_EPSILON        2.2204460492503131e-16
#endif

#define OMP_NUM_THREADS 4
#define Stringify(x) #x
#define OMP_FOR(x) _Pragma(Stringify(omp parallel for x))

#define ELEM_SWAP(a, b) {register Item t = a; a = b; b = t;}
#define PIX_SORT(a, b)  {if (p[a] > p[b]) ELEM_SWAP(p[a], p[b]);}
Item opt_med3(Item *p){
	PIX_SORT(0, 1); PIX_SORT(1, 2); PIX_SORT(0, 1);
	return(p[1]) ;
}
Item opt_med5(Item *p){
	PIX_SORT(0, 1); PIX_SORT(3, 4); PIX_SORT(0, 3);
	PIX_SORT(1, 4); PIX_SORT(1, 2); PIX_SORT(2, 3) ;
	PIX_SORT(1, 2);
	return(p[2]) ;
}
Item opt_med7(Item *p){
	PIX_SORT(0, 5); PIX_SORT(0, 3); PIX_SORT(1, 6);
	PIX_SORT(2, 4); PIX_SORT(0, 1); PIX_SORT(3, 5);
	PIX_SORT(2, 6); PIX_SORT(2, 3); PIX_SORT(3, 6);
	PIX_SORT(4, 5); PIX_SORT(1, 4); PIX_SORT(1, 3);
	PIX_SORT(3, 4); return (p[3]) ;
}
Item opt_med9(Item *p){
	PIX_SORT(1, 2); PIX_SORT(4, 5); PIX_SORT(7, 8);
	PIX_SORT(0, 1); PIX_SORT(3, 4); PIX_SORT(6, 7);
	PIX_SORT(1, 2); PIX_SORT(4, 5); PIX_SORT(7, 8);
	PIX_SORT(0, 3); PIX_SORT(5, 8); PIX_SORT(4, 7);
	PIX_SORT(3, 6); PIX_SORT(1, 4); PIX_SORT(2, 5);
	PIX_SORT(4, 7); PIX_SORT(4, 2); PIX_SORT(6, 4);
	PIX_SORT(4, 2); return(p[4]);
}
Item opt_med25(Item *p){
	PIX_SORT(0, 1)  ; PIX_SORT(3, 4)  ; PIX_SORT(2, 4) ;
	PIX_SORT(2, 3)  ; PIX_SORT(6, 7)  ; PIX_SORT(5, 7) ;
	PIX_SORT(5, 6)  ; PIX_SORT(9, 10) ; PIX_SORT(8, 10) ;
	PIX_SORT(8, 9)  ; PIX_SORT(12, 13); PIX_SORT(11, 13) ;
	PIX_SORT(11, 12); PIX_SORT(15, 16); PIX_SORT(14, 16) ;
	PIX_SORT(14, 15); PIX_SORT(18, 19); PIX_SORT(17, 19) ;
	PIX_SORT(17, 18); PIX_SORT(21, 22); PIX_SORT(20, 22) ;
	PIX_SORT(20, 21); PIX_SORT(23, 24); PIX_SORT(2, 5) ;
	PIX_SORT(3, 6)  ; PIX_SORT(0, 6)  ; PIX_SORT(0, 3) ;
	PIX_SORT(4, 7)  ; PIX_SORT(1, 7)  ; PIX_SORT(1, 4) ;
	PIX_SORT(11, 14); PIX_SORT(8, 14) ; PIX_SORT(8, 11) ;
	PIX_SORT(12, 15); PIX_SORT(9, 15) ; PIX_SORT(9, 12) ;
	PIX_SORT(13, 16); PIX_SORT(10, 16); PIX_SORT(10, 13) ;
	PIX_SORT(20, 23); PIX_SORT(17, 23); PIX_SORT(17, 20) ;
	PIX_SORT(21, 24); PIX_SORT(18, 24); PIX_SORT(18, 21) ;
	PIX_SORT(19, 22); PIX_SORT(8, 17) ; PIX_SORT(9, 18) ;
	PIX_SORT(0, 18) ; PIX_SORT(0, 9)  ; PIX_SORT(10, 19) ;
	PIX_SORT(1, 19) ; PIX_SORT(1, 10) ; PIX_SORT(11, 20) ;
	PIX_SORT(2, 20) ; PIX_SORT(2, 11) ; PIX_SORT(12, 21) ;
	PIX_SORT(3, 21) ; PIX_SORT(3, 12) ; PIX_SORT(13, 22) ;
	PIX_SORT(4, 22) ; PIX_SORT(4, 13) ; PIX_SORT(14, 23) ;
	PIX_SORT(5, 23) ; PIX_SORT(5, 14) ; PIX_SORT(15, 24) ;
	PIX_SORT(6, 24) ; PIX_SORT(6, 15) ; PIX_SORT(7, 16) ;
	PIX_SORT(7, 19) ; PIX_SORT(13, 21); PIX_SORT(15, 23) ;
	PIX_SORT(7, 13) ; PIX_SORT(7, 15) ; PIX_SORT(1, 9) ;
	PIX_SORT(3, 11) ; PIX_SORT(5, 17) ; PIX_SORT(11, 17) ;
	PIX_SORT(9, 17) ; PIX_SORT(4, 10) ; PIX_SORT(6, 12) ;
	PIX_SORT(7, 14) ; PIX_SORT(4, 6)  ; PIX_SORT(4, 7) ;
	PIX_SORT(12, 14); PIX_SORT(10, 14); PIX_SORT(6, 7) ;
	PIX_SORT(10, 12); PIX_SORT(6, 10) ; PIX_SORT(6, 17) ;
	PIX_SORT(12, 17); PIX_SORT(7, 17) ; PIX_SORT(7, 10) ;
	PIX_SORT(12, 18); PIX_SORT(7, 12) ; PIX_SORT(10, 18) ;
	PIX_SORT(12, 20); PIX_SORT(10, 20); PIX_SORT(10, 12) ;
	return (p[12]);
}
/*
Item quick_select(Item **arr, int n){
	int low, high;
	int median;
	int middle, ll, hh;
	float ret;
	low = 0 ; high = n-1 ; median = (low + high) / 2;
	for(;;){
		if(high <= low) // One element only
			break;
		if(high == low + 1){ // Two elements only
			PIX_SORT(arr[low], arr[high]) ;
			break;
		}
		// Find median of low, middle and high items; swap into position low
		middle = (low + high) / 2;
		PIX_SORT(arr[middle], arr[high]) ;
		PIX_SORT(arr[low], arr[high]) ;
		PIX_SORT(arr[middle], arr[low]) ;
		// Swap low item (now in position middle) into position (low+1)
		ELEM_SWAP(arr[middle], arr[low+1]) ;
		// Nibble from each end towards middle, swapping items when stuck
		ll = low + 1;
		hh = high;
		for(;;){
			do ll++; while (*arr[low] > *arr[ll]);
			do hh--; while (*arr[hh] > *arr[low]);
			if(hh < ll) break;
			ELEM_SWAP(arr[ll], arr[hh]) ;
		}
		// Swap middle item (in position low) back into correct position
		ELEM_SWAP(arr[low], arr[hh]) ;
		// Re-set active partition
		if (hh <= median) low = ll;
		if (hh >= median) high = hh - 1;
	}
	ret = *arr[median];
	return ret;
}
#undef PIX_SORT
#undef ELEM_SWAP
*/





#define ItemLess(a,b) ((a)<(b))
#define ItemMean(a,b) (((a)+(b))/2)

typedef struct Mediator_t{
	Item* data; // circular queue of values
	int* pos;   // index into `heap` for each value
	int* heap;  // max/median/min heap holding indexes into `data`.
	int N;      // allocated size.
	int idx;    // position in circular queue
	int ct;     // count of items in queue
} Mediator;

/*--- Helper Functions ---*/

#define minCt(m) (((m)->ct-1)/2) //count of items in minheap
#define maxCt(m) (((m)->ct)/2) //count of items in maxheap

//returns 1 if heap[i] < heap[j]
inline int mmless(Mediator* m, int i, int j){
	return ItemLess(m->data[m->heap[i]],m->data[m->heap[j]]);
}

//swaps items i&j in heap, maintains indexes
inline int mmexchange(Mediator* m, int i, int j){
	int t = m->heap[i];
	m->heap[i] = m->heap[j];
	m->heap[j] = t;
	m->pos[m->heap[i]] = i;
	m->pos[m->heap[j]] = j;
	return 1;
}

//swaps items i&j if i<j; returns true if swapped
inline int mmCmpExch(Mediator* m, int i, int j){
	return (mmless(m,i,j) && mmexchange(m,i,j));
}

//maintains minheap property for all items below i/2.
void minSortDown(Mediator* m, int i){
	for(; i <= minCt(m); i*=2){
		if(i>1 && i < minCt(m) && mmless(m, i+1, i)) ++i;
		if(!mmCmpExch(m,i,i/2)) break;
	}
}

//maintains maxheap property for all items below i/2. (negative indexes)
void maxSortDown(Mediator* m, int i){
	for(; i >= -maxCt(m); i*=2){
		if(i<-1 && i > -maxCt(m) && mmless(m, i, i-1)) --i;
	if(!mmCmpExch(m,i/2,i)) break;
	}
}

//maintains minheap property for all items above i, including median
//returns true if median changed
int minSortUp(Mediator* m, int i){
	while (i > 0 && mmCmpExch(m, i, i/2)) i /= 2;
	return (i == 0);
}

//maintains maxheap property for all items above i, including median
//returns true if median changed
int maxSortUp(Mediator* m, int i){
	while (i < 0 && mmCmpExch(m, i/2, i)) i /= 2;
	return (i == 0);
}

/*--- Public Interface ---*/


//creates new Mediator: to calculate `nItems` running median.
//mallocs single block of memory, caller must free.
Mediator* MediatorNew(int nItems){
	int size = sizeof(Mediator) + nItems*(sizeof(Item)+sizeof(int)*2);
	Mediator* m = malloc(size);
	m->data = (Item*)(m + 1);
	m->pos = (int*) (m->data + nItems);
	m->heap = m->pos + nItems + (nItems / 2); //points to middle of storage.
	m->N = nItems;
	m->ct = m->idx = 0;
	while (nItems--){ //set up initial heap fill pattern: median,max,min,max,...
		m->pos[nItems] = ((nItems+1)/2) * ((nItems&1)? -1 : 1);
		m->heap[m->pos[nItems]] = nItems;
	}
	return m;
}


//Inserts item, maintains median in O(lg nItems)
void MediatorInsert(Mediator* m, Item v){
	int isNew=(m->ct<m->N);
	int p = m->pos[m->idx];
	Item old = m->data[m->idx];
	m->data[m->idx]=v;
	m->idx = (m->idx+1) % m->N;
	m->ct+=isNew;
	if(p>0){ //new item is in minHeap
		if (!isNew && ItemLess(old,v)) minSortDown(m,p*2);
		else if (minSortUp(m,p)) maxSortDown(m,-1);
	}else if (p<0){ //new item is in maxheap
		if (!isNew && ItemLess(v,old)) maxSortDown(m,p*2);
		else if (maxSortUp(m,p)) minSortDown(m, 1);
	}else{ //new item is at median
		if (maxCt(m)) maxSortDown(m,-1);
		if (minCt(m)) minSortDown(m, 1);
	}
}

//returns median item (or average of 2 when item count is even)
Item MediatorMedian(Mediator* m){
	Item v = m->data[m->heap[0]];
	if ((m->ct&1) == 0) v = ItemMean(v, m->data[m->heap[-1]]);
	return v;
}

// median + min/max
Item MediatorStat(Mediator* m, Item *minval, Item *maxval){
	Item v= m->data[m->heap[0]];
	if ((m->ct&1) == 0) v = ItemMean(v,m->data[m->heap[-1]]);
	Item min = v, max = v;
	int i;
	for(i = -maxCt(m); i < 0; ++i){
		int v = m->data[m->heap[i]];
		if(v < min) min = v;
	}
	*minval = min;
	for(i = 1; i <= minCt(m); ++i){
		int v = m->data[m->heap[i]];
		if(v > max) max = v;
	}
	*maxval = max;
	return v;
}

/**
 * median by cross 3x3 pixels (5 pixels total)
 */
static void get_median_cross(IMAGE *img, IMAGE *out){
	size_t w = img->width, h = img->height;
	Item *med = out->data, *inputima = img->data;
#ifdef EBUG
	double t0 = dtime();
#endif
	OMP_FOR(shared(inputima, med))
	for(size_t x = 1; x < w - 1; ++x){
		size_t xx, xm = w + x + 2, y, ymax, xmin = xm - 3;
		Mediator* m = MediatorNew(5);
		// initial fill
		MediatorInsert(m, inputima[x]); // (x,0)
		for(xx = xmin; xx < xm; ++xx) // (line with y = 1)
			MediatorInsert(m, inputima[xx]);
		ymax = h - 1;
		xx = x + 2 * w;
		size_t medidx = x + w;
		for(y = 1; y < ymax; ++y, xx += w, medidx += w){
			MediatorInsert(m, inputima[xx]);
			med[medidx] = MediatorMedian(m);
		}
		free(m);
	}
	DBG("time for median filtering by cross 3x3 of image %zdx%zd: %gs", w, h,
		dtime() - t0);
}
/**
 * filter image by median (seed*2 + 1) x (seed*2 + 1)
 */
IMAGE *get_median(IMAGE *img, int seed){
	size_t w = img->width, h = img->height, siz = w*h, bufsiz = siz*sizeof(int);
	IMAGE *out = MALLOC(IMAGE, 1);
	Item *med = MALLOC(Item, bufsiz), *inputima = img->data;
	out->data = med;
	out->width = w;
	out->height = h;
	out->dtype = img->dtype;

	memcpy(med, inputima, bufsiz);
	if(seed == 0){
		get_median_cross(img, out);
		return out;
	}

	size_t blksz = seed * 2 + 1, fullsz = blksz * blksz;
#ifdef EBUG
	double t0 = dtime();
#endif
	OMP_FOR(shared(inputima, med))
	for(size_t x = seed; x < w - seed; ++x){
		size_t xx, yy, xm = x + seed + 1, y, ymax = blksz - 1, xmin = x - seed;
		Mediator* m = MediatorNew(fullsz);
		// initial fill
		for(yy = 0; yy < ymax; ++yy)
			for(xx = xmin; xx < xm; ++xx)
				MediatorInsert(m, inputima[xx + yy*w]);
		ymax = 2*seed*w;
		xmin += ymax;
		xm += ymax;
		ymax = h - seed;
		size_t medidx = x + seed * w;
		for(y = seed; y < ymax; ++y, xmin += w, xm += w, medidx += w){
			for(xx = xmin; xx < xm; ++xx)
				MediatorInsert(m, inputima[xx]);
			med[medidx] = MediatorMedian(m);
		}
		free(m);
	}
	DBG("time for median filtering %zdx%zd of image %zdx%zd: %gs", blksz, blksz, w, h,
		dtime() - t0);
	return out;
}

/**
 * procedure for finding median value in window 5x5
 * PROBLEM: bounds
 */
static Item adp_med_5by5(IMAGE *img, size_t x, size_t y){
	size_t w = img->width;
	Item arr[25], *arrptr = arr, *dataptr = &img->data[x + (y - 2) * w];
	for(int yy = 0; yy < 5; ++yy, dataptr += w, arrptr += 5)
		memcpy(arrptr, dataptr, 5*sizeof(Item));
	return opt_med25(arr);
}

static void get_adp_median_cross(IMAGE *img, IMAGE *out){
	size_t w = img->width, h = img->height;
	Item *med = out->data, *inputima = img->data;
#ifdef EBUG
	double t0 = dtime();
#endif
	OMP_FOR(shared(inputima, med))
	for(size_t x = 1; x < w - 1; ++x){
		size_t xx, xm = w + x + 2, y, ymax, xmin = xm - 3;
		Mediator* m = MediatorNew(5);
		// initial fill
		MediatorInsert(m, inputima[x]); // (x,0)
		for(xx = xmin; xx < xm; ++xx) // (line with y = 1)
			MediatorInsert(m, inputima[xx]);
		ymax = h - 1;
		xx = x + 2 * w;
		size_t medidx = x + w;
		for(y = 1; y < ymax; ++y, xx += w, medidx += w){
			MediatorInsert(m, inputima[xx]);
			Item s, l, md, I = inputima[medidx];
			md = MediatorStat(m, &s, &l);
			s += DBL_EPSILON, l -= DBL_EPSILON;
			if(s < md && md < l){
				if(s < I && I < l) med[medidx] = I;
				else med[medidx] = md;
			}else{
				med[medidx] = adp_med_5by5(img, x, y);
			}
		}
		free(m);
	}
	DBG("time for median filtering by cross 3x3 of image %zdx%zd: %gs", w, h,
		dtime() - t0);
}
/**
 * filter image by median (seed*2 + 1) x (seed*2 + 1)
 */
IMAGE *get_adaptive_median(IMAGE *img, int seed){
	size_t w = img->width, h = img->height, siz = w*h, bufsiz = siz*sizeof(int);
	IMAGE *out = MALLOC(IMAGE, 1);
	Item *med = MALLOC(Item, bufsiz), *inputima = img->data;
	out->data = med;
	out->width = w;
	out->height = h;
	out->dtype = img->dtype;

	memcpy(med, inputima, bufsiz);

	if(seed == 0){
		get_adp_median_cross(img, out);
		return out;
	}

	size_t blksz = seed * 2 + 1, fullsz = blksz * blksz;
#ifdef EBUG
	double t0 = dtime();
#endif
	OMP_FOR(shared(inputima, med))
	for(size_t x = seed; x < w - seed; ++x){
		size_t xx, yy, xm = x + seed + 1, y, ymax = blksz - 1, xmin = x - seed;
		Mediator* m = MediatorNew(fullsz);
		// initial fill
		for(yy = 0; yy < ymax; ++yy)
			for(xx = xmin; xx < xm; ++xx)
				MediatorInsert(m, inputima[xx + yy*w]);
		ymax = 2*seed*w;
		xmin += ymax;
		xm += ymax;
		ymax = h - seed;
		size_t curpos = x + seed * w;
		for(y = seed; y < ymax; ++y, xmin += w, xm += w, curpos += w){
			for(xx = xmin; xx < xm; ++xx)
				MediatorInsert(m, inputima[xx]);
			Item s, l, md, I = inputima[curpos];
			md = MediatorStat(m, &s, &l);
			s += DBL_EPSILON, l -= DBL_EPSILON;
			if(s < md && md < l){
				if(s < I && I < l) med[curpos] = I;
				else med[curpos] = md;
			}else{
				med[curpos] = adp_med_5by5(img, x, y);
			}
		}
		free(m);
	}
	DBG("time for adadptive median filtering %zdx%zd of image %zdx%zd: %gs", blksz, blksz, w, h,
		dtime() - t0);
	return out;
}


