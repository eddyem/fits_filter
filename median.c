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
#include <assert.h>

#include "median.h"
#include "usefull_macros.h"

// largest radius for adaptive median filter - 7x7 pix
#define LARGEST_ADPMED_RADIUS  (3)

#define ELEM_SWAP(a, b) {register Item t = a; a = b; b = t;}
#define PIX_SORT(a, b)  {if (p[a] > p[b]) ELEM_SWAP(p[a], p[b]);}

Item opt_med2(Item *p){
	return (p[0] + p[1]) * 0.5;
}
Item opt_med3(Item *p){
	PIX_SORT(0, 1); PIX_SORT(1, 2); PIX_SORT(0, 1);
	return(p[1]) ;
}
Item opt_med4(Item *p){
	PIX_SORT(0, 2); PIX_SORT(1, 3);
	PIX_SORT(0, 1); PIX_SORT(2, 3);
	return(p[1] + p[2]) * 0.5;
}
Item opt_med5(Item *p){
	PIX_SORT(0, 1); PIX_SORT(3, 4); PIX_SORT(0, 3);
	PIX_SORT(1, 4); PIX_SORT(1, 2); PIX_SORT(2, 3) ;
	PIX_SORT(1, 2);
	return(p[2]) ;
}
// even values are from "FAST, EFFICIENT MEDIAN FILTERS WITH EVEN LENGTH WINDOWS", J.P. HAVLICEK, K.A. SAKADY, G.R.KATZ
Item opt_med6(Item *p){
	PIX_SORT(1, 2); PIX_SORT(3, 4);
	PIX_SORT(0, 1); PIX_SORT(2, 3); PIX_SORT(4, 5);
	PIX_SORT(1, 2); PIX_SORT(3, 4);
	PIX_SORT(0, 1); PIX_SORT(2, 3); PIX_SORT(4, 5);
	PIX_SORT(1, 2); PIX_SORT(3, 4);
	return ( p[2] + p[3] ) * 0.5;
}
Item opt_med7(Item *p){
	PIX_SORT(0, 5); PIX_SORT(0, 3); PIX_SORT(1, 6);
	PIX_SORT(2, 4); PIX_SORT(0, 1); PIX_SORT(3, 5);
	PIX_SORT(2, 6); PIX_SORT(2, 3); PIX_SORT(3, 6);
	PIX_SORT(4, 5); PIX_SORT(1, 4); PIX_SORT(1, 3);
	PIX_SORT(3, 4); return (p[3]);
}
// optimal Batcher's sort for 8 elements (http://myopen.googlecode.com/svn/trunk/gtkclient_tdt/include/fast_median.h)
Item opt_med8(Item *p){
	PIX_SORT(0, 4); PIX_SORT(1, 5); PIX_SORT(2, 6);
	PIX_SORT(3, 7); PIX_SORT(0, 2); PIX_SORT(1, 3);
	PIX_SORT(4, 6); PIX_SORT(5, 7); PIX_SORT(2, 4);
	PIX_SORT(3, 5); PIX_SORT(0, 1); PIX_SORT(2, 3);
	PIX_SORT(4, 5); PIX_SORT(6, 7); PIX_SORT(1, 4);
	PIX_SORT(3, 6);
	return(p[3] + p[4]) * 0.5;
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
Item opt_med16(Item *p){
	PIX_SORT(0, 8); PIX_SORT(1, 9); PIX_SORT(2, 10); PIX_SORT(3, 11);
	PIX_SORT(4, 12); PIX_SORT(5, 13); PIX_SORT(6, 14); PIX_SORT(7, 15);
	PIX_SORT(0, 4); PIX_SORT(1, 5); PIX_SORT(2, 6); PIX_SORT(3, 7);
	PIX_SORT(8, 12); PIX_SORT(9, 13); PIX_SORT(10, 14); PIX_SORT(11, 15);
	PIX_SORT(4, 8); PIX_SORT(5, 9); PIX_SORT(6, 10); PIX_SORT(7, 11);
	PIX_SORT(0, 2); PIX_SORT(1, 3); PIX_SORT(4, 6); PIX_SORT(5, 7);
	PIX_SORT(8, 10); PIX_SORT(9, 11); PIX_SORT(12, 14); PIX_SORT(13, 15);
	PIX_SORT(2, 8); PIX_SORT(3, 9); PIX_SORT(6, 12); PIX_SORT(7, 13);
	PIX_SORT(2, 4); PIX_SORT(3, 5); PIX_SORT(6, 8); PIX_SORT(7, 9);
	PIX_SORT(10, 12); PIX_SORT(11, 13); PIX_SORT(0, 1); PIX_SORT(2, 3);
	PIX_SORT(4, 5); PIX_SORT(6, 7); PIX_SORT(8, 9); PIX_SORT(10, 11);
	PIX_SORT(12, 13); PIX_SORT(14, 15); PIX_SORT(1, 8); PIX_SORT(3, 10);
	PIX_SORT(5, 12); PIX_SORT(7, 14); PIX_SORT(5, 8); PIX_SORT(7, 10);
	return (p[7] + p[8]) * 0.5;
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
#undef PIX_SORT
#define PIX_SORT(a, b)  {if (a > b) ELEM_SWAP(a, b);}
/**
 * quick select - algo for approximate median calculation for array idata of size n
 */
Item quick_select(Item *idata, int n){
	int low, high;
	int median;
	int middle, ll, hh;
	Item *arr = MALLOC(Item, n);
	memcpy(arr, idata, n*sizeof(Item));
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
			do ll++; while (arr[low] > arr[ll]);
			do hh--; while (arr[hh] > arr[low]);
			if(hh < ll) break;
			ELEM_SWAP(arr[ll], arr[hh]) ;
		}
		// Swap middle item (in position low) back into correct position
		ELEM_SWAP(arr[low], arr[hh]) ;
		// Re-set active partition
		if (hh <= median) low = ll;
		if (hh >= median) high = hh - 1;
	}
	Item ret = arr[median];
	FREE(arr);
	return ret;
}
#undef PIX_SORT
#undef ELEM_SWAP

/**
 * calculate median of array idata with size n
 */
Item calc_median(Item *idata, int n){
	assert(idata); assert(n>0);
	typedef Item (*medfunc)(Item *p);
	medfunc fn = NULL;
	const medfunc fnarr[] = {opt_med2, opt_med3, opt_med4, opt_med5, opt_med6,
			opt_med7, opt_med8, opt_med9};
	if(n == 1) return *idata;
	if(n < 10) fn = fnarr[n - 1];
	else if(n == 16) fn = opt_med16;
	else if(n == 25) fn = opt_med25;
	if(fn){
		return fn(idata);
	}else{
		return quick_select(idata, n);
	}
}

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

static void get_adp_median_cross(IMAGE *img, IMAGE *out, int adp);

/**
 * filter image by median (seed*2 + 1) x (seed*2 + 1)
 */
IMAGE *get_median(IMAGE *img, Filter *f, _U_ Itmarray *i){
	int seed = f->w;
	size_t w = img->width, h = img->height;
	IMAGE *out = copyFITS(img);
	Item *med = out->data, *inputima = img->data;
	if(seed == 0){
		get_adp_median_cross(img, out, 0);
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
	size_t blocklen, w = img->width, h = img->height, yy, _2w = 2 * w;
	Item arr[25], *arrptr = arr, *dataptr, *currpix;
	int position = ((x < 1) ? 1 : 0)        // left columns
				 + ((x > w - 2) ? 2 : 0)    // right columns
				 + ((y < 1) ? 4 : 0)        // top rows
				 + ((y > w - 2) ? 8 : 0);   // bottom rows
	/* Now by value of "position" we know where is the point:
	 ***************************
	 * 5 *        4        * 6 *
	 ***************************
	 *   *                 *   *
	 *   *                 *   *
	 * 1 *        0        * 2 *
	 *   *                 *   *
	 *   *                 *   *
	 ***************************
	 * 9 *        8        *10 *
	 ***************************/
	currpix = &img->data[x + y * w]; // pointer to current pixel
	dataptr = currpix - _2w - 2;     // pointer to left upper corner of 5x5 square
	inline void copy5times(Item val){
		for(int i = 0; i < 5; ++i) *arrptr++ = val;
	}
	inline void copy9times(Item val){
		for(int i = 0; i < 9; ++i) *arrptr++ = val;
	}
	void copycolumn(Item *startpix){
		for(int i = 0; i < 5; ++i, startpix += w) *arrptr++ = *startpix;
	}
	inline void copyvertblock(size_t len){
		for(int i = 0; i < 5; ++i, dataptr += w, arrptr += len)
			memcpy(arrptr, dataptr, len * sizeof(Item));
	}
	inline void copyhorblock(size_t len){
		for(size_t i = 0; i < len; ++i, dataptr += w, arrptr += 5)
			memcpy(arrptr, dataptr, 5 * sizeof(Item));
	}
	inline void copyblock(){
		for(size_t i = 0; i < 4; ++i, dataptr += w, arrptr += 4)
			memcpy(arrptr, dataptr, 4 * sizeof(Item));
	}
	switch(position){
		case 1: // left
			copy5times(*currpix); // make 5 copies of current pixel
			if(x == 0){ // copy 1st column too
				dataptr += 2;
				copycolumn(dataptr);
				blocklen = 3;
			}else{ // 2nd column - no copy need
				++dataptr;
				blocklen = 4;
			}
			copyvertblock(blocklen);
		break;
		case 2: // right
			copy5times(*currpix);
			if(x == w - 1){ // copy last column too
				copycolumn(dataptr + 2);
				blocklen = 3;
			}else{ // 2nd column - no copy need
				blocklen = 4;
			}
			copyvertblock(blocklen);
		break;
		case 4: // top
			copy5times(*currpix);
			if(y == 0){
				dataptr += _2w;
				memcpy(arrptr, dataptr, 5 * sizeof(Item));
				blocklen = 3;
			}else{
				dataptr += w;
				blocklen = 4;
			}
			copyhorblock(blocklen);
		break;
		case 8: // bottom
			copy5times(*currpix);
			if(y == h - 1){
				memcpy(arrptr, dataptr + _2w, 5 * sizeof(Item));
				blocklen = 3;
			}else{
				blocklen = 4;
			}
			copyhorblock(blocklen);
		break;
		case 5: // top left corner: in all corners we just copy 4x4 square & 9 times this pixel
			copy9times(*currpix);
			dataptr = img->data;
			copyblock();
		break;
		case 6: // top right corner
			copy9times(*currpix);
			dataptr = &img->data[w - 4];
			copyblock();
		break;
		case 9: // bottom left cornet
			copy9times(*currpix);
			dataptr = &img->data[(y - 4) * w];
			copyblock();
		break;
		case 10: // bottom right cornet
			copy9times(*currpix);
			dataptr = &img->data[(y - 3) * w - 4];
			copyblock();
		break;
		default:  // 0
			for(yy = 0; yy < 5; ++yy, dataptr += w, arrptr += 5)
				memcpy(arrptr, dataptr, 5*sizeof(Item));
	}
	return opt_med25(arr);
}

/**
 * Adaptive median by cross 3x3
 * We have 5 datapoints and 4 inserts @ each step, so
 * better to use opt_med5 instead of Mediator
 * @param adp == 1 for adaptive filtering
 */
static void get_adp_median_cross(IMAGE *img, IMAGE *out, int adp){
	size_t w = img->width, h = img->height;
	Item *med = out->data, *inputima = img->data, *iptr;
#ifdef EBUG
	double t0 = dtime();
#endif
	OMP_FOR(shared(inputima, med))
	for(size_t x = 1; x < w - 1; ++x){
		Item buffer[5];
		size_t curpix = x + w, // index of current pixel image arrays
			y, ymax = h - 1;
		for(y = 1; y < ymax; ++y, curpix += w){
			Item md, *I = &inputima[curpix], Ival = *I;
			memcpy(buffer, I - 1, 3*sizeof(Item));
			buffer[3] = I[-w]; buffer[4] = I[w];
			md = opt_med5(buffer);
			if(adp){
				Item s, l;
				s = ITM_EPSILON + MIN(buffer[0], buffer[1]);
				l = MAX(buffer[3], buffer[4]) - ITM_EPSILON;
				if(s < md && md < l){
					if(s < Ival && Ival < l) med[curpix] = Ival;
					else med[curpix] = md;
				}else{
					med[curpix] = adp_med_5by5(img, x, y);
				}
			}else
				med[curpix] = md;
		}
	}
	// process corners (without adaptive)
	Item buf[5];
	// left top
	buf[0] = inputima[0]; buf[1] = inputima[0];
	buf[2] = inputima[1]; buf[3] = inputima[w];
	buf[4] = inputima[w + 1];
	med[0] = opt_med5(buf);
	// right top
	iptr = &inputima[w - 1];
	buf[0] = iptr[0]; buf[1] = iptr[0];
	buf[2] = iptr[-1]; buf[3] = iptr[w - 1];
	buf[4] = iptr[w];
	med[w - 1] = opt_med5(buf);
	// left bottom
	iptr = &inputima[(h - 1) * w];
	buf[0] = iptr[0]; buf[1] = iptr[0];
	buf[2] = iptr[-w]; buf[3] = iptr[1 - w];
	buf[4] = iptr[1];
	med[(h - 1) * w] = opt_med5(buf);
	// right bottom
	iptr = &inputima[h * w - 1];
	buf[0] = iptr[0]; buf[1] = iptr[0];
	buf[2] = iptr[-w-1]; buf[3] = iptr[-w];
	buf[4] = iptr[-1];
	med[h * w - 1] = opt_med5(buf);
	// process borders without corners
	// top
	OMP_FOR(shared(med))
	for(size_t x = 1; x < w - 1; ++x){
		Item *iptr = &inputima[x];
		buf[0] = buf[1] = *iptr;
		buf[2] = iptr[-1]; buf[3] = iptr[2];
		buf[4] = iptr[w];
		med[x] = opt_med5(buf);
	}
	// bottom
	size_t curidx = (h-2)*w;
	OMP_FOR(shared(curidx, med))
	for(size_t x = 1; x < w - 1; --x){
		Item *iptr = &inputima[curidx + x];
		buf[0] = buf[1] = *iptr;
		buf[2] = iptr[-w]; buf[3] = iptr[-1];
		buf[4] = iptr[1];
		med[curidx + x] = opt_med5(buf);
	}
	// left
	OMP_FOR(shared(med))
	for(size_t y = 1; y < h - 1; ++y){
		size_t cur = y * w;
		Item *iptr = &inputima[cur];
		buf[0] = buf[1] = *iptr;
		buf[2] = iptr[-w]; buf[3] = iptr[1];
		buf[4] = iptr[w];
		med[cur] = opt_med5(buf);
	}
	// right
	curidx = w - 1;
	OMP_FOR(shared(curidx, med))
	for(size_t y = 1; y < h - 1; ++y){
		size_t cur = curidx + y * w;
		Item *iptr = &inputima[cur];
		buf[0] = buf[1] = *iptr;
		buf[2] = iptr[-w]; buf[3] = iptr[-1];
		buf[4] = iptr[w];
		med[cur] = opt_med5(buf);
	}
	DBG("time for median filtering by cross 3x3 of image %zdx%zd: %gs", w, h,
		dtime() - t0);
}
/**
 * filter image by median (seed*2 + 1) x (seed*2 + 1)
 */
IMAGE *get_adaptive_median(IMAGE *img, Filter *f, _U_ Itmarray *i){
	int seed = f->w;
	size_t w = img->width, h = img->height, siz = w*h, bufsiz = siz*sizeof(Item);
	IMAGE *out = similarFITS(img, img->dtype);
	Item *med = out->data, *inputima = img->data;
	memcpy(med, inputima, bufsiz);
	if(seed == 0){
		get_adp_median_cross(img, out, 1);
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
			s += ITM_EPSILON, l -= ITM_EPSILON;
			if(s < md && md < l){
				if(s < I && I < l) med[curpos] = I;
				else med[curpos] = md;
			}else{
				if(seed > LARGEST_ADPMED_RADIUS)
					med[curpos] = I;
				else
					med[curpos] = adp_med_5by5(img, x, y);
			}
		}
		free(m);
	}
	DBG("time for adadptive median filtering %zdx%zd of image %zdx%zd: %gs", blksz, blksz, w, h,
		dtime() - t0);
	return out;
}


