/*
 * fits.c - cfitsio routines
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

#include <string.h>
#include <errno.h>

#include "fits.h"
#include "types.h"
#include "usefull_macros.h"
#include "cmdlnopts.h"

static int fitsstatus = 0;

/*
 * Macros for error processing when working with cfitsio functions
 */
// Try to run function f with arguments
#define TRYFITS(f, ...)							\
do{ fitsstatus = 0;								\
	f(__VA_ARGS__, &fitsstatus);				\
	if(fitsstatus){								\
		fits_report_error(stderr, fitsstatus);}	\
}while(0)
// the same as previous but with checking of funtion return
#define FITSFUN(f, ...)							\
do{ fitsstatus = 0;								\
	int ret = f(__VA_ARGS__, &fitsstatus);		\
	if(ret || fitsstatus){						\
		fits_report_error(stderr, fitsstatus);	\
		if(fitsstatus == 0) fitsstatus = ret;}	\
}while(0)
// write a keyword
#define WRITEKEY(...)							\
do{ fitsstatus = 0;								\
	fits_write_key(__VA_ARGS__, &fitsstatus);	\
	if(status) fits_report_error(stderr, status);\
}while(0)

KeyList *list_get_end(KeyList *list){
	if(!list) return NULL;
	KeyList *first = list;
	if(list->last) return list->last;
	while(list->next) list = list->next;
	//while(first->next){
		first->last = list;
	/*	first = first->next;
	}*/
	return list;
}

/**
 * add record to keylist
 * @param list (io) - pointer to root of list or NULL
 * 						if *root == NULL, just created node will be placed there
 * @param rec       - data inserted
 * @return pointer to created node
 */
KeyList *list_add_record(KeyList **list, char *rec){
	KeyList *node, *last;
	if((node = (KeyList*) MALLOC(KeyList, 1)) == 0)  return NULL; // allocation error
	node->record = strdup(rec); // insert data
	if(!node->record){
		/// "Не могу скопировать данные"
		WARNX(_("Can't copy data"));
		return NULL;
	}
	if(list){
		if(*list){ // there was root node - search last
			last = list_get_end(*list);
			last->next = node; // insert pointer to new node into last element in list
			(*list)->last = node;
		//	DBG("last node %s", (*list)->last->record);
		}else *list = node;
	}
	return node;
}

/**
 * return record with given key or NULL
 */
KeyList *list_find_key(KeyList *list, char *key){
	if(!list || !key) return NULL;
	size_t L = strlen(key);
	do{
		if(list->record){
			if(strncmp(list->record, key, L) == 0){ // key found
				return list;
			}
		}
		list = list->next;
	}while(list);
	return NULL;
}

/**
 * modify key value
 * return NULL if given key is absent
 */
KeyList *list_modify_key(KeyList *list, char *key, char *newval){
	char buf[FLEN_CARD];
	KeyList *rec = list_find_key(list, key);
	if(!rec) return NULL;
	char *comm = strchr(rec->record, '/');
	if(!comm) comm = "";
	snprintf(buf, FLEN_CARD, "%-8s=%21s %s", key, newval, comm);
	FREE(rec->record);
	rec->record = strdup(buf);
	return rec;
}

/**
 * remove record by key
 */
void list_remove_key(KeyList **keylist, char *key){
	if(!keylist || !key) return;
	size_t L = strlen(key);
	KeyList *prev = NULL, *list = *keylist, *last = list_get_end(list);
	do{
		if(list->record){
			if(strncmp(list->record, key, L) == 0){ // key found
				if(prev){ // not first record
					prev->next = list->next;
				}else{ // first or only record
					if(*keylist == last){
						*keylist = NULL; // the only record - erase it
					}else{ // first record - modyfy heading record
						*keylist = list->next;
						(*keylist)->last = last;
					}
				}
				DBG("remove record by key \"%s\":\n%s",key, list->record);
				FREE(list->record);
				FREE(list);
				return;
			}
		}
		prev = list;
		list = list->next;
	}while(list);
}
/**
 * remove records by any sample
 */
void list_remove_records(KeyList **keylist, char *sample){
	if(!keylist || !sample) return;
	KeyList *prev = NULL, *list = *keylist, *last = list_get_end(list);
	DBG("remove %s", sample);
	do{
		if(list->record){
			if(strstr(list->record, sample)){ // key found
				if(prev){
					prev->next = list->next;
				}else{
					if(*keylist == last){
						*keylist = NULL; // the only record - erase it
					}else{ // first record - modyfy heading record
						*keylist = list->next;
						(*keylist)->last = last;
					}
				}
				KeyList *tmp = list->next;
				FREE(list->record);
				FREE(list);
				list = tmp;
				continue;
			}
		}
		prev = list;
		list = list->next;
	}while(list);
}
/**
 * free list memory & set it to NULL
 */
void list_free(KeyList **list){
	KeyList *node = *list, *next;
	if(!list || !*list) return;
	do{
		next = node->next;
		FREE(node->record);
		free(node);
		node = next;
	}while(node);
	*list = NULL;
}
/**
 * make a full copy of given list
 */
KeyList *list_copy(KeyList *list){
	if(!list) return NULL;
	KeyList *newlist = NULL;
	#ifdef EBUG
	int n = 0;
	#endif
	do{
		list_add_record(&newlist, list->record);
		list = list->next;
		#ifdef EBUG
		++n;
		#endif
	}while(list);
	DBG("copy list of %d entries", n);
	return newlist;
}

void list_print(KeyList *list){
	while(list){
		printf("%s\n", list->record);
		list = list->next;
	}
}

void imfree(IMAGE **img){
	list_free(&(*img)->keylist);
	FREE((*img)->data);
	if((*img)->tables){
		size_t i, N = (*img)->tables->amount;
		FITStable **tbl = (*img)->tables->tables;
		for(i = 0; i < N; ++i, ++tbl){
			FITStable *cur = *tbl;
			tablefree(&cur);
		}
		FREE((*img)->tables->tables);
		FREE((*img)->tables);
	}
	FREE(*img);
}

void tablefree(FITStable **tbl){
	if(!tbl || !*tbl) return;
	FITStable *intab = *tbl;
	size_t i, N = intab->ncols;
	for(i = 0; i < N; ++i){
		table_column *col = &(intab->columns[i]);
		if(col->coltype == TSTRING && col->width){
			size_t r, R = col->repeat;
			void **cont = (void**) col->contents;
			for(r = 0; r < R; ++r) free(*(cont++));
		}
		FREE(col->contents);
		FREE(col);
	}
	FREE(*tbl);
}

/**
 * make full copy of FITStables structure & return pointer to it
 */
FITStables *table_copy(FITStables *intab){
	if(!intab || intab->amount == 0) return NULL;
	FITStables *tbl = MALLOC(FITStables, 1);
	size_t N = intab->amount, i;
	tbl->tables = MALLOC(FITStable*, N);
	for(i = 0; i < N; ++i){
		FITStable *cur = MALLOC(FITStable, 1);
		FITStable *in = intab->tables[i];
		memcpy(cur, in, sizeof(FITStable));
		size_t ncols = in->ncols, col;
		cur->columns = MALLOC(table_column, ncols);
		memcpy(cur->columns, in->columns, sizeof(table_column)*ncols);
		table_column *ocurcol = cur->columns, *icurcol = in->columns;
		for(col = 0; col < ncols; ++col, ++ocurcol, ++icurcol){
			if(ocurcol->coltype == TSTRING && ocurcol->width){ // string array - copy all
				size_t r, R = ocurcol->repeat;
				char **oarr = (char**)ocurcol->contents, **iarr = (char**)icurcol->contents;
				for(r = 0; r < R; ++r, ++oarr, ++iarr) *oarr = strdup(*iarr);
			}else memcpy(ocurcol->contents, icurcol->contents, icurcol->repeat * icurcol->width);
		}
		tbl->tables[i] = cur;
	}
	return tbl;
}

/**
 * add FITS table to image structure
 */
FITStable *table_read(IMAGE *img, fitsfile *fp){
	int ncols, i;
	long nrows;
	char extname[FLEN_VALUE];
	#define TRYRET(f, ...) do{TRYFITS(f, fp, __VA_ARGS__); if(fitsstatus) goto ret;}while(0)
	TRYRET(fits_get_num_rows, &nrows);
	TRYRET(fits_get_num_cols, &ncols);
	TRYRET(fits_read_key, TSTRING, "EXTNAME", extname, NULL);
	DBG("Table named %s with %ld rows and %d columns", extname, nrows, ncols);
	FITStable *tbl = table_new(img, extname);
	if(!tbl) return NULL;
	for(i = 1; i <= ncols; ++i){
		int typecode;
		long repeat, width;
		FITSFUN(fits_get_coltype, fp, i, &typecode, &repeat, &width);
		if(fitsstatus){
			WARNX(_("Can't read column %d!"), i);
			continue;
		}
		DBG("typecode=%d, repeat=%ld, width=%ld", typecode, repeat, width);
		table_column col = {.repeat = repeat, .width = width, .coltype = typecode};
		void *array = malloc(width*repeat);
		if(!array) ERRX("malloc");
		int anynul;
		int64_t nullval = 0;
		int j;
		for(j = 0; j < repeat; ++j){
			FITSFUN(fits_read_col, fp, typecode, i, j=1, 1, 1, (void*)nullval, array, &anynul);
			if(fitsstatus){
				WARNX(_("Can't read column %d row %d!"), i, j);
				continue;
			}
		}
		DBG("done");
		continue;
		col.contents = array;
		char keyword[FLEN_KEYWORD];
		int stat = 0;
		fits_make_keyn("TTYPE", i, keyword, &stat);
		if(stat){WARNX("???"); stat = 0;}
		fits_read_key(fp, TSTRING, keyword, col.colname, NULL, &stat);
		if(stat){ sprintf(col.colname, "noname"); stat = 0;}
		fits_make_keyn("TUNIT", i, keyword, &stat);
		if(stat){WARNX("???"); stat = 0;}
		fits_read_key(fp, TSTRING, keyword, col.unit, NULL, &stat);
		if(stat) *col.unit = 0;
		DBG("Column, cont[2]=%d, type=%d, w=%ld, r=%ld, nm=%s, u=%s", ((int*)col.contents)[2], col.coltype,
			col.width, col.repeat, col.colname, col.unit);
		//table_addcolumn(tbl, &col);
		FREE(array);
	}
//	DBG("fits_create_tbl nrows=%zd, ncols=%zd, colnms[0]=%s, formats[0]=%s, "
//			"units[0]=%s, name=%s", tbl->nrows, cols, tbl->colnames[0], tbl->formats[0], tbl->units[0], tbl->tabname);

	#undef TRYRET
ret:
	if(fitsstatus){
		tablefree(&tbl);
		--img->tables->amount;
	}
	return tbl;
}

/**
 * create empty FITS table for image 'img'
 */
FITStable *table_new(IMAGE *img, char *tabname){
	if(!img->tables) img->tables = MALLOC(FITStables, 1);
	size_t N = ++img->tables->amount;
	if(!(img->tables->tables = realloc(img->tables->tables, sizeof(FITStable*)*N))) ERR("realloc()");
	FITStable *tab = MALLOC(FITStable, 1);
	snprintf(tab->tabname, 80, "%s", tabname);
	DBG("add new table: %s", tabname);
	img->tables->tables[N-1] = tab;
	return tab;
}

/**
 * add to table 'tbl' column 'column'
 * Be carefull all fields of 'column' exept of 'format' should be filled
 * - if data is character array, 'width' should be equal 0
 * - all input data will be copied, so caller should run 'free' after this function!
 */
FITStable *table_addcolumn(FITStable *tbl, table_column *column){
	FNAME();
	if(!tbl || !column || !column->contents) return NULL;
	long nrows = column->repeat;
	int width = column->width;
	if(tbl->nrows < nrows) tbl->nrows = nrows;
	size_t datalen = nrows * width, cols = ++tbl->ncols;
	char *curformat = column->format;
	DBG("add column; width: %d, nrows: %ld, name: %s", width, nrows, column->colname);
	/*void convchar(){ // count maximum length of strings in array
		char **charr = (char**)column->contents, *dptr = charr;
		size_t n, N = column->repeat;
		for(n = 0; n < N; ++n){
			if(strlen(
			if(*dptr++ == 0){ --processed; if(len > maxlen) maxlen = len; len = 0; }
			else{ ++len; }
		}
	}*/
	#define CHKLEN(type) do{if(width != sizeof(type)) datalen = sizeof(type) * nrows;}while(0)
	switch(column->coltype){
		case TBIT:
			snprintf(curformat, FLEN_FORMAT, "%ldX", nrows);
			CHKLEN(int8_t);
		break;
		case TBYTE:
			snprintf(curformat, FLEN_FORMAT, "%ldB", nrows);
			CHKLEN(int8_t);
		break;
		case TLOGICAL:
			snprintf(curformat, FLEN_FORMAT, "%ldL", nrows);
			CHKLEN(int8_t);
		break;
		case TSTRING:
			if(width == 0){
				snprintf(curformat, FLEN_FORMAT, "%ldA", nrows);
				datalen = nrows;
			}else
				snprintf(curformat, FLEN_FORMAT, "%ldA%d", nrows, width);
		break;
		case TSHORT:
			snprintf(curformat, FLEN_FORMAT, "%ldI", nrows);
			CHKLEN(int16_t);
		break;
		case TLONG:
			snprintf(curformat, FLEN_FORMAT, "%ldJ", nrows);
			CHKLEN(int32_t);
		break;
		case TLONGLONG:
			snprintf(curformat, FLEN_FORMAT, "%ldK", nrows);
			CHKLEN(int64_t);
		break;
		case TFLOAT:
			snprintf(curformat, FLEN_FORMAT, "%ldE", nrows);
			CHKLEN(float);
		break;
		case TDOUBLE:
			snprintf(curformat, FLEN_FORMAT, "%ldD", nrows);
			CHKLEN(double);
		break;
		case TCOMPLEX:
			snprintf(curformat, FLEN_FORMAT, "%ldM", nrows);
			if(width != sizeof(float)*2) datalen = sizeof(float) * nrows * 2;
		break;
		case TDBLCOMPLEX:
			snprintf(curformat, FLEN_FORMAT, "%ldM", nrows);
			if(width != sizeof(double)*2) datalen = sizeof(double) * nrows * 2;
		break;
		case TINT:
			snprintf(curformat, FLEN_FORMAT, "%ldJ", nrows);
			CHKLEN(int32_t);
		break;
		case TSBYTE:
			snprintf(curformat, FLEN_FORMAT, "%ldS", nrows);
			CHKLEN(int8_t);
		break;
		case TUINT:
			snprintf(curformat, FLEN_FORMAT, "%ldV", nrows);
			CHKLEN(int32_t);
		break;
		case TUSHORT:
			snprintf(curformat, FLEN_FORMAT, "%ldU", nrows);
			CHKLEN(int16_t);
		break;
		default:
			WARNX(_("Unsupported column data type!"));
			return NULL;
	}

	DBG("new size: %ld, old: %ld", sizeof(table_column)*cols, sizeof(table_column)*(cols-1));
	if(!(tbl->columns = realloc(tbl->columns, sizeof(table_column)*cols))) ERRX("malloc");
	table_column *newcol = &(tbl->columns[cols-1]);
	memcpy(newcol, column, sizeof(table_column));
	newcol->contents = calloc(datalen, 1);
	if(!newcol->contents) ERRX("malloc");
	DBG("copy %zd bytes", datalen);
	if(column->coltype == TSTRING && width){
		long n;
		char **optr = (char**)newcol->contents, **iptr = (char**)column->contents;
		for(n = 0; n < nrows; ++n, ++optr, ++iptr) *optr = strdup(*iptr);
	}else
		memcpy(newcol->contents, column->contents, datalen);
	return tbl;
}

void table_print(FITStable *tbl){
	printf("\nTable name: %s\n", tbl->tabname);
	int c, cols = tbl->ncols;
	long r, rows = tbl->nrows;
	for(c = 0; c < cols; ++c){
		printf("%s", tbl->columns[c].colname);
		if(*tbl->columns[c].unit) printf(" (%s)", tbl->columns[c].unit);
		printf("\t");
	}
	printf("\n");
	for(r = 0; r < rows; ++r){
		for(c = 0; c < cols; ++c){
			double *dpair; float *fpair;
			table_column *col = &(tbl->columns[c]);
			if(col->repeat < r){ // table with columns of different length
				printf("(empty)\t");
				continue;
			}
			switch(col->coltype){
				case TBIT:
				case TBYTE:
					printf("%u\t", ((uint8_t*)col->contents)[r]);
				break;
				case TLOGICAL:
					printf("%s\t", ((int8_t*)col->contents)[r] == 0 ? "FALSE" : "TRUE");
				break;
				case TSTRING:
					if(col->width == 0) printf("%c\t", ((char*)col->contents)[r]);
					else printf("%s\t", ((char**)col->contents)[r]);
				break;
				case TSHORT:
					printf("%d\t", ((int16_t*)col->contents)[r]);
				break;
				case TLONG:
				case TINT:
					printf("%d\t", ((int32_t*)col->contents)[r]);
				break;
				case TLONGLONG:
					printf("%zd\t", ((int64_t*)col->contents)[r]);
				break;
				case TFLOAT:
					printf("%g\t", ((float*)col->contents)[r]);
				break;
				case TDOUBLE:
					printf("%g\t", ((double*)col->contents)[r]);
				break;
				case TCOMPLEX:
					fpair = (float*)col->contents + 2*r;
					printf("%g %s %g*i\t", fpair[0], fpair[1] > 0 ? "+" : "-", fpair[1]);
				break;
				case TDBLCOMPLEX:
					dpair = (double*)col->contents + 2*r;
					printf("%g %s %g*i\t", dpair[0], dpair[1] > 0 ? "+" : "-", dpair[1]);
				break;
				case TSBYTE:
					printf("%d\t", ((int8_t*)col->contents)[r]);
				break;
				case TUINT:
					printf("%d\t", ((uint32_t*)col->contents)[r]);
				break;
				case TUSHORT:
					printf("%d\t", ((uint16_t*)col->contents)[r]);
				break;
			}
		}
		printf("\n");
	}
}

void table_print_all(IMAGE *img){
	if(!img->tables || img->tables->amount < 1) return;
	size_t i, N = img->tables->amount;
	for(i = 0; i < N; ++i)
		table_print(img->tables->tables[i]);
}

/**
 * save all tables of given image into file
 */
void table_write(IMAGE *img, fitsfile *fp){
	FNAME();
	if(!img->tables || img->tables->amount == 0) return;
	size_t N = img->tables->amount, i;
	for(i = 0; i < N; ++i){
		FITStable *tbl = img->tables->tables[i];
		size_t c, cols = tbl->ncols;
		/*int hdutype;
		TRYFITS(fits_movabs_hdu, fp, ++img->lasthdu, &hdutype);
		if(fitsstatus){
			WARNX(_("Can't write table %s - cannot move to HDU #%d!"), tbl->tabname, img->lasthdu);
			continue;
		}*/
		char **columns = MALLOC(char*, cols);
		char **formats = MALLOC(char*, cols);
		char **units = MALLOC(char*, cols);
		table_column *col = tbl->columns;
		for(c = 0; c < cols; ++c, ++col){
			columns[c] = col->colname;
			formats[c] = col->format;
			units[c] = col->unit;
			DBG("col: %s, form: %s, unit: %s", columns[c], formats[c], units[c]);
		}
		FITSFUN(fits_create_tbl, fp, BINARY_TBL, tbl->nrows, cols,
			columns, formats, units, tbl->tabname);
		FREE(columns); FREE(formats); FREE(units);
		if(fitsstatus){
			WARNX(_("Can't write table %s!"), tbl->tabname);
			continue;
		}
		col = tbl->columns;
		for(c = 0; c < cols; ++c, ++col){
			DBG("write column %zd", c);
			TRYFITS(fits_write_col, fp, col->coltype, c+1, 1, 1, col->repeat, col->contents);
			if(fitsstatus){
				WARNX(_("Can't write column %s!"), col->colname);
			continue;
			}
		}
		/*int hdutype;
		TRYFITS(fits_movabs_hdu, fp, ++img->lasthdu, &hdutype);
		if(fitsstatus){
			WARNX(_("Can't write table %s - cannot move to HDU #%d!"), tbl->tabname, img->lasthdu);
			continue;
		}*/
	}
}

/**
 * read FITS file and fill 'IMAGE' structure (with headers and tables)
 * can't work with image stack - opens the first image met
 * works only with binary tables
 */
IMAGE* readFITS(char *filename, IMAGE **fits){
	FNAME();
	#define TRYRET(f, ...) do{TRYFITS(f, __VA_ARGS__); if(fitsstatus) goto returning;}while(0)
	fitsfile *fp;
	int i, j, hdunum = 0, hdutype, nkeys, keypos;
	int naxis;
	long naxes[2];
	char card[FLEN_CARD];
	IMAGE *img = MALLOC(IMAGE, 1);
	TRYRET(fits_open_file, &fp, filename, READONLY);
	FITSFUN(fits_get_num_hdus, fp, &hdunum);
	if(hdunum < 1){
		WARNX(_("Can't read HDU"));
		fitsstatus = 1;
		goto returning;
	}
	// get image dimensions
	TRYRET(fits_get_img_param, fp, 2, &img->dtype, &naxis, naxes);
	if(naxis > 2){
		WARNX(_("Images with > 2 dimensions are not supported"));
		fitsstatus = 1;
		goto returning;
	}
	img->width = naxes[0];
	img->height = naxes[1];
	DBG("got image %ldx%ld pix, bitpix=%d", naxes[0], naxes[1], img->dtype);
	// loop through all HDUs
	KeyList *list = img->keylist;
	int imghdu = -1;
	for(i = 1; !(fits_movabs_hdu(fp, i, &hdutype, &fitsstatus)); ++i){
		int hdutype;
		TRYFITS(fits_get_hdu_type, fp, &hdutype);
		// types: IMAGE_HDU , ASCII_TBL, BINARY_TBL
		if(fitsstatus) continue;
		DBG("HDU type %d", hdutype);
		//if(hdutype == ASCII_TBL || hdutype == BINARY_TBL){
		if(hdutype == BINARY_TBL){
			table_read(img, fp);
			continue;
		}
		if(imghdu < 1) imghdu = i;
		TRYFITS(fits_get_hdrpos, fp, &nkeys, &keypos);
		if(fitsstatus) continue;
		//DBG("HDU # %d of %d keys", i, nkeys);
		for(j = 1; j <= nkeys; ++j){
			FITSFUN(fits_read_record, fp, j, card);
			if(!fitsstatus){
				if(!list_add_record(&list, card)){
					/// "Не могу добавить запись в список"
					WARNX(_("Can't add record to list"));
				}
				//DBG("key %d: %s", j, card);
			}
		}
	}
	img->keylist = list;
	if(fitsstatus == END_OF_FILE){
		fitsstatus = 0;
	}else{
		fits_report_error(stderr, fitsstatus);
		goto returning;
	}
	if(fits_movabs_hdu(fp, imghdu, &hdutype, &fitsstatus)){
		WARNX(_("Can't open image HDU #%d"), imghdu);
		fitsstatus = 1;
		goto returning;
	}
	size_t sz = naxes[0] * naxes[1];
	img->data = MALLOC(double, sz);
	int stat = 0;
	TRYFITS(fits_read_img, fp, TDOUBLE, 1, sz, NULL, img->data, &stat);
	if(stat) WARNX(_("Found %d pixels with undefined value"), stat);
	DBG("ready");
	#undef TRYRET
returning:
	if(fitsstatus){
		imfree(&img);
	}
	FITSFUN(fits_close_file, fp);
	if(fits){
		FREE(*fits);
		*fits = img;
	}
	return img;
}

bool writeFITS(char *filename, IMAGE *fits){
	if(!filename || !fits) return FALSE;
	int w = fits->width, h = fits->height;
	long naxes[2] = {w, h};
	size_t sz = w * h;
	fitsfile *fp;
	TRYFITS(fits_create_file, &fp, filename);
	if(fitsstatus) return FALSE;
	TRYFITS(fits_create_img, fp, fits->dtype, 2, naxes);
	if(fitsstatus) return FALSE;
	if(fits->keylist){ // there's keys
		KeyList *records = fits->keylist;
		while(records){
			char *rec = records->record;
			records = records->next;
			if(strncmp(rec, "SIMPLE", 6) == 0 || strncmp(rec, "EXTEND", 6) == 0) // key "file does conform ..."
				continue;
				// comment of obligatory key in FITS head
			else if(strncmp(rec, "COMMENT   FITS", 14) == 0 || strncmp(rec, "COMMENT   and Astrophysics", 26) == 0)
				continue;
			else if(strncmp(rec, "NAXIS", 5) == 0 || strncmp(rec, "BITPIX", 6) == 0) // NAXIS, NAXISxxx, BITPIX
				continue;
			FITSFUN(fits_write_record, fp, rec);
		//	DBG("write key: %s", rec);
		}
	}
	//fits->lasthdu = 1;
	//FITSFUN(fits_write_record, fp, "COMMENT  modified by simple test routine");
	TRYFITS(fits_write_img, fp, TDOUBLE, 1, sz, fits->data);
	if(fitsstatus) return FALSE;
	if(fits->tables && !G.deltabs) table_write(fits, fp);
	TRYFITS(fits_close_file, fp);
	return TRUE;
}

/**
 * create an empty image without headers, assign data type to "dtype"
 */
IMAGE *newFITS(size_t h, size_t w, int dtype){
	size_t bufsiz = w*h;
	IMAGE *out = MALLOC(IMAGE, 1);
	out->data = MALLOC(Item, bufsiz);
	out->width = w;
	out->height = h;
	out->dtype = dtype;
	return out;
}

/**
 * build IMAGE image from data array indata
 */
IMAGE *buildFITSfromdat(size_t h, size_t w, int dtype, uint8_t *indata){
	size_t stride = 0;
	double (*fconv)(uint8_t *data) = NULL;
	double ubyteconv(uint8_t *data){return (double)*data;}
	double ushortconv(uint8_t *data){return (double)*(int16_t*)data;}
	double ulongconv(uint8_t *data){return (double)*(uint32_t*)data;}
	double ulonglongconv(uint8_t *data){return (double)*(uint64_t*)data;}
	double floatconv(uint8_t *data){return (double)*(float*)data;}
	IMAGE *out = newFITS(h, w, dtype);
	switch (dtype){
		case BYTE_IMG:
			stride = 1;
			fconv = ubyteconv;
		break;
		case SHORT_IMG:
			stride = 2;
			fconv = ushortconv;
		break;
		case LONG_IMG:
			stride = 4;
			fconv = ulongconv;
		break;
		case FLOAT_IMG:
			stride = 4;
			fconv = floatconv;
		break;
		case LONGLONG_IMG:
			fconv = ulonglongconv;
			stride = 8;
		break;
		case DOUBLE_IMG:
			memcpy(out->data, indata, sizeof(double)*w*h);
			return out;
		break;
		default:
		/// Неправильный тип данных
			ERRX(_("Wrong data type"));
	}
	size_t y, W = w*stride;
	Item *data = out->data;
	OMP_FOR(shared(data))
	for(y = 0; y < h; ++y){
		Item *dout = &data[y*w];
		uint8_t *din = &indata[y*W];
		size_t x;
		for(x = 0; x < w; ++x, din += stride)
			*dout++ = fconv(din);
	}
	return out;
}

/**
 * create an empty copy of image "in" without headers, assign data type to "dtype"
 */
IMAGE *similarFITS(IMAGE *img, int dtype){
	size_t w = img->width, h = img->height, bufsiz = w*h;
	IMAGE *out = MALLOC(IMAGE, 1);
	out->data = MALLOC(Item, bufsiz);
	out->width = w;
	out->height = h;
	out->dtype = dtype;
	return out;
}

/**
 * make full copy of image 'in'
 */
IMAGE *copyFITS(IMAGE *in){
	IMAGE *out = similarFITS(in, in->dtype);
	memcpy(out->data, in->data, sizeof(Item)*in->width*in->height);
	out->keylist = list_copy(in->keylist);
	out->tables = table_copy(in->tables);
	return out;
}

/*
 * Different file functions
 */
/**
 * Return TRUE if file _name_ not exists
 */
bool file_absent(char *name){
	struct stat filestat;
	if(!stat(name, &filestat)) return FALSE;
	if(errno == ENOENT) return TRUE;
	return FALSE;
}
/**
 * find the first non-exists filename like prefixXXXX.suffix & put it into buff
 */
char* make_filename(char *buff, size_t buflen, char *prefix, char *suffix){
	int num;
	for(num = 1; num < 10000; ++num){
		if(snprintf(buff, buflen, "%s_%04d.%s", prefix, num, suffix) < 1)
			return NULL;
		if(file_absent(buff)) // OK, file not exists
			return buff;
	}
	return NULL;
}
