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
#include <stdio.h>
#include <sys/stat.h>
#include <errno.h>

#include "fits.h"
#include "types.h"
#include "usefull_macros.h"

static int fitsstatus = 0;

/*
 * Macros for error processing when working with cfitsio functions
 */
#define TRYFITS(f, ...)							\
do{ fitsstatus = 0;								\
	f(__VA_ARGS__, &fitsstatus);				\
	if(fitsstatus){								\
		fits_report_error(stderr, fitsstatus);	\
		return FALSE;}							\
}while(0)
#define FITSFUN(f, ...)							\
do{ fitsstatus = 0;								\
	int ret = f(__VA_ARGS__, &fitsstatus);		\
	if(ret || fitsstatus)						\
		fits_report_error(stderr, fitsstatus);	\
}while(0)
#define WRITEKEY(...)							\
do{ fitsstatus = 0;								\
	fits_write_key(__VA_ARGS__, &fitsstatus);	\
	if(status) fits_report_error(stderr, status);\
}while(0)

KeyList *list_get_end(KeyList *list){
	if(!list) return NULL;
	while(list->next) list = list->next;
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
		}
		else *list = node;
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
	char buf[80];
	KeyList *rec = list_find_key(list, key);
	if(!rec) return NULL;
	char *comm = strchr(rec->record, '/');
	if(!comm) comm = "";
	FREE(rec->record);
	snprintf(buf, 80, "%-8s=%21s %s", key, newval, comm);
	rec->record = strdup(buf);
	DBG("modify: %s", buf);
	return rec;
}

/**
 * remove record by key
 */
void list_remove_key(KeyList **keylist, char *key){
	if(!keylist || !key) return;
	size_t L = strlen(key);
	KeyList *prev = NULL, *list = *keylist;
	do{
		if(list->record){
			if(strncmp(list->record, key, L) == 0){ // key found
				if(prev) prev->next = list->next;
				else *keylist = list->next; // first record
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
	KeyList *prev = NULL, *list = *keylist;
	DBG("remove %s", sample);
	do{
		if(list->record){
			if(strstr(list->record, sample)){ // key found
				if(prev) prev->next = list->next;
				else *keylist = list->next; // first record
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
	KeyList *newlist = NULL, *node = NULL, *nxt;
	do{
		nxt = list_add_record(&node, list->record);
		if(!newlist) newlist = node;
		node = nxt;
		list = list->next;
	}while(list);
	return newlist;
}

void imfree(IMAGE **img){
	list_free(&(*img)->keylist);
	FREE((*img)->data);
	FREE(*img);
}

IMAGE* readFITS(char *filename, IMAGE **fits){
	FNAME();
	bool ret = TRUE;
	fitsfile *fp;
//	float nullval = 0., imBits, bZero = 0., bScale = 1.;
	int i, j, hdunum, hdutype, nkeys, keypos;
	int naxis;
	long naxes[2];
	char card[FLEN_CARD];
	IMAGE *img = MALLOC(IMAGE, 1);

	TRYFITS(fits_open_file, &fp, filename, READONLY);
	FITSFUN(fits_get_num_hdus, fp, &hdunum);
	if(hdunum < 1){
		WARNX(_("Can't read HDU"));
		ret = FALSE;
		goto returning;
	}
	// get image dimensions
	TRYFITS(fits_get_img_param, fp, 2, &img->dtype, &naxis, naxes);
	if(naxis > 2){
		WARNX(_("Images with > 2 dimensions are not supported"));
		ret = FALSE;
		goto returning;
	}
	img->width = naxes[0];
	img->height = naxes[1];
	DBG("got image %ldx%ld pix, bitpix=%d", naxes[0], naxes[1], img->dtype);
	// loop through all HDUs
	KeyList *list = img->keylist;
	for(i = 1; !(fits_movabs_hdu(fp, i, &hdutype, &fitsstatus)); ++i){
		TRYFITS(fits_get_hdrpos, fp, &nkeys, &keypos);
	/*	if(!(img->keylist = realloc(img->keylist, sizeof(char*) * img->keynum))){
			ERR(_("Can't realloc"));
		}*/
	//	char **currec = &(img->keylist[oldnkeys]);
		DBG("HDU # %d of %d keys", i, nkeys);
		for(j = 1; j <= nkeys; ++j){
			FITSFUN(fits_read_record, fp, j, card);
			if(!fitsstatus){
				if(!list_add_record(&list, card)){
					/// "Не могу добавить запись в список"
					ERR(_("Can't add record to list"));
				}
				//*currec = MALLOC(char, FLEN_CARD);
				//memcpy(*currec, card, FLEN_CARD);
				DBG("key %d: %s", j, card);
				//++currec;
			}
		}
		img->keylist = list;
	}
	if(fitsstatus == END_OF_FILE){
		fitsstatus = 0;
	}else{
		fits_report_error(stderr, fitsstatus);
		ret = FALSE;
		goto returning;
	}
	size_t sz = naxes[0] * naxes[1];
	img->data = MALLOC(double, sz);
	int stat = 0;
	TRYFITS(fits_read_img, fp, TDOUBLE, 1, sz, NULL, img->data, &stat);
	if(stat) WARNX(_("Found %d pixels with undefined value"), stat);
	DBG("ready");

returning:
	FITSFUN(fits_close_file, fp);
	if(!ret){
		imfree(&img);
	}
	if(fits) *fits = img;
	return img;
}

bool writeFITS(char *filename, IMAGE *fits){
	int w = fits->width, h = fits->height;
	long naxes[2] = {w, h};
	size_t sz = w * h;
	fitsfile *fp;
	TRYFITS(fits_create_file, &fp, filename);
	TRYFITS(fits_create_img, fp, fits->dtype, 2, naxes);
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
			DBG("write key: %s", rec);
		}
	}
	FITSFUN(fits_write_record, fp, "COMMENT  modified by simple test routine");

	TRYFITS(fits_write_img, fp, TDOUBLE, 1, sz, fits->data);
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
 * make full copi of image 'in'
 */
IMAGE *copyFITS(IMAGE *in){
	IMAGE *out = similarFITS(in, in->dtype);
	memcpy(out->data, in->data, sizeof(Item)*in->width*in->height);
	out->keylist = list_copy(in->keylist);
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
