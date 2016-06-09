/*
 * fits.h
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
#ifndef __FITS_H__
#define __FITS_H__

#include <stdio.h>
#include <stdbool.h>
#include <fitsio.h>
#include <stdint.h>
#include <sys/stat.h>

typedef double Item;

/*
cfitsio.h BITPIX code values for FITS image types:
#define BYTE_IMG      8
#define SHORT_IMG    16
#define LONG_IMG     32
#define LONGLONG_IMG 64
#define FLOAT_IMG   -32
#define DOUBLE_IMG  -64
*/

typedef struct klist_{
	char *record;
	struct klist_ *next;
	struct klist_ *last;
} KeyList;

#define FLEN_FORMAT	(12)

typedef struct{
	void *contents;				// contents of table
	int coltype;				// type of columns
	long width;					// data width
	long repeat;				// amount of rows -> 'contents' size = width*repeat
	char colname[FLEN_KEYWORD];	// column name (arg ttype of fits_create_tbl)
	char format[FLEN_FORMAT];	// format codes (tform)
	char unit[FLEN_CARD];		// units (tunit)
}table_column;

typedef struct{
	int ncols;			// amount of columns
	long nrows;			// max amount of rows
	char tabname[80];	// table name
	table_column *columns;// array of structures 'table_column'
}FITStable;

typedef struct{
	size_t amount;		// amount of tables in file
	FITStable **tables;	// array of pointer to tables
} FITStables;

typedef struct{
	int width;			// width
	int height;			// height
	int dtype;			// data type
	//int lasthdu;		// last filled HDU number
	Item *data;			// picture data
	KeyList *keylist;	// list of options for each key
	FITStables *tables; // tables from FITS file
} IMAGE;


void list_free(KeyList **list);
KeyList *list_add_record(KeyList **list, char *rec);
KeyList *list_find_key(KeyList *list, char *key);
void list_remove_key(KeyList **list, char *key);
KeyList *list_modify_key(KeyList *list, char *key, char *newval);
void list_remove_records(KeyList **list, char *sample);
KeyList *list_copy(KeyList *list);
KeyList *list_get_end(KeyList *list);
void list_print(KeyList *list);

void tablefree(FITStable **tbl);
FITStable *table_new(IMAGE *img, char *tabname);
FITStable *table_read(IMAGE *img, fitsfile *fp);
FITStable *table_addcolumn(FITStable *tbl, table_column *column);
void table_write(IMAGE *img, fitsfile *fp);
void table_print(FITStable *tbl);
void table_print_all(IMAGE *img);

void imfree(IMAGE **ima);
IMAGE *readFITS(char *filename, IMAGE **fits);
bool writeFITS(char *filename, IMAGE *fits);
IMAGE *newFITS(size_t h, size_t w, int dtype);
IMAGE *similarFITS(IMAGE *in, int dtype);
IMAGE *copyFITS(IMAGE *in);
IMAGE *buildFITSfromdat(size_t h, size_t w, int dtype, uint8_t *indata);

extern struct stat filestat;
char* make_filename(char *buff, size_t buflen, char *prefix, char *suffix);
bool file_absent(char *name);


#endif // __FITS_H__
