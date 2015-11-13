/*
 * main.c
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
#include "usefull_macros.h"
#include "fits.h"
#include "median.h"

void signals(int signo){
	exit(signo);
}

int main(int argc, char **argv){
	IMAGE *fits;
	//size_t i, s;
	initial_setup();
	if(argc != 4) ERRX("Usage: %s infile outfile num", argv[0]);
	if(!readFITS(argv[1], &fits))
		ERR(_("Can't read input file!"));
	DBG("ima: %dx%d", fits->width, fits->height);
	//s = fits->width * fits->height;
	//double *img = fits->data;
	//for(i = 0; i < s; ++i) *img++ /= 2.;
	unlink(argv[2]);
	IMAGE *newfits = get_adaptive_median(fits, atoi(argv[3]));
	newfits->keylist = fits->keylist;
	newfits->keynum = fits->keynum;
	writeFITS(argv[2], newfits);
	imfree(&fits);
	FREE(newfits->data); FREE(newfits);
	return 0;
}
