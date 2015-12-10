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
#include "convfilter.h"
#include "linfilter.h"
#include "cmdlnopts.h"
#include "pipeline.h"

#ifndef BUFF_SIZ
#define BUFF_SIZ 4096
#endif

void signals(int signo){
	exit(signo);
}

int main(int argc, char **argv){
	IMAGE *fits = NULL, *newfit = NULL;
	bool pipe_needed = FALSE;
	char buff[BUFF_SIZ];
	//size_t i, s;
	initial_setup();
	parce_args(argc, argv);
	if(!G.infile){
		/// "Не задано имя входного файла"
		ERRX(_("Missed input file name!"));
	}
	if(!G.outfile){ // user didn't write out file name - check for prefix
		if(!G.rest_pars_num){
			/// "Задайте имя выходного файла (-o) или его префикс (без ключа)"
			ERRX(_("Set output file name (-o) or its prefix (without key)"));
		}
		if(!(G.outfile = make_filename(buff, BUFF_SIZ, G.rest_pars[0], "fits"))){
			/// "Все 9999 файлов вида %sXXXX.fits заняты!"
			ERRX(_("All 9999 files like %sXXXX.fits exists!"), G.rest_pars[0]);
		}
	}else{ // check whether file don't exists or there's a key '--rewrite'
		if(!file_absent(G.outfile)){
			if(rewrite_ifexists){
				/// "Не могу удалить файл %s"
				if(unlink(G.outfile)) ERR(_("Can't remove file %s"), G.outfile);
			}else{
				/// "Выходной файл существует"
				ERRX(_("The output file exists"));
			}
		}
	}
	// pre-check pipeline parameters
	if(G.conv){
		pipe_needed = get_pipeline_params();
	}
	if(!readFITS(G.infile, &fits)){
		// "Невозможно прочесть входной файл!"
		ERR(_("Can't read input file!"));
	}
	DBG("ima: %dx%d", fits->width, fits->height);

	if(show_stat){
		Item min, max, mean, std, med;
		get_statictics(fits, &min, &max, &mean, &std, &med);
		// "Статистика по входному изображению:\n"
		green(_("Input image statistics:\n"));
		printf("min = %g, max = %g, mean = %g, std = %g, median = %g\n",
				min, max, mean, std, med);
	}

	if(pipe_needed)
		newfit = process_pipeline(fits);

	if(show_stat && newfit){
		Item min, max, mean, std, med;
		get_statictics(newfit, &min, &max, &mean, &std, &med);
		// "Статистика по изображению после конвейера:\n"
		green(_("Image statistics after pipeline:\n"));
		printf("min = %g, max = %g, mean = %g, std = %g, median = %g\n",
				min, max, mean, std, med);
	}
	writeFITS(G.outfile, newfit);
/*
	IMAGE *newfit = get_median(fits, 2);
	Filter f = {LAPGAUSS, 200, 200, atoi(argv[3]), atoi(argv[3])};
	IMAGE *newf = DiffFilter(newfit, &f);
	Filter f1 = {STEP, 3, LOG, 0, 0};
	IMAGE *newfits = StepFilter(newf, &f1, NULL);

	//IMAGE *newfit = get_median(fits, atoi(argv[3]));
	//IMAGE *newfits = GradFilterSimple(newfit);

	//Filter f = {STEP, atoi(argv[3]), POW, 0, 0};
	//Item *levels;
	//IMAGE *newfits = StepFilter(fits, &f, &levels);
	//FREE(levels);
	newfits->keylist = fits->keylist;
	newfits->keynum = fits->keynum;
	writeFITS(G.outfile, newfits);
	imfree(&fits);
	FREE(newfits->data); FREE(newfits);
	*/
	return 0;
}
