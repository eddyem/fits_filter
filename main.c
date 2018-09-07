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
#include "group_operations.h"
#include "binmorph.h"

#ifndef BUFF_SIZ
#define BUFF_SIZ 4096
#endif

void signals(int signo){
    exit(signo);
}

#define swap_el(a, b)  {register Item t = a; a = b; b = t;}
/**
 * flip an image by X-axis (top <-> bottom)
 */
void flip_X(IMAGE *f){
    if(!f){
        WARNX("Empty arg!");
        return;
    }
    int w = f->width, h = f->height, h2 = h/2;
    OMP_FOR()
    for(int _col = 0; _col < w; ++_col){
        Item *pixa = &f->data[_col], *pixb = pixa + w*(h-1); // first & last pixels in column
        for(int _row = 0; _row < h2; ++_row, pixa += w, pixb -= w){
            swap_el(*pixa, *pixb);
        }
    }
}
/**
 * flip an image by Y-axis (left <-> right)
 */
void flip_Y(IMAGE *f){
    if(!f){
        WARNX("Empty arg!");
        return;
    }
    int w = f->width, h = f->height, w2 = w/2;
    OMP_FOR()
    for(int _row = 0; _row < h; ++_row){
        Item *pixa = &f->data[_row*w], *pixb = pixa + w - 1; // first & last pixels in row
        for(int _col= 0; _col < w2; ++_col, ++pixa, --pixb){
            swap_el(*pixa, *pixb);
        }
    }
}



int main(int argc, char **argv){
    IMAGE *fits = NULL, *newfit = NULL;
    bool pipe_need = FALSE;
    char buff[BUFF_SIZ];
    //size_t i, s;
    initial_setup();
    parse_args(argc, argv);
    // pre-check pipeline parameters
    if(G.conv){
        pipe_need = get_pipeline_params();
    }
    if(!G.infile && G.oper == MATH_NONE){
        /// "Не задано имя входного файла"
        ERRX(_("Missed input file name!"));
    }else{
        if(G.oper != MATH_NONE){
            if(G.infile){
                /// Параметр '-i' не используется при групповых операциях!
                ERRX(_("Parameter '-i' isn't used in group operations!"));
            }
            if(inplace){
                /// Не могу выполнить групповую операцию с сохранением 'в тот же файл'!
                ERRX(_("Can't make group operations 'in place'!"));
            }
        }
    }
    if((G.infile) && !readFITS(G.infile, &fits)){
        // "Невозможно прочесть входной файл!"
        ERR(_("Can't read input file!"));
    }
    if(!G.outfile){ // user didn't write out file name - check for prefix
        if(inplace){ // G.outfile is the same as G.infile
            G.outfile = G.infile;
        }else if(!show_stat && !G.listabs){
            if(!G.rest_pars_num){
                /// "Задайте имя выходного файла (-o) или его префикс (без ключа)"
                ERRX(_("Set output file name (-o) or its prefix (without key)"));
            }
            if(!(G.outfile = make_filename(buff, BUFF_SIZ, G.rest_pars[0], "fits"))){
                /// "Все 9999 файлов вида %sXXXX.fits заняты!"
                ERRX(_("All 9999 files like %sXXXX.fits exists!"), G.rest_pars[0]);
            }
            // move starting to next free parameter
            if(--G.rest_pars_num)
                G.rest_pars = &G.rest_pars[1];
            else
                G.rest_pars = NULL;
        }
    }else{ // check whether file don't exists or there's a key '--rewrite'
        if(!file_absent(G.outfile)){
            if(rewrite_ifexists || inplace){
                /// "Не могу удалить файл %s"
                if(unlink(G.outfile)) ERR(_("Can't remove file %s"), G.outfile);
            }else{
                /// "Выходной файл существует"
                ERRX(_("The output file exists"));
            }
        }
    }
    // Process pipeline only if there's
    if(G.infile){
        if(show_stat){
            Item min, max, mean, std, med;
            get_statictics(fits, &min, &max, &mean, &std, &med);
            // "Статистика по входному изображению:\n"
            green(_("Input image statistics:\n"));
            printf("min = %g, max = %g, mean = %g, std = %g, median = %g\n",
                    min, max, mean, std, med);
        }
        if(G.listabs) table_print_all(fits);
    }else{ // G.oper != MATH_NONE or some other (in future?)
        if(G.oper != MATH_NONE){ // process all files to make group operation
            if(G.rest_pars_num < 2){
                /// "Групповые операции требуют задания как минимум двух FITS-файлов
                ERRX(_("Group operations need at least two FITS-files"));
            }
            // now create variable "fits" with result of grouping operation
            fits = make_group_operation(G.rest_pars_num, G.rest_pars, G.oper);
            if(!fits){
                /// Ошибка при групповой обработке
                ERRX(_("Error in group operation"));
            }
            if(show_stat){
                Item min, max, mean, std, med;
                get_statictics(fits, &min, &max, &mean, &std, &med);
                // "Статистика по изображению после групповых операций:\n"
                green(_("Image statistics after group operations:\n"));
                printf("min = %g, max = %g, mean = %g, std = %g, median = %g\n",
                        min, max, mean, std, med);
            }
        }
    }
    // process pipeline both in case of single input file ('-i')
    if(pipe_need)
        newfit = process_pipeline(fits);

    if(show_stat && newfit){
        Item min, max, mean, std, med;
        get_statictics(newfit, &min, &max, &mean, &std, &med);
        // "Статистика по изображению после конвейера:\n"
        green(_("Image statistics after pipeline:\n"));
        printf("min = %g, max = %g, mean = %g, std = %g, median = %g\n",
                min, max, mean, std, med);
    }
    if(!newfit) newfit = fits;
    // process cuts & so on
    cut_bounds(newfit, G.low_bound, G.up_bound);
    if(!newfit) signals(-1);
    // here are operations that can't be done together
    if(G.binarize < DBL_MAX - 1.){
        IMAGE *tmp = get_binary(newfit, G.binarize);
        imfree(&newfit);
        newfit = tmp;
    }else if(G.conncomp4 < DBL_MAX - 1.){
        size_t n;
        IMAGE *tmp = cclabel4(newfit, G.conncomp4, &n);
        imfree(&newfit);
        newfit = tmp;
        green("Found %d 4-connected regions\n", n);
    }else if(G.conncomp8 < DBL_MAX - 1.){
        size_t n;
        IMAGE *tmp = cclabel8(newfit, G.conncomp8, &n);
        imfree(&newfit);
        newfit = tmp;
        green("Found %d 8-connected regions\n", n);
    }
    if(!newfit) return 0;
    /**************************************************************************************************************
     *
     * place here any other operations with [processed through pipeline] input file or result of group operations
     *    (newfit)
     *
     **************************************************************************************************************/
    // change keys in output file FITS-header
    if(keys2delete){
        DBG("delete keys");
        do{
            list_remove_key(&newfit->keylist, *keys2delete);
        }while(*(++keys2delete));
    }
    if(recs2delete){
        DBG("delete records");
        do{
            list_remove_records(&newfit->keylist, *recs2delete);
        }while(*(++recs2delete));
    }
    if(recs2add){
        DBG("add records");
        do{
            list_add_record(&newfit->keylist, *recs2add);
        }while(*(++recs2add));
    }

    // process flipping
    if(G.flip){
        char *flip = G.flip;
        if(strchr(flip, 'x') || strchr(flip, 'X')) flip_X(newfit);
        if(strchr(flip, 'y') || strchr(flip, 'Y')) flip_Y(newfit);
    }
    writeFITS(G.outfile, newfit);

    return 0;
}
