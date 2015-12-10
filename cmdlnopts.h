/*
 * cmdlnopts.h - comand line options for parceargs
 *
 * Copyright 2013 Edward V. Emelianoff <eddy@sao.ru>
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
#ifndef __CMDLNOPTS_H__
#define __CMDLNOPTS_H__

#include "parceargs.h"

/*
 * here are some typedef's for global data
 */
typedef struct{
	char *infile;					// input file name
	char *outfile;					// output file name
	int rest_pars_num;				// amount of "free" parameters
	char **rest_pars;				// array of those parameters
	char **conv;					// conversion pipeline parameters
} glob_pars;


// default & global parameters
extern glob_pars const Gdefault;
extern int rewrite_ifexists, verbose_level, show_stat;
glob_pars G;

glob_pars *parce_args(int argc, char **argv);

#endif // __CMDLNOPTS_H__
