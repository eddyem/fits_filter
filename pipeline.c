/*
 * pipeline.c - image conversion pipeline
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
#include "pipeline.h"
#include "cmdlnopts.h"
#include "convfilter.h"
#include "linfilter.h"
#include "usefull_macros.h"
#include "convfilter.h"
#include "linfilter.h"
#include "median.h"

static Filter **farray = NULL; // array of pipeline conversion types
static size_t farray_size = 0;

// intermediate structure for input data conversion into filter parameters
typedef struct{
	char *ftype;
	char *scale;
	int help;
	int xsz;
	int ysz;
	double xhw;
	double yhw;
	imfuncptr imfunc;
} pipepars;

typedef struct{
	FType FilterType;  // type
	char *parname;     // its name
	char *descr;       // short description
	char **arguments;  // list of arguments available
	imfuncptr imfunc;  // function called for this type of conversion
} ftypename;

/// "sx,sy\tпараметр сигма по осям x и y\nw,h\tширина и высота ненулевого окна фильтра"
char* lgargs = N_("sx,sy\tsigma by axes x & y\nw,h\tnon-zero window width & height");
/// "аргументы отсутствуют"
char* noneargs = N_("arguments are absent");
/// "r\tрадиус фильтра (беззнаковое целое или 0 для \"креста\" 3x3)"
char* medargs = N_("r\tradius of filter (uint, 0 for cross 3x3)");
/// "nsteps\tколичество градаций яркости\nscale\tфункция преобразования (uniform, log, exp, sqrt, pow)"
char* stepargs = N_("nsteps\tamount of steps\nscale\tscale type (uniform, log, exp, sqrt, pow)");

ftypename filter_names[] = {
	/// "медианный фильтр"
	{MEDIAN,    "median",    N_("median filter"), &medargs, get_median},
	/// "простой адаптивный медианный фильтр"
	{ADPT_MEDIAN,"adpmed",   N_("simple adaptive median filter"), &medargs, get_adaptive_median},
	/// "лапласиан гауссианы"
	{LAPGAUSS,  "lapgauss",  N_("laplasian of gaussian"), &lgargs, DiffFilter},
	/// "гауссов фильтр"
	{GAUSS,     "gauss",     N_("gaussian"), &lgargs, DiffFilter},
	/// "горизонтальный фильтр Собеля"
	{SOBELH,    "sobelh",    N_("horizontal Sobel"), &noneargs, DiffFilter},
	/// "вертикальный фильтр Собеля"
	{SOBELV,    "sobelv",    N_("vertical Sobel"), &noneargs, DiffFilter},
	/// "простой градиент (операторами Собеля)"
	{SIMPLEGRAD,"simplegrad",N_("simple gradient (by Sobel)"), &noneargs, GradFilterSimple},
	/// "горизонтальный фильтр Прюитта (простейшая производная)"
	{PREWITTH,  "prewitth",  N_("Prewitt horizontal - simple derivative"), &noneargs, DiffFilter},
	/// "вертикальный фильтр Прюитта"
	{PREWITTV,  "prewittv",  N_("Prewitt vertical"), &noneargs, DiffFilter},
	/// "горизонтальный фильтр Щарра (модифицированный Собель)"
	{SCHARRH,   "scharrh",   N_("Scharr (modified Sobel) horizontal"), &noneargs, DiffFilter},
	/// "вертикальный фильтр Щарра"
	{SCHARRV,   "scharrv",   N_("Scharr vertical"), &noneargs, DiffFilter},
	/// "\"постеризация\""
	{STEP,      "step",      N_("posterisation"), &stepargs, StepFilter},
	{FILTER_NONE, NULL, NULL, NULL, NULL}
};

void show_pipeline_pars(){
	int i = 0;
	/// "Параметры конвейера:\n"
	red(_("Pipeline parameters:\n"));
	while(filter_names[i].parname){
		printf("\t%-12s %s\n", filter_names[i].parname, filter_names[i].descr);
		++i;
	}
	signals(9);
}

void showparhelp(int idx){
	/// "Преобразование %s: %s\n"
	red(_("Conversion %s <%s> parameters:\n"), filter_names[idx].parname, _(filter_names[idx].descr));
	printf("%s\n", _(*filter_names[idx].arguments));
	signals(9);
}

Filter *parse_filter(char *pars){
	Filter *fltr;
	int idx = -1;
	pipepars popts;
	mysuboption pipeopts[] = {
		// help for given conversion
		{"help", NO_ARGS,  arg_none,   &popts.help},
		// type of conversion
		{"type", NEED_ARG, arg_string, &popts.ftype},
		// radius of median filter
		{"r",    NEED_ARG, arg_int,    &popts.xsz},
		// sigmax, sigmay for gauss/lapgauss
		{"sx",   NEED_ARG, arg_double, &popts.xhw},
		{"sy",   NEED_ARG, arg_double, &popts.yhw},
		// & window size for them
		{"w",    NEED_ARG, arg_int,    &popts.xsz},
		{"h",    NEED_ARG, arg_int,    &popts.ysz},
		// posterisation
		{"nsteps",NEED_ARG,arg_int,    &popts.xsz},
		{"scale",NEED_ARG, arg_string, &popts.scale},
		end_suboption
	};
	memset(&popts, 0, sizeof(pipepars));
	if(!get_suboption(pars, pipeopts)){
		return NULL;
	}else{
		if(!popts.ftype){ // no type given
			/// "Необходимо по крайней мере задать параметр 'type', например: '-c type=help'"
			ERRX(_("You should at least give parameter 'type', for example: '-c type=help'"));
		}
		DBG("type: %s", popts.ftype);
		if(strcmp(popts.ftype, "help") == 0){ // info about pipeline names
			show_pipeline_pars();
		}
		int i = 0;
		while(filter_names[i].parname){
			if(!strcmp(filter_names[i].parname, popts.ftype)){
				idx = i;
				break;
			}
			++i;
		}
		if(idx == -1){ // wrong type
			/// "Неправильный параметр 'type' конвейера"
			WARNX(_("Wrong pipeline 'type' parameter: %s"), popts.ftype);
			show_pipeline_pars();
		}
	}
	if(popts.help){ // help about current filter options
		showparhelp(idx);
	}
	fltr = MALLOC(Filter, 1);
	fltr->FilterType = filter_names[idx].FilterType;
	popts.imfunc = filter_names[idx].imfunc;
	DBG("idx: %d, ftype: %d", idx, fltr->FilterType);
	// check parameters & fill Filer fields
	if(popts.imfunc ==  get_median || popts.imfunc == get_adaptive_median){
		fltr->w = popts.xsz;
	}else if(popts.imfunc ==  DiffFilter &&
			(fltr->FilterType == LAPGAUSS || fltr->FilterType == GAUSS)){
		if(popts.xsz < 5){
			// "Ширина фильтра изменена на 5"
			WARNX(_("Filter window width changed to 5"));
			popts.xsz = 5;
		}
		if(popts.ysz < 5){
			// "Высота фильтра изменена на 5"
			WARNX(_("Filter window height changed to 5"));
			popts.ysz = 5;
		}
		if(popts.xhw < 1. || popts.yhw < 1.){
			/// "Полуширина фильтра должна быть не меньше 1."
			ERRX(_("Filter FWHM should be not less than 1."));
		}
		fltr->w = popts.xsz; fltr->h = popts.ysz;
		fltr->sx = popts.xhw; fltr->sy = popts.yhw;
	}else if(popts.imfunc ==  StepFilter){ // check levels & type
		if(popts.xsz < 2 || popts.xsz > 255){
			/// "Количество уровней градаций яркости должно быть от 2 до 255"
			ERRX(_("Brightness levels amount shoul be from 2 to 255"));
		}
		fltr->w = popts.xsz;
		int i = 0, idx = -1;
		DBG("name: %s", popts.scale);
		if(!popts.scale){
			/// "Не указан параметр scale"
			ERRX(_("You should set 'scale' parameter"));
		}
		while(scales[i].name){
			if(!strcmp(popts.scale, scales[i].name)){
				idx = i; break;
			}
			++i;
		}
		if(idx == -1){
			/// "Параметры фильтра постеризации должны быть такими:\n%s"
			ERRX(_("Posterisation filter parameters should be:\n%s"), stepargs);
		}
		DBG("idx: %d", idx);
		fltr->h = scales[idx].type;
	}
	fltr->imfunc = popts.imfunc;
	DBG("Got filter #%d: w=%d, h=%d, sx=%g, sy=%g\n", fltr->FilterType,
		fltr->w, fltr->h, fltr->sx, fltr->sy);
	return fltr;
}

/*
 * makes array of pipeline parameters
 * if any found, return TRUE
 * else return FALSE
 */
bool get_pipeline_params(){
	int i, N;
	char **p = G.conv;
	Filter *f;
	for(N = 0; *p; ++N, ++p);
	farray_size = N;
	if(N == 0) return FALSE;
	farray = MALLOC(Filter*, farray_size);
	p = G.conv;
	for(i = 0; i < N; ++i, ++p){
		DBG("filter: %s", *p);
		if(!(f = parse_filter(*p))){
			/// "Неправильно заданы параметры конвейера"
			ERRX(_("Wrong pipeline parameters!"));
		}
		farray[i] = f;
	}
	return TRUE;
}

IMAGE *process_pipeline(IMAGE *image){
	if(!image){
		/// "Не задано входное изображение"
		ERRX(_("No input image given"));
	}
	if(!farray || !farray_size){
		/// "Не заданы параметры конвейера"
		WARNX(_("No pipeline parameters given"));
	}
	size_t i;
	Filter **far = farray;
	IMAGE *in = copyFITS(image); // copy original image to leave it unchanged
	in->keylist = list_copy(image->keylist);
	IMAGE *processed = NULL;
	for(i = 0; i < farray_size; ++i, ++far){
		Filter *f = *far;
		DBG("Got filter #%d: w=%d, h=%d, sx=%g, sy=%g\n", f->FilterType,
			f->w, f->h, f->sx, f->sy);
		printf("try filter %zd\n", i);
		Itmarray oarg = {NULL, 0};
		processed = f->imfunc(in, f, &oarg);
		/// "Ошибка в обработке конвейера"
		if(!processed) ERRX(_("Error on pipeline processing!"));
		// TODO: what should I do with oarg???
		if(oarg.size){
			size_t i, l = oarg.size;
			green("got oarg: \n");
			for(i = 0; i < l; ++i){
				printf("%5zd: %g\n", i, oarg.data[i]);
			}
			FREE(oarg.data);
		}
		processed->keylist = in->keylist;
		// TODO: what's about writting changes into HISTORY?
		in->keylist = NULL; // prevent deleting global keylist
		imfree(&in);
		in = processed;
	}
	return processed;
}

