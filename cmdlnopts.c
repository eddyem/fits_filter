/*
 * cmdlnopts.c - the only function that parse cmdln args and returns glob parameters
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
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include "cmdlnopts.h"
#include "usefull_macros.h"

#define RAD 57.2957795130823
#define D2R(x) ((x) / RAD)
#define R2D(x) ((x) * RAD)

/*
 * here are global parameters initialisation
 */
int help;
glob_pars G;

int rewrite_ifexists = 0  // rewrite existing files == 0 or 1
    ,verbose_level = 0    // each -v increments this value, e.g. -vvv sets it to 3
    ,show_stat = 0        // show statistic parameters of input and output image
    ,inplace = 0          // save results into same file
;

char **keys2delete = NULL  // keylist for deletion
    ,**recs2delete = NULL  // records to delete by substring
    ,**recs2add    = NULL  // records to add
;

//            DEFAULTS
// default global parameters
glob_pars const Gdefault = {
     .infile = NULL
    ,.outfile = NULL
    ,.rest_pars_num = 0
    ,.rest_pars = NULL
    ,.conv = NULL
    ,.oper = MATH_NONE
    ,.low_bound = DBL_MAX
    ,.up_bound = DBL_MAX
    ,.binarize = DBL_MAX
    ,.conncomp4 = DBL_MAX
    ,.conncomp8 = DBL_MAX
    ,.flip = NULL
    ,.deltabs = 0
    ,.listabs = 0
};

/// "установить параметры конвейера, аргументы: type:[help]:...\n\t\ttype - тип преобразования (help для справки)\n\t\thelp - список доступных для данного 'type' опций"
const char FilPar[] = N_("set pipeline parameters, arg: type=type:[help]:...\n" \
        "\t\ttype - transformation type (help for list)\n" \
        "\t\thelp - list of available parameters for given 'type'");

/*
 * Define command line options by filling structure:
 *  name    has_arg flag    val     type        argptr          help
*/
myoption cmdlnopts[] = {
    /// "отобразить эту справку"
    {"help",    NO_ARGS,    NULL,   'h',    arg_int,    APTR(&help),        N_("show this help")},
    /// "входной файл"
    {"infile",  NEED_ARG,   NULL,   'i',    arg_string, APTR(&G.infile),    N_("input file")},
    /// "записать изменения в тот же файл"
    {"inplace", NO_ARGS,    &inplace,1,     arg_none,   NULL,               N_("save results into same file")},
    /// "выходной файл"
    {"outfile", NEED_ARG,   NULL,   'o',    arg_string, APTR(&G.outfile),   N_("output file")},
    /// "установить параметры конвейера, arg: type=type:[help]:...\n\t\ttype - тип преобразования (help для справки)\n\t\thelp - справка по параметрам типа"
    {"pipeline",MULT_PAR,   NULL,   'p',    arg_string, APTR(&G.conv),      FilPar},
    /// "перезаписать выходной файл, если он существует (только с опцией -i)"
    {"rewrite", NO_ARGS,    &rewrite_ifexists,1,arg_none,NULL,              N_("rewrite output file if exists (works only with option -i)")},
    /// "уровень подробностей вывода (каждый -v увеличивает его)"
    {"verbose", NO_ARGS,    NULL,   'v',    arg_none,   APTR(&verbose_level),N_("verbose level (each -v increase it)")},
    /// "отобразить статистические параметры входного и выходного изображения"
    {"stat",    NO_ARGS,    NULL,   's',    arg_int,    APTR(&show_stat),   N_("show statistic parameters of input and output image")},
    /// "удалить указанный ключ"
    {"del-key", MULT_PAR,   NULL,   'd',    arg_string, APTR(&keys2delete), N_("delete given key")},
    /// удалить все записи с указанной подстрокой
    {"del-rec", MULT_PAR,   NULL,   'D',    arg_string, APTR(&recs2delete), N_("delete all records with given substring")},
    /// добавить запись в FITS-шапку
    {"add-rec", MULT_PAR,   NULL,   'a',    arg_string, APTR(&recs2add),    N_("add record to FITS-header")},
    /// вычислить сумму всех перечисленных изображений
    {"sum",     NO_ARGS,    &G.oper,MATH_SUM,arg_none,  NULL,               N_("calculate sum of specified FITS-files")},
    /// вычислить разность первого и остальных перечисленных изображений
    {"diff",    NO_ARGS,    &G.oper,MATH_DIFF,arg_none, NULL,               N_("calculate difference of first and rest specified FITS-files")},
    /// вычислить среднее арифметическое всех перечисленных изображений
    {"mean",    NO_ARGS,    &G.oper,MATH_MEAN,arg_none, NULL,               N_("calculate mean of specified FITS-files")},
    /// вычислить медиану всех перечисленных изображений
    {"median",  NO_ARGS,    &G.oper,MATH_MEDIAN,arg_none,NULL,              N_("calculate median of specified FITS-files")},
    /// нижняя граница, все значения интенсивности ниже нее будут установлены в нее
    {"bottom",  NEED_ARG,   NULL,   'b',    arg_double, APTR(&G.low_bound), N_("the lowest value, all data less than 'bottom' vould be set to it")},
    /// верхняя граница, все значения интенсивности выше нее будут установлены в нее
    {"top",     NEED_ARG,   NULL,   't',    arg_double, APTR(&G.up_bound),  N_("the lowest value, all data more than 'top' vould be set to it")},
    /// бинаризация изображения по порогу (в %% от динамического диапазона)
    {"binarize",NEED_ARG,   NULL,   'B',    arg_double, APTR(&G.binarize),  N_("binarize image by threshold (in %%from dynamic range)")},
    /// маркировать 4-связные области по заданному порогу
    {"conn4",   NEED_ARG,   NULL,   '4',    arg_double, APTR(&G.conncomp4), N_("label 4-connected components with given threshold")},
    /// маркировать 8-связные области по заданному порогу
    {"conn8",   NEED_ARG,   NULL,   '8',    arg_double, APTR(&G.conncomp8), N_("label 8-connected components with given threshold")},
    /// зеркалировать изображение
    {"flip",    NEED_ARG,   NULL,   'f',        arg_string, APTR(&G.flip),      N_("flip image (arg = X, Y or XY)")},
    {"no-tabs", NO_ARGS,    &G.deltabs,1,   arg_none,   NULL,               N_("don't save any tables in output file")},
    {"list-tabs",NO_ARGS,   &G.listabs,1,   arg_none,   NULL,               N_("List all tables in input file")},
    end_option
};


// suboptions structure for get_mirpars
// in array last MUST BE {0,0,0}
typedef struct{
    double *val;    // pointer to result
    char *par;      // parameter name (CASE-INSENSITIVE!)
    bool isdegr;    // == TRUE if parameter is an angle in format "[+-][DDd][MMm][SS.S]"
} suboptions;

/**
 * Safely convert data from string to double
 *
 * @param num (o) - double number read from string
 * @param str (i) - input string
 * @return TRUE if success
 */
bool myatod(double *num, const char *str){
    double res;
    char *endptr;
    assert(str);
    res = strtod(str, &endptr);
    if(endptr == str || *str == '\0' || *endptr != '\0'){
        /// "Неправильный формат числа double!"
        WARNX(_("Wrong double number format!"));
        return FALSE;
    }
    *num = res;
    return TRUE;
}

/**
 * Convert string "[+-][DDd][MMm][SS.S]" into radians
 *
 * @param ang (o) - angle in radians or exit with help message
 * @param str (i) - string with angle
 * @return TRUE if OK
 */
bool get_radians(double *ret, char *str){
    double val = 0., ftmp, sign = 1.;
    char *ptr;
    assert(str);
    switch(*str){ // check sign
        case '-':
            sign = -1.;
        case '+':
            str++;
    }
    if((ptr = strchr(str, 'd'))){ // found DDD.DDd
        *ptr = 0; if(!myatod(&ftmp, str)) return FALSE;
        ftmp = fabs(ftmp);
        if(ftmp > 360.){
            /// "Угол должен быть меньше 360 градусов"
            WARNX(_("Degrees should be less than 360"));
            return FALSE;
        }
        val += ftmp;
        str = ptr + 1;
    }
    if((ptr = strchr(str, 'm'))){ // found DDD.DDm
        *ptr = 0; if(!myatod(&ftmp, str)) return FALSE;
        ftmp = fabs(ftmp);
        val += ftmp / 60.;
        str = ptr + 1;
    }
    if(strlen(str)){ // there is something more
        if(!myatod(&ftmp, str)) return FALSE;
        ftmp = fabs(ftmp);
        val += ftmp / 3600.;
    }
    *ret = D2R(val * sign); // convert degrees to radians
    return TRUE;
}

/**
 * Parse string of suboptions (--option=name1=var1:name2=var2... or -O name1=var1,name2=var2...)
 * Suboptions could be divided by colon or comma
 *
 * !!NAMES OF SUBOPTIONS ARE CASE-UNSENSITIVE!!!
 *
 * @param arg (i)    - string with description
 * @param V (io)     - pointer to suboptions array (result will be stored in sopts->val)
 * @return TRUE if success
 *
bool get_mirpars(void *arg, suboptions *V){
    char *tok, *val, *par;
    int i;
    tok = strtok(arg, ":,");
    do{
        if((val = strchr(tok, '=')) == NULL){ // wrong format
            printf("Wrong format: no value for keyword");
            return FALSE;
        }
        *val++ = '\0';
        par = tok;
        for(i = 0; V[i].val; i++){
            if(strcasecmp(par, V[i].par) == 0){ // found parameter
                if(V[i].isdegr){ // DMS
                    if(!get_radians(V[i].val, val)) // wrong angle
                        return FALSE;
                }else{ // simple double
                    if(!myatod(V[i].val, val)) // wrong number
                        return FALSE;
                }
                break;
            }
        }
        if(!V[i].val){ // nothing found - wrong format
            printf("Bad keyword!");
            return FALSE;
        }
    }while((tok = strtok(NULL, ":,")));
    return TRUE;
}*/

/**
 * functions of subargs parcing can looks as this
 */
/**
 * Parse string of mirror parameters (--mir-diam=...)
 *
 * @param arg (i) - string with description
 * @return TRUE if success
 *
bool get_mir_par(void *arg){
    suboptions V[] = { // array of mirror parameters and string keys for cmdln pars
        {&M.D,      "diam",     FALSE},
        {&M.F,      "foc",      FALSE},
        {&M.Zincl,  "zincl",    TRUE},
        {&M.Aincl,  "aincl",    TRUE},
        {0,0,0}
    };
    return get_mirpars(arg, V);
}*/

/**
 * Parse command line options and return dynamically allocated structure
 *      to global parameters
 * @param argc - copy of argc from main
 * @param argv - copy of argv from main
 * @return allocated structure with global parameters
 */
glob_pars *parse_args(int argc, char **argv){
    int i;
    void *ptr;
    ptr = memcpy(&G, &Gdefault, sizeof(G)); assert(ptr);
    /// "Использование: %s [аргументы] [префикс выходного файла]\n\n\tГде [аргументы]:\n"
    change_helpstring(_("Usage: %s [args] [outfile prefix] [file list for group operations]\n\n\tWhere args are:\n"));
    // parse arguments
    parseargs(&argc, &argv, cmdlnopts);
    if(help) showhelp(-1, cmdlnopts);
    if(argc > 0){
        G.rest_pars_num = argc;
        G.rest_pars = calloc(argc, sizeof(char*));
        for (i = 0; i < argc; i++)
            G.rest_pars[i] = strdup(argv[i]);
    }
    return &G;
}

