#include "nrm.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
/*
 * process configuration
 * Futing Fan
 * Kansas State University
 * Updates by Niema Moshiri (UC San Diego)
 * Last Modified: March 2018
 * Copyright (c) 2016, Futing Fan. All rights reserved.
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted
 */

/*
 * count item in some section
 *
 *input:  FILE* p_file     [ file pointer]
 *input:  char* target_section     [ name of target section]
 *return: int   [>=0: number of items in target section; <0: failure to locate section]
 */
int item_count( FILE* p_file, char* target_section);

/*
 * locate section without 0 item check
 *
 *input:  FILE* p_file     [ file pointer]
 *input:  char* target_section     [ name of target section]
 *return: int   [0: success; <0: failure to locate section]
 */
int locate_section_only( FILE* p_file, char* target_section);

/*
 *locate section with 0 itme check
 *
 *input:  FILE* p_file     [ file pointer]
 *input:  char* target_section     [ name of target section]
 *output: FILE* p_file     [ with pointer pointing to the target section first item]
 *return: int   [0: success; <0: failure to locate section or section has 0 item]
 */
int locate_section( FILE* p_file, char* target_section);

/*
 *count column number in the first row of target file
 *
 *input:  FILE* p_file     [ file pointer]
 *return: int   [number of columns in first row]
 */
int fcolumn_count( FILE* p_file);
/*
 *count column number in input string
 *
 *input:  char* string     [ file pointer]
 *return: size_t   [number of columns in first row]
 */
size_t column_count( char* string);

/*
 *analysis network metrics
 *
 *input:  FILE* fil_para   [ parameter file, with pointer pointing at the beginning of the first network data file name]
 *inout:  Graph* graph     [graph struct]
 *return: int   [0: success; <0: failure]
 */
int analysis_network(FILE* fil_para, Graph* graph);

/*
 *fgetline( skip blank line and trim space on both end)
 *
 *input:  FILE*  file       [ file pointer]
 *        size_t limit      [ buffer length]
 *output: char*  line       [ output buffer]
 *return: int    [0: reach end; >0: length of line read]
 */
int fgetline(FILE* file, char* line, size_t limit);

/*
 *skip_top_comment( skip comment lead by 'c' on top)
 *
 *input:  FILE* file       [ file pointer]
 *        char  c          [ comment mark]
 *return: int   [<0: error; >= 0: number of line skipped]
 */
int skip_top_comment(FILE* file, char c);

/*
 *fget_next_item( skip blank line and trim space on both end AND line start with # AND stop at line starts with '[')
 *
 *input:  FILE*  file       [ file pointer]
 *        size_t limit      [ buffer length]
 *output: char*  line       [ output buffer]
 *return: int    [0: stop; >0: length of item read]
 */
int fget_next_item(FILE* file, char* line, size_t limit);

/*
 *fcheck_config( check configure of the target section suites requirment
 *
 *input:  FILE*  file       [ file pointer]
 *        size_t line       [ check line number]
 *        size_t colun      [ check column number in each line]
 *return: int    [0: success; <0: failure]
 */
int fcheck_config(FILE* file, char* target_section, size_t row, size_t column);

/*
 *analysis network metrics
 *
 *input:  FILE* fil_para   [ parameter file, with pointer pointing at the beginning of the first network data file name]
 *        int   compartment_qtt [ number of compartment]
 *        int   node_qtt        [ number of node]
 *        int   BEGIN_NUM       [ beginning number of node]
 *        int   STS_BEGIN       [ beginning number of status]
 *output: int*  p_ini_sts_lst   [ status list, length is node_qtt]
 *        int*  ini_sts_cnt     [ status count list, length compartment_qtt]
 *return: int   [0: success; <0: failure]
 */
int initial_con( Graph* graph, Status* sts);

//sscanf with moving pointer
int _auto_sscanf(char** str, char* format, ...);
//get an integer value from section of fil
LONG getValInt( FILE* fil, char* section, int echo);
//get a size_t value from section of fil
size_t* getValSize_tLst( FILE* fil, char* section, size_t len, int echo);
//get a double value from section of fil
double getValDbl( FILE* fil, char* section, int echo);
//get a char array  value from section of fil
char *getValStr( FILE* fil, char* section, size_t limit, int echo);
//get a 3D matrix from section of fil, len X dim X dim
double*** getValMatrixLst( FILE* fil, char * section, size_t dim, size_t len, size_t skip, size_t line_limit, int echo);
//get a 2D matrix from section of fil, dim X dim
double** getValMatrix( FILE* fil, char * section, size_t dim, size_t skip, size_t line_limit, int echo);
