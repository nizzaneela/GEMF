#include "para.h"
#include "common.h"
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <limits.h>
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

int _auto_sscanf(char** str, char* format, ...){
    int i, ret;
    va_list args;
    char* _str;
    _str= *str;
    va_start( args, format);
    while( isspace(_str[0])){
        _str++;
    }
    ret= vsscanf( _str, format, args);
    va_end(args);
    i= ret;
    while(i> 0){
        while( !isspace(_str[0])){
            _str++;
            if( _str[0]== '\0'){
                *str= _str;
                return ret;
            }
        }
        while( isspace(_str[0])){
            if( _str[0]== '\0'){
                return ret;
            }
            _str++;
        }
        i --;
    }
    *str= _str;
    return ret;
}

/*
 * count item in some section
 *
 *input:  FILE* p_file     [ file pointer]
 *input:  char* target_section     [ name of target section]
 *return: int   [>=0: number of items in target section; <0: failure to locate section]
 */
int item_count( FILE* p_file, char* target_section){
    int ret= 0;
    LINE tmp_str;

    ret= locate_section_only( p_file, target_section);
    if( ret< 0) return ret;
    while( fgetline( p_file, tmp_str, MAX_LINE_LEN)){
        if( tmp_str[0]== '#') continue;
        if( tmp_str[0]== '[') break;
        ret++;
    }
    return ret;
}

/*
 * locate section without 0 item check
 *
 *input:  FILE* p_file     [ file pointer]
 *input:  char* target_section     [ name of target section]
 *return: int   [0: success; <0: failure to locate section]
 */
int locate_section_only( FILE* p_file, char* target_section){
    int find_flag= 0;
    LINE tmp_str;

    rewind(p_file);
    while( fgetline(p_file, tmp_str, MAX_LINE_LEN)){
        if( !strcmp(tmp_str, target_section)){
            find_flag= 1;
            break;
        }
    }
    if( find_flag!= 1){
        if( strcmp(target_section,"[RANDOM_SEED]") != 0){
            printf( "missing section %s\n", target_section);
        }
        return -1;
    }
    return 0;
}
/*
 *locate section with 0 itme check
 *
 *input:  FILE* p_file     [ file pointer]
 *input:  char* target_section     [ name of target section]
 *output: FILE* p_file     [ with pointer pointing to the target section first item]
 *return: int   [0: success; <0: failure to locate section or section has 0 item]
 */
int locate_section( FILE* p_file, char* target_section){
    if( item_count( p_file, target_section)<= 0){
        if( strcmp(target_section,"[RANDOM_SEED]") == 0){
            return -1; //allow for no specification of random seed (use time)
        }
        printf("missing config item in section[%s]\n", target_section);
        exit( - 1);
    }
    locate_section_only( p_file, target_section);
    return 0;
}
/*
 *count column number in the first row of target file
 *
 *input:  FILE* p_file     [ file pointer]
 *return: int   [number of columns in first row]
 */
int fcolumn_count( FILE* p_file){
    int s, c, i= 2;

    while( isspace(c= fgetc(p_file)));
    while( c!= '\n'){
        s= isspace(c);
        if( s^(isspace(c=fgetc(p_file)))){
            i++;
        }
    }
    return i/2;
}
/*
 *count column number in input string
 *
 *input:  char* string     [ file pointer]
 *return: size_t   [number of columns in first row]
 */
size_t column_count( char* string){
    int i= - 1;
    size_t c= 2;

    while( string[++i]!= '\0' && isspace(string[i]));
    if( string[i]== '\0') return 0;
    while( string[++i]!= '\0'){
        if( isspace( string[i-1])^isspace( string[i])){
            c++;
        }
    }
    return c/2;
}
/*
 *analysis network metrics
 *
 *input:  FILE* fil_para   [ parameter file, with pointer pointing at the beginning of the first network data file name]
 *inout:  Graph* graph     [graph struct]
 *return: int   [0: success; <0: failure]
 */
int analysis_network(FILE* fil_para, Graph* graph){
    int _weighted;
    size_t count, layer;
    LONG li, lj;
    int ret;
    NINT _end_num, _begin_num;
    double time0;
    double lw;
    LINE fil_nam;
    LINE tmp_str;
    FILE* fil_dat;

    _begin_num= INT_MAX;
    _end_num= 0;
    _weighted= -1;
    graph->E= (size_t*)malloc(sizeof(size_t)*graph->L);
    if( graph->E== NULL){
        printf("Memory allocation failure for edges number list, size[%zu]\n", sizeof(size_t)*graph->L);
        exit( - 1);
    }
    for( layer= 0; layer< graph->L; layer++){
        //fscanf( fil_para, "%s ", fil_nam);
        fget_next_item( fil_para, fil_nam, MAX_LINE_LEN);
        fil_dat= fopen( fil_nam, "r");
        if( fil_dat== NULL){
            printf("read layer[%zu] network file[%s] error\n", layer, fil_nam);
            exit( - 1);
        }
        fget_next_item( fil_dat, tmp_str, MAX_LINE_LEN);
        _weighted= (int)column_count( tmp_str);
        if( _weighted== 2){
            _weighted= 0;
        }
        else if( _weighted== 3){
            _weighted= 1;
        }
        else{
            printf("wrong column number[%d] in file[%s]\n", _weighted, fil_nam);
            return -1;
        }
        count= 0;
        skip_top_comment( fil_dat, '#');
        time0= gettimenow();
        while( 1){
            if( _weighted){
                ret= fscanf(fil_dat,"%lld %lld %lf", &li, &lj, &lw);
            }
            else{
                ret= fscanf(fil_dat, "%lld %lld", &li, &lj);
            }
            if( ret!= (2+ _weighted)) break;
            if( count&& !(count%5000000)){
                time_print( "[ ", gettimenow() - time0, " ]");
                kilobit_print(", [ ", (LONG)count, " ] edges detected\n");
            }

            check_int_range( li);
            check_int_range( lj);
            if( _begin_num> (NINT)li) _begin_num= (NINT)li;
            if( _begin_num> (NINT)lj) _begin_num= (NINT)lj;
            if( _end_num< (NINT)li) _end_num= (NINT)li;
            if( _end_num< (NINT)lj) _end_num= (NINT)lj;
            count++;
        }
        graph->E[layer]= count;
        printf("layer %zu edges", layer);
        kilobit_print("\t\t[ ", (LONG)count, " ]\n");
        fclose( fil_dat);
    }
    graph->_s= _begin_num;
    graph->_e= _end_num +1;
    graph->V= _end_num - _begin_num+ 1;
    graph->weighted= _weighted;
    return 0;
}
/*
 *initial status
 *
 *input:  FILE* fil_para   [ parameter file, with pointer pointing at the beginning of the first network data file name]
 *        int   node_qtt        [ number of node]
 *        int   STS_BEGIN       [ beginning number of status]
 *output: int*  p_ini_sts_lst   [ status list, length is node_qtt]
 *return: int   [0: success; <0: failure]
 */
int initial_con(FILE* fil_sts, Graph* graph, Status* sts){
    //MODE 1. full list mode
    /*
     *e.g. for 3 compartments, 10nodes
     *{
         0
         1
         1
         2
         1
         1
         1
         0
         1
         2
     *}
     *
     */
    //MODE 2. statistic mode
    /*
     *e.g. for 3 compartments, 10nodes
     *
     *{
         6 2 2
         or
         -1 2 2
     *}
     */
    //MODE 3. compensate list mode
    /*
     *e.g. for 3 compartments, 10nodes, node 1,2,9 is in compartment 1, node 4,8 is incompartment 2, all rest in the other compartmetn
     *{
         0
         1 1
         2 1
         9 1
         4 2
         8 2
     *}
     */
    size_t li, lj, dynamic, line_num, max_compartmet_value, ns;
    int j, k, val;
    NINT count, max_compartment, ni;
    LINE ch;
    sts->init_lst= (size_t*)malloc(sizeof(size_t)*graph->_e);
    sts->subsample_lst= (size_t*)malloc(sizeof(size_t)*graph->_e);
    sts->init_cnt= (NINT*)malloc(sizeof(NINT)*(sts->_s+sts->M));

    // line_num= 0;
    // while( fgetline( fil_sts, ch, MAX_LINE_LEN)) line_num++;
    // if( line_num== graph->V){
    //     //MODE 1.
    //     LOG(2, __FILE__, __LINE__, "Status file mode [1]\n");
    //     rewind( fil_sts);
    //     for (li= graph->_s; li< graph->_e; li++) {
    //         if( fscanf(fil_sts,"%zu", sts->init_lst+li)<= 0){
    //             printf("fatal error, missing initial status for node [%zu] and all above\n", li);
    //             return -1;
    //         }
    //         sts->init_cnt[*(sts->init_lst+li)]++;
    //     }
    // }
    // else if( line_num== 1){
    //     //MODE 2.
    //     LOG(2, __FILE__, __LINE__, "Status file mode [2]\n");
    //     srand((unsigned int)sts->random_seed);
    //     k= 0;
    //     count= 0;
    //     rewind( fil_sts);
    //     for( li= sts->_s; li< sts->M+ sts->_s; li++){
    //         if( fscanf( fil_sts, "%d", &val)> 0){
    //             if( val< 0){
    //                 k++;
    //                 dynamic= li;
    //             }
    //             else{
    //                 count+= (NINT)val;
    //             }
    //             sts->init_cnt[li]= (NINT)val;
    //         }
    //         else{
    //             for(; li< sts->M+ sts->_s; li++){
    //                 count+= (NINT)val;
    //                 sts->init_cnt[li]= (NINT)val;
    //             }
    //             break;
    //         }
    //     }
    //     if( k> 1){
    //         printf("Fatal error! wrong status config, only one or less compartment can be dynamic(negative)\n");
    //         exit( - 1);
    //     }
    //     else if( count> graph->V){
    //         printf("Fatal error! wrong status config, configure exceed node number\n");
    //         exit( - 1);
    //     }
    //     else if( k==0 && count< graph->V){
    //         printf("Fatal error! wrong status config, if total number is uncertain, spcify one compartment as -1\n");
    //         exit( - 1);
    //     }
    //     else if( k== 1){
    //         sts->init_cnt[dynamic]= graph->V - count;
    //     }
    //     //find max compartment
    //     max_compartment= 0;
    //     for( li=sts->_s; li< sts->M+ sts->_s; li++){
    //         if( sts->init_cnt[li]> max_compartment){
    //             max_compartment= sts->init_cnt[li];
    //             max_compartmet_value= li;
    //         }
    //     }
    //     for( li= graph->_s; li< graph->_e; li++){
    //         sts->init_lst[li]= max_compartmet_value;
    //     }
    //     for( li=sts->_s; li< sts->M+ sts->_s; li++){
    //         if( li!= max_compartmet_value){
    //             for( j= 0; j< sts->init_cnt[li];){
    //                 lj= (size_t)((rand()/(double)RAND_MAX)*graph->V);
    //                 if( sts->init_lst[lj]!= max_compartmet_value) continue;
    //                 sts->init_lst[graph->_s+ lj%graph->V]= li;
    //                 j++;
    //             }
    //         }
    //     }
    // }
    // else{
    //     //MODE 3.
    //     LOG(2, __FILE__, __LINE__, "Status file mode [3]\n");
    //     rewind( fil_sts);
    //     fscanf( fil_sts, "%s", ch);
    //     if( strcmp(ch, "default")){
    //         printf(" check status file, read manual\n");
    //         exit( - 1);
    //     }
    //     rewind( fil_sts);
    //     if( fcolumn_count( fil_sts)!= 2){
    //         printf(" check status file, read manual\n");
    //         exit( - 1);
    //     }
    //     rewind( fil_sts);
    //     fscanf( fil_sts, "%s %zu", ch, &max_compartmet_value);
    //     for( li= graph->_s; li< graph->_e; li++){
    //         sts->init_lst[li]= max_compartmet_value;
    //     }
    //     sts->init_cnt[max_compartmet_value]= graph->V;
    //     for( li= 0; li< line_num - 1; li++){
    //         fscanf( fil_sts, "%d %zu", &ni, &ns);
    //         sts->init_lst[ni]= ns;
    //         sts->init_cnt[ns]++;
    //         sts->init_cnt[max_compartmet_value] --;
    //     }
    // }
   return (int)line_num;
}
/*
 *fgetline( skip blank line and trim space on both end)
 *
 *input:  FILE*  file       [ file pointer]
 *        size_t limit      [ buffer length]
 *output: char*  line       [ output buffer]
 *return: int    [0: reach end; >0: length of line read]
 */
int fgetline(FILE* file, char* line, size_t limit){
    int c;
    size_t i;
    i=0;
    //skip blank line, trim prefix
    while( isspace(c=fgetc( file)));
    while( c!= EOF&& c!='\n'&& i<limit - 1){
        line[i++]= (char)c;
        c=fgetc( file);
    }
    if( limit== i+1){
        printf("line truncated, exceeding max length[%zu]\n", limit);
        line[i]='\0';
        while( (c= fgetc( file))!= EOF&& c!='\n');
        return -1;
    }
    //trim end
    while( i> 1 && isspace( line[i - 1])){
        i --;
    }
    line[i]='\0';
    return (int)i;
}
/*
 *fget_next_item( skip blank line and trim space on both end AND line start with # AND stop at line starts with '[')
 *
 *input:  FILE*  file       [ file pointer]
 *        size_t limit      [ buffer length]
 *output: char*  line       [ output buffer]
 *return: int    [0: stop; >0: length of item read]
 */
int fget_next_item(FILE* file, char* line, size_t limit){
    int ret;

    while( (ret= fgetline( file, line, limit)>0)){
        if( ret< 0) return ret;
        else if( line[0]== '[') return 0;
        else if( line[0]== '#') continue;
        return ret;
    }
    return ret;
}

/*
 *fcheck_config( check configure of the target section suites requirment
 *
 *input:  FILE*  file       [ file pointer]
 *        size_t line       [ check line number]
 *        size_t colun      [ check column number in each line]
 *return: int    [0: success; <0: failure]
 */
int fcheck_config(FILE* file, char* target_section, size_t row, size_t column){
    size_t _row= 0;
    LINE line;

    rewind( file);
    locate_section_only( file, target_section);

    while( fget_next_item( file, line, MAX_LINE_LEN)){
        _row++;
        if( column> 0&& column_count( line)!= column){
            printf("section%s line[%zu] column number [%zu] error, should be[%zu]\n", target_section, _row,
                        column_count( line), column);
            return -1;
        }
    }
    if( row> 0&&_row!= row){
        printf("section%s line number [%zu] error, should be[%zu]\n", target_section, _row, row);
        return -1;
    }
    return 0;
}
/*
 *skip_top_comment( skip comment lead by 'c' on top)
 *
 *input:  FILE* file       [ file pointer]
 *        char  c          [ comment mark]
 *return: int   [<0: error; >= 0: number of line skipped]
 */
int skip_top_comment(FILE* file, char c){
    int i, comment_line= 0;
    LINE tmp_str;

    rewind( file);
    while( fgetline( file, tmp_str, MAX_LINE_LEN)){
        if( tmp_str[0]== c) comment_line++;
        else break;
    }
    rewind( file);
    for( i= 0; i< comment_line; i++){
        fgetline( file, tmp_str, MAX_LINE_LEN);
    }
    return comment_line;
}
LONG getValInt( FILE* fil, char* section, int echo){
    LONG ret;
    if( locate_section( fil, section)){
        ret = (LONG)time(NULL);
    } else {
        fscanf(fil, " %lld", &ret);
    }
    if( echo){
        printf("%s\t\t[%lld]\n", section, ret);
    }
    return ret;
}
double getValDbl( FILE* fil, char* section, int echo){
    double ret;
    locate_section( fil, section);
    fscanf(fil, " %lf", &ret);
    if( echo){
        printf("%s\t\t[%.4g]\n", section, ret);
    }
    return ret;
}
char *getValStr( FILE* fil, char* section, size_t limit, int echo){
    char* str;
    int ret;
    str= (char*)malloc(sizeof(char)*(limit+1));
    locate_section( fil, section);
    ret= fget_next_item( fil, str, limit);
    if( ret<= 0){
        printf("error: get %s failed\n", section);
        exit(-1);
    }
    if( echo){
        printf("%s\t\t[%s]\n", section, str);
    }
    return str;
}
size_t* getValSize_tLst( FILE* fil, char* section, size_t len, int echo){
    if( fcheck_config( fil, section, 1, len)< 0){
        printf("error: expecting 1 by %zu array for %s\n", len, section);
        exit( -1);
    }
    size_t* ret= (size_t*)malloc(sizeof(size_t)*len);
    locate_section( fil, section);
    if( echo){
        printf("%s\n", section);
    }
    if( ret== NULL){
        printf("error: Memory allocation failure for inducer list, size[%zu]\n", sizeof(size_t)*len);
        exit( - 1);
    }
    size_t l;
    for( l= 0; l< len; l++){
        fscanf(fil, " %zu", &ret[l]);
        if(echo){
            printf("  %zu", ret[l]);
        }
    }
    printf("\n");
    return ret;
}
double** getValMatrix( FILE* fil, char * section, size_t dim, size_t skip, size_t limit, int echo){
    double** mtx= (double **)malloc(sizeof(double*)*(dim+ skip));
    size_t i, j;
    if( mtx== NULL){
        printf("Memory allocation failure for nodal transisition matrix , size[%zu]\n", sizeof(double*)*(dim+ skip));
        exit( - 1);
    }
    for( j= skip; j< dim+ skip; j++){
        mtx[j]= (double*)malloc(sizeof(double)*(dim+ skip+ 1));
        if( mtx[j]== NULL){
            printf("Memory allocation failure for nodal transisition matrix row[%zu], size[%zu]\n", j, sizeof(double)*(dim+ skip+ 1));
            exit( - 1);
        }
    }
    //read in nodal transition rate matrix
    if( fcheck_config( fil, section, dim, dim)< 0){
        printf("read nodal transition rate matrix failed.\n");
        exit(-1);
    }
    locate_section( fil, section);
    double tmp_rat;
    char *tmp_str= (char*)malloc(sizeof(char)*limit);
    char* p_char;
    for( i= skip; i< dim+ skip; i++){
        tmp_rat= 0.0;
        fget_next_item( fil, tmp_str, limit);
        p_char= tmp_str;
        for( j= skip; j< dim+ skip; j++){
            _auto_sscanf( &p_char, "%lf ", &mtx[i][j]);
            tmp_rat+= mtx[i][j];
        }
        mtx[i][j]= tmp_rat;
    }
    if(echo){
        printf("%s\n", section);
        for( i= skip; i< dim+ skip; i++){
            for( j= skip; j< dim+ skip; j++){
                printf("%.2g\t", mtx[i][j]);
            }
            printf("\n");
        }
    }
    free(tmp_str);
    return mtx;
}
double*** getValMatrixLst( FILE* fil, char * section, size_t dim, size_t len, size_t skip, size_t limit, int echo){
    double*** mtxLst= (double ***)malloc(sizeof(double**)*(len));
    if( mtxLst== NULL){
        printf("Memory allocation failure for edge based transisition, size[%zu]\n", sizeof(double*)*(len));
        exit( - 1);
    }
    size_t i, j, k;
    for( i= 0; i< len; i++){
        mtxLst[i]= (double**)malloc(sizeof(double*)*(dim+ skip));
        if( mtxLst[i]== NULL){
            printf("Memory allocation failure for edge based transisition matrix layer[%zu], size[%zu]\n", i, sizeof(double*)*(dim+ skip));
            exit( - 1);
        }
        for( j= skip; j< dim+ skip; j++){
            mtxLst[i][j]= (double*)malloc(sizeof(double)*(dim+ skip+ 1));
            if( mtxLst[i][j]== NULL){
                printf("Memory allocation failure for edge based transisition matrix layer[%zu], row[%zu], size[%zu]\n", i, j,
                                                                                sizeof(double*)*(dim+ skip));
                exit( - 1);
            }
        }
    }
    //read in edge based transition rate matrix
    if( fcheck_config( fil, section, len*dim, dim)< 0){
        printf("read edge based transition rate matrix failed.\n");
        exit(-1);
    }
    locate_section( fil, section);
    double tmp_rat;
    char *tmp_str= (char*)malloc(sizeof(char)*limit);
    char* p_char;
    for( i= 0; i< len; i++){
        for( j= skip; j< dim+ skip; j++){
            tmp_rat= 0.0;
            fget_next_item( fil, tmp_str, limit);
            p_char= tmp_str;
            for( k= skip; k< dim+ skip; k++){
                _auto_sscanf( &p_char,  "%lf ", &mtxLst[i][j][k]);
                tmp_rat+= mtxLst[i][j][k];
            }
            mtxLst[i][j][k]= tmp_rat;
        }
    }
    if(echo){
        printf("%s\n", section);
        for( i= 0; i< len; i++){
            for( j= skip; j< dim+ skip; j++){
                for( k= skip; k< dim+ skip; k++){
                    printf("%.2g\t", mtxLst[i][j][k]);
                }
                printf("\n");
            }
            printf("\n");
        }
    }
    free(tmp_str);
    return mtxLst;
}
