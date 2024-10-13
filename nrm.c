#include "nrm.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include "pcg_basic.h"

#define SORT_NAME Edge
#define SORT_TYPE Edge
#define SORT_CMP(x, y)  ((x).i < (y).i ? -1 : ((y).i < (x).i ? 1 : 0))
#include "sort.h"

/*
 * nrm.c of GEMF in C language
 * Futing Fan
 * Kansas State University
 * Updates by Niema Moshiri (UC San Diego)
 * Last Modified: March 2018
 * Copyright (c) 2016, Futing Fan. All rights reserved.
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted
 */

void* malloc1( size_t l, size_t s);
double ** malloc2Dbl( size_t m, size_t n);
int ** malloc2Int( size_t  m, size_t  n);
NINT ** malloc2NINT( size_t m, size_t n);
double** init_inducer(Graph* graph, Status* sts, Transition* tran);
NINT** init_index(Graph* graph);
double get_rat_lst(Graph* graph, Transition* tran, Status* sts, double** p_raw_rat_lst, double** p_inducer_cal_lst);
void heart_beat( Heart_beat *hb);
void heap_init(Heap* heap, Graph* graph);
void heap_sort(Heap* heap);
double cal_new_tau(double r_old, double r_new, double t_old, double t);
double get_tau( Heap* heap, NINT n);
void heap_update( Heap* heap, Reaction *reaction);
void dump_heap( Heap* heap);
void print_inducer( Graph* graph, Transition* tran, Status *sts, Event* evt, FILE* fil_out);
void generate_graph(Graph* graph);
int nrm(Graph* graph, Transition* tran, Status* sts, Run* run){
    FILE* fil_out;
    size_t j, layer, compartment, section;
    size_t count= 0;
    int k, primary_case;
    NINT beg_num, end_num, cur_nod, i;
    int** p_nsim_avg_lst= NULL;
    //double T= 0.0;
    double tmp_double, elapse_tim;
    double *p_raw_rat_lst;
    double **p_inducer_cal_lst;
    double timer0, timer1;
    double R= 0.0;
    Heart_beat hb;
    Event evt;
    struct{
        double R;
        size_t* init_lst;
        double *p_raw_rat_lst;
        double **p_inducer_cal_lst;
    } restore;
    Heap heap;
    Reaction reaction;

    if(_LOGLVL_> 1){
        dump_transition(tran);
        dump_graph(graph);
        dump_status(sts);
    }
    pcg32_srandom((unsigned int)sts->random_seed, (unsigned int)sts->random_seed);
    //start timer
    timer0= gettimenow();

    heap.reaction= NULL;
    heap.idx= NULL;
    heap._s= 0;
    heap._e= 0;
    heap.V= 0;
    //generate and sort graph
    generate_graph(graph);

    // set all nodes to susceptible and out of subsample
    for( i= graph->_s; i< graph->_e; i++){
        if( sts->init_lst[i]){
            sts->init_lst[i] = 0;
        }
        if( sts->subsample_lst[i]){
            sts->subsample_lst[i] = 0;
        }
    }
    // set a new primary case (at compartment 2, as in Pekar)
    primary_case = pcg32_boundedrand(graph->_e);
    sts->init_lst[primary_case] = 2;
    // flag it as part of the subsample
    sts->subsample_lst[primary_case] = 1;
    // set compartment populations
    sts->init_cnt[0] = 4999999;
    sts->init_cnt[1] = 0;
    sts->init_cnt[2] = 1;
    for( compartment= 3; compartment< sts->M+ sts->_s; compartment++){
        sts->init_cnt[compartment] = 0;
    }
    printf("[initial population]\t[");
    for( i= sts->_s; i< sts->M+ sts->_s; i++){
        kilobit_print(" ", (LONG)sts->init_cnt[i], " ");
    }
    printf("]\n");

    //init inducer list
    p_inducer_cal_lst=  init_inducer( graph, sts, tran);
    //init adjacency index list
    init_index(graph);


    //calculate initial rate Ri for i in N
    R= get_rat_lst( graph, tran, sts, &p_raw_rat_lst, p_inducer_cal_lst);

    //open output file
    fil_out= fopen( run->out_file, "w");
    if( fil_out== NULL){
        printf("open output file[%s] faild\n", run->out_file);
        return -1;
    }

    // ***********************events happen***************************************
    if( run->sim_rounds> 1){
        //save initial status
        p_nsim_avg_lst= malloc2Int(sts->M, run->interval_num+ 1);
        restore.init_lst= (size_t*)malloc1(graph->_e, sizeof(size_t));
        restore.p_raw_rat_lst= (double*)malloc1(graph->_e, sizeof(double));
        restore.p_inducer_cal_lst= malloc2Dbl(graph->L, (size_t)graph->_e);

        restore.R= R;
        memcpy( restore.init_lst, sts->init_lst, sizeof(size_t)*(graph->_e));
        memcpy( restore.p_raw_rat_lst, p_raw_rat_lst, sizeof(double)*(graph->_e));
        for( layer= 0; layer< graph->L; layer++){
            memcpy( restore.p_inducer_cal_lst[layer], p_inducer_cal_lst[layer],  sizeof(double)*(graph->_e));
        }
        LOG(1, __FILE__, __LINE__, "save initial status success\n");
    }

    //alert&reset timer
    timer1= gettimenow() - timer0;
    time_print("preprocess time cost[ ", timer1, "]\n");
    hb.count= &count;
    hb.timer0= gettimenow();

    heap_init(&heap, graph);
    //repeat N times
    size_t round= 1;
    while(1){
        LOG(1, __FILE__, __LINE__, "Start simulation round [%zu/%zu]\n", round, run->sim_rounds);
        //reset count
        count= 0;
        //initial tau for all i
        heap.V= graph->_e- graph->_s;
        for( i= graph->_s; i< graph->_e; i++){
            heap.reaction[i].n= i;
            if( p_raw_rat_lst[i]> FLT_EPSILON){
                heap.reaction[i].t= - log(ldexp(pcg32_random(), -32))/(p_raw_rat_lst[i]);
            }
            else{
                heap.reaction[i].t= DBL_MAX;
            }
        }

        //make heap
        heap_sort(&heap);
        /*
            printf("before update[%d]\n", __LINE__);
        dump_heap( &heap);
            printf("after update[%d]\n", __LINE__);
            */

        while( 1){
            elapse_tim= heap.reaction[heap._s].t;
            if (run->max_time < elapse_tim){
                printf("T [%.6g] \treach limit [%6g], stop at [%zu] events.\t", elapse_tim, run->max_time, count);
                break;
            }
            else if(count>= run->max_events){
                printf("N [%zu] \treach limit [%zu], stop.\t", count, run->max_events);
                break;
            }
            //get a weighted radom node, ns-- active node, ni-- past_status, nj-- present_status
            evt.ns= heap.reaction[heap._s].n;

            get_next_evt(p_raw_rat_lst, p_inducer_cal_lst, graph, tran, sts, &evt, &heap);
            count++;
            sts->init_lst[evt.ns]= evt.nj;
            // track subsample
            if( (sts->init_cnt[0]>4950000) && (evt.ni==0)){
                sts->subsample_lst[evt.ns] = 1;
            }
            LOG(2, __FILE__, __LINE__, "event[%d], time[%.4g]\n", count, elapse_tim);
            //if run only once, output events details, else calculate intervals
            if( run->sim_rounds<=1){
                sts->init_cnt[evt.ni] --;
                sts->init_cnt[evt.nj] ++;
                // Only print subsample
                if( sts->subsample_lst[evt.ns]==1 &&
                    // and only if relevant to sampling, i.e. primary case
                    (evt.ns == primary_case ||
                        // ascertainment events
                        evt.nj == 3 ||
                        // recoveries of ascertained cases
                        (evt.ni == 3 && evt.nj == 7) ||
                        // hospitalizations
                        evt.nj == 6 ||
                        // recoveries of hospitalized cases
                        (evt.ni == 6 && evt.nj == 7)
                    )
                ){
                    fprintf( fil_out, "%lf %lf "fmt_n" %zu %zu", elapse_tim, R, evt.ns, evt.ni, evt.nj);
                    for( compartment= sts->_s; compartment< sts->M+ sts->_s; compartment++){
                        fprintf( fil_out, " %d", sts->init_cnt[compartment]);
                    }
                    if(run->show_inducer){
                        print_inducer( graph, tran, sts, &evt, fil_out);
                    }
                    fprintf( fil_out, "\n");
                }
            }
            else{
                //calculate intervals
                section= (size_t)((double)run->interval_num*(elapse_tim/ run->max_time));
                if( section>= run->interval_num){
                    printf("fatal error, wrong interval point value[%zu], max[%zu]\n", section, run->interval_num);
                    return -1;
                }
                p_nsim_avg_lst[evt.ni - sts->_s][section] --;
                p_nsim_avg_lst[evt.nj - sts->_s][section] ++;
            }

            //update rates
            //1. ni->nj
            //nodal transition rate
            tmp_double= tran->nodal_trn[evt.nj][sts->M+ sts->_s];
            //edge based transition rate
            for( layer= 0; layer< graph->L; layer++){
                tmp_double+= tran->edge_trn[layer][evt.nj][sts->M+ sts->_s]* p_inducer_cal_lst[layer][evt.ns];
            }
            if( tmp_double> FLT_EPSILON){
                reaction.t= - log(ldexp(pcg32_random(), -32))/(tmp_double)+ elapse_tim;
            }
            else{
                reaction.t= DBL_MAX;
            }
            reaction.n= evt.ns;
            /*
            printf("before update[%d]\n", __LINE__);
            dump_heap(&heap);
            printf("reaction[%d][%lf]LINE[%d]\n", reaction.n, reaction.t, __LINE__);
            printf("after update[%d]\n", __LINE__);
            */
            heap_update(&heap, &reaction);
            /*
            dump_heap(&heap);
            printf("after update[%d]\n", __LINE__);
            */
            R= R+ tmp_double - p_raw_rat_lst[evt.ns];
            p_raw_rat_lst[evt.ns]= tmp_double;
            //2. inducer_neighbour++/--
            for( layer= 0; layer< graph->L; layer++){
                k= 0;
                if( evt.ni== tran->inducer_lst[layer]){
                    k= - 1;
                }
                else if( evt.nj== tran->inducer_lst[layer]){
                    k= 1;
                }
                if( k != 0){
                    if( evt.ns== graph->_s){
                       beg_num= 0;
                    }
                    else{
                        beg_num= graph->index[layer][evt.ns];
                    }
                    end_num= graph->index[layer][evt.ns+1];
                    while( beg_num< end_num){
                        double change;
                        if( graph->weighted){
                            cur_nod= graph->edge_w[layer][beg_num].j;
                            change= k*graph->edge_w[layer][beg_num].w;
                        }
                        else{
                            cur_nod= graph->edge[layer][beg_num].j;
                            change= (double)k;
                        }
                        p_inducer_cal_lst[layer][cur_nod]+= change;
                        //adjust neighbour rate& total rate
                        tmp_double= change* tran->edge_trn[layer][sts->init_lst[cur_nod]][sts->M+ sts->_s];
                        R+= tmp_double;
                        //update affected rates and time
                        reaction.n= cur_nod;
                        reaction.t= cal_new_tau(p_raw_rat_lst[cur_nod], p_raw_rat_lst[cur_nod]+tmp_double, get_tau(&heap, cur_nod), elapse_tim);
                        /*
            printf("before update[%d]\n", __LINE__);
            dump_heap(&heap);
            printf("reaction[%d][%lf]LINE[%d]\n", reaction.n, reaction.t, __LINE__);
            printf("after update[%d]\n", __LINE__);
            */
                        heap_update(&heap, &reaction);
                        /*
            dump_heap(&heap);
            printf("after update[%d]\n", __LINE__);
            */
                        p_raw_rat_lst[cur_nod]+= tmp_double;
                        beg_num++;
                    }
                }
            }
            heart_beat(&hb);
        }
        printf("stop simulation round [%zu/%zu]\n", round, run->sim_rounds);
        // stop if the simulation was successful (at least 400 total infections and one active infection at end)
        if(++round> run->sim_rounds &&
                    (sts->init_cnt[0] <= 4999600 &&
                        (sts->init_cnt[1] ||
                            sts->init_cnt[2] ||
                            sts->init_cnt[3] ||
                            sts->init_cnt[4] ||
                            sts->init_cnt[5] ||
                            sts->init_cnt[6]
                        )
                    )
            ){
            break;
        }
        // otherwise, reset and try again
        // reset all nodes to susceptible and out of subsample
        for( i= graph->_s; i< graph->_e; i++){
            if( sts->init_lst[i]){
                sts->init_lst[i] = 0;
            }
            if( sts->subsample_lst[i]){
                sts->subsample_lst[i] = 0;
            }
        }
        // set a new primary case (to 2, as in Pekar)
        primary_case = (int)pcg32_boundedrand(graph->_e);
        sts->init_lst[primary_case] = 2;
        // flag it as part of the subsample
        sts->subsample_lst[primary_case] = 1;
        // set compartment populations
        sts->init_cnt[0] = 4999999;
        sts->init_cnt[1] = 0;
        sts->init_cnt[2] = 1;
        for( compartment= 3; compartment< sts->M+ sts->_s; compartment++){
            sts->init_cnt[compartment] = 0;
        }
        // regenerate graph
        generate_graph(graph);
        //re init adjacency index list
        init_index(graph);
        // free rate and inducer lists
        for( layer= 0; layer< graph->L; layer++){
            free( p_inducer_cal_lst[layer]);
        }
        free( p_inducer_cal_lst);
        free( p_raw_rat_lst);
        // intialise new rate and inducer lists
        p_inducer_cal_lst=  init_inducer( graph, sts, tran);
        R= get_rat_lst( graph, tran, sts, &p_raw_rat_lst, p_inducer_cal_lst);
        // close output and open a new one
        fclose( fil_out);
        fil_out= fopen( run->out_file, "w");
        if( fil_out== NULL){
            printf("open output file[%s] faild\n", run->out_file);
            return -1;
        }
        // Olde code for rerunning with the same intitial conditions
        //restore original status and run again
        // R= restore.R;
        // memcpy( sts->init_lst, restore.init_lst, sizeof(size_t)*(graph->_e));
        // memcpy( p_raw_rat_lst, restore.p_raw_rat_lst, sizeof(double)*(graph->_e));
        // for(layer= 0; layer< graph->L; layer++){
        //     memcpy( p_inducer_cal_lst[layer], restore.p_inducer_cal_lst[layer],  sizeof(double)*(graph->_e));
        // }
        // LOG(1, __FILE__, __LINE__, "End simulation round [%zu/%zu]\n", round, run->sim_rounds);
    }
    //post population
    if( run->sim_rounds<=1){
        printf("last moment population[ ");
        for( compartment= sts->_s; compartment< sts->M+ sts->_s; compartment++){
            kilobit_print(" ", sts->init_cnt[compartment], "");
        }
        printf(" ]\n");
    }
    else{
        //save results
        //output initial status count
        fprintf( fil_out, "0.0");
        tmp_double= elapse_tim/ run->interval_num;
        for( j= sts->_s; j< sts->M+ sts->_s; j++){
            fprintf( fil_out, " %lf", (double)sts->init_cnt[j]);
        }
        fprintf( fil_out, "\n%lf ", tmp_double);
        //output first interval
        for( j= 0; j< sts->M; j++){
            p_nsim_avg_lst[j][0]+= run->sim_rounds*sts->init_cnt[j+ sts->_s];
            fprintf( fil_out, "%lf ", p_nsim_avg_lst[j][0]/ (double)run->sim_rounds);
        }
        //output rest interval
        for( section= 1; section< run->interval_num; section++){
            fprintf( fil_out, "\n%lf ", (section+1)* tmp_double);
            for( j= 0; j< sts->M; j++){
                p_nsim_avg_lst[j][section]+= p_nsim_avg_lst[j][section - 1];
                fprintf( fil_out, "%lf ", p_nsim_avg_lst[j][section]/ (double)run->sim_rounds);
            }
        }
    }
    // ***********************events end***************************************
    //alert&end timer
    kilobit_print("events number[ ", (LONG)count, " ]\n");
    time_print("preprocess time cost[ ", timer1, "]\n");
    time_print("simulation time cost[ ", gettimenow() - hb.timer0, "]\n");

    //clean up
    LOG(1, __FILE__, __LINE__, "Begin clean up\n");
    for( layer= 0; layer< graph->L; layer++){
        free( p_inducer_cal_lst[layer]);
    }
    if( run->sim_rounds> 1){
        for( j= 0; j< sts->M; j++){
            free( p_nsim_avg_lst[j]);
        }
        free( p_nsim_avg_lst);
        free( restore.init_lst);
        free( restore.p_raw_rat_lst);
        for( layer= 0; layer< graph->L; layer++){
            free( restore.p_inducer_cal_lst[layer]);
        }
        free( restore.p_inducer_cal_lst);
    }
    free( p_inducer_cal_lst);
    free( p_raw_rat_lst);
    if( heap.reaction!= NULL){
        free( heap.reaction);
    }
    if( heap.idx!= NULL){
        free( heap.idx);
    }

    fclose( fil_out);
    LOG(1, __FILE__, __LINE__, "End clean up\n");
    return 0;
}

//chose a weighted random node from the network
size_t weighed_rat_rand(double* rat_lst, size_t len){
    double* tmp_rat_sum;
    double key;
    size_t i, left, right, mid, ret;
    tmp_rat_sum = (double*)malloc(sizeof(double)* (len+1))+ 1;
    if( tmp_rat_sum== NULL){
        printf("Memory allocation failure for temp rate sumation list, size[%zu]\n", sizeof(double)*(len+ 1));
        exit( - 1);
    }
    tmp_rat_sum[- 1]= 0.0;

    //accumulate weight list
    for( i= 0; i< len; i++){
        tmp_rat_sum[i]= rat_lst[i]+ tmp_rat_sum[i - 1];
    }
    left= 0;
    right= len - 1;
    key= (ldexp(pcg32_random(), -32))*tmp_rat_sum[right];

    //binary search target section
    while(1){
        //stop condition
        if( left== right){
            ret= left;
            break;
        }
        else if( left== right - 1){
            if( key< tmp_rat_sum[left]){
                ret= left;
                break;
            }
            else{
                ret= right;
                break;
            }
        }

        //binary
        mid= (left+right)/2;
        if( key < tmp_rat_sum[mid]){
            right= mid;
        }
        else{
            left= mid;
        }
    }
    while(1){
        if( ret== 0) break;
        if( tmp_rat_sum[ret] < tmp_rat_sum[ret - 1]){
            ret --;
        }
        else{
            break;
        }
    }
    free(tmp_rat_sum - 1);
    return ret;
}

int get_next_evt(double* p_raw_rat_lst, double** p_inducer_cal_lst, Graph* graph, Transition* tran, Status* sts, Event* evt, Heap* heap){
    double* nodal_tmp_rat_lst, *edgeb_tmp_rat_lst;
    double nodal_rat_ttl, edgeb_rat_ttl;
    size_t layer, i, j;
    //pick out one node randomly by weight
    //evt->ns= weighed_rat_rand( p_raw_rat_lst+ graph->_s, graph->V)+ graph->_s;
    evt->ni= sts->init_lst[evt->ns];
    nodal_tmp_rat_lst= (double*)malloc(sizeof(double)*(sts->M+ sts->_s+ 1))+ 1;
    if( nodal_tmp_rat_lst== NULL){
        printf("Memory allocation failure for temp nodal rate sumation list, size[%zu]\n",
                                        sizeof(double)*(sts->M+ sts->_s+ 1));
        exit( - 1);
    }
    edgeb_tmp_rat_lst= (double*)malloc(sizeof(double)*(graph->L*(sts->M+sts->_s+ 1)))+ 1;
    if( edgeb_tmp_rat_lst== NULL){
        printf("Memory allocation failure for temp edge based rate sumation list, size[%zu]\n",
                                        sizeof(double)*(sts->M*graph->L+ sts->_s+ 1));
        exit( - 1);
    }
    nodal_tmp_rat_lst[sts->_s - 1]= 0.0;
    edgeb_tmp_rat_lst[sts->_s - 1]= 0.0;
    nodal_rat_ttl= 0.0;
    edgeb_rat_ttl= 0.0;
    for( j= sts->_s; j< sts->M+ sts->_s; j++){
        nodal_tmp_rat_lst[j]= tran->nodal_trn[evt->ni][j];
        nodal_rat_ttl+= nodal_tmp_rat_lst[j];
    }
    for( layer= 0; layer< graph->L; layer++){
        for( j= sts->_s; j< sts->M+ sts->_s; j++){
            edgeb_tmp_rat_lst[layer* (sts->M)+ j]= tran->edge_trn[layer][evt->ni][j]*p_inducer_cal_lst[layer][evt->ns];
            edgeb_rat_ttl+= edgeb_tmp_rat_lst[layer* (sts->M)+ j];
        }
    }
    if(ldexp(pcg32_random(), -32)< nodal_rat_ttl/(nodal_rat_ttl+ edgeb_rat_ttl)){
        //nodal
        i= weighed_rat_rand(nodal_tmp_rat_lst+ sts->_s, sts->M);
    }
    else{
        //edgebased
        i= weighed_rat_rand(edgeb_tmp_rat_lst+ sts->_s, sts->M*graph->L);
    }
    evt->nj= i%sts->M+ sts->_s;
    free( nodal_tmp_rat_lst - 1);
    free( edgeb_tmp_rat_lst - 1);
    return 0;
}
void* malloc1( size_t l, size_t s){
    void* ret = malloc(s*l);
    if( ret== NULL){
        printf("Memory allocation failure, size[%zu]\n", s*l);
        exit( - 1);
    }
    memset(ret, 0, s*l);
    return ret;
}
double** init_inducer(Graph* graph, Status* sts, Transition* tran){
    LOG(1, __FILE__, __LINE__, " initial inducer\n");
    double** mtx=  malloc2Dbl( graph->L, (size_t)graph->_e);
    size_t l;
    size_t li;
    if( graph->weighted== 1){
        Edge_w * p_ew;
        for( l=0; l< graph->L; l++){
            p_ew= graph->edge_w[l];
            for( li= 0; li< graph->E[l]; li++){
                if( sts->init_lst[p_ew->i] == tran->inducer_lst[l]){
                    mtx[l][p_ew->j]+= p_ew->w;
                }
                p_ew++;
            }
        }
    }
    else{
        Edge * p_e;
        for( l=0; l< graph->L; l++){
            p_e= graph->edge[l];
            for( li= 0; li< graph->E[l]; li++){
                if( sts->init_lst[p_e->i] == tran->inducer_lst[l]){
                    mtx[l][p_e->j]++;
                }
                p_e++;
            }
        }
    }
    LOG(1, __FILE__, __LINE__, " initial inducer success\n");
    return mtx;
}
NINT** init_index(Graph* graph){
    LOG(1, __FILE__, __LINE__, " initial index\n");
    graph->index= malloc2NINT( graph->L, (size_t)graph->_e+1);
    NINT cur_nod, next_nod, count;
    int init_flag;
    size_t layer, li;
    for( layer= 0; layer< graph->L; layer++){
        //cur_nod= -1;
        init_flag= 1;
        for( li= 0; li< graph->E[layer]; li++){
            //next_nod= next_node
            if( graph->weighted){
                next_nod= graph->edge_w[layer][li].i;
            }
            else{
                next_nod= graph->edge[layer][li].i;
            }
            //if( cur_nod== -1){
            if( init_flag){
                init_flag= 0;
                cur_nod= next_nod;
                count= 1;
            }
            else if( cur_nod== next_nod){
                count++;
            }
            else if( cur_nod> next_nod){
                printf("cur_nod "fmt_n" next_nod "fmt_n"\n", cur_nod, next_nod);
                printf("fatal error:[data matrix sort failed]\n");
                exit( -1);
            }
            else{
                while(cur_nod< next_nod){
                    graph->index[layer][cur_nod+1]= count;
                    cur_nod++;
                }
                count++;
            }
        }
        graph->index[layer][cur_nod+1]= count;
    }
    LOG(1, __FILE__, __LINE__, " initial index success\n");
    return graph->index;
}
double get_rat_lst(Graph* graph, Transition* tran, Status* sts, double** p_raw_rat_lst, double** p_inducer_cal_lst){
    LOG(1, __FILE__, __LINE__, "calculate initial Ri\n");
    double ret= 0.0, tmp_double;
    NINT i;
    size_t layer;
    *p_raw_rat_lst= (double*)malloc1((size_t)graph->_e, sizeof(double));
    for( i= graph->_s; i< graph->_e; i++){
        tmp_double= 0;
        //add nodal tran rate
        tmp_double+= tran->nodal_trn[sts->init_lst[i]][sts->M+ sts->_s];
        //add edge based tran rate
        for( layer= 0; layer< graph->L; layer++){
            if( p_inducer_cal_lst[layer][i] > 0){
                tmp_double+= tran->edge_trn[layer][sts->init_lst[i]][sts->M+ sts->_s]* p_inducer_cal_lst[layer][i];
            }
        }
        (*p_raw_rat_lst)[i]= tmp_double;
        ret+= tmp_double;
    }
    LOG(1, __FILE__, __LINE__, "calculate initial Ri success\n");
    return ret;
}
int Edge_cmp( const void *a, const void *b){
    LONG ret= (LONG)((Edge*)a)->i - (LONG)((Edge*)b)->i;
    if( ret== 0){ return 0; }
    else if( ret > 0){ return 1; }
    else return -1;
}
void heart_beat( Heart_beat *hb){
    if( -- hb->count_down <= 0){
        if( hb->count_down< 0){
            hb->last_report_time= gettimenow();
            hb->timer2= gettimenow() - hb->timer0;
            hb->count_down= 1;
        }
        else{
            //first 1 minute, report every 5 seconds
            if( hb->timer2< 60){
                hb->frequency= 5;
            }
            //first 10 minute, report every 30 seconds
            else if( hb->timer2< 10*60){
                hb->frequency= 30;
            }
            //first 60 minute, report every 2 minutes
            else if( hb->timer2< 60*60){
                hb->frequency= 2*60;
            }
            //report every 10 minutes
            else{
                hb->frequency= 10*60;
            }
            hb->count_down= (int)(hb->frequency/((gettimenow() - hb->last_report_time)/(*(hb->count) - hb->last_report)));
            hb->last_report_time= gettimenow();
            hb->timer2= gettimenow() - hb->timer0;
            time_print("elapse time[", hb->timer2, "]");
            kilobit_print(", [", (LONG)*(hb->count), "]events generated\n");
            hb->last_report= *(hb->count);
        }
    }
}

double ** malloc2Dbl( size_t m, size_t n){
    double ** ret = (double**)malloc(sizeof(double*)*(m));
    size_t l;
    if( ret== NULL){
        printf("2-D Memory allocation failure, size[%zu]\n", sizeof(double*)*m);
        exit( - 1);
    }
    for( l= 0; l< m; l++){
        ret[l]= (double*)malloc(sizeof(double)*(n));
        if( ret[l]== NULL){
            printf("2-D Memory allocation failure for layer[%zu/%zu], size[%zu]\n", l+ 1, m, sizeof(double)*(n));
            exit( - 1);
        }
        memset(ret[l], 0, sizeof(double)*(n));
    }
    return ret;
}
NINT ** malloc2NINT( size_t m, size_t n){
    NINT ** ret = (NINT**)malloc(sizeof(NINT*)*(m));
    size_t l;
    if( ret== NULL){
        printf("2-D Memory allocation failure, size[%zu]\n", sizeof(NINT*)*m);
        exit( - 1);
    }
    for( l= 0; l< m; l++){
        ret[l]= (NINT*)malloc(sizeof(NINT)*(n));
        if( ret[l]== NULL){
            printf("2-D Memory allocation failure for layer[%zu/%zu], size[%zu]\n", l+ 1, m, sizeof(NINT)*(n));
            exit( - 1);
        }
        memset(ret[l], 0, sizeof(NINT)*(n));
    }
    return ret;
}
int ** malloc2Int( size_t m, size_t n){
    int ** ret = (int**)malloc(sizeof(int*)*(m));
    size_t l;
    if( ret== NULL){
        printf("2-D Memory allocation failure, size[%zu]\n", sizeof(int*)*m);
        exit( - 1);
    }
    for( l= 0; l< m; l++){
        ret[l]= (int*)malloc(sizeof(int)*(n));
        if( ret[l]== NULL){
            printf("2-D Memory allocation failure for layer[%zu/%zu], size[%zu]\n", l+ 1, m, sizeof(int)*(n));
            exit( - 1);
        }
        memset(ret[l], 0, sizeof(int)*(n));
    }
    return ret;
}
int Reaction_cmp(const void* a, const void* b){
    double d = ((Reaction*)a)->t - ((Reaction*)b)->t;
    if( d > 0) return 1;
    if( d < 0) return -1;
    return 0;
}
void heap_init(Heap* heap, Graph* graph){
    heap->_s= graph->_s;
    heap->_e= graph->_e;
    heap->reaction= (Reaction*)malloc(sizeof(Reaction)*heap->_e);
    if( heap->reaction== NULL){
        printf("malloc reaction list failed.\n");
        exit(-1);
    }
    heap->idx= (NINT*)malloc1((size_t)heap->_e, sizeof(NINT));
}
void heap_sort(Heap* heap){
    qsort( heap->reaction+ heap->_s, heap->V, sizeof( Reaction), Reaction_cmp);
    NINT i;
    for( i= heap->_s; i< heap->_e; i++){
        heap->idx[heap->reaction[i].n]= i;
    }
}
void heap_swap(Heap* heap, NINT a, NINT b){
    Reaction tr;
    tr.n= heap->reaction[a].n;
    tr.t= heap->reaction[a].t;
    heap->reaction[a].n= heap->reaction[b].n;
    heap->reaction[a].t= heap->reaction[b].t;
    heap->reaction[b].n= tr.n;
    heap->reaction[b].t= tr.t;
    heap->idx[heap->reaction[a].n]= a;
    heap->idx[heap->reaction[b].n]= b;
}
void heap_update_aux( Heap* heap, NINT n){
    NINT parent= (n+heap->_s-1)/2;
    NINT lchild= 2*n+ 1-heap->_s;
    NINT rchild= 2*n+ 2-heap->_s;
    if( n> heap->_s&& heap->reaction[parent].t > heap->reaction[n].t){
        heap_swap(heap, n, parent);
        heap_update_aux(heap, parent);
    }
    else if( lchild< (heap->_s+ heap->V)){
        if( rchild< (heap->_s+ heap->V)&& heap->reaction[rchild].t< heap->reaction[lchild].t && heap->reaction[rchild].t< heap->reaction[n].t){
            heap_swap(heap, rchild, n);
            heap_update_aux(heap, rchild);
        }
        else if( heap->reaction[lchild].t< heap->reaction[n].t){
            heap_swap(heap, lchild, n);
            heap_update_aux(heap, lchild);
        }
    }
}
void heap_update( Heap* heap, Reaction *reaction){
    NINT n= heap->idx[reaction->n];
    heap->reaction[n].t= reaction->t;;
    heap_update_aux(heap, n);
}
/*void heap_update( Heap* heap, Reaction *reaction){
    int n= heap->idx[reaction->n];
    //swap until ok
    int parent= (n+ heap->_s- 1)/2;
    while( n> heap->_s&& heap->reaction[parent].t > reaction->t){
        heap->reaction[n].n= heap->reaction[parent].n;
        heap->reaction[n].t = heap->reaction[parent].t;
        heap->idx[heap->reaction[n].n]= n;
        n= parent;
        parent= (n+ heap->_s- 1)/2;
    }
    heap->reaction[n].n= reaction->n;
    heap->reaction[n].t = reaction->t;
    heap->idx[reaction->n]= n;
}
*/
double cal_new_tau(double r_old, double r_new, double t_old, double t){
    if( r_new< FLT_EPSILON) return DBL_MAX;
    if( r_old< FLT_EPSILON) return (- log(ldexp(pcg32_random(), -32))/(r_new)+ t);
    return (r_old/r_new)*(t_old- t)+ t;
}
double get_tau( Heap* heap, NINT n){
    return heap->reaction[heap->idx[n]].t;
}
void dump_heap( Heap* heap){
    printf("begin dump heap\n");
    for( int i= 0; i< heap->_e; i++){
        if( heap->idx[i]> 10000){ exit(-1);}
        printf("[%d][%d][%.5g] index[" fmt_n "]\n", i, heap->reaction[i].n, heap->reaction[i].t, heap->idx[i]);
    }
    printf("end dump heap\n");
    fflush(stdout);
}
void print_inducer( Graph* graph, Transition* tran, Status* sts, Event* evt, FILE* fil_out){
    fprintf( fil_out, " [");
    if( tran->nodal_trn[evt->ni][evt->nj]> 0){
        fprintf( fil_out, "%d", evt->ns);
    }
    for( int layer= 0; layer< graph->L; layer++){
        fprintf( fil_out, "],[");
        if( tran->edge_trn[layer][evt->ni][evt->nj]> 0){
            int flag= 0;
            if( graph->weighted){
                for( int i= graph->index[layer][evt->ns]; i< graph->index[layer][evt->ns+1]; i++){
                    if( sts->init_lst[graph->edge_w[layer][i].j] == tran->inducer_lst[layer]){
                        if( flag){
                            fprintf( fil_out, ",");
                        }
                        else{
                            flag= 1;
                        }
                        fprintf( fil_out, "%d", graph->edge_w[layer][i].j);
                    }
                }
            }
            else{
                for( int i= graph->index[layer][evt->ns]; i< graph->index[layer][evt->ns+1]; i++){
                    if( sts->init_lst[graph->edge[layer][i].j] == tran->inducer_lst[layer]){
                        if( flag){
                            fprintf( fil_out, ",");
                        }
                        else{
                            flag= 1;
                        }
                        fprintf( fil_out, "%d", graph->edge[layer][i].j);
                    }
                }
            }
        }
    }
    fprintf( fil_out, "]");
}

void generate_graph(Graph* graph){
    printf("Generating network...\n");

    const int m = 8;         // Initial number of nodes
    const int n = 5000000;    // Total number of nodes to be generated
    int i, j, k;
    int new_nodes_edge_targets[m];
    int *edge_ends = malloc((m + 2 * m * (n - m)) * sizeof(int));
    size_t layer;
    double t0 = gettimenow();

    if (edge_ends == NULL) {
        perror("Error allocating memory");
        exit(EXIT_FAILURE);
    }

    // Initialize the edge end list with eight unconnected nodes
    for (i = 0; i < m; ++i) {
        edge_ends[i] = i;
    }

    // Add nodes until we have n total
    for (i = m; i < n; ++i) {
        // For each new node, connect it with edges to m different nodes from the edge end list
        for (j = 0; j < m; ++j) {
            int is_duplicate;
            do {
                // draw random node from the edge end list (this provides preferential attachment)
                new_nodes_edge_targets[j] = edge_ends[(int)pcg32_boundedrand(m + 2 * m * (i - m))];
                // Check each new node against the previous ones and draw again if duplicate
                is_duplicate = 0;
                for (k = 0; k < j; ++k) {
                    if (new_nodes_edge_targets[k] == new_nodes_edge_targets[j]) {
                        is_duplicate = 1;
                        break;
                    }
                }
            } while (is_duplicate);
        }
        // for each edge of the new node, add the nodes to edge ends and the edges to the edge list
        for (j = 0; j < m; ++j) {
            // Add the nodes to edge ends
            edge_ends[m + 2 * m * (i - m) + 2 * j] = i;
            edge_ends[m + 2 * m * (i - m) + 2 * j + 1] = new_nodes_edge_targets[j];
            // Add the edges to the edge list
            graph->edge[0][m*(i - m) + j].i= i;
            graph->edge[0][m*(i - m) + j].j= new_nodes_edge_targets[j];
            graph->edge[0][graph->E[0]/2+ m*(i - m) + j].i= new_nodes_edge_targets[j];
            graph->edge[0][graph->E[0]/2+ m*(i - m) + j].j= i;
        }
    }
    free(edge_ends);
    // sort the ede list and then copy to the other layers
    Edge_tim_sort( graph->edge[0], graph->E[0]);
    for( layer= 1; layer< graph->L; layer++){
        memcpy(graph->edge[layer], graph->edge[0], graph->E[0] * sizeof(Edge));
    }
    time_print( "generation time cost[ ", gettimenow() - t0, " ]\n");
}
