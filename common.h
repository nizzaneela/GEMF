#ifndef COMMONH
#define COMMONH

#define MAX_INT_LEN 20
#define MAX_LINE_LEN 16*1024

#ifdef WIN_X64
#define _CRT_SECURE_NO_DEPRECATE
//#define _CRT_SECURE_NO_WARNINGS
#include <windows.h>
#include <time.h>
#else
#include <sys/time.h>
#endif
#include <stddef.h>

typedef char LINE[MAX_LINE_LEN];
typedef long long LONG;

//node serial type
typedef unsigned int NINT;
#define fmt_n "%d"

typedef struct
{
    //adjacency list i->j
    NINT i;
    NINT j;
} Edge;
typedef struct
{
    //weighted adjacency list i->j and weight w
    NINT i;
    NINT j;
    double w;
} Edge_w;
typedef struct
{
    Edge **edge;
    Edge_w **edge_w;
    //number of nodes
    NINT V;
    //nodes start from _s, end at _e - 1
    NINT _s;
    NINT _e;
    //edge num list for each layer
    size_t *E;
    //weighted flag, 0 for unweighted, otherwise weighted
    int weighted;
    //directed flag, 0 for undirected, otherwise directed
    int directed;
    //number of layers
    size_t L;
    //if network is sorted, index[i] is the position of the last edge oriented from node i;
    NINT** index;
} Graph;
typedef struct
{
    //number of compartments
    size_t M;
    //number of layers
    size_t L;
    //compartments start from _s, end at _s+M -1
    size_t _s;
    //2D array, (_s+M+1) by (_s+M)
    //1. if nodal_trn[i][j], value is invalid if i< _s or j< _s
    //2. if (i>= _s && j== _s+M)(meaning the last col for each valid row), nodal_trn[i][j]= sum(nodal_trn[i][_s]...nodal_trn[i][_s+M-1])
    double **nodal_trn;
    //3D array, L by (_s+M+1) by (_s+M)
    //1. if nodal_trn[*][i][j], value is invalid if i< _s or j< _s
    //2. if (i>= _s && j== _s+M)(meaning the last col for each valid row), nodal_trn[*][i][j]= sum(nodal_trn[*][i][_s]...nodal_trn[*][i][_s+M-1])
    double ***edge_trn;
    //1 by L array, inducer for each layer
    size_t *inducer_lst;
} Transition;
typedef struct
{
    //number of compartments
    size_t M;
    //compartments start from _s, end at _s+M -1
    size_t _s;
    //number of nodes
    NINT _node_V;
    //nodes start from _node_s, end at _node_e - 1
    NINT _node_s;
    NINT _node_e;
    //1 by _node_e list, initial status for each node
    size_t *init_lst;
    //1 by _node_e list, flag if subsample
    size_t *subsample_lst;
    //1 by (_s+M) array
    //the population of each compartment
    NINT *init_cnt;
    //random number seed
    int random_seed;
} Status;
typedef struct
{
    //arbitrary stop time
    double max_time;
    //maximum events number
    size_t  max_events;
    //number of rounds
    size_t sim_rounds;
    //sampling interval number for multiple simulation( histogram like)
    size_t interval_num;
    char *out_file;
    int show_inducer;
} Run;
typedef struct{
    //node of the event
    NINT ns;
    //status change from ni to nj
    size_t ni;
    size_t nj;
} Event;

typedef struct{
    //time t
    double t;
    //node
    NINT n;
} Reaction;
typedef struct{
    //reactoin list
    Reaction* reaction;
    //index fo locate the position of each node
    NINT* idx;
    //nodes start from _s, end at _e-1
    NINT _s;
    NINT _e;
    //nodes number V
    NINT V;
} Heap;

double gettimenow();

extern int _LOGLVL_;
int LOG(int loglvl, const char* file, int line, char* format, ...);
//dump graph
void dump_graph(Graph* graph);
//print graph size info
void print_graph_size(Graph* graph);
int check_int_range( LONG li);
//exponential
LONG exp10( int pow);
//print with kilobit
int kilobit_print( char* prefix, LONG val, char* suffix);
//print time Hour:min:sec
int time_print( char* prefix, double val, char* suffix);
//dump transition matrices
void dump_transition( Transition* tran);
//dump status info
void dump_status( Status* sts);

#endif
