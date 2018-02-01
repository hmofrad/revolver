/*
 * libgraph.h Graph processing header file
 * (c) Mohammad H. Mofrad, 2017
 * (e) mohammad.hmofrad@pitt.edu
 */
#ifndef LIBGRAPH_H
#define LIBGRAPH_H
 
#include "macros.h"

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>
#include <pthread.h>
#include <sys/types.h>
#include <math.h>
#include <libgen.h>
#include <sys/mman.h>
#include <fcntl.h>
//#include <sys/types.h>
//#include <sys/stat.h>


#include "libmath.h"
#include "libla.h"


struct Adjacency_list_node
{
    unsigned int destination;
    int count;  // Outgoing count
    int icount; // Ingoing count
    struct Adjacency_list_node *next;
};

struct Adjacency_list
{
	float score;
    unsigned int source;
    int degree;  // Outgoing degree
    int idegree; // Ingoing degree
    int label;
	int temp_label;
	char allocated;
    struct Adjacency_list_node *head;
	/* Revolver */
	int signal;
	int action;
	float *probability;
    float *signal_;
    float probability_pr; // pagerank probability
};

struct List_of_Adjacency_list_node
{
    struct Adjacency_list *nodes;
    struct List_of_Adjacency_list_node *next;
};

struct List_of_Adjacency_list
{
    char allocated;
    struct List_of_Adjacency_list_node *head;
//    struct Adjacency_list *nodes;
//    struct List_of_Adjacency_list *next;
};

struct Algorithm
{
    char name[50];
	char dir_name[100];
    char base_name[1000];
    int num_partitions;
    int file_offset;
	float l; // lambda
    float c; // load imbalance
    float o; // omega
    int    e; // epsilon
    unsigned int *label;
    unsigned int *capacity;
    unsigned int *load;
    float *load_weight;
    signed int *remaining_capacity;
	unsigned int *candidate_vertices;
	float *migration_probability;
	/* Revolver */
	float alpha;
	float beta;
};


struct Graph
{
    long long num_edges;
    long long num_alloced_edges; // Undirected edges
    long long num_nodes;
    long long num_alloced_nodes;
    unsigned int *mapped_nodes;
    //unsigned int *mapped_nodes_;
    struct Adjacency_list *nodes;
    struct List_of_Adjacency_list *nodes_list;
    struct Algorithm *algorithm;
};

struct pthread_args_struct{
	int iteration;
	int index;
	unsigned int low;
	unsigned int high;
    FILE *fd;
    char *fmt;
    int cols;
};

extern struct Graph *graph;

extern pthread_barrier_t barrier;
extern pthread_mutex_t mutex;
extern pthread_mutex_t mutex1;
extern float score;
extern char  sentinel;


extern float score_array[MAX_THREADS_NUM];
extern long long cut_array[MAX_THREADS_NUM];
extern long long local_edges_array[MAX_THREADS_NUM];
extern unsigned int *external_edges_array_per_partition[MAX_THREADS_NUM];
extern unsigned int *internal_edges_array_per_partition[MAX_THREADS_NUM];
extern unsigned int *ingoing_cut_array_per_partition[MAX_THREADS_NUM];
extern unsigned int *candidate_vertices_array[MAX_THREADS_NUM];
extern unsigned int max_source_array[MAX_THREADS_NUM];
extern unsigned int max_destination_array[MAX_THREADS_NUM];
extern unsigned int max_array[MAX_THREADS_NUM];
extern unsigned int min_source_array[MAX_THREADS_NUM];
extern unsigned int min_destination_array[MAX_THREADS_NUM];
extern unsigned int min_array[MAX_THREADS_NUM];
extern long long    num_nodes_array[MAX_THREADS_NUM];
extern unsigned int num_alloced_nodes_array[MAX_THREADS_NUM];
extern unsigned int num_edges_array[MAX_THREADS_NUM];
extern unsigned int num_alloced_edges_array[MAX_THREADS_NUM];
extern struct Adjacency_list *nodes_array[MAX_THREADS_NUM];
extern unsigned int *mapped_nodes_array[MAX_THREADS_NUM];
extern unsigned int thread_low_index_array[MAX_THREADS_NUM];
extern unsigned int thread_high_index_array[MAX_THREADS_NUM];
extern unsigned int probability_array[MAX_THREADS_NUM];
extern long long external_messages_array[MAX_THREADS_NUM];
extern long long  internal_messages_array[MAX_THREADS_NUM];

int  edge_open_file(char *file_name, char *comment_style, char *delimiter, int columns);
char *edge_format_specifier(FILE* file_descriptor, char *comment_style, char *delimiter);
int  vertex_open_file(char *file_name, char *delimiter);
char *vertex_format_specifier(char *delimiter);
int  vertex_map();
int  add_edge(struct Adjacency_list **nodes, long long *num_nodes, unsigned int source, unsigned int destination, char type);
struct Adjacency_list_node *malloc_edge(int destination);
int  algorithm_initialization();
int  label_initialization();
int  load_initialization();
int  score_initialization();
int  migration_initialization();
int  temporary_array_initialization();
int  learning_automata_initialization();
int  computation_is_halted(float score, float score_old, int *ep);

void print_graph();
void print_edges_per_vertex();

void *revolver(void *arguments);
void *spinner(void *arguments);
void statistics(int iteration);
int write_partition();
int write_partition_e(char *file_name, char *comment_style, char *delimiter, int columns);

int  return_label(float *array, int count, int label);
void score_function(float *array, int count, int index);
void objective_function(float *array, int count, int index);

void enumerate_nodes();
void enumerate_edges();
void enumerate_candidate_vertices();
void calculate_remaining_capacity();

long long phi(unsigned int low, unsigned int high); // ratio of local edges
long long cut_size(unsigned int low, unsigned int high);
unsigned int enumerate_partition_edges(unsigned int low, unsigned int high, int label, int type);
unsigned int enumerate_partition_cut(unsigned int low, unsigned int high, int label, int type);

void compute_score(unsigned int low, unsigned int high);
int  migrate_vertex_spinner(int index);
void calculate_migration_probability();
int  migrate_vertex_revolver(int index);
void calculate_migration_probability_();
void calculate_first_migration_probability();

int  para_edge_open_file(char *file_name, char *comment_style, char *delimiter, int columns);
void multiprocessing_pool(unsigned int total_size, int num_chunks, int *k1, unsigned int *k1_num, int *k2, unsigned int *k2_num);
int para_add_edge(struct Adjacency_list **nodes, unsigned int **mapped_nodes, long long *num_nodes, unsigned int source, unsigned int destination, char edge_type);
void *random_partitioning(void *arguments);
void *hash_partitioning(void *arguments);
void *range_partitioning(void *arguments);
void *read_chunk(void *arguments);
void *pagerank(void *arguments);
long long int refine_position(FILE *file_descriptor, long long position);
int edge_open_file_mmapped(char *file_name, char *comment_style, char *delimiter, int columns);

void free_all();

#endif

