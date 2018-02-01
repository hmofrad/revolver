/*
 * libstat.c Statistics library
 * (c) Mohammad H. Mofrad, 2017
 * (e) mohammad.hmofrad@pitt.edu
 */

#include "libgraph.h"

struct Graph *graph;
float score_array[MAX_THREADS_NUM];
long long cut_array[MAX_THREADS_NUM];
long long local_edges_array[MAX_THREADS_NUM];
unsigned int *external_edges_array_per_partition[MAX_THREADS_NUM];
unsigned int *internal_edges_array_per_partition[MAX_THREADS_NUM];
unsigned int *ingoing_cut_array_per_partition[MAX_THREADS_NUM];
unsigned int *candidate_vertices_array[MAX_THREADS_NUM];

int main(int argc, char *argv[])
{
    long long i = 0; // num_nodes/thread index
    long long j = 0; // mapped_nodes index
	int k = 0; // num_partitions index   
    
    if((argc != 13))
    {
        fprintf(stderr, "USAGE: %s -n <num_partitions> -a <revolver|spinner> -f [<edge_file_name> <comment_style> <delimiter> <columns>] -v [<vetex_file_name> <delimiter>] \n", argv[0]);
        for(i = 0; i < argc; i++)
            fprintf(stderr, "argv[%2lli]: %s\n", i, argv[i]);
        
        exit(EXIT_FAILURE);
    }
    
	graph = malloc(sizeof(struct Graph));
    if(!graph)
        exit(EXIT_FAILURE);
    memset(graph, '\0', sizeof(struct Graph));
	
    char comment[strlen(argv[7])];
    memset(comment, '\0', strlen(argv[7]));
    memcpy(comment, argv[7], strlen(argv[7]));
    
    char delimiter[strlen(argv[8])];
    memset(delimiter, '\0', strlen(argv[8]));
    memcpy(delimiter, argv[8], strlen(argv[8]));
    
    int columns = atoi(argv[9]);
    (void) columns;
/*    
    int status = edge_open_file(argv[6], argv[7], argv[8], atoi(argv[9]));
	if(status == EXIT_ERROR)
	{
		fprintf(stderr, "Error: open_edge_file()\n");
        exit(EXIT_FAILURE);
	}
*/
    int status = para_edge_open_file(argv[6], argv[7], argv[8], atoi(argv[9]));
    if(status == EXIT_ERROR)
	{
		fprintf(stderr, "Error: para_edge_open_file()\n");
        exit(EXIT_FAILURE);
	}
    
    // Refine node indicces
    status = vertex_map();
	if(status == EXIT_ERROR)
	{
		fprintf(stderr, "Error: vertex_map()\n");
        exit(EXIT_FAILURE);
	}
	
	printf("Sort the vertex file before calling open_vertex_file(%s)\n", argv[11]);
	
	printf("    sort -k 1 -n  %s\n", argv[11]);
	status = vertex_open_file(argv[11], argv[12]);
	if(status == EXIT_ERROR)
	{
		fprintf(stderr, "Error: open_file()\n");
        exit(EXIT_FAILURE);
	}
	
	graph->algorithm = malloc(sizeof(struct Algorithm)); 
    if(!graph->algorithm)
        exit(EXIT_FAILURE);
    memset(graph->algorithm, 0, sizeof(struct Algorithm));
    struct Algorithm *algorithm = graph->algorithm;
	algorithm->num_partitions = atoi(argv[2]);    
    memset(algorithm->name, 0, sizeof(algorithm->name));
    strcpy(algorithm->name, argv[4]);
    memset(algorithm->dir_name, 0, sizeof(algorithm->dir_name));
    strcpy(algorithm->dir_name, dirname(argv[0]));
    memset(algorithm->base_name, 0, sizeof(algorithm->base_name));
    strcpy(algorithm->base_name, basename(argv[6]));
    
    algorithm->l =     1; // lambda
    algorithm->c =  0.05; // load imbalance
	algorithm->o = 0.001; // omega
    algorithm->e =     5; // epsilon
    
    
    // Initialize load and capacity
    status = load_initialization();
    if(status == EXIT_ERROR)
	{
		fprintf(stderr, "Error: load_initialization()\n");
        return(EXIT_ERROR);
	}
    
    // Initialize temporary arrays
    status = temporary_array_initialization();
    if(status == EXIT_ERROR)
	{
		fprintf(stderr, "Error: temporary_array_initialization()\n");
        return(EXIT_ERROR);
	}
  
	float temp_score_array[algorithm->num_partitions];
	memset(temp_score_array, 0, sizeof(temp_score_array));
	for(i = 0; i < graph->num_nodes; i++)
    {
		j = graph->mapped_nodes[i];
		score_function(temp_score_array, algorithm->num_partitions, j);
		graph->nodes[j].score = temp_score_array[graph->nodes[j].label];
		score_array[0] += graph->nodes[j].score;
    }
    
	for(k = 0; k < algorithm->num_partitions; k++)
	{
		external_edges_array_per_partition[0][k] = enumerate_partition_edges(0, graph->num_nodes - 1, k, 0);
		internal_edges_array_per_partition[0][k] = enumerate_partition_edges(0, graph->num_nodes - 1, k, 1);
        ingoing_cut_array_per_partition[0][k] = enumerate_partition_cut(0, graph->num_nodes - 1, k, 1);
	}
	
	cut_array[0] = cut_size(0, graph->num_nodes - 1);
	local_edges_array[0] = phi(0, graph->num_nodes - 1);
    
    statistics(0);
  /*  
    status = write_partition_e(argv[6], comment, delimiter, columns);
    if(status == EXIT_ERROR)
	{
		fprintf(stderr, "Error: write_partition_e()\n");
        exit(EXIT_FAILURE);
	}
*/
	free_all();
    
	return(EXIT_OK);
}

