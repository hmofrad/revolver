/*
 * libgraph.c Graph processing library
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
unsigned int max_source_array[MAX_THREADS_NUM];
unsigned int max_destination_array[MAX_THREADS_NUM];
unsigned int max_array[MAX_THREADS_NUM];
unsigned int min_source_array[MAX_THREADS_NUM];
unsigned int min_destination_array[MAX_THREADS_NUM];
unsigned int min_array[MAX_THREADS_NUM];
long long num_nodes_array[MAX_THREADS_NUM];
unsigned int num_alloced_nodes_array[MAX_THREADS_NUM];
unsigned int num_edges_array[MAX_THREADS_NUM];
unsigned int num_alloced_edges_array[MAX_THREADS_NUM];
struct Adjacency_list *nodes_array[MAX_THREADS_NUM];
unsigned int *mapped_nodes_array[MAX_THREADS_NUM];
unsigned int thread_low_index_array[MAX_THREADS_NUM];
unsigned int thread_high_index_array[MAX_THREADS_NUM];
unsigned int probability_array[MAX_THREADS_NUM];
long long external_messages_array[MAX_THREADS_NUM];
long long  internal_messages_array[MAX_THREADS_NUM];

pthread_barrier_t barrier;
pthread_mutex_t mutex;
pthread_mutex_t mutex1;

static const char LOG_DIR[] = "../outputs";

char sentinel;
float score;

void multiprocessing_pool(unsigned int total_size, int num_chunks, int *k1, unsigned int *k1_num, int *k2, unsigned int *k2_num)
{
    // Distributing nodes accross threads
	*k1 = total_size % num_chunks; // mod
	*k1_num = 1 + ((total_size - 1) / num_chunks); //ceiling division
	*k2 = num_chunks - *k1; // diff
	*k2_num = (unsigned int)((double) (total_size/num_chunks)); //floor division
	
	#ifdef VERBOSE
        printf("\n");
        printf("Multiprocesing chunks (Starting from 0)\n");
        printf("   Number of chunks  %3d\n", num_chunks);
        printf("   Chunks x size %3d x %3u\n",  *k1, *k1_num);
        printf("   Chunks x size %3d x %3u\n",  *k2, *k2_num);
        printf("   Total size %u == Computed size %u\n", total_size, (unsigned int)((*k1 * *k1_num) + (*k2 * *k2_num)));
	#endif
}


void *read_chunk(void *arguments)
{
    struct pthread_args_struct *args = arguments;
	int index = args->index;
	unsigned int low = args->low;
	unsigned int high = args->high;
    FILE *file_descriptor = args->fd;
    char *format = args->fmt;
    int columns = args->cols;
    
    unsigned int source = 0;
	unsigned int max_source = 0;
    unsigned int min_source = 0;
    unsigned int destination = 0;
    unsigned int max_destination = 0;
    unsigned int min_destination = 0;
    long long num_edges = 0;
    long long num_alloced_edges = 0;
    long long num_alloced_nodes = 0;
	num_edges_array[index] = 0;
    num_alloced_edges_array[index] = 0;
    
    (void) source;
    (void) destination;
    (void) columns;
    (void) format;
    (void) file_descriptor;
    (void) max_source;
    (void) max_destination;
    (void) high;
    (void) low;
    (void) index;
    (void) num_edges;
    (void) min_source;
    (void) min_destination;
    (void) num_alloced_edges;
    (void) num_alloced_nodes;
    
    
 
    fseek(file_descriptor, low, SEEK_SET);
    long long position = ftell(file_descriptor);
    char epsilon = 3; // ' ' + '\n' + 1 extra charachter
    min_source_array[index] = MAX_UNISIGNED_INT;
    min_destination_array[index] = MAX_UNISIGNED_INT;
    while((fscanf(file_descriptor, format, &source, &destination) == columns))
    {
        num_edges_array[index]++;
        num_alloced_edges_array[index]++;
        if(source > max_source_array[index])
            max_source_array[index] = source;
        
        if(destination > max_destination_array[index])
            max_destination_array[index] = destination;
        
        if(source < min_source_array[index])
            min_source_array[index] = source;
        
        if(destination < min_destination_array[index])
            min_destination_array[index] = destination;
        
        position = ftell(file_descriptor);
        if (position + epsilon > high)
        {
            //printf("[%d] %u %u\n", index, high, position);
            break;
        }   
    }
    
    if(min_source_array[index] < min_destination_array[index])
        min_array[index] = min_source_array[index];
    else
        min_array[index] = min_destination_array[index];
    
    if(max_source_array[index] > max_destination_array[index])
        max_array[index] = max_source_array[index];
    else
        max_array[index] = max_destination_array[index];
    num_alloced_nodes_array[index] = max_array[index] - min_array[index] + 1; // Number of nodes for this thread

    nodes_array[index] = malloc(num_alloced_nodes_array[index] * sizeof(struct Adjacency_list));
    if(!nodes_array[index])
        exit(EXIT_FAILURE);
    memset(nodes_array[index], 0, num_alloced_nodes_array[index] * sizeof(struct Adjacency_list));    

    num_nodes_array[index] = 0;    

    pthread_barrier_wait(&barrier);
    pthread_mutex_lock(&mutex);
    if(!sentinel)
    {
        sentinel = 1;
        #if GCC_VERSION > GCC_GENERIC_FEATURE_SUPPORT
            max_source = ffmax(max_source_array, MAX_THREADS_NUM);
            max_destination = ffmax(max_destination_array, MAX_THREADS_NUM);
            num_edges = ffsum(num_edges_array, MAX_THREADS_NUM);
            //num_alloced_edges = ffsum(num_alloced_edges_array, MAX_THREADS_NUM);
            graph->num_edges = ffsum(num_edges_array, MAX_THREADS_NUM);
            //graph->num_alloced_edges = ffsum(num_alloced_edges_array, MAX_THREADS_NUM);
        #else
            max_source = max_value_uint(max_source_array, MAX_THREADS_NUM);
            max_destination = max_value_uint(max_destination_array, MAX_THREADS_NUM);
            num_edges = sum_uint(num_edges_array, MAX_THREADS_NUM);
            //num_alloced_edges = sum_uint(num_alloced_edges_array, MAX_THREADS_NUM);
            graph->num_edges = sum_uint(num_edges_array, MAX_THREADS_NUM);
            //graph->num_alloced_edges = sum_uint(num_alloced_edges_array, MAX_THREADS_NUM);
        #endif
        if(max_source > max_destination)
        {
            num_alloced_nodes = max_source;
            graph->num_alloced_nodes = max_source;
        }
        else
        {
            num_alloced_nodes = max_destination;
            graph->num_alloced_nodes = max_destination;
        }
        num_alloced_nodes++; // Take account of 0
        graph->num_alloced_nodes++;  // Take account of 0
        
        //printf("%u %u %u %u %u\n", max_source, max_destination, num_edges, num_alloced_edges, num_alloced_nodes);
        //printf("%u %u %u\n", graph->num_edges, graph->num_alloced_edges, graph->num_alloced_nodes);
        
        graph->nodes = malloc(graph->num_alloced_nodes * sizeof(struct Adjacency_list));
        if(!graph->nodes)
            exit(EXIT_FAILURE);
        memset(graph->nodes, 0, graph->num_alloced_nodes * sizeof(struct Adjacency_list));
        
    }
    pthread_mutex_unlock(&mutex);
    
    mapped_nodes_array[index] = malloc(graph->num_alloced_nodes * sizeof(unsigned int));	
    if(!mapped_nodes_array[index])
        exit(EXIT_ERROR);
    memset(mapped_nodes_array[index], 0, graph->num_alloced_nodes * sizeof(unsigned int));
   // printf("%d %u %u %u %u %u %u %u\n", index, min_source_array[index], min_destination_array[index], max_source_array[index], max_destination_array[index], min_array[index], max_array[index], num_alloced_nodes_array[index]);
    
    unsigned int i = 0;
    for(i = 0; i < num_alloced_nodes_array[index]; i++)
        mapped_nodes_array[index][min_array[index] + i] =  i;
    
    int status = 0;
    fseek(file_descriptor, low, SEEK_SET);
    position = ftell(file_descriptor);
    
    while(fscanf(file_descriptor, format, &source, &destination) == columns)
    {
        status = para_add_edge(&(nodes_array[index]), &(mapped_nodes_array[index]), &(num_nodes_array[index]), source, destination, 0);
        if(status == EXIT_ERROR)
        {
            fprintf(stderr, "ERROR: add_edge(): Adding the new edge failed\n");
            return(NULL);
        }
        
        // We treat all graphs as undirected graphs for graph operations 
        // and directed for graph statistics
        status = para_add_edge(&(nodes_array[index]), &(mapped_nodes_array[index]), &(num_nodes_array[index]), destination, source, 1);
        if(status == EXIT_ERROR)
        {
            fprintf(stderr, "ERROR: add_edge(): Adding the new edge failed\n");
            return(NULL);
        }
        num_alloced_edges_array[index]++;
        position = ftell(file_descriptor);
        if (position + epsilon > high)
        {
            //printf("[%d] %u %u\n", index, high, position);
            break;
        }
    }
    
    pthread_barrier_wait(&barrier);    
    pthread_mutex_lock(&mutex);
    if(sentinel == 1)
    {
        sentinel = 0;
        #if GCC_VERSION > GCC_GENERIC_FEATURE_SUPPORT
            graph->num_alloced_edges = ffsum(num_alloced_edges_array, MAX_THREADS_NUM);
        #else
            graph->num_alloced_edges = sum_uint(num_alloced_edges_array, MAX_THREADS_NUM);
        #endif
        
    }        
    pthread_mutex_unlock(&mutex);
    
    unsigned int j = 0;;
    unsigned int k = 0;
    unsigned int l = 0;
    char node_not_found = 0;
    struct Adjacency_list_node *src_node = NULL;
    struct Adjacency_list_node *new_node = NULL;
    struct Adjacency_list_node *tmp_node = NULL;
    struct Adjacency_list_node *old_node = NULL;
    
    
    graph->nodes_list = malloc(graph->num_alloced_nodes * sizeof(struct List_of_Adjacency_list));
    if(!graph->nodes_list)
            exit(EXIT_FAILURE);
    memset(graph->nodes_list, 0, graph->num_alloced_nodes * sizeof(struct List_of_Adjacency_list));
    struct List_of_Adjacency_list_node *src_list = NULL;
    struct List_of_Adjacency_list_node *old_list = NULL;
    (void) old_list;
    (void) old_node;
    (void) tmp_node;
    (void) src_node;
    (void) node_not_found;
    (void) l;
    (void) new_node;
    
    pthread_barrier_wait(&barrier);
    pthread_mutex_lock(&mutex);
    if(sentinel == 0)
    {
        sentinel = 1;
        for(i = 0; i < graph->num_alloced_nodes; i++)
        {
            for(j = 0; j < MAX_THREADS_NUM; j++)
            {
                k = mapped_nodes_array[j][i];
                if((i >= min_array[j]) && (i <= max_array[j]) && nodes_array[j][k].allocated)
                {
                    if(!graph->nodes_list[i].allocated)
                    {
                        graph->nodes_list[i].head = malloc(sizeof(struct List_of_Adjacency_list_node));
                        if(!graph->nodes_list[i].head)
                            exit(EXIT_FAILURE);
                        memset(graph->nodes_list[i].head, 0, sizeof(struct List_of_Adjacency_list_node));
                        graph->nodes_list[i].allocated = 1;
                        src_list = graph->nodes_list[i].head;
                        
                        src_list->nodes = malloc(sizeof(struct Adjacency_list));
                        memcpy(src_list->nodes, &nodes_array[j][k], sizeof(struct Adjacency_list));
                        src_list->next = NULL;
                        graph->num_nodes++;
                    }
                    else
                    {
                        src_list = graph->nodes_list[i].head;
                        while(src_list->next)
                        {
                            src_list = src_list->next;
                        }
                        
                        src_list->next = malloc(sizeof(struct List_of_Adjacency_list_node));
                        if(!src_list->next)
                            exit(EXIT_FAILURE);
                        memset(src_list->next, 0, sizeof(struct List_of_Adjacency_list_node));
                        src_list = src_list->next;
                        
                        src_list->nodes = malloc(sizeof(struct Adjacency_list));
                        memcpy(src_list->nodes, &nodes_array[j][k], sizeof(struct Adjacency_list));
                        
                        src_list->next = NULL;
                    }
                    
                }
            }
        }
        
        
        int k1 = 0;
        unsigned int k1_num = 0;
        int k2 = 0;
        unsigned int k2_num = 0;
        multiprocessing_pool(graph->num_alloced_nodes, MAX_THREADS_NUM, &k1, &k1_num, &k2, &k2_num);
        for(i = 0; i < MAX_THREADS_NUM; i++)
        {
            if(i < k1)
            {
                low = i * k1_num;
                high = ((i + 1)* k1_num) - 1;
            }
            else if((i >= k1) && (i < k1 + k2))
            {
                low = (k1 * k1_num) + ((i - k1) * k2_num);
                high = (k1 * k1_num) + (((i + 1) - k1) * k2_num) - 1;
            }
            thread_low_index_array[i] = low;
            thread_high_index_array[i] = high;
        }
    }
    pthread_mutex_unlock(&mutex);
    free(nodes_array[index]);

    low = thread_low_index_array[index];
    high = thread_high_index_array[index];
    for(i = low; i <= high; i++)
    {
         if(graph->nodes_list[i].allocated)
         {
            src_list = graph->nodes_list[i].head;
            while(src_list)
            {      
                if(src_list->nodes->source != i)
                {
                    printf("++++++++++++++ %u %d %d %u %p\n", src_list->nodes->source, src_list->nodes->degree, src_list->nodes->idegree, i, src_list->nodes);
                    exit(0);
                }
                if(!graph->nodes[i].allocated)
                {
                    graph->nodes[i].allocated = 1;
                    graph->nodes[i].source = i;
                }
                graph->nodes[i].degree += src_list->nodes->degree;
                graph->nodes[i].idegree += src_list->nodes->idegree;                
                new_node = src_list->nodes->head;
                
                while(new_node)
                {
                    node_not_found = 0;
                    l = new_node->destination;   
                    if(graph->nodes[i].head)
                    {
                        src_node = graph->nodes[i].head; 
                        while(src_node)
                        {          
                            if(l == src_node->destination)
                            {
                                node_not_found = 1;
                                src_node->count  += new_node->count;
                                src_node->icount += new_node->icount;
                            }
                            old_node = src_node;
                            src_node = src_node->next; 
                        }
                    }
                    if(!node_not_found)
                    {
                        tmp_node = malloc_edge(l);
                        tmp_node->count  = new_node->count;
                        tmp_node->icount = new_node->icount;
                        tmp_node->next = NULL;
                        if(!graph->nodes[i].head)
                            graph->nodes[i].head = tmp_node;
                        else
                            old_node->next = tmp_node;
                    }
                    new_node = new_node->next;
                }
                
                // clean
                src_node = src_list->nodes->head;
                while(src_node)
                {
                    old_node = src_node->next;
                    free(src_node);
                    src_node = old_node;
                }
                
                old_list = src_list->next;
                free(src_list);
                src_list = old_list;
            }
         }
    }
        
    if(!fclose(file_descriptor))
         ; // NoOp, we're good to go!
    else
    {
        fprintf(stderr, "Error on closing the input file\n");
        exit(EXIT_FAILURE);
    }

    free(args);
    return(NULL);
    
    
}


long long refine_position(FILE *file_descriptor, long long position) {
    char ch = 0;
    fseek(file_descriptor, position, SEEK_SET);
    do
    {
        ch = fgetc(file_descriptor);
    } while(ch != '\n');
    position = ftell(file_descriptor);
    return(position);
}


int para_edge_open_file(char *file_name, char *comment_style, char *delimiter, int columns)
{	
	FILE *file_descriptor = fopen(file_name,"r");
    if(!file_descriptor)
    {
        fprintf(stderr, "Error on opening the input file\n");
        return(EXIT_ERROR);
    }
	
    char  *format = edge_format_specifier(file_descriptor, comment_style, delimiter);
	if(!format)
	{
		fprintf(stderr, "Error: edge_format_specifier()\n");
		return(EXIT_ERROR);
	}
	int offset = ftell(file_descriptor);
    printf("offset: %d\n", offset);
    
    fseek(file_descriptor, 0L, SEEK_END);
    long long file_length = ftell(file_descriptor);
    rewind(file_descriptor);
    printf("File length %llu\n", file_length - offset);
    int k1 = 0;
    unsigned int k1_num = 0;
    int k2 = 0;
    unsigned int k2_num = 0;
    multiprocessing_pool(file_length, MAX_THREADS_NUM, &k1, &k1_num, &k2, &k2_num);
    
    int i = 0;
    unsigned int left_position = 0;
    unsigned int right_position = 0;;

    pthread_t tid[MAX_THREADS_NUM];
	memset(tid, 0, sizeof(tid));
	struct pthread_args_struct *args = NULL;
    int status = 0;
    sentinel = 0;

    FILE *file_descriptor_per_thread;
    pthread_barrier_init(&barrier, NULL, MAX_THREADS_NUM);
    for(i = 0; i < MAX_THREADS_NUM; i++)
    {	
        
        args = malloc(sizeof(struct pthread_args_struct) * MAX_THREADS_NUM);
        if(!args)
            return(EXIT_FAILURE);
        memset(args, 0, sizeof(struct pthread_args_struct));

        file_descriptor_per_thread =  fopen(file_name,"r");;
        if(!file_descriptor_per_thread)
        {
            fprintf(stderr, "Error on opening the input file\n");
            return(EXIT_ERROR);
        }

        args->index = i;
        args->fd = file_descriptor_per_thread;
        args->fmt = format;
        args->cols = columns;
        
        if(i < k1)
        {
            if(!left_position && !right_position)
            {
                left_position = offset;
            }
            else 
            {
                left_position = right_position;
            }

            right_position = ((i + 1)* k1_num) - 1; 

            right_position = refine_position(file_descriptor, right_position);

            args->low = left_position;
            args->high = right_position;
        }
        else if((i >= k1) && (i < k1 + k2))
        {
            if(!left_position && !right_position)
            {                
                left_position = offset;
            }
            else 
            {
                left_position = right_position;
            }
           
           
            right_position = ((k1 * k1_num) + (((i + 1) - k1) * k2_num) - 1);
            right_position = refine_position(file_descriptor, right_position);
            args->low = left_position;
            args->high = right_position;
        }
        
        status = pthread_create(&tid[i], NULL, (void*)read_chunk, (void *)args);
            
        if (status != EXIT_SUCCESS)
        {
            fprintf(stderr, "Error: pthread_create()\n");
            exit(EXIT_FAILURE);
        }
        //printf("[%d] %u %u\n", i, left_position, right_position);
    }
    for(i = 0; i < MAX_THREADS_NUM; i++)
    {
        pthread_join(tid[i], NULL);
    }
	
    #ifdef VERBOSE
	    printf("#num_alloced_nodes: %lli\n", graph->num_alloced_nodes);
	    printf("#num_nodes:         %lli\n", graph->num_nodes);
		printf("#num_edges:         %lli\n", graph->num_edges);
        printf("#num_alloced_edges: %lli\n", graph->num_alloced_edges);
	#endif
    
	if(!fclose(file_descriptor))
         ; // NoOp, we're good to go!
    else
    {
        fprintf(stderr, "Error on closing the input file\n");
        return(EXIT_ERROR);
    }
    
    free(format);
	return(EXIT_SUCCESS);
}

int edge_open_file(char *file_name, char *comment_style, char *delimiter, int columns)
{	
	FILE *file_descriptor = fopen(file_name,"r");
    if(!file_descriptor)
    {
        fprintf(stderr, "Error on opening the input file\n");
        return(EXIT_ERROR);
    }
	
    char  *format = edge_format_specifier(file_descriptor, comment_style, delimiter);
	if(!format)
	{
		fprintf(stderr, "Error: edge_format_specifier()\n");
		return(EXIT_ERROR);
	}
	int offset = ftell(file_descriptor);
	
	int status = 0;
    unsigned int source = 0;
	unsigned int max_source = 0;
    unsigned int destination = 0;
    unsigned int max_destination = 0;
	
	while(fscanf(file_descriptor, format, &source, &destination) == columns)
    {
        graph->num_edges++;
        graph->num_alloced_edges++;
        if(source > max_source)
            max_source = source;
        if(destination > max_destination)
            max_destination = destination;
    }

	if(max_source > max_destination)
        graph->num_alloced_nodes = max_source;
    else
        graph->num_alloced_nodes = max_destination;
    graph->num_alloced_nodes++; // Take account of 0
	 
	fseek(file_descriptor, offset, SEEK_SET);
					
	graph->nodes = malloc(graph->num_alloced_nodes * sizeof(struct Adjacency_list));
    if(!graph->nodes)
        exit(EXIT_FAILURE);
    memset(graph->nodes, 0, graph->num_alloced_nodes * sizeof(struct Adjacency_list));
    
    while(fscanf(file_descriptor, format, &source, &destination) == columns)
    {
        status = add_edge(&(graph->nodes), &(graph->num_nodes), source, destination, 0);
        if(status == EXIT_ERROR)
        {
            fprintf(stderr, "ERROR: add_edge(): Adding the new edge failed\n");
            return(EXIT_ERROR);
        }
        
        // We treat all graphs as undirected graphs for graph operations 
        // and directed for graph statistics
        status = add_edge(&(graph->nodes), &(graph->num_nodes), destination, source, 1);
        if(status == EXIT_ERROR)
        {
            fprintf(stderr, "ERROR: add_edge(): Adding the new edge failed\n");
            return(EXIT_ERROR);
        }
        graph->num_alloced_edges++;
        //printf("(%u %u) %lli\n", source, destination, graph->num_alloced_edges);
    }
	
	#ifdef VERBOSE
	    printf("#num_alloced_nodes: %lli\n", graph->num_alloced_nodes);
	    printf("#num_nodes:         %lli\n", graph->num_nodes);
		printf("#num_edges:         %lli\n", graph->num_edges);
        printf("#num_alloced_edges: %lli\n", graph->num_alloced_edges);
	#endif
	free(format);
	
	if(!fclose(file_descriptor))
         ; // NoOp, we're good to go!
    else
    {
        fprintf(stderr, "Error on closing the input file\n");
        return(EXIT_ERROR);
    }
	
	return(EXIT_SUCCESS);
}

int edge_open_file_mmapped(char *file_name, char *comment_style, char *delimiter, int columns)
{	
    int file_descriptor = open(file_name, O_RDONLY);
    if(file_descriptor == -1)
    {
        fprintf(stderr, "Error on opening the input file %s\n", file_name);
        return(EXIT_ERROR);
    }
    long long file_length = lseek(file_descriptor, 0L, SEEK_END);
    //printf("%lli\n", file_length);
    lseek(file_descriptor, 0L, SEEK_SET);
    //printf("%lli\n", file_length);
    printf("mmap\n");
    char *file = mmap(0, (size_t) file_length, PROT_READ, MAP_SHARED, file_descriptor, 0);

    //printf("%c %c\n", file[0], file[1]);
    if((void *) file == MAP_FAILED)
    {
        fprintf(stderr, "Error on mapping memory\n");
        close(file_descriptor);
        return(EXIT_ERROR);
    }
    
    FILE *file_descriptor_ = fopen(file_name,"r");
    if(!file_descriptor_)
    {
        fprintf(stderr, "Error on opening the input file\n");
        return(EXIT_ERROR);
    }
    
    char *format = edge_format_specifier(file_descriptor_, comment_style, delimiter);
	if(!format)
	{
		fprintf(stderr, "Error: edge_format_specifier()\n");
		return(EXIT_ERROR);
	}
	int offset = ftell(file_descriptor_);
    char *file_pointer = file + (char) offset;
    
 	if(!fclose(file_descriptor_))
         ; // NoOp, we're good to go!
    else
    {
        fprintf(stderr, "Error on closing the input file\n");
        return(EXIT_ERROR);
    }   
    
    printf("max_source\n");
    int status = 0;
    unsigned int source = 0;
	unsigned int max_source = 0;
    unsigned int destination = 0;
    unsigned int max_destination = 0;
    printf("%s\n", format);
	while(sscanf(file_pointer, format, &source, &destination) == columns)
    {
        graph->num_edges++;
        graph->num_alloced_edges++;
        if(source > max_source)
            max_source = source;
        if(destination > max_destination)
            max_destination = destination;
        printf("%u %u\n", source, destination);

        do
        {
            file_pointer++;
        }
        while(*file_pointer != '\n');
        file_pointer++;
    }

	if(max_source > max_destination)
        graph->num_alloced_nodes = max_source;
    else
        graph->num_alloced_nodes = max_destination;
    graph->num_alloced_nodes++; // Take account of 0
	 
    file_pointer = file + (char) offset;
    
    printf("reading ...\n");
    graph->nodes = malloc(graph->num_alloced_nodes * sizeof(struct Adjacency_list));
    if(!graph->nodes)
        exit(EXIT_FAILURE);
    memset(graph->nodes, 0, graph->num_alloced_nodes * sizeof(struct Adjacency_list));
	
	while(sscanf(file_pointer, format, &source, &destination) == columns)
    {
        status = add_edge(&(graph->nodes), &(graph->num_nodes), source, destination, 0);
        if(status == EXIT_ERROR)
        {
            fprintf(stderr, "ERROR: add_edge(): Adding the new edge failed\n");
            return(EXIT_ERROR);
        }
        
        // We treat all graphs as undirected graphs for graph operations 
        // and directed for graph statistics
        status = add_edge(&(graph->nodes), &(graph->num_nodes), destination, source, 1);
        if(status == EXIT_ERROR)
        {
            fprintf(stderr, "ERROR: add_edge(): Adding the new edge failed\n");
            return(EXIT_ERROR);
        }
        graph->num_alloced_edges++;
        
        do
        {
            file_pointer++;
        }
        while(*file_pointer != '\n');
        file_pointer++;
    }
	
	#ifdef VERBOSE
	    printf("#num_alloced_nodes: %lli\n", graph->num_alloced_nodes);
	    printf("#num_nodes:         %lli\n", graph->num_nodes);
		printf("#num_edges:         %lli\n", graph->num_edges);
        printf("#num_alloced_edges: %lli\n", graph->num_alloced_edges);
	#endif
	
    free(format);
    if(munmap(file, file_length) == -1) {
        fprintf(stderr, "Error unmapping memory\n");
        return(EXIT_ERROR);
    }

	return(EXIT_SUCCESS);
}


char *edge_format_specifier(FILE* file_descriptor, char *comment_style, char *delimiter)
{
    char comment[2];
    memset(comment, '\0',sizeof(comment));

    char * line = NULL;
    size_t len = 0;
    ssize_t read;

    do
    {
        read = getline(&line, &len, file_descriptor);
        if(read == -1)
        { 
            fprintf(stderr, "Error: getline()\n");
            return(NULL);
        }
        strncpy(comment, &line[0], 1);
    }
    while(!strcmp(comment,comment_style));
    fseek(file_descriptor, -read, SEEK_CUR);
	
    char *buffer = malloc(strlen(line)+1);
	memset(buffer, '\0', strlen(line)+1);
	memcpy(buffer, line, strlen(line));
	free(line);
	
    int count = 0;    
    char *token = NULL;
	char tab[3];
	memset(tab, '\0', sizeof(tab));
	tab[0] = '\\';
	tab[1] = 't';
	tab[2] = '\0';
	
	char ctab = '\t';
	if((delimiter[0] == tab[0]) && (delimiter[1] == tab[1]))
		strcpy(delimiter, &ctab);
	
    token = strtok(buffer, delimiter);
	
    while(token)
    {
        token = strtok(NULL, delimiter);
        count++;
    }
	
    free(buffer);
	
	if(count == 1)
	{
		fprintf(stderr, "Error: strtok(line, delimiter = \"%s\"), use $\"%s\" instead\n", delimiter, delimiter);
		return(NULL);
	}
    int format_len =  10;
    //int format_len =  16;
    char *format = malloc(format_len * sizeof(char));
	if(!format)
		return(NULL);
    memset(format, '\0', format_len * sizeof(char));
    switch(count)
    {
        case 2:
            strncpy(format, "%u %u", 5);
		    //strncpy(format, "%d %d", 5);
            //strncpy(format, "%lli %lli", 9);
            break;    
        case 3:
            strncpy(format, "%u %u %*u", 9);
            //strncpy(format, "%d %d %*d", 9);
            //strncpy(format, "%lli %lli %*lli", 15);
            break;    
    }
	return(format);
}

int vertex_open_file(char *file_name, char *delimiter)
{
	struct Adjacency_list *nodes = graph->nodes;
	
	FILE *file_descriptor = fopen(file_name,"r");
    if(!file_descriptor)
    {
        fprintf(stderr, "Error on opening the input file\n");
        return(EXIT_ERROR);
    }
	
    char *format = vertex_format_specifier(delimiter);
	if(!format)
	{
		fprintf(stderr, "Error: vertex_format_specifier()\n");
		return(EXIT_ERROR);
	}
    
    unsigned int i = 0;
	unsigned int j = 0;
    
	unsigned int node = 0;
    int label = 0;
	int columns = 2;
    (void) i ; 

	while(fscanf(file_descriptor, format, &node, &label) == columns)
    {

        if((node > 0) || (node < graph->num_alloced_nodes))
        {
            nodes[node].label = label;
        } else {
            fprintf(stderr, "Error on reading the input file\n");
	    fprintf(stderr, "Tuple (Index=%u Node=%u)\n", j, node);
	    printf("Make sure %s is sorted\n", file_name);
	    exit(EXIT_ERROR);
        }


    /*
		j = graph->mapped_nodes[i];
		i++;
		if(j == node)
			nodes[j].label = label;
		else
		{
		    fprintf(stderr, "Error on reading the input file\n");
			fprintf(stderr, "Tuple (Index=%u Node=%u)\n", j, node);
			printf("Make sure %s is sorted\n", file_name);
			//exit(EXIT_ERROR);
                        i--;
		}
    */
    }
	free(format);
	
	if(!fclose(file_descriptor))
         ; // NoOp, we're good to go!
    else
    {
        fprintf(stderr, "Error on closing the input file\n");
        return(EXIT_ERROR);
    }
	
    return(EXIT_OK);
}

char *vertex_format_specifier(char *delimiter)
{
	int format_len =  6;
    char *format = malloc(format_len * sizeof(char));
	if(!format)
		return((char *) EXIT_ERROR);
    memset(format, '\0', format_len * sizeof(char));
	
	memcpy(format, "%d", 2);
	memcpy(format+2, delimiter, 1);
	memcpy(format+3, "%d", 2);
	
	return(format);
}

int add_edge(struct Adjacency_list **nodes, long long *num_nodes, unsigned int source, unsigned int destination, char edge_type)
{
    struct Adjacency_list_node *temp_node = NULL;
    struct Adjacency_list_node *new_node = malloc_edge(destination);

    if(!new_node)
        goto error;
    if(!(*nodes)[source].head)
    {
		if(!(*nodes)[source].allocated)
		{
		    (*num_nodes)++;
		    (*nodes)[source].allocated = 1;
            (*nodes)[source].source = source;
            (*nodes)[source].head = new_node;
            new_node->count = 1;
            if(!edge_type)
            {
                (*nodes)[source].degree = 1;
                new_node->count = 1;
            }
            else
            {
                (*nodes)[source].idegree = 1;
                new_node->icount = 1;
            }
            new_node->next = NULL;
		}
		else
		{
            (*nodes)[source].head = new_node;
            if(!edge_type)
            {
                (*nodes)[source].degree = 1;
                new_node->count = 1;
            }
            else
            {
                (*nodes)[source].idegree = 1;
                new_node->icount = 1;
            }
            new_node->next = NULL;
		}
    }
    else
    {
        temp_node = (*nodes)[source].head;
        while(temp_node)
        {
            if(temp_node->destination == destination)
            {

                if(!edge_type)
                {
                    (*nodes)[source].degree++;
                    temp_node->count++;
                }
                else
                {
                    (*nodes)[source].idegree++;
                    temp_node->icount++;
                }
                free(new_node);
                break;    
            }
            else if(!temp_node->next)
            {
                temp_node->next = new_node;
                if(!edge_type)
                {
                    (*nodes)[source].degree++;
                    new_node->count = 1;
                }
                else
                {
                    (*nodes)[source].idegree++;
                    new_node->icount = 1;
                }
                new_node->next = NULL;
                break;
            }        
            else
                temp_node = temp_node->next;
        }
    }
	
    if(!(*nodes)[destination].head && !(*nodes)[destination].allocated)
    {
		(*num_nodes)++;
		(*nodes)[destination].allocated = 1;
        (*nodes)[destination].source = destination;        
    }

    return(EXIT_OK);

error:
    (*num_nodes)--;
    return(EXIT_ERROR);
}


int para_add_edge(struct Adjacency_list **nodes, unsigned int **mapped_nodes, long long *num_nodes, unsigned int source, unsigned int destination, char edge_type)
{

    struct Adjacency_list_node *temp_node = NULL;
    struct Adjacency_list_node *new_node = malloc_edge(destination);

    unsigned int mapped_source = (*mapped_nodes)[source];
    unsigned int mapped_destination = (*mapped_nodes)[destination];
    
    if(!new_node)
        goto error;
    if(!(*nodes)[mapped_source].head)
    {
		if(!(*nodes)[mapped_source].allocated)
		{
		    (*num_nodes)++;
		    (*nodes)[mapped_source].allocated = 1;
            (*nodes)[mapped_source].source = source;
            (*nodes)[mapped_source].head = new_node;
            new_node->count = 1;
            if(!edge_type)
            {
                (*nodes)[mapped_source].degree = 1;
                new_node->count = 1;
            }
            else
            {
                (*nodes)[mapped_source].idegree = 1;
                new_node->icount = 1;
            }
            new_node->next = NULL;
		}
		else
		{
            (*nodes)[mapped_source].head = new_node;
            if(!edge_type)
            {
                (*nodes)[mapped_source].degree = 1;
                new_node->count = 1;
            }
            else
            {
                (*nodes)[mapped_source].idegree = 1;
                new_node->icount = 1;
            }
            new_node->next = NULL;
		}
    }
    else
    {
        temp_node = (*nodes)[mapped_source].head;
        while(temp_node)
        {
            if(temp_node->destination == destination)
            {

                if(!edge_type)
                {
                    (*nodes)[mapped_source].degree++;
                    temp_node->count++;
                }
                else
                {
                    (*nodes)[mapped_source].idegree++;
                    temp_node->icount++;
                }
                free(new_node);
                break;    
            }
            else if(!temp_node->next)
            {
                temp_node->next = new_node;
                if(!edge_type)
                {
                    (*nodes)[mapped_source].degree++;
                    new_node->count = 1;
                }
                else
                {
                    (*nodes)[mapped_source].idegree++;
                    new_node->icount = 1;
                }
                new_node->next = NULL;
                break;
            }        
            else
                temp_node = temp_node->next;
        }
    }
	
    if(!(*nodes)[mapped_destination].head && !(*nodes)[mapped_destination].allocated)
    {
		(*num_nodes)++;
		(*nodes)[mapped_destination].allocated = 1;
        (*nodes)[mapped_destination].source = destination;        
    }

    return(EXIT_OK);

error:
    (*num_nodes)--;
    return(EXIT_ERROR);
}

struct Adjacency_list_node *malloc_edge(int destination)
{
    struct Adjacency_list_node *edge = malloc(sizeof(struct Adjacency_list_node));
    if(!edge)
    {
        fprintf(stderr, "Error: malloc_edge()\n");
        return(NULL);
    }
    memset(edge, '\0', sizeof(struct Adjacency_list_node));

    edge->destination = destination;
    edge->count = 0;   // Redundant
    edge->next = NULL; // Redundant
    return(edge);
}

int vertex_map()
{
   struct Adjacency_list *nodes = graph->nodes;
   
    unsigned int i = 0;
    unsigned int j = 0;

	graph->mapped_nodes = malloc(graph->num_nodes * sizeof(unsigned int));	
    if(!graph->mapped_nodes)
        return(EXIT_ERROR);
    memset(graph->mapped_nodes, 0, graph->num_nodes * sizeof(unsigned int));

    for(i = 0; i < graph->num_alloced_nodes; i++)
    { 
        if(nodes[i].allocated)
        {
            graph->mapped_nodes[j] = i;
            j++;
        }
    } 

    // Uncommect if you want to store partitioning results
    /*
    graph->mapped_nodes_ = malloc(graph->num_alloced_nodes * sizeof(unsigned int));	
    if(!graph->mapped_nodes_)
        return(EXIT_ERROR);
    memset(graph->mapped_nodes_, 0, graph->num_alloced_nodes * sizeof(unsigned int));
    for(i = 0; i < graph->num_alloced_nodes; i++)
    { 
        if(nodes[i].allocated)
        {            
            graph->mapped_nodes_[i] = j;
            j++;
        }
    } 
    */
    return(EXIT_OK);
}

int algorithm_initialization()
{
    int status;
    
    // Randomly assign labels for the first time
    status = label_initialization();
    if(status == EXIT_ERROR)
	{
		fprintf(stderr, "Error: label_initialization()\n");
        return(EXIT_ERROR);
	}
   
    // Initialize load and capacity
    status = load_initialization();
    if(status == EXIT_ERROR)
	{
		fprintf(stderr, "Error: load_initialization()\n");
        return(EXIT_ERROR);
	}
   
    // Initialize LA
    status = learning_automata_initialization();
    if(status == EXIT_ERROR)
	{
		fprintf(stderr, "Error: learning_automata_initialization()\n");
        return(EXIT_ERROR);
	}
   
    // Compute initial scores
    status = score_initialization();
    if(status == EXIT_ERROR)
	{
		fprintf(stderr, "Error: score_initialization()\n");
        return(EXIT_ERROR);
	}
    
    // Compute initial migration probability
    status = migration_initialization();
    if(status == EXIT_ERROR)
	{
		fprintf(stderr, "Error: migration_initialization()\n");
        return(EXIT_ERROR);
	}
    
    // Initialize temporary arrays
    status = temporary_array_initialization();
    if(status == EXIT_ERROR)
	{
		fprintf(stderr, "Error: temporary_array_initialization()\n");
        return(EXIT_ERROR);
	}
    

    
    return(EXIT_OK);
}



int label_initialization()
{
    struct Adjacency_list *nodes = graph->nodes;
    struct Algorithm *algorithm = graph->algorithm;
    
    unsigned int i = 0;
    unsigned int j = 0;
    for(i = 0; i < graph->num_nodes; i++)
    {
		j = graph->mapped_nodes[i];
		#if GCC_VERSION > GCC_GENERIC_FEATURE_SUPPORT
            nodes[j].label = ffrandom(algorithm->num_partitions - 1);
		#else
            nodes[j].label = random_at_most(algorithm->num_partitions - 1);
        #endif
    }
    return(EXIT_OK);
}

int load_initialization()
{
    struct Algorithm *algorithm = graph->algorithm;
    
    int k = 0;
    
    algorithm->load = malloc(algorithm->num_partitions * sizeof(unsigned int));
    if(!algorithm->load)
	{
        return(EXIT_ERROR);
	}
	memset(algorithm->load, 0, algorithm->num_partitions * sizeof(unsigned int));
	enumerate_edges();
    
    algorithm->load_weight = malloc(algorithm->num_partitions * sizeof(float));
    if(!algorithm->load_weight)
	{
        return(EXIT_ERROR);
	}
	memset(algorithm->load_weight, 0, algorithm->num_partitions * sizeof(float));
    
    algorithm->capacity = malloc(algorithm->num_partitions * sizeof(unsigned int));
    if(!algorithm->capacity)
	{
        return(EXIT_ERROR);
	}
	memset(algorithm->capacity, 0, algorithm->num_partitions * sizeof(unsigned int));
    //int lambda = 32;
    for(k = 0; k < algorithm->num_partitions; k++)
    {
        //algorithm->load_weight[k] = 0.01 + 0.09 * (exp(((float) (lambda * k)/(algorithm->num_partitions - 1)) - lambda));
      //  algorithm->load_weight[k] = 0.01 + 0.09 * (exp((float) (lambda * k)/(algorithm->num_partitions - 1)) * exp(-lambda));
        
        algorithm->load_weight[k] = algorithm->c;
        algorithm->capacity[k] = (unsigned int) ((algorithm->l + algorithm->load_weight[k]) * graph->num_edges) / algorithm->num_partitions;
    }
    /*
    for(k = 0; k < algorithm->num_partitions; k++)
    {
        printf("load [%d] %f %u\n", k, algorithm->l + algorithm->load_weight[k], algorithm->capacity[k]);
    }
    printf("load = %f %u\n", algorithm->l + algorithm->c, (unsigned int) ((algorithm->l + algorithm->c) * graph->num_edges)/algorithm->num_partitions);
    */
    algorithm->label = malloc(algorithm->num_partitions * sizeof(unsigned int));
    if(!algorithm->label)
	{
        return(EXIT_ERROR);
	}
	memset(algorithm->label, 0, algorithm->num_partitions * sizeof(unsigned int));
	enumerate_nodes();
    
    return(EXIT_OK);
}

int score_initialization()
{
    compute_score(0, graph->num_nodes - 1);
    return(EXIT_OK);
}

int migration_initialization()
{
    struct Algorithm *algorithm = graph->algorithm;
    
    algorithm->candidate_vertices = malloc(algorithm->num_partitions * sizeof(unsigned int));
    if(!algorithm->candidate_vertices)
	{
        return(EXIT_ERROR);
	}
	memset(algorithm->candidate_vertices, 0, algorithm->num_partitions * sizeof(unsigned int));
    enumerate_candidate_vertices();
    
	algorithm->remaining_capacity = malloc(algorithm->num_partitions * sizeof(signed int));
    if(!algorithm->remaining_capacity)
	{
        return(EXIT_ERROR);
	}
	memset(algorithm->remaining_capacity, 0, algorithm->num_partitions * sizeof(signed int));
	calculate_remaining_capacity();

	algorithm->migration_probability = malloc(algorithm->num_partitions * sizeof(float));
    if(!algorithm->migration_probability)
	{
        return(EXIT_ERROR);
	}
	memset(algorithm->migration_probability, 0, algorithm->num_partitions * sizeof(float));
    calculate_first_migration_probability();

    return(EXIT_OK);
}

int temporary_array_initialization()
{
    struct Algorithm *algorithm = graph->algorithm;
    
    int i = 0;
    
    for(i = 0; i < MAX_THREADS_NUM; i++)
	{
		external_edges_array_per_partition[i] = malloc(sizeof(unsigned int) * (algorithm->num_partitions));
        if(!external_edges_array_per_partition[i])
        {
            return(EXIT_ERROR);
        }
		memset(external_edges_array_per_partition[i], 0, sizeof(unsigned int) * (algorithm->num_partitions));
		
		internal_edges_array_per_partition[i] = malloc(sizeof(unsigned int) * (algorithm->num_partitions));
        if(!internal_edges_array_per_partition[i])
        {
            return(EXIT_ERROR);
        }
		memset(internal_edges_array_per_partition[i], 0, sizeof(unsigned int) * (algorithm->num_partitions));
       
        ingoing_cut_array_per_partition[i] = malloc(sizeof(unsigned int) * (algorithm->num_partitions));
        if(!ingoing_cut_array_per_partition[i])
        {
            return(EXIT_ERROR);
        }
		memset(ingoing_cut_array_per_partition[i], 0, sizeof(unsigned int) * (algorithm->num_partitions));
        
		candidate_vertices_array[i] = malloc(sizeof(unsigned int) * (algorithm->num_partitions));
        if(!candidate_vertices_array[i])
        {
            return(EXIT_ERROR);
        }
		memset(candidate_vertices_array[i], 0, sizeof(unsigned int) * (algorithm->num_partitions));
        
	}
    return(EXIT_OK);
}

int learning_automata_initialization()
{
    struct Adjacency_list *nodes = graph->nodes;
    struct Algorithm *algorithm = graph->algorithm;
    
    unsigned int i = 0;
    unsigned int j = 0;
    int k = 0;
    
    float temp_score_array[algorithm->num_partitions];
	memset(temp_score_array, 0, sizeof(temp_score_array));

    for(i = 0; i < graph->num_nodes; i++)
	{
		j = graph->mapped_nodes[i];
		graph->nodes[j].probability = malloc(algorithm->num_partitions * sizeof(float));
        if(!graph->nodes[j].probability)
        {
            return(EXIT_ERROR);
        }  
		memset(nodes[j].probability, 0, algorithm->num_partitions * sizeof(float));
        
        //graph->nodes[j].signal_ = NULL;
        graph->nodes[j].signal_ = malloc(algorithm->num_partitions * sizeof(float));
        if(!graph->nodes[j].signal_)
        {
            return(EXIT_ERROR);
        }  
		memset(graph->nodes[j].signal_, 0, algorithm->num_partitions * sizeof(float));
        
        
        score_function(temp_score_array, algorithm->num_partitions, j);
        
        //memcpy(nodes[j].probability, temp_score_array, sizeof(temp_score_array));
        
        for(k = 0; k < algorithm->num_partitions; k++)
            graph->nodes[j].probability[k] = (float) 1/algorithm->num_partitions;

        nodes[j].action = trimmer_tool(nodes[j].probability, algorithm->num_partitions);
        nodes[j].label = nodes[j].action;
        nodes[j].score = temp_score_array[nodes[j].label];
        nodes[j].signal = 1;
	}    
    enumerate_nodes();
    enumerate_edges();
    return(EXIT_OK);
}

int computation_is_halted(float score, float score_old, int *ep)
{
    struct Algorithm *algorithm = graph->algorithm;
    
    int halt = 1;
    
    if(score - score_old < algorithm->o)
    {
        *ep = *ep - 1;
        if(*ep < 0)
            halt = 0;
    }
    else
        *ep = algorithm->e;
    return(halt);
    
}

void enumerate_nodes()
{
    struct Adjacency_list *nodes = graph->nodes;
    struct Algorithm *algorithm = graph->algorithm;
    
    unsigned int i = 0;
	unsigned int j = 0;
    
	memset(algorithm->label, 0, algorithm->num_partitions * sizeof(unsigned int));
	
    for(i = 0; i < graph->num_nodes; i++)
    {
        j = graph->mapped_nodes[i];
    	algorithm->label[nodes[j].label]++;
    }
}

void enumerate_edges()
{
    struct Adjacency_list *nodes = graph->nodes;
    struct Algorithm *algorithm = graph->algorithm;
    
    unsigned int i = 0;
	unsigned int j = 0;

    memset(algorithm->load, 0, algorithm->num_partitions * sizeof(unsigned int));
	
    for(i = 0; i < graph->num_nodes; i++)
	{
		j = graph->mapped_nodes[i];
        algorithm->load[nodes[j].label] += nodes[j].degree;
        //algorithm->load[nodes[j].label] += nodes[j].degree + nodes[j].idegree;
	}
}

void enumerate_candidate_vertices()
{
    struct Adjacency_list *nodes = graph->nodes;
    struct Algorithm *algorithm = graph->algorithm;
    
    unsigned int i = 0;
    unsigned int j = 0;
    int k = 0;
    
    for(k = 0; k < algorithm->num_partitions; k++)
        algorithm->candidate_vertices[k] = 0;
    for(i = 0; i < graph->num_nodes; i++)
	{
		j = graph->mapped_nodes[i];
		if(nodes[j].label != graph->nodes[j].temp_label)
        {
		algorithm->candidate_vertices[nodes[j].temp_label] += nodes[j].degree;
        	//algorithm->candidate_vertices[nodes[j].temp_label] += nodes[j].degree + nodes[j].idegree;
        }
        
	}
}


void calculate_remaining_capacity()
{
	struct Algorithm *algorithm = graph->algorithm;
	int k = 0;
	for(k = 0; k < algorithm->num_partitions; k++)
    {
        algorithm->remaining_capacity[k] = algorithm->capacity[k] - algorithm->load[k];
    }
}

void *revolver(void *arguments)
{
	struct pthread_args_struct *args = arguments;
	int iteration = args->iteration;
	int index = args->index;
	unsigned int low = args->low;
	unsigned int high = args->high;
	
	struct Adjacency_list *nodes = graph->nodes;
    struct Algorithm *algorithm = graph->algorithm;

    unsigned int i = 0;
	unsigned int j = 0;
	int k = 0;
    action_selection(low, high);   
    probability_array[index] = 0;
    for(k = 0; k < algorithm->num_partitions; k++)
	{
		candidate_vertices_array[index][k] = 0;
	}
    
	for(i = low; i <= high; i++)
	{
		j = graph->mapped_nodes[i];
		if(nodes[j].label != nodes[j].action)
        {
			//candidate_vertices_array[index][nodes[j].action] += (nodes[j].degree + nodes[j].idegree);
            candidate_vertices_array[index][nodes[j].action] += nodes[j].degree;
        }
	}
   
	pthread_barrier_wait(&barrier);
	
	pthread_mutex_lock(&mutex);	
	if(!sentinel)
	{
		sentinel = 1;
		calculate_migration_probability();
                
        printf("Migration probability :");
        for(k = 0; k < algorithm->num_partitions; k++)
            printf("(%d %0.3f)", k, algorithm->migration_probability[k]);
        printf("\n");
        
        printf("Partition capacity    :");
        for(k = 0; k < algorithm->num_partitions; k++)
            printf("(%d %u)", k, algorithm->capacity[k]);
        printf(" %u\n", sum_uint(algorithm->capacity, algorithm->num_partitions));
        
        
        
        printf("Partition edges (load):");
        for(k = 0; k < algorithm->num_partitions; k++)
            printf("(%d %u)", k, algorithm->load[k]);
        printf(" %u %lli\n", sum_uint(algorithm->load, algorithm->num_partitions), graph->num_edges);
        
        printf("Partition vertices    :");
        for(k = 0; k < algorithm->num_partitions; k++)
            printf("(%d %u)", k, algorithm->label[k]);
        printf(" \n");
        
        printf("Remaining capacity    :");
        for(k = 0; k < algorithm->num_partitions; k++)
            printf("(%d %d)", k, algorithm->remaining_capacity[k]);
        printf(" \n");
        
        printf("Candidates Vertices   :");
        for(k = 0; k < algorithm->num_partitions; k++)
            printf("(%d %u)", k, algorithm->candidate_vertices[k]);
        printf(" \n");
	}
    pthread_mutex_unlock(&mutex);
    pthread_barrier_wait(&barrier);
    
    score_array[index] = .0;
	float temp_score_array[algorithm->num_partitions];
    

    for(i = low; i <= high; i++)
	{
		j = graph->mapped_nodes[i];        
        score_function(temp_score_array, algorithm->num_partitions, j);        
        if((nodes[j].label != nodes[j].action) && !migrate_vertex_revolver(j))
        //if((nodes[j].label != nodes[j].action))
        {
            pthread_mutex_lock(&mutex);    
                
                algorithm->load[nodes[j].label] -=  nodes[j].degree;
                algorithm->load[nodes[j].action]+= nodes[j].degree;
                algorithm->label[nodes[j].label]--;
                algorithm->label[nodes[j].action]++;
                
                /*
                algorithm->load[nodes[j].label]  -=  (nodes[j].degree + nodes[j].idegree);
                algorithm->load[nodes[j].action] += (nodes[j].degree + nodes[j].idegree);
                algorithm->label[nodes[j].label]--;
                algorithm->label[nodes[j].action]++;
                */
                
                nodes[j].label = nodes[j].action;
            pthread_mutex_unlock(&mutex);    
        }
        // Update score
        nodes[j].score = temp_score_array[nodes[j].label];
        score_array[index] += nodes[j].score;
        // Calculate reinforcement signal 
        objective_function(temp_score_array, algorithm->num_partitions, j);
	}
    pthread_barrier_wait(&barrier);

	// Update probability
	weighted_probability_update(low, high, iteration);
    for(i = low; i <= high; i++)
	{
		j = graph->mapped_nodes[i];    
        if(nodes[j].probability[nodes[j].label] > MAX_PROBABILITY)
            probability_array[index]++;
        //printf("%u %f %d\n", j, nodes[j].probability[nodes[j].label], nodes[j].probability[nodes[j].label] > MAX_PROBABILITY);
    }
    
    
    
    //probability_update_0(low, high);
    //pthread_barrier_wait(&barrier);
	// Calculate cut size
	cut_array[index] = 0;
	cut_array[index] += cut_size(low, high);
	
	for(k = 0; k < algorithm->num_partitions; k++)
	{
        external_edges_array_per_partition[index][k] = 0;
		internal_edges_array_per_partition[index][k] = 0;
        ingoing_cut_array_per_partition[index][k] = 0;
		external_edges_array_per_partition[index][k] = enumerate_partition_edges(low, high, k, 0);
		internal_edges_array_per_partition[index][k] = enumerate_partition_edges(low, high, k, 1);
        ingoing_cut_array_per_partition[index][k] = enumerate_partition_cut(low, high, k, 1);
	}
   
	// Calcuate Local edge ratio 
    local_edges_array[index] = 0;
	local_edges_array[index] += phi(low, high);
	
	pthread_barrier_wait(&barrier);
    
	pthread_mutex_lock(&mutex);
	if(sentinel == 1)
	{

        sentinel = -1;
        calculate_remaining_capacity();	
        
		statistics(iteration);
	}
    pthread_mutex_unlock(&mutex);
    //printf("%d %d\n", iteration, index);

	free(args);
	return(NULL);
}

void *spinner(void *arguments)
{
	struct pthread_args_struct *args = arguments;
	int iteration = args->iteration;
	int index = args->index;
	unsigned int low = args->low;
	unsigned int high = args->high;
	
	struct Adjacency_list *nodes = graph->nodes;
    struct Algorithm *algorithm = graph->algorithm;

	unsigned int i = 0;
	unsigned int j = 0;
	int k = 0;

	compute_score(low, high);

	for(k = 0; k < algorithm->num_partitions; k++)
	{
		candidate_vertices_array[index][k] = 0;
	}

	for(i = low; i <= high; i++)
	{
		j = graph->mapped_nodes[i];
		if(nodes[j].label != nodes[j].temp_label)
        {
	//		candidate_vertices_array[index][nodes[j].temp_label] += nodes[j].degree + nodes[j].idegree;
            candidate_vertices_array[index][nodes[j].temp_label] += nodes[j].degree;
        }
	}
    pthread_barrier_wait(&barrier);
	
	// First thread do the work
	// others will have to wait
    pthread_mutex_lock(&mutex);
	if(!sentinel)
	{
        sentinel = 1;
		calculate_migration_probability();
	}
    pthread_mutex_unlock(&mutex);

    score_array[index] = .0;
    float temp_score_array[algorithm->num_partitions];

    for(i = low; i <= high; i++)
	{
		//l = i - low;
		j = graph->mapped_nodes[i];
		if((nodes[j].label != nodes[j].temp_label) && !migrate_vertex_spinner(j))  // If you wanna migrate and you can
        {      
            pthread_mutex_lock(&mutex);          
                
                algorithm->load[nodes[j].label] -= nodes[j].degree;
                algorithm->label[nodes[j].label]--;
                algorithm->load[nodes[j].temp_label] += nodes[j].degree;
                algorithm->label[nodes[j].temp_label]++;
                
                /*
                algorithm->load[nodes[j].label] -= nodes[j].degree - nodes[j].idegree;
                algorithm->label[nodes[j].label]--;
                algorithm->load[nodes[j].temp_label] += nodes[j].degree + nodes[j].idegree;
                algorithm->label[nodes[j].temp_label]++;
                 */   
                nodes[j].label = nodes[j].temp_label;
            pthread_mutex_unlock(&mutex);  
        }
        score_function(temp_score_array, algorithm->num_partitions, j);		
        nodes[j].score = temp_score_array[nodes[j].label];
        score_array[index] += nodes[j].score;
	}

	pthread_barrier_wait(&barrier);

	// Calculate cut size
	cut_array[index] = 0;
	cut_array[index] += cut_size(low, high);
	
	for(k = 0; k < algorithm->num_partitions; k++)
	{
        external_edges_array_per_partition[index][k] = 0;
		internal_edges_array_per_partition[index][k] = 0;
        ingoing_cut_array_per_partition[index][k] = 0;
		external_edges_array_per_partition[index][k] = enumerate_partition_edges(low, high, k, 0);
		internal_edges_array_per_partition[index][k] = enumerate_partition_edges(low, high, k, 1);
        ingoing_cut_array_per_partition[index][k] = enumerate_partition_cut(low, high, k, 1);
	}
	
	// Calcuate Local edge ratio 
    local_edges_array[index] = 0;
	local_edges_array[index] += phi(low, high);
	
	pthread_barrier_wait(&barrier);
	pthread_mutex_lock(&mutex);	
	if(sentinel == 1)
	{
        sentinel = -1;
        
        printf("load ");
        for(k = 0; k < algorithm->num_partitions; k++)
            printf("%u ", algorithm->load[k]);
        printf("\n");
        calculate_remaining_capacity();
		statistics(iteration);
	}
	pthread_mutex_unlock(&mutex);
    
	free(args);
	return(NULL);
}

void *random_partitioning(void *arguments) {
    
    struct pthread_args_struct *args = arguments;
	int iteration = args->iteration;
	int index = args->index;
	unsigned int low = args->low;
	unsigned int high = args->high;
    
    struct Adjacency_list *nodes = graph->nodes;
    struct Algorithm *algorithm = graph->algorithm;

    unsigned int i = 0;
	unsigned int j = 0;
	int k = 0;

    
    for(i = low; i <= high; i++)
    {
		j = graph->mapped_nodes[i];
		#if GCC_VERSION > GCC_GENERIC_FEATURE_SUPPORT
            nodes[j].temp_label = ffrandom(algorithm->num_partitions - 1);
		#else
            nodes[j].temp_label = random_at_most(algorithm->num_partitions - 1);
        #endif
    }

	for(k = 0; k < algorithm->num_partitions; k++)
	{
		candidate_vertices_array[index][k] = 0;
	}
 
	for(i = low; i <= high; i++)
	{
		j = graph->mapped_nodes[i];
		if(nodes[j].label != nodes[j].temp_label)
        {
			candidate_vertices_array[index][nodes[j].temp_label] += nodes[j].degree;
        }
	}
    pthread_barrier_wait(&barrier);

    score_array[index] = .0;
    float temp_score_array[algorithm->num_partitions];


    for(i = low; i <= high; i++)
	{
		//l = i - low;
		j = graph->mapped_nodes[i];
		if(nodes[j].label != nodes[j].temp_label)  // If you wanna migrate and you can
        {    
            pthread_mutex_lock(&mutex);          
                algorithm->load[nodes[j].label] -= nodes[j].degree;
                algorithm->label[nodes[j].label]--;
                algorithm->load[nodes[j].temp_label] += nodes[j].degree;
                algorithm->label[nodes[j].temp_label]++;
                    
                nodes[j].label = nodes[j].temp_label;
            pthread_mutex_unlock(&mutex);  
        }
        score_function(temp_score_array, algorithm->num_partitions, j);		
        nodes[j].score = temp_score_array[nodes[j].label];
        score_array[index] += nodes[j].score;
	}

	pthread_barrier_wait(&barrier);

	// Calculate cut size
	cut_array[index] = 0;
	cut_array[index] += cut_size(low, high);
	
	for(k = 0; k < algorithm->num_partitions; k++)
	{
        external_edges_array_per_partition[index][k] = 0;
		internal_edges_array_per_partition[index][k] = 0;
        ingoing_cut_array_per_partition[index][k] = 0;
		external_edges_array_per_partition[index][k] = enumerate_partition_edges(low, high, k, 0);
		internal_edges_array_per_partition[index][k] = enumerate_partition_edges(low, high, k, 1);
        ingoing_cut_array_per_partition[index][k] = enumerate_partition_cut(low, high, k, 1);
	}
	
	// Calcuate Local edge ratio 
    local_edges_array[index] = 0;
	local_edges_array[index] += phi(low, high);
	
	pthread_barrier_wait(&barrier);
	pthread_mutex_lock(&mutex);	
	if(sentinel == 0)
	{
        sentinel = 1;
        
        printf("load ");
        for(k = 0; k < algorithm->num_partitions; k++)
            printf("%u ", algorithm->load[k]);
        printf("\n");
		statistics(iteration);
	}
	pthread_mutex_unlock(&mutex);
    
	free(args);
	return(NULL);    
}


void *hash_partitioning(void *arguments) {
    
    struct pthread_args_struct *args = arguments;
	int iteration = args->iteration;
	int index = args->index;
	unsigned int low = args->low;
	unsigned int high = args->high;
    
    struct Adjacency_list *nodes = graph->nodes;
    struct Algorithm *algorithm = graph->algorithm;

    unsigned int i = 0;
	unsigned int j = 0;
	int k = 0;

    
    for(i = low; i <= high; i++)
    {
		j = graph->mapped_nodes[i];
        nodes[j].temp_label = j % algorithm->num_partitions;
    }

	for(k = 0; k < algorithm->num_partitions; k++)
	{
		candidate_vertices_array[index][k] = 0;
	}
 
	for(i = low; i <= high; i++)
	{
		j = graph->mapped_nodes[i];
		if(nodes[j].label != nodes[j].temp_label)
        {
			candidate_vertices_array[index][nodes[j].temp_label] += nodes[j].degree;
        }
	}
    pthread_barrier_wait(&barrier);

    score_array[index] = .0;
    float temp_score_array[algorithm->num_partitions];


    for(i = low; i <= high; i++)
	{
		//l = i - low;
		j = graph->mapped_nodes[i];
		if(nodes[j].label != nodes[j].temp_label)  // If you wanna migrate and you can
        {    
            pthread_mutex_lock(&mutex);          
                algorithm->load[nodes[j].label] -= nodes[j].degree;
                algorithm->label[nodes[j].label]--;
                algorithm->load[nodes[j].temp_label] += nodes[j].degree;
                algorithm->label[nodes[j].temp_label]++;
                    
                nodes[j].label = nodes[j].temp_label;
            pthread_mutex_unlock(&mutex);  
        }
        score_function(temp_score_array, algorithm->num_partitions, j);		
        nodes[j].score = temp_score_array[nodes[j].label];
        score_array[index] += nodes[j].score;
	}

	pthread_barrier_wait(&barrier);

	// Calculate cut size
	cut_array[index] = 0;
	cut_array[index] += cut_size(low, high);
	
	for(k = 0; k < algorithm->num_partitions; k++)
	{
        external_edges_array_per_partition[index][k] = 0;
		internal_edges_array_per_partition[index][k] = 0;
        ingoing_cut_array_per_partition[index][k] = 0;
		external_edges_array_per_partition[index][k] = enumerate_partition_edges(low, high, k, 0);
		internal_edges_array_per_partition[index][k] = enumerate_partition_edges(low, high, k, 1);
        ingoing_cut_array_per_partition[index][k] = enumerate_partition_cut(low, high, k, 1);
	}
	
	// Calcuate Local edge ratio 
    local_edges_array[index] = 0;
	local_edges_array[index] += phi(low, high);
	
	pthread_barrier_wait(&barrier);
	pthread_mutex_lock(&mutex);	
	if(sentinel == 0)
	{
        sentinel = 1;
        
        printf("load ");
        for(k = 0; k < algorithm->num_partitions; k++)
            printf("%u ", algorithm->load[k]);
        printf("\n");
		statistics(iteration);
	}
	pthread_mutex_unlock(&mutex);
    
	free(args);
	return(NULL);    
}

void *range_partitioning(void *arguments) {
    
    struct pthread_args_struct *args = arguments;
	int iteration = args->iteration;
	int index = args->index;
	unsigned int low = args->low;
	unsigned int high = args->high;
    
    struct Adjacency_list *nodes = graph->nodes;
    struct Algorithm *algorithm = graph->algorithm;

    unsigned int i = 0;
	unsigned int j = 0;
	int k = 0;

    
    for(i = low; i <= high; i++)
    {
		j = graph->mapped_nodes[i];
//        printf("%u %u %u \n", j, graph->num_alloced_nodes / algorithm->num_partitions, (j * al//gorithm->num_partitions)/graph->num_alloced_nodes);
        
        nodes[j].temp_label = (nodes[j].source * algorithm->num_partitions)/graph->num_alloced_nodes;
        if(nodes[j].temp_label > algorithm->num_partitions)
        {
            fprintf(stderr, "Error: range_partitioning()\n");
            exit(EXIT_FAILURE);
        }
        
    }

	for(k = 0; k < algorithm->num_partitions; k++)
	{
		candidate_vertices_array[index][k] = 0;
	}
 
	for(i = low; i <= high; i++)
	{
		j = graph->mapped_nodes[i];
		if(nodes[j].label != nodes[j].temp_label)
        {
			candidate_vertices_array[index][nodes[j].temp_label] += nodes[j].degree;
        }
	}
    pthread_barrier_wait(&barrier);

    score_array[index] = .0;
    float temp_score_array[algorithm->num_partitions];


    for(i = low; i <= high; i++)
	{
		//l = i - low;
		j = graph->mapped_nodes[i];
		if(nodes[j].label != nodes[j].temp_label)  // If you wanna migrate and you can
        {    
            pthread_mutex_lock(&mutex);          
                algorithm->load[nodes[j].label] -= nodes[j].degree;
                algorithm->label[nodes[j].label]--;
                algorithm->load[nodes[j].temp_label] += nodes[j].degree;
                algorithm->label[nodes[j].temp_label]++;
                    
                nodes[j].label = nodes[j].temp_label;
            pthread_mutex_unlock(&mutex);  
        }
        score_function(temp_score_array, algorithm->num_partitions, j);		
        nodes[j].score = temp_score_array[nodes[j].label];
        score_array[index] += nodes[j].score;
	}

	pthread_barrier_wait(&barrier);

	// Calculate cut size
	cut_array[index] = 0;
	cut_array[index] += cut_size(low, high);
	
	for(k = 0; k < algorithm->num_partitions; k++)
	{
        external_edges_array_per_partition[index][k] = 0;
		internal_edges_array_per_partition[index][k] = 0;
        ingoing_cut_array_per_partition[index][k] = 0;
		external_edges_array_per_partition[index][k] = enumerate_partition_edges(low, high, k, 0);
		internal_edges_array_per_partition[index][k] = enumerate_partition_edges(low, high, k, 1);
        ingoing_cut_array_per_partition[index][k] = enumerate_partition_cut(low, high, k, 1);
	}
	
	// Calcuate Local edge ratio 
    local_edges_array[index] = 0;
	local_edges_array[index] += phi(low, high);
	
	pthread_barrier_wait(&barrier);
	pthread_mutex_lock(&mutex);	
	if(sentinel == 0)
	{
        sentinel = 1;
        
        printf("load ");
        for(k = 0; k < algorithm->num_partitions; k++)
            printf("%u ", algorithm->load[k]);
        printf("\n");
		statistics(iteration);
	}
	pthread_mutex_unlock(&mutex);
    
	free(args);
	return(NULL);    
}



void *pagerank(void *arguments) {
    struct pthread_args_struct *args = arguments;
	int iteration = args->iteration;
	int index = args->index;
	unsigned int low = args->low;
	unsigned int high = args->high;
    
    (void) iteration;
    (void) index;
    (void) low;
    (void) high;
    //printf("%d\n", iteration);
    
    struct Adjacency_list *nodes = graph->nodes;
    struct Algorithm *algorithm = graph->algorithm;
    struct Adjacency_list_node *src_node = NULL;
    external_messages_array[index] = 0;
    internal_messages_array[index] = 0;
    (void) nodes;
    (void) algorithm;
    
    unsigned int i = 0;
	unsigned int j = 0;
	int k = 0;
    unsigned int l = 0;
    (void) i;
    (void) j;
    (void) k;
    float temp_probability = 0;
    
    long long external_messages = 0;
    long long internal_messages = 0;
    
    if(iteration == 1) {
        for(i = low; i <= high; i++)
        {
            j = graph->mapped_nodes[i];
            //nodes[j].probability_pr = (float) 1/graph->num_nodes;
            nodes[j].probability_pr = (float) 1;
            //printf("%u %d\n", j, nodes[j].label);
        }
    } else {
        //nodes[j].probability_pr
        for(i = low; i <= high; i++)
        {
            j = graph->mapped_nodes[i];
            
            if(nodes[j].idegree)
            {
                src_node = nodes[j].head;
                temp_probability = 0;
                while(src_node)
                {   
                    if(src_node->icount)
                    {
                        l = src_node->destination;    
                        temp_probability += (float) nodes[l].probability_pr / nodes[l].degree;
                        if(nodes[j].label != nodes[l].label)
                        {
                            //pthread_mutex_lock(&mutex);	
                                //printf("  %u %d %u %d\n", j, nodes[j].label, l, nodes[l].label);
                            external_messages_array[index] += src_node->icount;
                            //pthread_mutex_unlock(&mutex);	
                        } else {
                             internal_messages_array[index] += src_node->icount;
                        }
                    }
                    src_node = src_node->next;  
                }
                nodes[j].probability_pr = temp_probability;
            }
        }
    }
    
    pthread_barrier_wait(&barrier);
	pthread_mutex_lock(&mutex);	
	if((iteration > 1) && (sentinel == 0))
	{
        sentinel = 1;
        for(i = 0; i < MAX_THREADS_NUM; i++) 
        {
            external_messages += external_messages_array[i];
            internal_messages += internal_messages_array[i];
        }
        
        printf("Iteration=%d, Num Edges=%lli,  ExternalMessages=%lli   InternalMessages=%lli\n", iteration, graph->num_edges, external_messages, internal_messages);
        
        #ifdef FILEIO
            char file_name[MAX_PATH_LEN];
            memset(file_name, '\0', sizeof(file_name));	
            snprintf(file_name, MAX_PATH_LEN, "%s/%s/%d_%s_stats_%s", algorithm->dir_name, LOG_DIR, algorithm->num_partitions, algorithm->name, algorithm->base_name);

            FILE *file_descriptor = fopen(file_name, "a");
            if(!file_descriptor)
            {
                fprintf(stderr, "Error on opening the input file %s\n", file_name);
                exit(EXIT_FAILURE);
            }
            // Add text file header
            if(iteration == 2)
            {
                file_descriptor = freopen(file_name, "w", file_descriptor);
                fprintf(file_descriptor, "Step NumEdges InterPartitionMessages IntraPartitionMessages \n");
            }
            fprintf(file_descriptor, "%d %lli %lli\n", iteration, external_messages, internal_messages);
            
            if(!fclose(file_descriptor))
            ; // NoOp, we're good to go!
            else
            {
                fprintf(stderr, "Error on closing the input file\n");
                exit(EXIT_FAILURE);
            }
        #endif
        
        
        
        
        
		//statistics(iteration);
	}
	pthread_mutex_unlock(&mutex);
    
    
	free(args);
	return(NULL);    
}


void calculate_first_migration_probability()
{
	struct Algorithm *algorithm = graph->algorithm;
    
	int k = 0;

	for(k = 0; k < algorithm->num_partitions; k++)
	{
		if(algorithm->remaining_capacity[k] <= 0 || !algorithm->candidate_vertices[k])
			algorithm->migration_probability[k] = 0.0;
		else
		{
			algorithm->migration_probability[k] = algorithm->remaining_capacity[k] / algorithm->candidate_vertices[k];
			if(algorithm->migration_probability[k] > 1)
				algorithm->migration_probability[k] = 1;
		}
	}
}


void calculate_migration_probability()
{
	struct Algorithm *algorithm = graph->algorithm;
    
	int i = 0;
	int k = 0;

	for(k = 0; k < algorithm->num_partitions; k++)
	{
		algorithm->candidate_vertices[k] = 0;
		for(i = 0; i < MAX_THREADS_NUM; i++)
		{
			algorithm->candidate_vertices[k]  += candidate_vertices_array[i][k];
		}
	}
//        unsigned int threshold = (unsigned int) ((double) graph->num_edges / (algorithm->num_partitions * (algorithm->num_partitions - 1)));
	for(k = 0; k < algorithm->num_partitions; k++)
	{
		if(algorithm->remaining_capacity[k] <= 0 || !algorithm->candidate_vertices[k])
			algorithm->migration_probability[k] = 0.0;
		else
		{
			algorithm->migration_probability[k] = (float) algorithm->remaining_capacity[k] / algorithm->candidate_vertices[k];
			if(algorithm->migration_probability[k] > 1)
				algorithm->migration_probability[k] = 1;
		}
	}
}

int migrate_vertex_spinner(int index)
{
	struct Adjacency_list *nodes = graph->nodes;
	struct Algorithm *algorithm = graph->algorithm;
	int migrate = 1;
	float random_number = 0.0;
	random_number = (float) random()/RAND_MAX;
	if((random_number < algorithm->migration_probability[nodes[index].temp_label]))
		migrate = 0;
	return(migrate);
}

int migrate_vertex_revolver(int index)
{
	struct Adjacency_list *nodes = graph->nodes;
	struct Algorithm *algorithm = graph->algorithm;
	int migrate = 1;
	float random_number = 0.0;
	random_number = (float) random()/RAND_MAX;
	if(random_number < algorithm->migration_probability[nodes[index].action])
		migrate = 0;
	return(migrate);
}

void compute_score(unsigned int low, unsigned int high)
{
	struct Adjacency_list *nodes = graph->nodes;
	struct Algorithm *algorithm = graph->algorithm;
	float temp_score_array[algorithm->num_partitions];
	memset(temp_score_array, 0, sizeof(temp_score_array));
	unsigned int i = 0;
	unsigned int j = 0;
	
    if(high >= graph->num_nodes)
    {
        fprintf(stderr, "Error: Exceeding graph size.\n");
    }
    else
    {
        for(i = low; i <= high; i++)
        {
            j = graph->mapped_nodes[i];
            score_function(temp_score_array, algorithm->num_partitions, j);
            nodes[j].temp_label = return_label(temp_score_array, algorithm->num_partitions, nodes[j].label);
        }
    }
}

 // Global ratio of local edges (phi)
long long phi(unsigned int low, unsigned int high)
{
	struct Adjacency_list *nodes = graph->nodes;
	struct Adjacency_list_node *src_node = NULL;
	unsigned int i = 0;
	unsigned int j = 0;
    unsigned int l = 0;
	long long local_edges = 0;
    if(high >= graph->num_nodes)
    {
        fprintf(stderr, "Error: Exceeding graph size.\n");
    }
    else
    {    
        for(i = low; i <= high; i++)
        {
            j = graph->mapped_nodes[i];
            src_node = nodes[j].head;
            while(src_node)
            {
                //dst_index = graph->mapped_nodes[src_node->destination];
                //l = graph->mapped_nodes_[src_node->destination];
                l = src_node->destination;
                if(nodes[j].label == nodes[l].label)
                //if(nodes[j].label == nodes[src_node->destination].label)
                {
                    //local_edges += (src_node->count + src_node->icount);
                    local_edges += src_node->count;
                }
                src_node = src_node->next;
            }
        }
    }    
	return(local_edges);
}

long long cut_size(unsigned int low, unsigned int high)
{
	struct Adjacency_list *nodes = graph->nodes;
    struct Adjacency_list_node *src_node = NULL;
    
	unsigned int i = 0;
	unsigned int j = 0;
    unsigned int l = 0;
	long long cut_edges = 0;
    
    if(high >= graph->num_nodes)
    {
        fprintf(stderr, "Error: Exceeding graph size.\n");
    }
    else
    {
        for(i = low; i <= high; i++)
        {
            j = graph->mapped_nodes[i];
            src_node = nodes[j].head;
            while(src_node)
            {
                //dst_index = graph->mapped_nodes[src_node->destination];
                //l = graph->mapped_nodes_[src_node->destination];
                l = src_node->destination;
                if(nodes[j].label != nodes[l].label)
                //if(nodes[j].label != nodes[src_node->destination].label)
                {
                    //cut_edges += (src_node->count + src_node->icount);
                    cut_edges += src_node->count;
                }
                src_node = src_node->next;
            }    
        }
    }
	return(cut_edges);
}

// Outgoing/Ingoing edges per partition 0/1
unsigned int enumerate_partition_cut(unsigned int low, unsigned int high, int label, int type)
{
    
	struct Adjacency_list *nodes = graph->nodes;
    struct Adjacency_list_node *src_node = NULL;
    //struct Adjacency_list_node *dst_node = NULL;
	unsigned int i = 0;
	unsigned int j = 0;
    unsigned int l = 0;
    
	unsigned int cut_edges = 0;
    
    if(high >= graph->num_nodes)
    {
        fprintf(stderr, "Error: Exceeding graph size.\n");
    }
    else
    {
        for(i = low; i <= high; i++)
        {
            j = graph->mapped_nodes[i];
            src_node = nodes[j].head;
            if(nodes[j].label == label)
            {
                while(src_node)
                {
                    l = src_node->destination;
                    if((type == 0) && nodes[j].label != nodes[l].label)
                    {
                        cut_edges += src_node->count;
                    }
                    if((type == 1) && nodes[j].label != nodes[l].label)
                    {
                        cut_edges += src_node->icount;
                    }
                    src_node = src_node->next;
                }
            }            
        }
    }
	return(cut_edges);
}

// Outgoing/Ingoing edges per partition 0/1
unsigned int enumerate_partition_edges(unsigned int low, unsigned int high, int label, int type)
{
	struct Adjacency_list *nodes = graph->nodes;
    struct Adjacency_list_node *src_node = NULL;
	unsigned int i = 0;
	unsigned int j = 0;
    unsigned int l = 0;
	unsigned int edges = 0;
    if(high >= graph->num_nodes)
    {
        fprintf(stderr, "Error: Exceeding graph size.\n");
    }
    else
    {    
        for(i = low; i <= high; i++)
        {
            j = graph->mapped_nodes[i];
            src_node = nodes[j].head;
            if(nodes[j].label == label)
            {
                while(src_node)
                {
                    //l = graph->mapped_nodes_[src_node->destination];
                    l = src_node->destination;
                    if((type == 0) && (nodes[j].label != nodes[l].label))
                    {
                        edges += src_node->count;
                        //edges += (src_node->count + src_node->icount);
                    }
                    else if((type == 1) && (nodes[j].label == nodes[l].label))
                    {
                        edges += src_node->count;
                        //edges += (src_node->count + src_node->icount);
                    }
                    src_node = src_node->next;
                }    
            }
        }
    }
	return(edges);
}

void objective_function(float *array, int count, int index)
{
	int k = 0; // Number of partitions
    (void) k;
    unsigned int l = 0; // destination vertex index
    float temp_weigth = 0;
    (void) temp_weigth;

    struct Adjacency_list_node *src_node = NULL;
    (void) src_node;
    struct Adjacency_list *nodes = graph->nodes;
    struct Algorithm *algorithm = graph->algorithm;
    (void) algorithm;
    
    #if GCC_VERSION > GCC_GENERIC_FEATURE_SUPPORT
        int max_i =  ffmaxi(array, count, 0);
    #else
        int max_i = max_index_float(array, count, 0);
    #endif    
    
    src_node = nodes[index].head;
    while(src_node)
	{
        l = src_node->destination;
        temp_weigth = 0;
        temp_weigth = src_node->count + src_node->icount;
        if(nodes[l].action == max_i)
            nodes[l].signal_[max_i] += temp_weigth;
        else if(algorithm->migration_probability[max_i])
            nodes[l].signal_[max_i] += temp_weigth;
        src_node = src_node->next;      
    }
}

void score_function(float *array, int count, int index)
{
	int k = 0; // Number of partitions
    unsigned int l = 0; // destination vertex index
    float sum_weight = .0;
	//float penalty_function = .0;
    //float label_propagation_algorithm = 0.0;
    float temp_weigth = 0;
    
    struct Adjacency_list *nodes = graph->nodes;
    struct Algorithm *algorithm = graph->algorithm;
    
    struct Adjacency_list_node *src_node = NULL;
    struct Adjacency_list_node *dst_node = NULL;

    memset(array, 0, sizeof(float) * count);
    
    (void) dst_node;
    
    float lpa[algorithm->num_partitions];
    memset(lpa, 0, sizeof(lpa));
    (void) lpa;
    float pf[algorithm->num_partitions];
    memset(pf, 0, sizeof(pf));
    (void) pf;
    float min_pf = 0.0;    
    float max_pf = 0.0;
    float range_pf = 0.0;
    float sum_pf = 0.0;
    char contained_negative = 1;
    
    //int energy[algorithm->num_partitions];
    //memset(energy, 0, sizeof(energy));
    
    (void) min_pf;
    (void) max_pf;
    (void) range_pf;
    (void) sum_pf;
    (void) contained_negative;
    int dst_index = 0;
    (void) dst_index;


    
    src_node = nodes[index].head;
    while(src_node)
    {   
        l = src_node->destination;
        temp_weigth = 0;
        temp_weigth = src_node->count + src_node->icount;
        array[nodes[l].label] += temp_weigth;
        sum_weight += temp_weigth;  
        src_node = src_node->next;  
        //energy[nodes[l].label]++;
    }
    
    //for(k = 0; k < count; k++)
    //    printf("%f ", array[k]);
    //printf("\n");
    
    for(k = 0; k < count; k++)
	{
        pf[k] = algorithm->l - (float) algorithm->load[k]/algorithm->capacity[k];
        //if(pf[k]  < 0.9 * (algorithm->l - algorithm->c))
            //pf[k] = -pf[k];
        lpa[k] = 0;        
		if(array[k])
			lpa[k] = array[k]/sum_weight;
        if(pf[k] < 0)
            contained_negative = 0;
        
       array[k] = lpa[k] + pf[k];
	}
    
    
    if(!contained_negative)
    {
        #if GCC_VERSION > GCC_GENERIC_FEATURE_SUPPORT
            max_pf = ffmax(pf, count);
            min_pf = ffmin(pf, count);
        #else
            max_pf = max_value_float(pf, count);
            min_pf = min_value_float(pf, count);
        #endif
    
        range_pf = max_pf - min_pf;
        for(k = 0; k < count; k++)
            pf[k] = (pf[k] - min_pf) / range_pf;
    }
    
    #if GCC_VERSION > GCC_GENERIC_FEATURE_SUPPORT
        sum_pf = ffsum(pf, count);
    #else
        sum_pf = sum_float(pf, count);
    #endif
    
    for(k = 0; k < count; k++)
        pf[k] = pf[k] / sum_pf;
    
    
    for(k = 0; k < count; k++)
        array[k] = (lpa[k] + pf[k]) / 2;
    
    //for(k = 0; k < count; k++)
    //    printf("%f ", array[k]);
    //printf(" %f\n", sum_float(array, count));
    //printf("%d %d %d %f\n", nodes[index].label, energy[nodes[index].label], sum_int(energy, count), my_energy);
    
    

    //exit(0);
      //  array[k] = (lpa[k] + pf[k]);
}


int return_label(float *array, int count, int label)
{	
	int k = 0;
	int num_labels = 0;
	int returned_label = -1;
	
	#if GCC_VERSION > GCC_GENERIC_FEATURE_SUPPORT
	    int *returned_labels = ffmaxi(array, count, &num_labels);
	#else
		int *returned_labels = max_indices_float(array, count, &num_labels);
	#endif
	
	if(num_labels > 1)
	{
		for(k = 0; k < num_labels; k++)
		{
			if(returned_labels[k] == label)
			{
				returned_label = returned_labels[k];
				break;
			}
		}
		if(returned_label == -1)
		{
			#if GCC_VERSION > GCC_GENERIC_FEATURE_SUPPORT
                returned_label = returned_labels[ffrandom(num_labels - 1)];
		    #else
                returned_label = returned_labels[random_at_most(num_labels - 1)];
            #endif			
		}
	}
	else 
		returned_label = returned_labels[0];	
	free(returned_labels);
	return(returned_label);
}

void free_all()
{
    struct Adjacency_list *nodes = graph->nodes;
    struct Algorithm *algorithm = graph->algorithm;
    
    unsigned int i = 0;
    unsigned int j = 0;
    
    free(algorithm->load);
    free(algorithm->capacity);
    free(algorithm->label);	
    free(algorithm->candidate_vertices);
	free(algorithm->remaining_capacity);
	free(algorithm->migration_probability);
    free(algorithm);
	
	for(i = 0; i < MAX_THREADS_NUM; i++)
	{
		free(external_edges_array_per_partition[i]);
		free(internal_edges_array_per_partition[i]);
        free(ingoing_cut_array_per_partition[i]);        
		free(candidate_vertices_array[i]);
	}
    
    for(i = 0; i < graph->num_nodes; i++)
    {
		j = graph->mapped_nodes[i];
        free(nodes[j].probability);
    }
	
	free(graph->nodes);
	free(graph->mapped_nodes);
	free(graph);
}

void statistics(int iteration)
{
	struct Algorithm *algorithm = graph->algorithm;
    
	int i = 0;
	int k = 0;
	
	score         = .0;
    long long cut_edges =  0;
    float cut_edges_ratio  = .0;	
    
    unsigned int external_edges_per_partition[algorithm->num_partitions];
	memset(external_edges_per_partition, 0, sizeof(external_edges_per_partition));
	unsigned int internal_edges_per_partition[algorithm->num_partitions];
	memset(internal_edges_per_partition, 0, sizeof(internal_edges_per_partition));
	unsigned int ingoing_cut_per_partition[algorithm->num_partitions];
	memset(ingoing_cut_per_partition, 0, sizeof(ingoing_cut_per_partition));
    
    float max_load   = .0;
    float min_load   = .0;
    float expected_load = .0;   
    float max_norm_load = .0;    
    float max_min_load = .0;
    float local_edges_ratio  = .0;	
    long long local_edges  = 0;	
    
    float probability_ratio = 0.;
    
	for(i = 0; i < MAX_THREADS_NUM; i++)
	{
		score += score_array[i];
		cut_edges += cut_array[i];
        cut_edges_ratio += cut_array[i];
		local_edges_ratio += local_edges_array[i];
        local_edges += local_edges_array[i];
        probability_ratio += probability_array[i];
	}		
    probability_ratio = probability_ratio/graph->num_nodes;
    
    //score = (float) score/graph->num_nodes;
    //printf("%lli\n", graph->num_nodes);
    local_edges_ratio = local_edges_ratio/graph->num_edges;
    cut_edges_ratio = cut_edges_ratio/graph->num_edges;
	
	for(k = 0; k < algorithm->num_partitions; k++)
	{
		for(i = 0; i < MAX_THREADS_NUM; i++)
		{
			external_edges_per_partition[k] += external_edges_array_per_partition[i][k];
			internal_edges_per_partition[k] += internal_edges_array_per_partition[i][k];
            ingoing_cut_per_partition[k] += ingoing_cut_array_per_partition[i][k];
		}
	}
    
    //expected_load = (float) graph->num_alloced_edges/algorithm->num_partitions;
    expected_load = (float) graph->num_edges/algorithm->num_partitions;
    #if GCC_VERSION > GCC_GENERIC_FEATURE_SUPPORT
        max_load = (float) ffmax(algorithm->load, algorithm->num_partitions);
        min_load = (float) ffmin(algorithm->load, algorithm->num_partitions);
    #else
        max_load = (float) max_value_uint(algorithm->load, algorithm->num_partitions);
        min_load = (float) min_value_uint(algorithm->load, algorithm->num_partitions);
    #endif
    max_norm_load = max_load/expected_load;
    max_min_load  = max_load/min_load;

	int unique = 0;
	for(k = 0; k < algorithm->num_partitions; k++)
	{
		if(!algorithm->label[k])
		{
			fprintf(stderr, "Error: Algorithm eliminated one or more partitions.\n");
            exit(0);
		}
		else
			unique++;
	}

    //int min_idx = 0;
    //int max_idx = algorithm->num_partitions - 1;
    //float min_cap = (float) algorithm->load[min_idx]/algorithm->capacity[min_idx];
    //float max_cap = (float) algorithm->load[max_idx]/algorithm->capacity[max_idx];
    
	#ifdef ALL
        printf("\n");
	    printf("Step   : %d\n", iteration);
		printf("Score  : %f\n", score);
		printf("Cut Edg: %f (%lli)\n", cut_edges_ratio, cut_edges);
    	printf("Loc Edg: %f (%lli)\n", local_edges_ratio, local_edges);
    	printf("Rho    : %f\n", max_norm_load);
		printf("MML    : %f\n", max_min_load);
        printf("Pratio : %f\n", probability_ratio);
        //printf("MinMaxC: %f\n", (float) (max_cap + min_cap)/2);

    	printf("Load per partitions\n");
    	#if GCC_VERSION > GCC_GENERIC_FEATURE_SUPPORT
	        printf("    Max: %u\n",   ffmax(algorithm->load, algorithm->num_partitions));
	        printf("    Min: %u\n",   ffmin(algorithm->load, algorithm->num_partitions));
    	    printf("   Mean: %f\n",   ffmean(algorithm->load, algorithm->num_partitions));
            printf("    Sum: %u\n", ffsum(algorithm->load, algorithm->num_partitions));
            printf("Max/E/k: %f\n", (float) ffmax(algorithm->load, algorithm->num_partitions)/(graph->num_edges/algorithm->num_partitions));
            printf("Max idx: %u\n", ffmaxi(algorithm->load, algorithm->num_partitions, 0));
        #else
    	    printf("    Max: %u\n", max_value_uint(algorithm->load, algorithm->num_partitions));
	        printf("    Min: %u\n", min_value_uint(algorithm->load, algorithm->num_partitions));
	        printf("   Mean: %u\n", (unsigned int) mean_uint(algorithm->load, algorithm->num_partitions));
            printf("    Sum: %u\n", sum_uint(algorithm->load, algorithm->num_partitions));
            printf("Max/E/k: %f\n", (float) max_value_uint(algorithm->load, algorithm->num_partitions)/((float) graph->num_edges/algorithm->num_partitions));
            printf("Max idx: %d\n", max_index_uint(algorithm->load, algorithm->num_partitions, 0));
        #endif
	
	    printf("Nodes per partitions\n");
	    #if GCC_VERSION > GCC_GENERIC_FEATURE_SUPPORT
    	    printf("    Max: %u\n",      ffmax(algorithm->label, algorithm->num_partitions));
	        printf("    Min: %u\n",      ffmin(algorithm->label, algorithm->num_partitions));
	        printf("   Mean: %u\n",      (unsigned int) ffmean(algorithm->label, algorithm->num_partitions));
            printf("    Sum: %u\n",      ffsum(algorithm->label, algorithm->num_partitions));
    	#else
	        printf("    Max: %u\n",      max_value_uint(algorithm->label, algorithm->num_partitions));
	        printf("    Min: %u\n",      min_value_uint(algorithm->label, algorithm->num_partitions));
    	    printf("   Mean: %u\n",      (unsigned int) mean_uint(algorithm->label, algorithm->num_partitions));
            printf("    Sum: %u\n",      sum_uint(algorithm->label, algorithm->num_partitions));
        #endif	

	    printf("External edges per partitions (cut size)\n");
    	#if GCC_VERSION > GCC_GENERIC_FEATURE_SUPPORT
	        printf("    Max: %u\n",      ffmax(external_edges_per_partition, algorithm->num_partitions));
	        printf("    Min: %u\n",      ffmin(external_edges_per_partition, algorithm->num_partitions));
    	    printf("   Mean: %u\n",      (unsigned int) ffmean(external_edges_per_partition, algorithm->num_partitions));
            printf("    Sum: %u\n",      ffsum(external_edges_per_partition, algorithm->num_partitions));
            printf("Max idx: %d\n",      ffmaxi(external_edges_per_partition, algorithm->num_partitions, 0));
	    #else
    	    printf("    Max: %u\n",      max_value_uint(external_edges_per_partition, algorithm->num_partitions));
	        printf("    Min: %u\n",      min_value_uint(external_edges_per_partition, algorithm->num_partitions));
	        printf("   Mean: %u\n",      (unsigned int) mean_uint(external_edges_per_partition, algorithm->num_partitions));
            printf("    Sum: %u\n",      sum_uint(external_edges_per_partition, algorithm->num_partitions));
            printf("Max idx: %d\n",      max_index_uint(external_edges_per_partition, algorithm->num_partitions, 0));
        #endif	

    	printf("Internal edges per partitions (local edges)\n");
	    #if GCC_VERSION > GCC_GENERIC_FEATURE_SUPPORT
	        printf("    Max: %u\n",      ffmax(internal_edges_per_partition, algorithm->num_partitions));
    	    printf("    Min: %u\n",      ffmin(internal_edges_per_partition, algorithm->num_partitions));
	        printf("   Mean: %u\n",      (unsigned int) ffmean(internal_edges_per_partition, algorithm->num_partitions));
            printf("    Sum: %u\n",      ffsum(internal_edges_per_partition, algorithm->num_partitions));
            printf("  Sum/E: %f\n",      (float) ffsum(internal_edges_per_partition, algorithm->num_partitions)/graph->num_edges);
            printf("Max idx: %d\n",      ffmaxi(internal_edges_per_partition, algorithm->num_partitions, 0));
    	#else
	        printf("    Max: %u\n",      max_value_uint(internal_edges_per_partition, algorithm->num_partitions));
	        printf("    Min: %u\n",      min_value_uint(internal_edges_per_partition, algorithm->num_partitions));
    	    printf("   Mean: %u\n",      (unsigned int) mean_uint(internal_edges_per_partition, algorithm->num_partitions));
            printf("    Sum: %u\n",      sum_uint(internal_edges_per_partition, algorithm->num_partitions));
            printf("  Sum/E: %f\n",      (float) sum_uint(internal_edges_per_partition, algorithm->num_partitions)/graph->num_edges);
            printf("Max idx: %d\n",      max_index_uint(internal_edges_per_partition, algorithm->num_partitions, 0));
        #endif	
        
        printf("Ingoing cut per partitions (ingoing cut size)\n");
	    #if GCC_VERSION > GCC_GENERIC_FEATURE_SUPPORT
	        printf("    Max: %u\n",      ffmax(ingoing_cut_per_partition, algorithm->num_partitions));
    	    printf("    Min: %u\n",      ffmin(ingoing_cut_per_partition, algorithm->num_partitions));
	        printf("   Mean: %u\n",      (unsigned int) ffmean(ingoing_cut_per_partition, algorithm->num_partitions));
            printf("    Sum: %u\n",      ffsum(ingoing_cut_per_partition, algorithm->num_partitions));
            printf("Max idx: %d\n",      ffmaxi(ingoing_cut_per_partition, algorithm->num_partitions, 0));
    	#else
	        printf("    Max: %u\n",      max_value_uint(ingoing_cut_per_partition, algorithm->num_partitions));
	        printf("    Min: %u\n",      min_value_uint(ingoing_cut_per_partition, algorithm->num_partitions));
    	    printf("   Mean: %u\n",      (unsigned int) mean_uint(ingoing_cut_per_partition, algorithm->num_partitions));
            printf("    Sum: %u\n",      sum_uint(ingoing_cut_per_partition, algorithm->num_partitions));
            printf("Max idx: %d\n",      max_index_uint(ingoing_cut_per_partition, algorithm->num_partitions, 0));
        #endif
	#endif

	#ifdef FILEIO
        char file_name[MAX_PATH_LEN];
        memset(file_name, '\0', sizeof(file_name));	
        snprintf(file_name, MAX_PATH_LEN, "%s/%s/%d_%s_stats_%s", algorithm->dir_name, LOG_DIR, algorithm->num_partitions, algorithm->name, algorithm->base_name);

        FILE *file_descriptor = fopen(file_name, "a");
        if(!file_descriptor)
        {
            fprintf(stderr, "Error on opening the input file %s\n", file_name);
            exit(EXIT_FAILURE);
        }
        // Add text file header
        if(iteration <= 1)
        {
            file_descriptor = freopen(file_name, "w", file_descriptor);
            fprintf(file_descriptor, "                                                      Load                        Node                       ExEdges                    INEdges\n");
            fprintf(file_descriptor, "Step         Score NCutEdges  Phi   Rho MML       max      min     mean  sum  max/E/K     max      min     mean sum      max      min     mean sum      max      min     mean sum sum/E      max      min     mean sum p_ratio min_cap max_cap \n");
        }
        #if GCC_VERSION > GCC_GENERIC_FEATURE_SUPPORT
            fprintf(file_descriptor, "%d %f %lli %f %f %f \
            %u %u %u %u %f \
            %u %u %u %u \
            %u %u %u %u \
            %u %u %u %u %f \
            %u %u %u %u \
            %f\n",
            iteration, score, cut_edges, local_edges_ratio, max_norm_load, max_min_load,
            
            ffmax(algorithm->load, algorithm->num_partitions), ffmin(algorithm->load, algorithm->num_partitions), (unsigned int) ffmean(algorithm->load, algorithm->num_partitions), ffsum(algorithm->load, algorithm->num_partitions), (float) ffmax(algorithm->load, algorithm->num_partitions)/((float) graph->num_edges/algorithm->num_partitions),
            
            ffmax(algorithm->label, algorithm->num_partitions), ffmin(algorithm->label, algorithm->num_partitions), (unsigned int) ffmean(algorithm->label, algorithm->num_partitions), ffsum(algorithm->label, algorithm->num_partitions),
            
            ffmax(external_edges_per_partition, algorithm->num_partitions), ffmin(external_edges_per_partition, algorithm->num_partitions), (unsigned int) ffmean(external_edges_per_partition, algorithm->num_partitions), ffsum(external_edges_per_partition, algorithm->num_partitions),
            
            ffmax(internal_edges_per_partition, algorithm->num_partitions), ffmin(internal_edges_per_partition, algorithm->num_partitions), (unsigned int) ffmean(internal_edges_per_partition, algorithm->num_partitions), ffsum(internal_edges_per_partition, algorithm->num_partitions), (float) ffsum(internal_edges_per_partition, algorithm->num_partitions)/graph->num_edges,
            
            ffmax(ingoing_cut_per_partition, algorithm->num_partitions), ffmin(ingoing_cut_per_partition, algorithm->num_partitions), (unsigned int) ffmean(ingoing_cut_per_partition, algorithm->num_partitions), ffsum(ingoing_cut_per_partition, algorithm->num_partitions),
            
            probability_ratio);
            
        #else
            //fprintf(file_descriptor, "%d %f %d %f %f %f %d %d %d %d %f %d %d %f %d %d %d %f %d %d %d %f %d %f %d %d %f %d\n",
            fprintf(file_descriptor, "%d %f %lli %f %f %f \
            %u %u %u %u %f \
            %u %u %u %u \
            %u %u %u %u \
            %u %u %u %u %f \
            %u %u %u %u \
            %f\n",
            iteration, score, cut_edges, local_edges_ratio, max_norm_load, max_min_load,
            
            max_value_uint(algorithm->load, algorithm->num_partitions), min_value_uint(algorithm->load, algorithm->num_partitions), (unsigned int) mean_uint(algorithm->load, algorithm->num_partitions), sum_uint(algorithm->load, algorithm->num_partitions), (float) max_value_uint(algorithm->load, algorithm->num_partitions)/((float) graph->num_edges/algorithm->num_partitions),
            
            max_value_uint(algorithm->label, algorithm->num_partitions), min_value_uint(algorithm->label, algorithm->num_partitions), (unsigned int) mean_uint(algorithm->label, algorithm->num_partitions), sum_uint(algorithm->label, algorithm->num_partitions),
            
            max_value_uint(external_edges_per_partition, algorithm->num_partitions), min_value_uint(external_edges_per_partition, algorithm->num_partitions), (unsigned int) mean_uint(external_edges_per_partition, algorithm->num_partitions), sum_uint(external_edges_per_partition, algorithm->num_partitions),
            
            max_value_uint(internal_edges_per_partition, algorithm->num_partitions), min_value_uint(internal_edges_per_partition, algorithm->num_partitions), (unsigned int) mean_uint(internal_edges_per_partition, algorithm->num_partitions), sum_uint(internal_edges_per_partition, algorithm->num_partitions), (float) sum_uint(internal_edges_per_partition, algorithm->num_partitions)/graph->num_edges,
    
            max_value_uint(ingoing_cut_per_partition, algorithm->num_partitions), min_value_uint(ingoing_cut_per_partition, algorithm->num_partitions), (unsigned int) mean_uint(ingoing_cut_per_partition, algorithm->num_partitions), sum_uint(ingoing_cut_per_partition, algorithm->num_partitions),
            
            probability_ratio);
        
        #endif
        
        if(!fclose(file_descriptor))
            ; // NoOp, we're good to go!
        else
        {
            fprintf(stderr, "Error on closing the input file\n");
            exit(EXIT_FAILURE);
        }

	#endif
}

int write_partition()
{
	struct Adjacency_list *nodes = graph->nodes;
	struct Algorithm *algorithm = graph->algorithm;

	unsigned int i = 0;
	unsigned int j = 0;
    
    char file_name[MAX_PATH_LEN];
    memset(file_name, '\0', sizeof(file_name));	
    snprintf(file_name, MAX_PATH_LEN, "%s/%s/%d_%s_parts_%s", algorithm->dir_name, LOG_DIR, algorithm->num_partitions, algorithm->name, algorithm->base_name);

	
	FILE *file_descriptor = fopen(file_name, "w");
    if(!file_descriptor)
    {
        fprintf(stderr, "Error on opening the input file\n");
        return(EXIT_ERROR);
    }
	// Write partitioning results
	for(i = 0; i < graph->num_nodes; i++)
	{
		j = graph->mapped_nodes[i];
		fprintf(file_descriptor,"%u %d\n", j, nodes[j].label);
	}
	
	if(!fclose(file_descriptor))
		; // NoOp, we're good to go!
	else
	{
		fprintf(stderr, "Error on closing the input file\n");
        return(EXIT_ERROR);
	}
    return(EXIT_OK);
}

void print_graph()
{
    struct Adjacency_list *nodes = graph->nodes;
    struct Adjacency_list_node *src_node = NULL;
    

    unsigned int i = 0;
    unsigned int j = 0;
    unsigned int l = 0;
    printf("Source ID (outdegree, indegree): Destination ID (outgoing, ingoing)\n");
    for(i = 0; i < graph->num_nodes; i++)
    {
        j = graph->mapped_nodes[i];
        printf("%u (%d %d)--> ",j, nodes[j].degree, nodes[j].idegree);
        src_node = nodes[j].head;
        while(src_node)
        {
            //l = graph->mapped_nodes_[src_node->destination];
            l = src_node->destination;
            printf(" %u (%d %d)", l, src_node->count, src_node->icount);
            src_node = src_node->next;
        }
        printf("\n");
    }
}


void print_edges_per_vertex() 
{
    struct Adjacency_list *nodes = graph->nodes;
    struct Algorithm *algorithm = graph->algorithm;

    unsigned int i = 0;
    unsigned int j = 0;
    /*
    unsigned int min_number_edges = 0;
    unsigned int max_number_edges = 0;
    

    for(i = 0; i < graph->num_nodes; i++)
    {
        j = graph->mapped_nodes[i];
        //printf("%u %d\n",j, nodes[j].degree);
        if(max_number_edges < nodes[j].degree)
            max_number_edges = nodes[j].degree;
    }
    
    printf("Min=%d, Max=%d\n", min_number_edges, max_number_edges);
    unsigned int *edges_degree = malloc(sizeof(unsigned int) * (max_number_edges + 1));
    memset(edges_degree, 0, sizeof(unsigned int) * (max_number_edges + 1));
    
    for(i = 0; i < graph->num_nodes; i++)
    {
       j = graph->mapped_nodes[i];
       edges_degree[nodes[j].degree]++;
    }
    */
    
    char file_name[MAX_PATH_LEN];
    memset(file_name, '\0', sizeof(file_name));	
    snprintf(file_name, MAX_PATH_LEN, "%s/%s/%d_%s_distEdges_%s", algorithm->dir_name, LOG_DIR, algorithm->num_partitions, algorithm->name, algorithm->base_name);

    FILE *file_descriptor = fopen(file_name, "w");
    if(!file_descriptor)
    {
        fprintf(stderr, "Error on opening the input file %s\n", file_name);
        exit(EXIT_FAILURE);
    }

    for(i = 0; i < graph->num_nodes; i++)
    {
        j = graph->mapped_nodes[i];
        fprintf(file_descriptor, "%u %d\n", j, nodes[j].degree);
        //printf("%u %d\n",i, edges_degree[i]);
    }
    
    /*
    for(i = 0; i < graph->num_nodes; i++)
    {
       j = graph->mapped_nodes[i];
       edges_degree[nodes[j].degree]++;
    }
    */
    
    if(!fclose(file_descriptor))
        ; // NoOp, we're good to go!
    else
    {
        fprintf(stderr, "Error on closing the input file\n");
        exit(EXIT_FAILURE);
    }
    
    
}

int write_partition_e(char *file_name, char *comment_style, char *delimiter, int columns)
{
    struct Adjacency_list *nodes = graph->nodes;
    struct Algorithm *algorithm = graph->algorithm;
    
    FILE *file_descriptor = fopen(file_name,"r");
    if(!file_descriptor)
    {
        fprintf(stderr, "Error on opening the input file %s\n", file_name);
        return(EXIT_ERROR);
    }

    char  *format = edge_format_specifier(file_descriptor, comment_style, delimiter);
	if(!format)
	{
		fprintf(stderr, "Error: edge_format_specifier()\n");
		return(EXIT_ERROR);
	}
    
    
    char file_name_e[MAX_PATH_LEN];
    memset(file_name_e, '\0', sizeof(file_name_e));	
    snprintf(file_name_e, MAX_PATH_LEN, "%s/%s/%d_%s_edges_%s", algorithm->dir_name, LOG_DIR, algorithm->num_partitions, algorithm->name, algorithm->base_name);
    
	FILE *file_descriptor_e = fopen(file_name_e, "w");
    if(!file_descriptor_e)
    {
        fprintf(stderr, "Error on opening the input file %s\n", file_name_e);
        return(EXIT_ERROR);
    }

    unsigned int source = 0;
    unsigned int destination = 0;
	unsigned int source_idx = 0;
    unsigned int old_source = 0;
    unsigned int destination_idx = 0;
    char flag = 0;
    
    int offset = ftell(file_descriptor);
    if(fscanf(file_descriptor, format, &source, &destination) == columns)
    {
        fseek(file_descriptor, offset, SEEK_SET);
        old_source = source;
    }
    else
    {
        fprintf(stderr, "Error on fscanf(%s)\n", file_name_e);
        return(EXIT_ERROR);
    }
    
	while(fscanf(file_descriptor, format, &source, &destination) == columns)
    { 
        source_idx = graph->mapped_nodes[source]; 
        //destination_idx = graph->mapped_nodes_[destination]; // Uncomment last part of vertex_map()
        
        if(old_source != source)
        {
            fprintf(file_descriptor_e, "\n");
            flag = 0;
        }
        
        if(!flag)
        {
            fprintf(file_descriptor_e, "%u_%d %u_%d ", source, nodes[source_idx].label, destination, nodes[destination_idx].label);  
            flag = 1;
        }
        else
        {
            fprintf(file_descriptor_e, "%u_%d ", destination, nodes[destination_idx].label);  
        }
    
        old_source = source;
        //fprintf(file_descriptor_e, "%d_%d %d_%d\n", source, nodes[source_idx].label, destination, nodes[destination_idx].label);  // This is old
    }
    free(format);

	if(!fclose(file_descriptor_e))
		; // NoOp, we're good to go!
	else
	{
		fprintf(stderr, "Error on closing the input file\n");
		return(EXIT_ERROR);
	}
    

	if(!fclose(file_descriptor))
         ; // NoOp, we're good to go!
    else
    {
        fprintf(stderr, "Error on closing the input file\n");
        return(EXIT_ERROR);
    }

	return(EXIT_OK);
}

