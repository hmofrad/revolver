/*
 * distrograph.c Base implementation of Revolver/Spinner graph partitioning algorithms
 * (c) Mohammad H. Mofrad, 2017
 * (e) mohammad.hmofrad@pitt.edu
 */
 
 #include "libgraph.h"
 
 int main(int argc, char *argv[])
{
    int i = 0; // num_nodes/thread index
    //int j = 0; // mapped_nodes index
    //int k = 0; // num_partitions index
    
    srand(time(NULL));
    struct timeval t1, t2;
	double elapsedTime;
    gettimeofday(&t1, NULL);

    if(((argc != 10) && (argc != 13)) || !(!strcmp(argv[4],"spinner") || !strcmp(argv[4],"revolver") || !strcmp(argv[4],"random") || !strcmp(argv[4],"hash") || !strcmp(argv[4],"range") || !strcmp(argv[4],"pagerank")))
    {
        fprintf(stderr, "USAGE: %s -n <num_partitions> -a <revolver|spinner|random|hash|range|pagerank> -f [<edge_file_name> <comment_style> <delimiter> <columns>] -v [<vetex_file_name> <delimiter>]\n", argv[0]);
        for(i = 0; i < argc; i++)
            fprintf(stderr, "argv[%2d]: %s\n", i, argv[i]);
        
        exit(EXIT_FAILURE);
    }

    // Allocate graph data structure
	graph = malloc(sizeof(struct Graph));
    if(!graph)
        exit(EXIT_FAILURE);
    memset(graph, 0, sizeof(struct Graph));
	
    char comment[strlen(argv[7])];
    memset(comment, '\0', strlen(argv[7]));
    memcpy(comment, argv[7], strlen(argv[7]));
    
    char delimiter[strlen(argv[8])];
    memset(delimiter, '\0', strlen(argv[8]));
    memcpy(delimiter, argv[8], strlen(argv[8]));
    
    int columns = atoi(argv[9]);
    (void) columns;

/*
    int status = edge_open_file_mmapped(argv[6], argv[7], argv[8], atoi(argv[9]));
    if(status == EXIT_ERROR)
	{
		fprintf(stderr, "Error: edge_open_file_mmapped()\n");
        exit(EXIT_FAILURE);
	}   
*/    
/*    
    int status = para_edge_open_file(argv[6], argv[7], argv[8], atoi(argv[9]));
    if(status == EXIT_ERROR)
	{
		fprintf(stderr, "Error: para_edge_open_file()\n");
        exit(EXIT_FAILURE);
	}
 */  

    
    // Read input file
    int status = edge_open_file(argv[6], argv[7], argv[8], atoi(argv[9]));
	if(status == EXIT_ERROR)
	{
		fprintf(stderr, "Error: open_edge_file()\n");
        exit(EXIT_FAILURE);
	}
    
    // Refine node indicces
    status = vertex_map();
	if(status == EXIT_ERROR)
	{
		fprintf(stderr, "Error: vertex_map()\n");
        exit(EXIT_FAILURE);
	}
    
    //print_graph();
    //exit(0);

    if(!strcmp(argv[4], "pagerank")) 
    {
        printf("Sort the vertex file before calling open_vertex_file(%s)\n", argv[11]);	
        printf("    sort -k 1 -n  %s\n", argv[11]);
        status = vertex_open_file(argv[11], argv[12]);
        if(status == EXIT_ERROR)
        {
            fprintf(stderr, "Error: open_file()\n");
            exit(EXIT_FAILURE);
        }
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
    
    
    //print_edges_per_vertex();
    //exit(0);
    
    // LA alpha and beta parameters
    algorithm->alpha = 1;  // LA reward signal
	algorithm->beta  = .1;  // LA penalty signal

    if(strcmp(argv[4], "pagerank")) 
    {
        status = algorithm_initialization();
        if(status == EXIT_ERROR)
        {
            fprintf(stderr, "Error: algorithm_initialization()\n");
            exit(EXIT_FAILURE);
        }
    }    
    
    int k1;
	unsigned int k1_num;
	int k2;
	unsigned int k2_num;
    multiprocessing_pool(graph->num_nodes, MAX_THREADS_NUM, &k1, &k1_num, &k2, &k2_num);
    
	pthread_t tid[MAX_THREADS_NUM];
	memset(tid, 0, sizeof(tid));
	struct pthread_args_struct *args = NULL;
    
	int iteration = 0;
    int ep = algorithm->e;
	float score_old = .0;

    while(iteration < MAX_ITER)
    {	
        iteration++;
        sentinel = 0;
		score_old = score;
		score = .0;
		pthread_barrier_init(&barrier, NULL, MAX_THREADS_NUM);
    	for(i = 0; i < MAX_THREADS_NUM; i++)
    	{	
            args = malloc(sizeof(struct pthread_args_struct) * MAX_THREADS_NUM);
            if(!args)
                return(EXIT_FAILURE);
            memset(args, 0, sizeof(struct pthread_args_struct));
            args->iteration = iteration;
            args->index = i;
    
		    if(i < k1)
		    {
	    		args->low = i * k1_num;
		    	args->high = ((i + 1)* k1_num) - 1;
	    	}
    		else if((i >= k1) && (i < k1 + k2))
	    	{
			    args->low = (k1 * k1_num) + ((i - k1) * k2_num);
    			args->high = (k1 * k1_num) + (((i + 1) - k1) * k2_num) - 1;
	    	}

            if(!strcmp(algorithm->name, "revolver"))
                status = pthread_create(&tid[i], NULL, (void*)revolver, (void *)args);	
            else if(!strcmp(algorithm->name, "spinner"))
                status = pthread_create(&tid[i], NULL, (void*)spinner, (void *)args);
            else if(!strcmp(algorithm->name, "random"))
                status = pthread_create(&tid[i], NULL, (void*)random_partitioning, (void *)args);
            else if(!strcmp(algorithm->name, "hash"))
                status = pthread_create(&tid[i], NULL, (void*)hash_partitioning, (void *)args);
            else if(!strcmp(algorithm->name, "range"))
                status = pthread_create(&tid[i], NULL, (void*)range_partitioning, (void *)args);
            else if(!strcmp(algorithm->name, "pagerank"))
                status = pthread_create(&tid[i], NULL, (void*)pagerank, (void *)args);
            
            
            if (status != EXIT_SUCCESS)
            {
                fprintf(stderr, "Error: pthread_create()\n");
                exit(EXIT_FAILURE);
            }
            
	    }

    	for(i = 0; i < MAX_THREADS_NUM; i++)
		{
	        pthread_join(tid[i], NULL);
		}
    
        status = computation_is_halted(score, score_old, &ep);
        if(status == 0)
        {
            printf("Computation has been converged.\n");
            break;
        }
        if((!strcmp(argv[4], "pagerank")) && iteration == 3)
        {
            printf("Pagerank is done\n");
            break;
        }
        

    }
    
    /*
    long long j = 0;
    for(i = 0; i < graph->num_nodes; i++){
        j = graph->mapped_nodes[i];
        printf("%lli %d\n", j, graph->nodes[j].label);
    }
    */
    
	status = write_partition();
    if(status == EXIT_ERROR)
	{
		fprintf(stderr, "Error: write_partition()\n");
        exit(EXIT_FAILURE);
	}
  
/*
    status = write_partition_e(argv[6], comment, delimiter, columns);
    if(status == EXIT_ERROR)
	{
		fprintf(stderr, "Error: write_partition_e()\n");
        exit(EXIT_FAILURE);
	}
	printf("[x] %s (%d)\n", algorithm->name, status);
*/
    //print_graph();
    free_all();

    // compute and print the elapsed time in milliseconds and seconds
    gettimeofday(&t2, NULL);
    elapsedTime  = (t2.tv_sec - t1.tv_sec) * 1000;      // sec to ms
    elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000;   // us to ms
    printf("%.0f milliseconds | %.2lf seconds\n", elapsedTime, elapsedTime / 1000);
    
    return(EXIT_SUCCESS);
}
