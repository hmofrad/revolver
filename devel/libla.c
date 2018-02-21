/*
 * libla.c Learning Automata (LA) implementation
 * (c) Mohammad H. Mofrad, 2017
 * (e) mohammad.hmofrad@pitt.edu
 */
 
#include "libla.h"
struct Graph *graph;
void action_selection(unsigned int low, unsigned int high)
{
	struct Adjacency_list *nodes = graph->nodes;
	struct Algorithm *algorithm = graph->algorithm;
	unsigned int i = 0;
	unsigned int j = 0;
    //int k = 0;
    
    // roulette_wheel, ruler_tool, scissors_tool, trimmer_tool 
    // are available for selecting the next action
	for(i = low; i <= high; i++)
	{
		j = graph->mapped_nodes[i];
		nodes[j].action = trimmer_tool(nodes[j].probability, algorithm->num_partitions); // ruler_tool, roulette_wheel
	}
}

void probability_update(unsigned int low, unsigned int high)
{
    struct Adjacency_list *nodes = graph->nodes;
	struct Algorithm *algorithm = graph->algorithm;
	unsigned int i = 0;
	unsigned int j = 0;
	int k = 0;
	for(i = low; i <= high; i++)
	{
		j = graph->mapped_nodes[i];

		if(!nodes[j].signal)
		{
			for(k = 0; k < algorithm->num_partitions; k++)
			{
			    if(k == nodes[j].action)
				    nodes[j].probability[k] = nodes[j].probability[k] + algorithm->alpha * (1 - nodes[j].probability[k]);
				else
					nodes[j].probability[k] = (1 - algorithm->alpha) * nodes[j].probability[k];
			}
		}
		else
		{
			for(k = 0; k < algorithm->num_partitions; k++)
			{
			    if(k == nodes[j].action)
				{
		        	nodes[j].probability[k] = (1 - algorithm->beta) * nodes[j].probability[k];
				}
			    else
				{
					nodes[j].probability[k] = (algorithm->beta / (algorithm->num_partitions - 1)) + (1 - algorithm->beta) * nodes[j].probability[k];
				}
			}
		}
        nodes[j].signal  = 1;
	}
}

 
void weighted_probability_update(unsigned int low, unsigned int high, int iteration)
{
    struct Adjacency_list *nodes = graph->nodes;
	struct Algorithm *algorithm = graph->algorithm;
	unsigned int i = 0;
	unsigned int j = 0;
	int k = 0;
    int k_ = 0;
    int kk = 0;
    int kkk = 0;
    (void) kkk;
    
    float alpha_[algorithm->num_partitions];
    memset(alpha_, 0, sizeof(alpha_));
    float beta_[algorithm->num_partitions];
    memset(beta_, 0, sizeof(beta_));
    
    float sig_[algorithm->num_partitions];
    memset(sig_, 0, sizeof(sig_));
    
    int indices[algorithm->num_partitions];
    memset(indices, 0, sizeof(indices));
    
    float sum_signal = 0;

    char sign[algorithm->num_partitions];
    memset(sign, 0, sizeof(sign));
    int positive_num = 0;
    int negative_num = 0;        
    
    //int negative_num = 0;
    //int negative_ind = 0;        
    float seprator = 0.0;
    float error_sum = 0.0;
    float threshold = (float) 1/(algorithm->num_partitions * sqrt(algorithm->num_partitions));
    (void) threshold;
    
    float w0 = .9;
    float w1 = .4;
    float weight = ((w0 - w1) * iteration * sqrt(algorithm->num_partitions)) / MAX_ITER;
    int weight_index = 0;
    //printf("%f\n", weight);
    
    int idx = -1;
	for(i = low; i <= high; i++)
	{
		j = graph->mapped_nodes[i];
        positive_num = 0;
        negative_num = 0;
            
            #if GCC_VERSION > GCC_GENERIC_FEATURE_SUPPORT
                seprator = ffsum(nodes[j].signal_, algorithm->num_partitions)/algorithm->num_partitions;
            #else
                seprator = sum_float(nodes[j].signal_, algorithm->num_partitions)/algorithm->num_partitions;
            #endif  
            
            #if GCC_VERSION > GCC_GENERIC_FEATURE_SUPPORT
                weight_index = ffmaxi(nodes[j].signal_, algorithm->num_partitions, 0);
            #else
                weight_index = max_index_float(nodes[j].signal_, algorithm->num_partitions, 0);
            #endif  
            
            if(j == 10) {
                printf("Original Probability: ");
                for(k = 0; k < algorithm->num_partitions; k++) {
                    printf("%f ", nodes[j].probability[k]);
                }
                printf("\n");
                
                printf("Original Signal: ");
                for(k = 0; k < algorithm->num_partitions; k++) {
                    printf("%f ", nodes[j].signal_[k]);
                }
                printf("\n");
            }
            
            
            nodes[j].signal_[weight_index] = nodes[j].signal_[weight_index] + nodes[j].signal_[weight_index] * weight;
            //(void) weight;
            //(void) weight_index;
            // > seprator for count only
            // >= seprator for count + icount
            for(k = 0; k < algorithm->num_partitions; k++)
            {
                if(nodes[j].signal_[k] >= seprator)
                {
                    sign[k] = 0;
                    positive_num++;
                }
                else
                {
                    sign[k] = 1;
                    negative_num++;
                }
            }

            float temp_prob[algorithm->num_partitions];
            memset(temp_prob, 0, sizeof(temp_prob));
            if(j == idx)
                memcpy(temp_prob, nodes[j].probability, algorithm->num_partitions * sizeof(float));
            
            float positive_part[positive_num];
            memset(positive_part, 0, sizeof(positive_part));
            
            int positive_ind[positive_num];
            memset(positive_ind, 0, sizeof(positive_ind));

            float negative_part[negative_num];
            memset(negative_part, 0, sizeof(negative_part));
            
            int negative_ind[negative_num];
            memset(negative_ind, 0, sizeof(negative_ind));
           
            positive_num = 0;
            negative_num = 0;
            
            for(k = 0; k < algorithm->num_partitions; k++)
            {
                if(nodes[j].signal_[k] >= seprator)
                {
                    positive_part[positive_num] = nodes[j].signal_[k];
                    positive_ind[positive_num] = k;
                    positive_num++;                    
                }
                else 
                {
                    negative_part[negative_num] = nodes[j].signal_[k];
                    negative_ind[negative_num] = k;
                    negative_num++;
                }   
            }
            
            #if GCC_VERSION > GCC_GENERIC_FEATURE_SUPPORT
                sum_signal = ffsum(positive_part, positive_num);
            #else
                sum_signal = sum_float(positive_part, positive_num);
            #endif  
            
            if(sum_signal > 0)
            {
                for(k = 0; k < positive_num; k++)
                {
                    if(positive_part[k])
                        positive_part[k] = positive_part[k] / sum_signal;
                }
            }
            
            
            for(k = 0; k < negative_num; k++)
            {
                if(negative_part[k] < 0)
                    negative_part[k] = -negative_part[k];
            }
            
            #if GCC_VERSION > GCC_GENERIC_FEATURE_SUPPORT
                sum_signal = ffsum(negative_part, negative_num);
            #else
                sum_signal = sum_float(negative_part, negative_num);
            #endif  

            if(sum_signal > 0)
            {
                for(k = 0; k < negative_num; k++)
                    negative_part[k] = negative_part[k] / sum_signal;
            }
            else
            {
                for(k = 0; k < negative_num; k++)
                    negative_part[k] = (float) 1 / negative_num;
            }
            
            positive_num = 0;            
            negative_num = 0;

            for(k = 0; k < algorithm->num_partitions; k++)
            {
                if(!sign[k])
                {
                    sig_[k] = positive_part[positive_num];
                    positive_num++;
                }
                else
                {
                    sig_[k] = negative_part[negative_num];
                    negative_num++;
                }
            }

            bubble_sort(positive_part, positive_ind, positive_num);

            bubble_sort(negative_part, negative_ind, negative_num);
            
            for(k = 0; k < negative_num; k++)
            {
                sig_[negative_ind[k]] = 0;
            }           
            


            memcpy(indices, negative_ind, sizeof(negative_ind));
            memcpy(indices+negative_num, positive_ind, sizeof(positive_ind));
            
            if(j == 10) {
                printf("Modified Signal: ");
                for(k = 0; k < algorithm->num_partitions; k++) {
                    printf("%f ", nodes[j].signal_[k]);
                }
                printf("\n");
                printf("Normalized Signal: ");
                for(k = 0; k < algorithm->num_partitions; k++) {
                    printf("%f ", sig_[k]);
                }
                printf("\n");
                
                printf("Signal Indices: ");
                for(k = 0; k < algorithm->num_partitions; k++) {
                    printf("%d ", k);
                }
                printf("\n");
                
                printf("Positive Values: ");
                for(k = 0; k < positive_num; k++)
                {
                    printf("%f ", positive_part[k]);
                }
                printf("\n");
                
                printf("Positive Indices: ");
                for(k = 0; k < positive_num; k++)
                {
                    printf("%d ", positive_ind[k]);
                }
                printf("\n");
                
                printf("Negative Values: ");
                for(k = 0; k < negative_num; k++)
                {
                    printf("%f ", negative_part[k]);
                }
                printf("\n");
                
                printf("Negative Indices: ");
                for(k = 0; k < negative_num; k++)
                {
                    printf("%d ", negative_ind[k]);
                }
                printf("\n");
            }

            
            for(k = 0; k < algorithm->num_partitions; k++)
            {
                //printf("%d %f %f %d\n", kk, nodes[j].probability[kk], sig_[kk], sign[kk]);
                kk = indices[k];
                if(!sign[kk] && sig_[kk])
                {
                    
                    alpha_[kk] = sig_[kk] * algorithm->alpha;
                    nodes[j].probability[kk] = nodes[j].probability[kk] + alpha_[kk] * (1 - nodes[j].probability[kk]);
                    for(k_ = 0; k_ < algorithm->num_partitions; k_++)
                    {
                        if(k_ != kk)
                            nodes[j].probability[k_] = (1 - alpha_[kk]) * nodes[j].probability[k_];
                    }
                    
                    /*
                    if(kk == weight_index)
                        nodes[j].probability[kk] = nodes[j].probability[kk] + algorithm->alpha * (1 - nodes[j].probability[kk]);
                    else
                        nodes[j].probability[kk] = (1 -  algorithm->alpha) * nodes[j].probability[kk];
                    */
                }
                else if(sign[kk] && sig_[kk])
                {
                    beta_[kk] = sig_[kk] * algorithm->beta;
                    nodes[j].probability[kk] = (1 - beta_[kk]) * nodes[j].probability[kk];
                    for(k_ = 0; k_ < algorithm->num_partitions; k_++)
                    {
                        if(k_ != kk)
                            nodes[j].probability[k_] = (beta_[kk] / (algorithm->num_partitions - 1)) + (1 - beta_[kk]) * nodes[j].probability[k_];
                    }
                }
                //printf("%d %f\n", kk, nodes[j].probability[kk]);
                
            }
            

            kk  = algorithm->num_partitions - 1;
            for(k = algorithm->num_partitions - 1; k <= 0; k++)
            {
                kkk = indices[k];
                if(nodes[j].probability[kkk] < threshold)
                {
                    if(nodes[j].probability[kkk] > 0)
                    {
                    //    printf("%d %f %f %f\n", j, nodes[j].probability[k], threshold, ffsum(nodes[j].probability, algorithm->num_partitions));
                        nodes[j].probability[kk] += nodes[j].probability[kkk];
                        nodes[j].probability[kkk] = 0;
                        #if GCC_VERSION > GCC_GENERIC_FEATURE_SUPPORT
                            error_sum = ffsum(nodes[j].probability, algorithm->num_partitions);
                        #else
                            error_sum = sum_float(nodes[j].probability, algorithm->num_partitions);
                        #endif  
                        if((1 - error_sum) > 0.0001)
                        {
                            nodes[j].probability[kk] += 1 - error_sum;
                        }
                        kk--;
                    }
                    else if(nodes[j].probability[kkk] == 0)
                        break;
                }
            }
            if(j == 10) {
            
            printf("New Probability: ");
                for(k = 0; k < algorithm->num_partitions; k++) {
                    printf("%f ", nodes[j].probability[k]);
                }
                printf("\n");
            }
            
            
            if(j == idx)
            {
                printf("%d  %d %f %d\n", nodes[j].action, nodes[j].temp_label, seprator, indices[algorithm->num_partitions - 1]);
                printf(" k  i New sig.   New sig. Old prob.   New prob.\n");
                for(k = 0; k < algorithm->num_partitions; k++)
                {
                    printf("%2d %2d %6.3f %6.3f %f %f\n", k, indices[k], nodes[j].signal_[indices[k]] , sig_[indices[k]], temp_prob[indices[k]], nodes[j].probability[indices[k]]);
                }
            }
            
            memset(nodes[j].signal_, 0, algorithm->num_partitions * sizeof(float));
            
            #if GCC_VERSION > GCC_GENERIC_FEATURE_SUPPORT
                error_sum = ffsum(nodes[j].probability, algorithm->num_partitions);
            #else
                error_sum = sum_float(nodes[j].probability, algorithm->num_partitions);
            #endif  

            if((1 - error_sum) > 0.0001)
            {
                
                printf("==>%d  %d %f %d\n", nodes[j].action, nodes[j].temp_label, seprator, indices[algorithm->num_partitions - 1]);
                printf(" k  i sig_     sig      o_prob   n_prob\n");
                for(k = 0; k < algorithm->num_partitions; k++)
                {
                    printf("%2d %2d %f %f %f %f\n", k, indices[k], nodes[j].signal_[indices[k]] , sig_[indices[k]], temp_prob[indices[k]], nodes[j].probability[indices[k]]);
                }
                
                printf("ERROR %d %f\n", j, error_sum); 
                exit(0);
            }

            nodes[j].signal  = 1;
	}
}

void reinforcement_signal(unsigned int low, unsigned int high)
{
	struct Adjacency_list *nodes = graph->nodes;
    struct Algorithm *algorithm = graph->algorithm;
    (void) algorithm;
	unsigned int i = 0;
	unsigned int j = 0;
    
	for(i = low; i <= high; i++)
	{
		j = graph->mapped_nodes[i];
        //if((nodes[j].action == nodes[j].temp_label) || (nodes[j].local_edges <= nodes[j].local_edges_old))
        //if(nodes[j].action == nodes[j].temp_label)
        if(nodes[j].action == nodes[j].temp_label)
        {
        //    if(nodes[j].probability[nodes[j].action] - nodes[j].temp_score < 0)
                nodes[j].signal = 0;
          //  else
            //    nodes[j].signal = 1;
            //printf("%f %f %f\n", nodes[j].probability[nodes[j].action], nodes[j].temp_score, nodes[j].probability[nodes[j].action] - nodes[j].temp_score);
            //if((fabs(nodes[j].probability[nodes[j].action] - nodes[j].temp_score) < 0.1) && (iteration > 50))
            //    nodes[j].signal = 0;
            //else
              //  nodes[j].signal = 1;
        }
        else
            nodes[j].signal = 1;
        
        /*
        if(iteration < 90)
        {
        //if(nodes[j].local_edges >= nodes[j].local_edges_old)
            if(nodes[j].action == nodes[j].temp_label)
		        nodes[j].signal = 0;
            else
	    		nodes[j].signal = 1;
        }
        else
        {
            if(nodes[j].local_edges <= nodes[j].local_edges_old)
		        nodes[j].signal = 0;
            else
	    		nodes[j].signal = 1;
        }
		*/
	}
}

int ruler_tool(float *probabilities, int num_events)
{
	float sum = .0;
	int k = 0;
	
	float random_number = (float)rand() / (float)RAND_MAX;

	do
	{
		sum += probabilities[k];
		k++;
	}
	while(sum < random_number && k < num_events);
	return(k-1);
}

/* Better implementation of roulette_wheel */
int roulette_wheel(float *probabilities, int num_events)
{
	int RAND_NUM = num_events * 2;
	
	int random_number;
	do
	{
		random_number = (int) (RAND_NUM * (float) random()/RAND_MAX);
	}while((random_number > RAND_NUM) || (random_number < 0));

	if((random_number > RAND_NUM) || (random_number < 0))
	{
		fprintf(stderr, "Error on roulette_wheel(): random_number = %d, RAND_NUM = %d\n", random_number, RAND_NUM);
		exit(EXIT_ERROR);
	}
	int i = 0;
	int k = 0;

	int sequence_length[num_events];
	memset(sequence_length, 0, sizeof(sequence_length));

	for(k = 0; k < num_events; k++)
	{
		sequence_length[k] = (int) (RAND_NUM * probabilities[k]);
	}
	
	#if GCC_VERSION > GCC_GENERIC_FEATURE_SUPPORT
	    int diff = RAND_NUM - ffsum(sequence_length, num_events);
	#else
	    int diff = RAND_NUM - sum_int(sequence_length, num_events);
	#endif
    
	while(diff)
	{
		#if GCC_VERSION > GCC_GENERIC_FEATURE_SUPPORT
            sequence_length[ffrandom(num_events - 1)]++;
		#else
            sequence_length[random_at_most(num_events - 1)]++;
		#endif			
		diff--;
	}

	int random_sequence[RAND_NUM];
	memset(random_sequence, 0, sizeof(random_sequence));
	int sequence_begin = 0;
	int sequence_end = 0;
	for(k = 0; k < num_events; k++)
	{
		sequence_end += sequence_length[k];
		for(i = sequence_begin; i < sequence_end; i++)
		    random_sequence[i] = k;
		sequence_begin += sequence_length[k];
	}
	
    #if GCC_VERSION > GCC_GENERIC_FEATURE_SUPPORT
        ffshuffle(random_sequence, RAND_NUM);
    #else
        ishuffle(random_sequence, RAND_NUM);
    #endif
    
    return(random_sequence[random_number]);
}

int scissors_tool(float *probabilities, int num_events)
{
	int left_num = 0;
	int right_num = 0;
	int *left_indices = NULL;
	float *left_probabilities = NULL;
	int *right_indices = NULL;
	float *right_probabilities = NULL;
    float *probability = malloc(num_events * sizeof(float));
	memcpy(probability, probabilities, num_events * sizeof(float));
	float sum_probability = 1.0;
	int factor = 2;
	float seprator = (float) sum_probability/factor;
	int remanining_events = num_events;
	int min_idx = 0;
	
	float random_number = 0.0;
	int returned_index = 0;
	
	while(remanining_events > factor)
	{
		random_number = (float) random()/RAND_MAX;
	    split_probability(probability, remanining_events, min_idx, &left_num, &left_probabilities, &left_indices, &right_num, &right_probabilities, &right_indices);
		free(probability);
		if(random_number < seprator)
		{
			probability = malloc(left_num * sizeof(float));
			memcpy(probability, left_probabilities, left_num * sizeof(float));
		    remanining_events = left_num;
			min_idx = left_indices[0];
		}
		else
		{
	        probability = malloc(right_num * sizeof(float));
			memcpy(probability, right_probabilities, right_num * sizeof(float));
		    remanining_events = right_num;
			min_idx = right_indices[0];
		}
		free(left_indices);
		free(left_probabilities);		
		free(right_indices);
		free(right_probabilities);
	}

    random_number = (float) random()/RAND_MAX;
	returned_index = min_idx + 0;
	if(random_number > probability[0])
		returned_index = min_idx + 1;
	free(probability);
	return(returned_index);
}

void split_probability(float *probabilities, int num_events, int min_idx, int *left_num, float **left_probabilities, int **left_indices, int *right_num, float **right_probabilities, int **right_indices)
{
	float	probability = 0.0;
	float sum_probability = 1.0;
	int factor = 2;
	float seprator = (float) sum_probability/factor;
	float epsilon = 1e-6;
	int left_max = 0;
	*left_num = 0;
	int right_min = 0;
    *right_num = 0;
	
	int j = 0;
	int k = 0;
	do
	{
		probability += probabilities[k];
		k++;
	}
	#if GCC_VERSION > GCC_GENERIC_FEATURE_SUPPORT
	    while(ffabs(probability - seprator) > epsilon && (probability < seprator));
	#else
		while(abs_float(probability - seprator) > epsilon && (probability < seprator));
	#endif
	
	right_min = k;
	left_max = right_min - 1;
	*left_num = left_max + 1;
	*right_num = num_events - *left_num;
	#if GCC_VERSION > GCC_GENERIC_FEATURE_SUPPORT
	    if(ffabs(probability - seprator) > epsilon)
	#else
		if(abs_float(probability - seprator) > epsilon)
	#endif
	{
		right_min--;
		*right_num = *right_num + 1;
	}	
	
	*left_indices = malloc(*left_num * sizeof(int));
	if(!*left_indices)
    {
        fprintf(stderr, "Error: split_probability()\n");
        exit(EXIT_FAILURE);
    }
	*left_probabilities = malloc(*left_num * sizeof(float));
	if(!*left_probabilities)
    {
        fprintf(stderr, "Error: split_probability()\n");
        exit(EXIT_FAILURE);
    }
	
	for(k = left_max; k >= *left_num - left_max - 1; k--)
	{
		j = *left_num - k - 1;
		(*left_indices)[j] = j + min_idx;
		(*left_probabilities)[j] = probabilities[j];
	}
	
	*right_indices = malloc(*right_num * sizeof(int));
	if(!*right_indices)
    {
        fprintf(stderr, "Error: split_probability()\n");
        exit(EXIT_FAILURE);
    }
	
	*right_probabilities = malloc(*right_num * sizeof(float));
	if(!*right_probabilities)
    {
        fprintf(stderr, "Error: split_probability()\n");
        exit(EXIT_FAILURE);
    }
	
	for(k = right_min; k < *right_num + right_min; k++)
	{
		j = k - right_min;
		(*right_indices)[j] = k + min_idx;
		(*right_probabilities)[j] = probabilities[k];
	}
	#if GCC_VERSION > GCC_GENERIC_FEATURE_SUPPORT
    	if(ffabs(probability - seprator) > epsilon)
	#else
		if(abs_float(probability - seprator) > epsilon)
	#endif
	    {
		#if GCC_VERSION > GCC_GENERIC_FEATURE_SUPPORT
	        (*left_probabilities)[*left_num - 1] = (*left_probabilities)[*left_num - 1] - (ffsum((*left_probabilities), *left_num) - seprator);
	        (*right_probabilities)[0] = (*right_probabilities)[0] - (ffsum((*right_probabilities), *right_num) - seprator);
		#else
			(*left_probabilities)[*left_num - 1] = (*left_probabilities)[*left_num - 1] - (sum_float((*left_probabilities), *left_num) - seprator);
	        (*right_probabilities)[0] = (*right_probabilities)[0] - (sum_float((*right_probabilities), *right_num) - seprator);
		#endif
	    }
	
	for(k = 0; k < *left_num; k++)
		(*left_probabilities)[k] *= factor;

	for(k = 0; k < *right_num; k++)
		(*right_probabilities)[k] *= factor;	
}

int trimmer_tool(float *probability, int num_event)
{
    int returned_index = -1;
    
    #if GCC_VERSION > GCC_GENERIC_FEATURE_SUPPORT
        float max_probability = ffmax(probability, num_event);
    #else
        float max_probability = max_value_float(probability, num_event);
    #endif
    
    if((1 - max_probability) < 0.001)
    {
        #if GCC_VERSION > GCC_GENERIC_FEATURE_SUPPORT
            returned_index = ffmaxi(probability, num_event, 0);     
        #else
            returned_index = max_index_float(probability, num_event, 0);     
        #endif
        return(returned_index);
    }
    
	int factor = 2;
	int k = 0;
	int num_events = num_event;
	float *probabilities = malloc(num_events * sizeof(float));
	if(!probabilities)
	{
		fprintf(stderr, "Error: trimmer_tool()\n");
		exit(EXIT_FAILURE);
	}
	memset(probabilities, 0, num_events * sizeof(float));
	memcpy(probabilities, probability, num_events * sizeof(float));
	
	int *indices = malloc(num_events * sizeof(int));
	if(!indices)
	{
		fprintf(stderr, "Error: trimmer_tool()\n");
		exit(EXIT_FAILURE);
	}
	memset(indices, 0, num_events * sizeof(int));
	for(k = 0; k < num_events; k++)
		indices[k] = k;
	
	float random_number = 0.0;
    
    while(num_events > factor)
	{
		random_number = (float) random()/RAND_MAX;
	    trim_probability(&probabilities, &indices, &num_events, random_number);
	}

	if(num_events == 1)
		returned_index = indices[0];
	else if(num_events == factor)
	{
		random_number = (float) random()/RAND_MAX;
		returned_index = (random_number < probabilities[0]) ? indices[0] : indices[1];
	}
	else
	{
		fprintf(stderr, "Error: trimmer_tool()\n");
		exit(EXIT_FAILURE);
	}
	free(probabilities);
	free(indices);
	return(returned_index);
	
}

void trim_probability(float **probabilities, int **indices, int *num_events, float random_number)
{
	
	float	probability = 0.0;	
	float sum_probability = 1.0;
	int factor = 2;
	float seprator = (float) sum_probability/factor;
	float epsilon = 1e-6;
	
	int left_max = 0;
	int left_num = 0;
	int right_min = 0;
    int right_num = 0;
	int min_idx = (*indices)[0];
	
	int j = 0;
	int k = 0;
	
	do
	{
		probability += (*probabilities)[k];
		k++;
	}
	#if GCC_VERSION > GCC_GENERIC_FEATURE_SUPPORT
	    while(ffabs(probability - seprator) > epsilon && (probability < seprator));
	#else
		while(abs_float(probability - seprator) > epsilon && (probability < seprator));
	#endif
	
	
	right_min = k;
	left_max = right_min - 1;
	left_num = left_max + 1;
	right_num = *num_events - left_num;
	#if GCC_VERSION > GCC_GENERIC_FEATURE_SUPPORT
	    if(ffabs(probability - seprator) > epsilon)
	#else
		if(abs_float(probability - seprator) > epsilon)
	#endif
	{
		right_min--;
		right_num = right_num + 1;
	}	
	
	if(random_number < seprator)
	{
		int *left_indices = malloc(left_num * sizeof(int));
	    if(!left_indices)
        {
            fprintf(stderr, "Error: trim_probability()\n");
            exit(EXIT_FAILURE);
        }
    	float *left_probabilities = malloc(left_num * sizeof(float));
	    if(!left_probabilities)
        {
            fprintf(stderr, "Error: trim_probability()\n");
            exit(EXIT_FAILURE);
        }
	
	    for(k = left_max; k >= left_num - left_max - 1; k--)
    	{
	    	j = left_num - k - 1;
		    left_indices[j] = j + min_idx;
    		left_probabilities[j] = (*probabilities)[j];
	    }
		
		#if GCC_VERSION > GCC_GENERIC_FEATURE_SUPPORT
    	    if(ffabs(probability - seprator) > epsilon)
	    #else
		    if(abs_float(probability - seprator) > epsilon)
	    #endif
	        {
		#if GCC_VERSION > GCC_GENERIC_FEATURE_SUPPORT
	            left_probabilities[left_num - 1] = left_probabilities[left_num - 1] - (ffsum(left_probabilities, left_num) - seprator);
		#else
			    left_probabilities[left_num - 1] = left_probabilities[left_num - 1] - (sum_float(left_probabilities, left_num) - seprator);
		#endif
	        }
	
	    for(k = 0; k < left_num; k++)
	    {
    		left_probabilities[k] *= factor;
	    }
		
		*num_events  = left_num;
		free(*probabilities);
	    free(*indices);		
	    *probabilities = malloc(*num_events * sizeof(float));
	    *indices = malloc(*num_events * sizeof(int));
		
		memcpy(*probabilities, left_probabilities, *num_events * sizeof(float));
		memcpy(*indices, left_indices, *num_events * sizeof(int));
		free(left_probabilities);
		free(left_indices);
	}
	else
	{	
	    int *right_indices = malloc(right_num * sizeof(int));
	    if(!right_indices)
        {
            fprintf(stderr, "Error: trim_probability()\n");
            exit(EXIT_FAILURE);
        }
	
	    float *right_probabilities = malloc(right_num * sizeof(float));
	    if(!right_probabilities)
        {
            fprintf(stderr, "Error: trim_probability()\n");
            exit(EXIT_FAILURE);
        }
	
	    for(k = right_min; k < right_num + right_min; k++)
	    {
		    j = k - right_min;
		    right_indices[j] = k + min_idx;
		    right_probabilities[j] = (*probabilities)[k];
	    }
	    #if GCC_VERSION > GCC_GENERIC_FEATURE_SUPPORT
    	    if(ffabs(probability - seprator) > epsilon)
	    #else
		    if(abs_float(probability - seprator) > epsilon)
	    #endif
            {
		#if GCC_VERSION > GCC_GENERIC_FEATURE_SUPPORT
	            right_probabilities[0] = right_probabilities[0] - (ffsum(right_probabilities, right_num) - seprator);
		#else
	            right_probabilities[0] = right_probabilities[0] - (sum_float(right_probabilities, right_num) - seprator);
		#endif
            }

	    for(k = 0; k < right_num; k++)
		    right_probabilities[k] *= factor;	
		
		*num_events = right_num;
		free(*probabilities);
	    free(*indices);		
	    *probabilities = malloc(*num_events * sizeof(float));
	    *indices = malloc(*num_events * sizeof(int));
		
		memcpy(*probabilities, right_probabilities, *num_events * sizeof(float));
		memcpy(*indices, right_indices, *num_events * sizeof(int));
		free(right_probabilities);
		free(right_indices);
	}
}
