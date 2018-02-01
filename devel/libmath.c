/*
 * libmath.c Basic math library
 * (c) Mohammad H. Mofrad, 2017
 * (e) mohammad.hmofrad@pitt.edu
 */
 
#include "libmath.h"
/* int */
int sum_int(int *array, int count)
{
	int i = 0;
	int sum = 0;
	for(i = 0; i < count; i++)
		sum += array[i];
	return(sum);
}

float mean_int(int *array, int count)
{
	return(sum_int(array, count)/(1.0 * count));
}

int max_index_int(int *array, int count, int dummy)
{
	int i = 0;
	int idx = 0;
	for(i = 1; i < count; i++)
	{
		if(array[i] > array[idx])
			idx = i;
	}
	return(idx);
}

int max_value_int(int *array, int count)
{
    return(array[max_index_int(array, count, 0)]);
}

int min_index_int(int *array, int count, int dummy)
{
	int i = 0;
	int idx = 0;
	for(i = 1; i < count; i++)
	{
		if(array[i] < array[idx])
			idx = i;
	}
	return(idx);
}

int min_value_int(int *array, int count)
{
    return(array[min_index_int(array, count, 0)]);
}

int abs_int(int number)
{
	return((int) abs_ldouble((long double) number));
}

// unsigned int
unsigned int sum_uint(unsigned int *array, int count)
{
	int i = 0;
	unsigned int sum = 0;
	for(i = 0; i < count; i++)
		sum += array[i];
	return(sum);
}

float mean_uint(unsigned int *array, int count)
{
	return(sum_uint(array, count)/(1.0 * count));
}

int max_index_uint(unsigned int *array, int count, int dummy)
{
	int i = 0;
	int idx = 0;
	for(i = 1; i < count; i++)
	{
		if(array[i] > array[idx])
			idx = i;
	}
	return(idx);
}

unsigned int max_value_uint(unsigned int *array, int count)
{
    return(array[max_index_uint(array, count, 0)]);
}

int min_index_uint(unsigned int *array, int count, int dummy)
{
	int i = 0;
	int idx = 0;
	for(i = 1; i < count; i++)
	{
		if(array[i] < array[idx])
			idx = i;
	}
	return(idx);
}

unsigned int min_value_uint(unsigned int *array, int count)
{
    return(array[min_index_uint(array, count, 0)]);
}

unsigned int abs_uint(unsigned int number) // It doesn't make sense, yet I'm adding it :-/
{
    unsigned int absolute_value = number;
    if (absolute_value < 0)
		exit(0);
    else
        ;
	return(absolute_value);
}

// float
float sum_float(float *array, int count)
{
	int i = 0;
	float sum = .0;
	for(i = 0; i < count; i++)
		sum += array[i];
	return(sum);
}

float mean_float(float *array, int count)
{
	return(sum_float(array, count)/(1.0 * count));
}

int max_index_float(float *array, int count, int dummy)
{
	int i = 0;
	int idx = 0;
	for(i = 1; i < count; i++)
	{
		if(array[i] > array[idx])
			idx = i;
	}
	return(idx);
}

float max_value_float(float *array, int count)
{
	return(array[max_index_float(array, count, 0)]);
}

int *max_indices_float(float *array, int count, int *num_indices)
{
	int i = 0;
	int j = 0;
	float max_value = max_value_float(array, count);
	int indices[count];
	memset(indices, 0, sizeof(indices));

	for(j = 0; j < count; j++)
	{
		if(max_value == array[j])
		{
			indices[j] = 1;
			*num_indices = *num_indices + 1;
		}
	}

	int *indices_array = malloc(*num_indices * sizeof(int));
	for(j = 0; j < count; j++)
	{
		if(indices[j])
		{
		    indices_array[i] = j;
			i++;
		}
	}
	return(indices_array);
}

int min_index_float(float *array, int count, int dummy)
{
	int i = 0;
	int idx = 0;
	for(i = 1; i < count; i++)
	{
		if(array[i] < array[idx])
			idx = i;
	}
	return(idx);
}

float min_value_float(float *array, int count)
{
	return(array[min_index_float(array, count, 0)]);
}

float abs_float(float number)
{
    float absolute_value = number;
    if (absolute_value < 0)
		absolute_value = -absolute_value;
	return(absolute_value);
}

// double
double sum_double(double *array, int count)
{
	int i = 0;
	double sum = .0;
	for(i = 0; i < count; i++)
		sum += array[i];
	return(sum);
}

double mean_double(double *array, int count)
{
	return(sum_double(array, count)/(1.0 * count));
}

int max_index_double(double *array, int count, int dummy)
{
	int i = 0;
	int idx = 0;
	for(i = 1; i < count; i++)
	{
		if(array[i] > array[idx])
			idx = i;
	}
	return(idx);
}

double max_value_double(double *array, int count)
{
	return(array[max_index_double(array, count, 0)]);
}

int *max_indices_double(double *array, int count, int *num_indices)
{
	int i = 0;
	int j = 0;
	double max_value = max_value_double(array, count);
	int indices[count];
	memset(indices, 0, sizeof(indices));

	for(j = 0; j < count; j++)
	{
		if(max_value == array[j])
		{
			indices[j] = 1;
			*num_indices = *num_indices + 1;
		}
	}

	int *indices_array = malloc(*num_indices * sizeof(int));
	for(j = 0; j < count; j++)
	{
		if(indices[j])
		{
		    indices_array[i] = j;
			i++;
		}
	}
	return(indices_array);
}

int min_index_double(double *array, int count, int dummy)
{
	int i = 0;
	int idx = 0;
	for(i = 1; i < count; i++)
	{
		if(array[i] < array[idx])
			idx = i;
	}
	return(idx);
}

double min_value_double(double *array, int count)
{
	return(array[min_index_double(array, count, 0)]);
}

double abs_double(double number)
{
	return((double) abs_ldouble((long double) number));
}

// long long
long long sum_llong(long long *array, int count)
{
	int i = 0;
	long long sum = 0;
	for(i = 0; i < count; i++)
		sum += array[i];
	return(sum);
}

double mean_llong(long long *array, int count)
{
	return(sum_llong(array, count)/(1.0 * count));
}

int max_index_llong(long long *array, int count, int dummy)
{
	int i = 0;
	int idx = 0;
	for(i = 1; i < count; i++)
	{
		if(array[i] > array[idx])
			idx = i;
	}
	return(idx);
}

long long max_value_llong(long long *array, int count)
{
	return(array[max_index_llong(array, count, 0)]);
}


int min_index_llong(long long *array, int count, int dummy)
{
	int i = 0;
	int idx = 0;
	for(i = 1; i < count; i++)
	{
		if(array[i] < array[idx])
			idx = i;
	}
	return(idx);
}

long long min_value_llong(long long *array, int count)
{
    return(array[min_index_llong(array, count, 0)]);
}

long double abs_ldouble(long double number)
{
	long double absolute_value = number;
    if (absolute_value < 0)
		absolute_value = -number;
	return(absolute_value);
}

/*
long double sum_ldouble(long double *array, int count)
{
	int i = 0;
	long double sum = .0;
	for(i = 0; i < count; i++)
		sum += array[i];
	return(sum);
}
*/

void ishuffle(int *array, int count)
{
	int i = 0;
	int j = 0;
	int t = 0;
	
    if (count > 1) 
    {
        for(i = 0; i < count - 1; i++) 
        {
          j = i + random() / (RAND_MAX / (count - i) + 1);
          t = array[j];
          array[j] = array[i];
          array[i] = t;
        }
    }
}

int random_gen(int min_num, int max_num)
{
	int random_number = min_num + random() / (RAND_MAX / (max_num - min_num + 1) + 1);
    if((random_number < min_num) && (random_number >= max_num))
    {
        fprintf(stderr, "Error: random_gen():\n");
        fprintf(stderr, "                     random_number = %d\n", random_number);
        fprintf(stderr, "                     min_num = %d\n", min_num);
        fprintf(stderr, "                     max_num = %d\n", max_num);
        exit(EXIT_ERROR);
    }

    return(random_number);
}

// https://stackoverflow.com/a/822361/5412470
// Assumes 0 <= max <= RAND_MAX
// Returns in the closed interval [0, max]
int random_at_most(int max)
{
    // max <= RAND_MAX < ULONG_MAX, so this is okay.
    unsigned int num_bins = (unsigned long) max + 1;
    unsigned int num_rand = (unsigned long) RAND_MAX + 1;
    unsigned int bin_size = num_rand / num_bins;
    unsigned int defect   = num_rand % num_bins;

    int x;
    do
    {
        x = random();
    }
    // This is carefully written not to overflow
    while(num_rand - defect <= (unsigned int)x);

    // Truncated division is intentional
    return(x/bin_size);
}

// https://stackoverflow.com/a/213322/5412470
int power_int(int base, int exp)
{
    if (exp == 0)
        return 1;
    else if (exp % 2)
        return base * power_int(base, exp - 1);
    else
    {
        int temp = power_int(base, exp / 2);
        return temp * temp;
    }
}

void bubble_sort(float *values, int *indices, int count)
{
    int k = 0;
    int l = 0;
    float tmp_value = 0.0;
    int    tmp_index = 0;
    for(k = 0; k < count; k++)
    {
        for(l = k + 1; l < count; l++)
        {
            if(values[k] > values[l])
            {
                tmp_value = values[k];
                values[k] = values[l];
                values[l] = tmp_value;
                tmp_index = indices[k];
                indices[k] = indices[l];
                indices[l] = tmp_index;
            }
        }
    }
}

