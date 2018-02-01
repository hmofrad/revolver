/*
 * libmath.h Basic math header file
 * (c) Mohammad H. Mofrad, 2017
 * (e) mohammad.hmofrad@pitt.edu
 */

#ifndef LIBMATH_H
#define LIBMATH_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//#include <math.h>
#include "macros.h"

int    sum_int(int *array, int count);
unsigned int sum_uint(unsigned int *array, int count);
float sum_float(float *array, int count);
double sum_double(double *array, int count);
long long sum_llong(long long *array, int count);

		   
int    max_index_int(int *array, int count, int dummy);
int    max_index_uint(unsigned int *array, int count, int dummy);
int    max_index_float(float *array, int count, int dummy);
int    max_index_double(double *array, int count, int dummy);
int    max_index_llong(long long *array, int count, int dummy);
int   *max_indices_float(float *array, int count, int *num_indices);
int   *max_indices_double(double *array, int count, int *num_indices);
							  
int    max_value_int(int *array, int count);
unsigned int max_value_uint(unsigned int *array, int count);
float max_value_float(float *array, int count);
double max_value_double(double *array, int count);
long long max_value_llong(long long *array, int count);

int min_index_int(int *array, int count, int dummy);
int min_index_uint(unsigned int *array, int count, int dummy);
int min_index_float(float *array, int count, int dummy);
int min_index_double(double *array, int count, int dummy);
int min_index_llong(long long *array, int count, int dummy);

int min_value_int(int *array, int count);
unsigned int min_value_uint(unsigned int *array, int count);
float min_value_float(float *array, int count);
double min_value_double(double *array, int count);
long long min_value_llong(long long *array, int count);

float mean_int(int *array, int count);
float mean_uint(unsigned int *array, int count);
float mean_float(float *array, int count);
double mean_double(double *array, int count);
double mean_llong(long long *array, int count);

int abs_int(int number);
unsigned int abs_uint(unsigned int number);
float abs_float(float number);
double abs_double(double number);
long double abs_ldouble(long double number);

void   ishuffle(int *array, int count);
							   
int    random_gen(int min_num, int max_num);

int    random_at_most(int max);

int power_int(int base, int exp);
double power_double(double base, double exp);

void bubble_sort(float *values, int *indices, int count);

/* Test for GCC > 3.2.0 */
    #if GCC_VERSION > GCC_GENERIC_FEATURE_SUPPORT
        #define ffsum(_1, _2) _Generic((_1),                  \
                                           int *: sum_int,     \
                                  unsigned int *: sum_uint,     \
                                         float *: sum_float,     \
                                        double *: sum_double,     \
                                     long long *: sum_llong)(_1, _2)

       #define ffmaxi(_1, _2, _3) _Generic((_3),                    \
                                         int *: *max_indices_double, \
                                         int  : _Generic((_1),        \
			                             int *: max_index_int,         \
                                unsigned int *: max_index_uint,         \
                                       float *: max_index_float,         \
			                          double *: max_index_double,         \
                                   long long *: max_index_llong)(_1, _2, _3)   

        #define ffmax(_1, _2) _Generic((_1),       \
                       int    *: max_value_int,     \
                 unsigned int *: max_value_uint,     \
                        float *: max_value_float,     \
                       double *: max_value_double,     \
                    long long *: max_value_llong)(_1, _2)									 

        #define ffmini(_1, _2, _3) _Generic((_1),          \
                              int *: min_index_int,         \
                     unsigned int *: min_index_uint,         \
                            float *: min_index_float,         \
                           double *: min_index_double,         \
                        long long *: min_index_llong)(_1, _2, _3)
                    
        #define ffmin(_1, _2) _Generic((_1),           \
                              int *: min_value_int,     \
                     unsigned int *: min_value_uint,     \
                            float *: min_value_float,     \
                           double *: min_value_double,     \
                        long long *: min_value_llong)(_1, _2)

        #define ffmean(_1, _2) _Generic((_1),     \
                              int *: mean_int,     \
                     unsigned int *: mean_uint,     \
                            float *: mean_float,     \
                           double *: mean_double,     \
                        long long *: mean_llong)(_1, _2)	

	   #define ffabs(_1) _Generic((_1),                  \
                                         int : abs_int,   \
                                unsigned int : abs_uint,   \
                                       float : abs_float,   \
                                      double : abs_double,   \
						         long double : abs_ldouble)(_1)		
						   
				     
        #define ffshuffle(_1, _2) _Generic((_1),                \
                                       int *   : ishuffle)(_1, _2)

        #define ffrandom(_1) _Generic((_1),                \
                                    int : random_at_most)(_1)
                                     
        #define ffrandom_(_1, _2) _Generic((_1),           \
                                    int : random_gen)(_1, _2)
                                     
        #define ffexp(_1, _2) _Generic((_1),               \
                                     int : power_int)(_1, _2)                                     
                                           
    #endif								   
				   
#endif


