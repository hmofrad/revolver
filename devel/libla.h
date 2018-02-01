/*
 * libla.h Learning Automata (LA) header file
 * (c) Mohammad H. Mofrad, 2017
 * (e) mohammad.hmofrad@pitt.edu
 */
 
#ifndef LIBLA_H
#define LIBLA_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "libmath.h"
#include "libgraph.h"


void action_selection(unsigned int low, unsigned int high);
void probability_update(unsigned int low, unsigned int high);
void weighted_probability_update(unsigned int low, unsigned int high, int iteration);
void reinforcement_signal(unsigned int low, unsigned int high);
int roulette_wheel(float *probabilities, int num_events);
int ruler_tool(float *probabilities, int num_events);
int scissors_tool(float *probabilities, int num_events);
void split_probability(float *probabilities, int num_events, int min_idx, int *left_num, float **left_probabilities, int **left_indices, int *right_num, float **right_probabilities, int **right_indices);
int trimmer_tool(float *probabilities, int num_events);
void trim_probability(float **probabilities, int **indices, int *num_events, float random_number);

#endif


