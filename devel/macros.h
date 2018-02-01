/*
 * macros.h Macros header file
 * (c) Mohammad H. Mofrad, 2017
 * (e) mohammad.hmofrad@pitt.edu
 */
 
#ifndef MACROS_H
#define MACROS_H

#define _GNU_SOURCE
 
#define VERBOSE
#define ALL
#define FILEIO

#define NUM_CORES (int)(sysconf(_SC_NPROCESSORS_ONLN))
#define MAX_THREADS_NUM (2)
#define MAX_ITER (290)
#define MAX_PATH_LEN (1000)
#define MAX_LONG_LONG (9223372036854775807)
#define MAX_UNISIGNED_INT (4294967295)

#define MAX_PROBABILITY (.9)
//#include <unistd.h>
//static inline int num_cores() { return((int) sysconf(_SC_NPROCESSORS_ONLN)); }

#define EXIT_ERROR (-1)
#define EXIT_OK    (0)

#define GCC_VERSION (__GNUC__ * 10000      \
                     + __GNUC_MINOR__ * 100 \
                     + __GNUC_PATCHLEVEL__   )
#define GCC_GENERIC_FEATURE_SUPPORT 60900
#endif

