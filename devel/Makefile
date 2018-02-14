#!/usr/bin/make
# Makefile: makefile for Spinner and Revolver algorithms
# (c) Mohammad H. Mofrad, 2017
# (e) mohammad.hmofrad@pitt.edu

CC     := gcc
CFLAGS := -g -Wall -Werror
PFLAGS := -pthread
LDFLAGS:= -lm

ANSI_RED   := \033[0;31m
ANSI_GREEN := \033[0;32m
ANSI_RESET := \033[0m

PRINT_STATS_OBJ := print_stats.o
MAIN_OBJ        := main.o

PRINT_STATS := print_stats
MAIN        := main
COMMON_LIB  := libgraph.a

LIB_INCLUDES :=   \
		libgraph.h \
		libla.h     \
		libmath.h    \

LIB_OBJECTS :=    \
		libgraph.o \
		libla.o     \
		libmath.o    \

ALG_OBJECTS :=         \
		print_stats.o   \
		main.o           \

.PHONY: all clean caption
			
all: caption $(ALG_OBJECTS) $(COMMON_LIB)
	@echo "[building]"
	$(CC) $(CFLAGS) $(PFLAGS) -o $(PRINT_STATS) $(PRINT_STATS_OBJ) $(COMMON_LIB) $(LDFLAGS)
	$(CC) $(CFLAGS) $(PFLAGS) -o $(PRINT_STATS) $(PRINT_STATS_OBJ) $(COMMON_LIB) $(LDFLAGS)
	$(CC) $(CFLAGS) $(PFLAGS) -o $(MAIN)        $(MAIN_OBJ)        $(COMMON_LIB) $(LDFLAGS)
	@echo "$(ANSI_GREEN)Revolver compiled successfully$(ANSI_RESET)"

caption:
	@echo "[Compiling]"

$(COMMON_LIB): $(LIB_OBJECTS)
	@echo "[Archiving]"
	ar rcs $(COMMON_LIB) $(LIB_OBJECTS)	
	
$(LIB_OBJECTS): %.o: %.c
	$(CC) $(CFLAGS) $(PFLAGS) -o $@ -c $< $(LDFLAGS)
	
$(ALG_OBJECTS): %.o: %.c
	$(CC) $(CFLAGS) $(PFLAGS) -o $@ -c $< $(LDFLAGS)	
	
clean:
	@echo "[Deleting]"
	$(RM) $(LIB_OBJECTS)
	$(RM) $(COMMON_LIB)
	$(RM) $(MAIN_OBJ)
	$(RM) $(MAIN)
	$(RM) $(PRINT_STATS_OBJ)
	$(RM) $(PRINT_STATS)
