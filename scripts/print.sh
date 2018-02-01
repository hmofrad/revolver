#!/bin/bash
# print.sh: Print script
# (c) Mohammad HMofrad, 2017
# (e) mohammad.hmofrad@pitt.edu

# echo on
#set -x
#set -v

LOCAL_PATH="/tmp";
SLAVES_PATH="${LOCAL_PATH}/hadoop";
SLAVES_CONFIG="${SLAVES_PATH}/conf/slaves";
LOG_PATH="${SLAVES_PATH}/logs/userlogs";
USERNAME="moh18";
declare -a SLAVES
let i=0
while IFS=$'\n' read -r line
do
    SLAVES[i]="${line}"
    ((++i))
done < ${SLAVES_CONFIG}

for ((i=0; i < ${#SLAVES[*]}; i++))
do
    ssh -l ${USERNAME} ${SLAVES[i]} "hostname; tail -n 100  \$(\ls -1dt ${LOG_PATH}/*/ | head -n 1)/attemp*/stdout; echo; exit;";
done