#!/bin/bash
# spinner.sh: Spinner run script
# (c) Mohammad HMofrad, 2017
# (e) mohammad.hmofrad@pitt.edu

# echo on
#set -x
#set -v

if [ "$#" -ne 1 ]
then
    echo "$0  HD_OUTPUT_PREFIX";
    exit;
fi

#MASTER="delorean.cs.pitt.edu";
#SLAVE="delorean.cs.pitt.edu";
#NODES="delorean.cs.pitt.edu";
#BASE=/home/moh18
#HADOOP_LOCAL_PATH="${BASE}/distrograph/hadoop/hadoop-1.2.1";
#OKAPI_JAR_FILE="${BASE}/distrograph/okapi/okapi-master/target/okapi-0.3.5-SNAPSHOT-jar-with-dependencies.jar";

#MASTER="opa-n2.sam.pitt.edu";
#SLAVES=("opa-n2.sam.pitt.edu" "opa-n3.sam.pitt.edu" "opa-n4.sam.pitt.edu" "opa-n5.sam.pitt.edu");
#NODES="opa-n2.sam.pitt.edu opa-n3.sam.pitt.edu opa-n4.sam.pitt.edu opa-n5.sam.pitt.edu";


BASE="/ihome/rmelhem/moh18";
OKAPI_JAR_FILE="${BASE}/distrograph/okapi/okapi-master/target/okapi-0.3.5-SNAPSHOT-jar-with-dependencies.jar";
LOCAL_PATH="/tmp";
SLAVES_PATH="${LOCAL_PATH}/hadoop";
SLAVES_CONFIG="${SLAVES_PATH}/conf/slaves";
NUM_WORKERS=`cat ${SLAVES_CONFIG} | wc -l`; # Master is an slave too

#GRAPH="soc-LiveJournal1.txt"
#GRAPH="hollywood-2011.txt"
#GRAPH="eu-2015-host.txt"
#GRAPH="com-orkut.ungraph.txt";
#GRAPHS="facebook_combined.txt";
GRAPHS="usa-roads.txt";
HD_BASE="/graphs";
HD_INPUT_PATH="${HD_BASE}/${GRAPH}";
DISK_BASE=${BASE}/distrograph;
DISK_OUTPUT_PATH="${DISK_BASE}/outputs";
N=$1;
PARTITIONS="8";
#PARTITIONS="2 4 8 16 32 64 128 192 256";
for G in ${GRAPHS}    
do
    for P in ${PARTITIONS}
    do
        HD_OUTPUT_PATH="/output/$N";
        LOG_FILE="${DISK_OUTPUT_PATH}/${P}_logs_${G}";
        PART_FILE="${DISK_OUTPUT_PATH}/${P}_raw_${G}";
        SORTED_FILE="${DISK_OUTPUT_PATH}/${P}_parts_${G}";
            { time hadoop jar ${OKAPI_JAR_FILE} org.apache.giraph.GiraphRunner ml.grafos.okapi.spinner.Spinner\$ConverterPropagate \
                -mc ml.grafos.okapi.spinner.Spinner\$PartitionerMasterCompute \
                -eip ${HD_INPUT_PATH} \
                -eif ml.grafos.okapi.spinner.Spinner\$SpinnerEdgeInputFormat \
                -op ${HD_OUTPUT_PATH} \
                -vof ml.grafos.okapi.spinner.Spinner\$SpinnerVertexValueOutputFormat \
                -w ${NUM_WORKERS} \
                -ca spinner.numberOfPartitions=$P \
                -ca giraph.outEdgesClass=ml.grafos.okapi.spinner.OpenHashMapEdges; } 2>&1 | tee $LOG_FILE
                #{ time cmd } 2> $LOG_FILE
        wait;
    
        for (( W=1; W<=${NUM_WORKERS}; W++ ))
            do
            if  [ "$W" -lt 10 ]
            then
                hadoop dfs -copyToLocal ${HD_OUTPUT_PATH}/"part-m-0000${W}" ${DISK_OUTPUT_PATH}/part-m-0000${W};
            else
                hadoop dfs -copyToLocal ${HD_OUTPUT_PATH}/"part-m-000${W}" ${DISK_OUTPUT_PATH}/part-m-0000${W};
            fi
        done
    
        wait;
        cat ${DISK_OUTPUT_PATH}/part-m-* > ${PART_FILE};
        rm -rf ${DISK_OUTPUT_PATH}/part-m-*;
        sort -k 1 -n ${PART_FILE} > ${SORTED_FILE};
        wait;
        N=$((N+1))
    done
done

