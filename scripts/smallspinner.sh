#!/bin/bash
# smallspinner.sh: Single node spinner
# (c) Mohammad HMofrad, 2017
# (e) mohammad.hmofrad@pitt.edu

# echo on
# set -x
# set -v


MASTER="delorean.cs.pitt.edu";
SLAVE="delorean.cs.pitt.edu";
LOCAL_PATH="/home/moh18";
OKAPI_JAR_FILE="${LOCAL_PATH}/distrograph/okapi/okapi-master/target/okapi-0.3.5-SNAPSHOT-jar-with-dependencies.jar";

MASTER_PATH="${LOCAL_PATH}/distrograph/hadoop/hadoop-1.2.1/bin";
LOG_PATH="${LOCAL_PATH}/distrograph/hadoop/hadoop-1.2.1/logs/userlogs";

GRAPH="facebook_combined.txt";
#GRAPH="soc-LiveJournal1.txt"
#GRAPH="hollywood-2011.txt"
#GRAPH="eu-2015-host.txt"
#GRAPH="com-orkut.ungraph.txt";

DISK_BASE="/home/moh18/distrograph";
DISK_OUTPUT_PATH="${DISK_BASE}/outputs";
HD_BASE="/graphs";
HD_INPUT_PATH="${HD_BASE}/${GRAPH}";

if [ "$#" -ne 5 ]
then
    echo "./smallspinner -[P|R] NP NR INPUT_PARTITIONING HD_OUTPUT_PREFIX"
    echo "-[P|R]: Partition|Repartition"
    echo "NP: NumPartitions"
    echo "NR: NumRepartitions"
    echo "E.g: ./smallspinner.sh -P 8 0 0 9"
    echo "E.g: ./smallspinner.sh -R 8 16 /output/22/part-m-00001 9"
    exit;    
fi

P=$2;
R=$3;
INPUT_PARTITIONING=$4;
N=$5;

HD_OUTPUT_PATH="/output/${N}";
LOG_FILE="${DISK_OUTPUT_PATH}/${P}_spinner_logs_${GRAPH}";
PART_FILE="${DISK_OUTPUT_PATH}/${P}_spinner_raw_${GRAPH}";
SORTED_FILE="${DISK_OUTPUT_PATH}/${P}_spinner_parts_${GRAPH}";

if [ $1 = "-P" ]; then
        hadoop jar ${OKAPI_JAR_FILE} org.apache.giraph.GiraphRunner ml.grafos.okapi.spinner.Spinner\$ConverterPropagate \
            -mc ml.grafos.okapi.spinner.Spinner\$PartitionerMasterCompute \
            -eip ${HD_INPUT_PATH} \
            -eif ml.grafos.okapi.spinner.Spinner\$SpinnerEdgeInputFormat \
            -op ${HD_OUTPUT_PATH} \
            -vof ml.grafos.okapi.spinner.Spinner\$SpinnerVertexValueOutputFormat \
            -w 1 \
            -ca spinner.numberOfPartitions=$P \
            -ca giraph.outEdgesClass=ml.grafos.okapi.spinner.OpenHashMapEdges;
    wait;
    hadoop dfs -copyToLocal ${HD_OUTPUT_PATH}/"part-m-00001" ${DISK_OUTPUT_PATH}/part_${P};
elif [ $1 = "-R" ]; then
#echo "$0 $1 $2 $3 $4 $5"
    hadoop jar ${OKAPI_JAR_FILE} org.apache.giraph.GiraphRunner ml.grafos.okapi.spinner.Spinner\$ConverterPropagate \
        -mc ml.grafos.okapi.spinner.Spinner\$PartitionerMasterCompute \
        -eip ${HD_INPUT_PATH} \
        -eif ml.grafos.okapi.spinner.Spinner\$SpinnerEdgeInputFormat \
        -op ${HD_OUTPUT_PATH} \
        -vof ml.grafos.okapi.spinner.Spinner\$SpinnerVertexValueOutputFormat \
        -w 1 \
        -vif ml.grafos.okapi.spinner.Spinner\$SpinnerVertexValueInputFormat \
        -vip ${INPUT_PARTITIONING} \
        -ca spinner.numberOfPartitions=$P \
        -ca spinner.repartition=$R \
        -ca giraph.outEdgesClass=ml.grafos.okapi.spinner.OpenHashMapEdges;
    wait;
    hadoop dfs -copyToLocal ${HD_OUTPUT_PATH}/"part-m-00001" ${DISK_OUTPUT_PATH}/part_${P}_repart_${R};        
fi
tail -n 25 $(\ls -1dt ${LOG_PATH}/*/ | head -n 1)/attemp*/stdout;
