#!/bin/bash
# okapi.sh: Okapi compile script
# ./okapi.sh
# (c) Mohammad HMofrad, 2017
# (e) mohammad.hmofrad@pitt.edu

#set -x
#set -v

if [ "$#" -ne 1 ]; then
    echo "$0 <revolver|spinner|hash|test>"
    exit;
fi

#module load maven
# unzip BASE/okapi/okapi-master.zip
#BASE="/ihome/cs3580_2017F/moh18";
#BASE="/afs/cs.pitt.edu/usr0/hasanzadeh/private";
BASE="/home/moh18";
BASE="/ihome/rmelhem/moh18";
OKAPI="${BASE}/distrograph/okapi/okapi-master";
OKAPI_JAR_FILE="${OKAPI}/target/okapi-0.3.5-SNAPSHOT-jar-with-dependencies.jar"



REVOLVER_SOURCE_FILE="${BASE}/distrograph/spinner/Revolver.java";
SPINNER_SOURCE_FILE="${BASE}/distrograph/spinner/Spinner.java";
HASH_SOURCE_FILE="${BASE}/distrograph/spinner/Hash.java";
TEST_SOURCE_FILE="${BASE}/distrograph/spinner/Test.java";
PAGERANK_SOURCE_FILE="${BASE}/distrograph/spinner/PageRank.java";

SPINNER_DESTINATION_FILE="${OKAPI}/src/main/java/ml/grafos/okapi/spinner/Spinner.java";

if   [ $1 = "revolver" ]; then
    cp ${REVOLVER_SOURCE_FILE} ${SPINNER_DESTINATION_FILE}
elif [ $1 = "spinner" ]; then 
    cp ${SPINNER_SOURCE_FILE} ${SPINNER_DESTINATION_FILE}
elif [ $1 = "hash" ]; then 
    cp ${HASH_SOURCE_FILE} ${SPINNER_DESTINATION_FILE}
elif [ $1 = "pagerank" ]; then 
    cp ${PAGERANK_SOURCE_FILE} ${SPINNER_DESTINATION_FILE}
elif [ $1 = "test" ]; then 
    cp ${TEST_SOURCE_FILE} ${SPINNER_DESTINATION_FILE}    
else
    exit;
fi

cd ${OKAPI};
mvn clean;
mvn package -DskipTests;

if [ ! -d ${OKAPI_JAR_FILE} ]
    then
        echo "${OKAPI_JAR_FILE} successfully created.";
    else
        echo "Build failed.";
        exit;
fi

USERNAME="moh18";
LOCAL_PATH="/tmp";
SLAVES_PATH="${LOCAL_PATH}/hadoop";
SLAVES_CONFIG="${SLAVES_PATH}/conf/slaves";

declare -a SLAVES
let i=0
while IFS=$'\n' read -r line
do
    SLAVES[i]="${line}"
    ((++i))
done < ${SLAVES_CONFIG}

for ((i=0; i < ${#SLAVES[*]}; i++))
do
    echo "${SLAVES[i]}: scping ${SLAVES_PATH} ... ";
    (scp ${OKAPI_JAR_FILE} ${USERNAME}@${SLAVES[i]}:${SLAVES_PATH};) & # Master is an slave too
    wait;
done
