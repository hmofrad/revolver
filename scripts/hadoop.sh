#!/bin/bash
#   hadoop.sh: Configure a hadoop cluster
#                  given masters and slaves
#   (c) Mohammad HMofrad, 2017
#   (e) mohammad.hmofrad@pitt.edu

#module load java;
#echo on
#set -x
#set -v

if [ "$#" -ne 1 ]
then
    echo "$0  --[format-all|deploy-all]";
    exit;
fi

if [ `jps | wc -l` != 1 ]
then
    echo "Stopping $MASTER_PATH ... ";
    stop-all.sh;
    wait;
else
    echo "No running hadoop process ...";
fi

NFS_STORAGE="/ihome/rmelhem/moh18"
ZFX_STORAGE="/zfs1/cs3580_2017F/moh18";
STORAGE=${NFS_STORAGE};
#tar xzvf distrograph/hadoop/hadoop-1.2.1.tar.gz
HADOOP_PATH=${STORAGE}/distrograph/hadoop/hadoop-1.2.1;

SRC_CONFIG="${STORAGE}/distrograph/hadoop/conf"
MASTERS_SRC_CONFIG="${SRC_CONFIG}/masters";
SLAVES_SRC_CONFIG="${SRC_CONFIG}/slaves";

declare -a SLAVES
let i=0
while IFS=$'\n' read -r line
do
    #SLAVES[i]="${line}"
    SLAVES[i]="${line::-1}" # if \r
    ((++i))
done < ${SLAVES_SRC_CONFIG}

IFS=$'\n' read -d '' -r -a lines < ${MASTERS_SRC_CONFIG}
MASTER="${lines}"

for ((i=0; i < ${#SLAVES[*]}; i++))
do
    echo "${SLAVES[i]}"
done
echo "${MASTER}"

#MASTERR="opa-n2.sam.pitt.edu";
#SLAVESS="opa-n2.sam.pitt.edu opa-n3.sam.pitt.edu opa-n4.sam.pitt.edu opa-n5.sam.pitt.edu";
LOCAL_PATH="/tmp"
CONFIG_PATH="${NFS_STORAGE}/distrograph/hadoop/conf";
MASTER_PATH="${LOCAL_PATH}/hadoop";
SLAVE_PATH="${LOCAL_PATH}/hadoop";

USERNAME="moh18"

#cd ~
#ssh-keygen -t rsa -P ""
#cat $HOME/.ssh/id_rsa.pub >> $HOME/.ssh/authorized_keys
if [ "$1" = "--ssh-copy-id" ]
then
    for ((i=0; i < ${#SLAVES[*]}; i++))
    do
        echo ${USERNAME}@${SLAVES[i]};
        if [[ ("${SLAVES[i]}" != ${MASTER}) && (! -f $HOME/.ssh/authorized_keys) ]]
        then
             echo "${MASTER}: Generating ssh key...";
             ssh-keygen -t rsa -P ""
             cat $HOME/.ssh/id_rsa.pub >> $HOME/.ssh/authorized_keys
        else
            echo "${SLAVES[i]}: Copying ssh key...";
            ssh-copy-id ${USERNAME}@${SLAVES[i]};
            #ssh-copy-id -i $HOME/.ssh/id_rsa.pub ${USERNAME}@${SLAVES[i]};
        fi
        wait;
    done
    exit;
fi

if [ "$1" = "--format-all" ]
then
    echo "Master-${MASTER}: Formatting $MASTER_PATH ... ";
    rm -rf ${MASTER_PATH}
    wait;
    for ((i=0; i < ${#SLAVES[*]}; i++))
    do
        if [ "${SLAVES[i]}" != ${MASTER} ]
        then   
            echo "Slave-${SLAVES[i]}: Formatting $SLAVE_PATH ... ";
            (ssh -l ${USERNAME} ${SLAVES[i]} "rm -rf $SLAVE_PATH;") &
        fi
        wait;
    done
    exit;
fi

if [ "$1" = "--deploy-all" ]
then
    if [ -d ${MASTER_PATH} ]
    then
        echo "$MASTER_PATH exists";
    else
        echo "Copying $MASTER_PATH ... ";
        cp -r ${HADOOP_PATH} ${MASTER_PATH};
        mkdir ${MASTER_PATH}/tmp
        MASTERS_CONFIG="${MASTER_PATH}/conf/masters";
        SLAVES_CONFIG="${MASTER_PATH}/conf/slaves";
        MAPRED_SITE="${MASTER_PATH}/conf/mapred-site.xml";
        CORE_SITE="${MASTER_PATH}/conf/core-site.xml";
        HDFS_SITE="${MASTER_PATH}/conf/hdfs-site.xml";
        HADOOP_ENV="${MASTER_PATH}/conf/hadoop-env.sh";

        cp ${CONFIG_PATH}/masters ${MASTERS_CONFIG}
        cp ${CONFIG_PATH}/slaves ${SLAVES_CONFIG}
        cp ${CONFIG_PATH}/mapred-site.xml ${MAPRED_SITE}
        cp ${CONFIG_PATH}/core-site.xml ${CORE_SITE}
        cp ${CONFIG_PATH}/hdfs-site.xml ${HDFS_SITE}
        cp ${CONFIG_PATH}/hadoop-env.sh ${HADOOP_ENV}
        
        echo "${MASTER}" > ${MASTERS_CONFIG};
        cat ${MASTERS_CONFIG};
        #echo ${SLAVESS} | tr " " "\n" > ${SLAVES_CONFIG};
        > ${SLAVES_CONFIG};
        for ((i=0; i < ${#SLAVES[*]}; i++))
        do
            echo "${SLAVES[i]}" >> ${SLAVES_CONFIG};
        done
        cat ${SLAVES_CONFIG};
        
        sed -i '/<name>\mapred\.job\.tracker<\/name>/!b;n;c<value>'${MASTER}:54311'</value>' ${MAPRED_SITE};
        sed -i '/<name>\mapred\.local\.dir<\/name>/!b;n;c<value>'${MASTER_PATH}/tmp'</value>' ${MAPRED_SITE};

        sed -i '/<name>\hadoop\.tmp\.dir<\/name>/!b;n;c<value>'${MASTER_PATH}/tmp'</value>' ${CORE_SITE};
        sed -i '/<name>\dfs\.data\.dir<\/name>/!b;n;c<value>'${MASTER_PATH}/tmp'</value>' ${CORE_SITE};
        sed -i '/<name>\dfs\.name\.dir<\/name>/!b;n;c<value>'${MASTER_PATH}/tmp'</value>' ${CORE_SITE};
        sed -i '/<name>fs\.default\.name<\/name>/!b;n;c<value>'hdfs://${MASTER}:54310'</value>' ${CORE_SITE};
        
        sed -i '/export JAVA_HOME/c\export JAVA_HOME='${JAVA_HOME}'' ${HADOOP_ENV};
        
        echo "Y" | hadoop namenode -format;
        
        for ((i=0; i < ${#SLAVES[*]}; i++))
        do
            if [ ${SLAVES[i]} != ${MASTER} ]
            then 
                # echo "export JAVA_HOME=/usr/lib/jvm/java" >> ~/.bashrc; cat ~/.bashrc; $JAVA_HOME/bin/java -version; exit;
                EXPORT_JAVA_HOME="sed -i '/export JAVA_HOME/c\export JAVA_HOME='\$JAVA_HOME'' $SLAVE_PATH/conf/hadoop-env.sh";
                
                #SCP="scp -r ${USERNAME}@${MASTER}:$MASTER_PATH $LOCAL_PATH"
                (scp -r  ${MASTER_PATH} ${USERNAME}@${SLAVES[i]}:${LOCAL_PATH};) &
                #(ssh -l ${USERNAME} ${SLAVES[i]} "hostname; ${SCP}; ${EXPORT_JAVA_HOME};") &
                (ssh -l ${USERNAME} ${SLAVES[i]} "hostname; ${EXPORT_JAVA_HOME};") &
            fi
            wait;
        done    
    fi
    echo "Starting $MASTER_PATH ... ";
    start-all.sh
    wait;
fi
T=15;
echo "sleeping ${T}";
sleep ${T};
GRAPHS="soc-LiveJournal1.txt facebook_combined.txt"
HD_BASE="/graphs";
DISK_BASE=${STORAGE}/distrograph;

for G in ${GRAPHS}
do
    DISK_INPUT_PATH="${DISK_BASE}/inputs/${G}";
    HD_INPUT_PATH="${HD_BASE}/${G}";
    hadoop dfs -copyFromLocal ${DISK_INPUT_PATH} ${HD_INPUT_PATH};
    wait;
done
