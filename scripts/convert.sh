#!/bin/bash
#   convert.sh: Combines graph edges with graph partitions to
#               create an output graph adjacency list with partition indices
#   (c) Mohammad HMofrad, 2017
#   (e) mohammad.hmofrad@pitt.edu

#BASE="/afs/cs.pitt.edu/usr0/hasanzadeh/private";
BASE="/ihome/cs3580_2017F/moh18";
BASE1="/zfs1/cs3580_2017F/moh18";
FILES="${BASE1}/spinner";
PARTITIONS="8";
ALGORITHMS="spinner";
GRAPHS="soc-LiveJournal1.txt";
#GRAPHS="com-orkut.ungraph.txt";
#GRAPHS="eu-2015-host.txt";
#GRAPHS="hollywood-2011.txt";
#GRAPHS="facebook_combined.txt";
PRINT="${BASE}/revolver/devel/print_stats";

for P in ${PARTITIONS}
do
    for A in ${ALGORITHMS}
    do
        for G in ${GRAPHS}
        do
            #echo "${DIR}/${P}_${A}_parts_${G}";
            ls ${BASE1}/graphs/${G};
            ls ${FILES}/${P}_${A}_parts_${G};
            #echo ${EDGE_PATH}/${P}_${A}_parts_${G}
            ${PRINT} -n ${P} -a ${A} -f ${BASE1}/graphs/${G} "#" "\t" 2  -v ${FILES}/${P}_${A}_parts_${G} " ";
            #ls ${FILES}/${P}_${A}_stats_${G};
            #echo ${FILES}/${P}_${A}_${G}.tar.gz
            #wait;
            #tar zcvf ${FILES}/${P}_${A}_${G}.tar.gz  ${FILES}/${P}_${A}_edges_${G};
            #wait;
        done
    done
done

exit;
