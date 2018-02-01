#!/bin/bash
#   revolver.sh: Run revolver
#   (c) Mohammad HMofrad, 2017
#   (e) mohammad.hmofrad@pitt.edu

#BASE="/afs/cs.pitt.edu/usr0/hasanzadeh/private";
#BASE="/ihome/cs3580_2017F/moh18";
#BASE="/zfs1/cs3580_2017F/moh18/";
BASE="/home/moh18/";
SOURCE="${BASE}/distrograph/devel/distrograph";
ALGORITHMS="range";
#FILES="${BASE}/facebook_combined.txt";
#FILES="${BASE}/graphs/soc-LiveJournal1.txt";
#FILES="${BASE}/graphs/com-orkut.ungraph.txt";
#FILES="${BASE}/graphs/hollywood-2011.txt";
#FILES="${BASE}/graphs/eu-2015-host.txt";
#FILES="${BASE}/graphs/higgs-twitter.txt";
#FILES="${BASE}/graphs/higgs-twitter.txt ${BASE}/graphs/uk-2007-05@1000000.txt ${BASE}/graphs/enwiki-2013.txt";
FILES="${BASE}/graphs/uk-2007-05@1000000.txt";

PARTITIONS="64 128 192 256";
COMMENT="#";
NUM_ITEMS=2;
DELIMITER=" ";
for A in ${ALGORITHMS}
do
    for F in ${FILES}
    do
        for P in ${PARTITIONS}
        do 
            ${SOURCE} -n ${P} -a ${A} -f ${F} ${COMMENT} " " ${NUM_ITEMS}
        done
    done
done
