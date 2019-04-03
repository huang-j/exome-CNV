#!/bin/bash
## Run Collect read counts on multiple samples
## Just an alternate version of the other one because I used a different folder for this.
usage() {
    echo "Usage: $0 [-d <string>] [-l <string>] [-b <string>] "
    echo "        -d location directory"
    echo "        -l file containing names of samples to be used"
    echo "                e.g. MK29-T_S4"
    echo "        -b bam qualifier, (sorted, gc_corrected etc.)"
    echo "                e.g. MK29-T_S4.sorted.bam"
    exit 1;
}
while getopts ":d:l:b" o; do
    case "${o}" in
        d)
            d=${OPTARG}
            ;;
        l)
            l=${OPTARG}
            ;;
        b)
            b=${OPTARG}
            ;;
        *)
            usage
            ;;
    esac
done
if [ -z "${l}" ]; then
    usage
fi
if [ -z "${d}" ]; then
    d="Agilent/firstSet"
fi
if [ -z "${b}" ]; then
    b="sorted"
fi

for s in $(cat ${l}); do
    gatk CollectReadCounts \
        -I ${d}/NB_Bams/${s} \
        -L ref/crev2.preprocessed.interval_list \
        --interval-merging-rule OVERLAPPING_ONLY \
        --read-filter MappingQualityReadFilter --minimum-mapping-quality 10 \
        -DF NotDuplicateReadFilter \
        -O CNV/${s}.counts.hdf5
done
