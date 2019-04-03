#!/bin/bash
## Run Collect read counts on multiple samples
usage() {
    echo "Usage: $0 [-d <string>] [-l <string>] [-b <string>] [-i <string>] [-e <string>]"
    echo "        -d location directory"
    echo "        -l file containing names of samples to be used"
    echo "                e.g. MK29-T_S4"
    echo "        -b bam qualifier, (sorted, gc_corrected etc.)"
    echo "                e.g. MK29-T_S4.sorted.bam"
    echo "        -i preprocessed interval list"
    echo "        -e extra stuff to name"
    exit 1;
}
while getopts ":d:l:b:i:e:" o; do
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
        i)
            i=${OPTARG}
            ;;

        e)
            e=${OPTARG}
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
    usage   
#  b="sorted"
fi
if [ -z "${i}" ]; then
    usage
fi
if [ -z ${e} ]; then
    e=""
fi
for s in $(cat ${l}); do
    gatk CollectReadCounts \
        -I ${d}/Bams/${s}/${s}.${b}.bam \
        -L ${i} \
        --interval-merging-rule OVERLAPPING_ONLY \
        --read-filter MappingQualityReadFilter --minimum-mapping-quality 10 \
        -DF NotDuplicateReadFilter \
        -O CNV/${s}${e}.counts.hdf5
done
