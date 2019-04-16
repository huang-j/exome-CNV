#!/bin/bash
## Run Collect read counts on multiple samples
usage() {
    echo "Usage: $0 [-f <string>] [-i <string>] [-e <string>] [-o <string>]"
    echo "        -f file"
    echo "        -i preprocessed interval list"
    echo "        -e extra stuff to name"
    echo "        -o output"
    exit 1;
}
while getopts ":f:i:e:o:" x; do
    case "${x}" in
        f)
            f=${OPTARG}
            ;;
        i)
            i=${OPTARG}
            ;;
        e)
            e=${OPTARG}
            ;;
        o)
            o=${OPTARG}
            ;;
        *)
            usage
            ;;
    esac
done
if [ -z "${f}" ]; then
    usage
fi
if [ -z "${i}" ]; then
    usage
fi
if [ -z ${e} ]; then
    e=""
fi
if [ -z ${o} ]; then
    o="${f%.bam}.counts.hdf5"
fi
s=
gatk CollectReadCounts \
    -I ${f} \
    -L ${i} \
    --interval-merging-rule OVERLAPPING_ONLY \
    --read-filter MappingQualityReadFilter --minimum-mapping-quality 10 \
    -DF NotDuplicateReadFilter \
    -O ${o}
