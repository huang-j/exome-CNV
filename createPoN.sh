#!/bin/bash
## pass in a file full of file names to create a panel of normals.
usage() {
    echo "Usage: $0 [-l <string>] [-o <string>] "
    echo "        -l file containing names of samples to be used"
    echo "                e.g. MK29-T_S4.counts.hdf5"
    echo "        -o output file name"
    exit 1;
}
while getopts ":l:o:" x; do
    case "${x}" in
        l)
            l=${OPTARG}
            ;;
        o)
            o=${OPTARG}
            ;;
        *)
            usage
            ;;
    esac
done
if [ -z "${l}" ]; then
    usage
fi
if [ -z "${o}" ]; then
    ${o}="cnvPoN"
fi

vars=""
for s in $(cat ${l}); do
    vars=${vars}"-I ${s} "
done

echo ${vars}

gatk --java-options "-Xmx32G" CreateReadCountPanelOfNormals ${vars} \
     --minimum-interval-median-percentile 5.0 \
     -O ${o}.pon.hdf5
