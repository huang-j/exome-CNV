#!/bin/bash
## pass in a file full of file names to create a panel of normals.
usage() {
    echo "Usage: $0 [-l <string>] [-o <string>] "
    echo "        -l file containing names of samples to be used"
    echo "                e.g. MK29-T_S4"
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
    vars=${vars}"-I CNV/${s}.counts.hdf5 "
done

echo ${vars}

gatk --java-options "-Xmx64G" CreateReadCountPanelOfNormals ${vars} \
     --minimum-interval-median-percentile 5.0 \
     -O CNV/${o}.pon.hdf5
