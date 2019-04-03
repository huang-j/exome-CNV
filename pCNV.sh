#!/bin/bash
usage() {
    echo "Usage: $0 [-i <string>] [-p <string>] [-a <string>]"
    echo "        -i input file "
    echo "        -p PoN"
    echo "        -a annotated intervals"
    exit 1;
}
while getopts ":i:p:a:" x; do
    case "${x}" in
        i)
            i=${OPTARG}
            ;;
        p)
            p=${OPTARG}
            ;;
        a)
            a=${OPTARG}
            ;;

        *)
            usage
            ;;
    esac
done
if [ -z "${i}" ]; then
    usage
fi
if [ -z "${p}" ]; then
    ${p}="cnvPoN"
fi
if [ -z ${a} ]; then
    usage
fi

mkdir ${i}
mkdir ${i}/plots

gatk --java-options "-Xmx24G" DenoiseReadCounts \
    -I ${i}.counts.hdf5 \
    --count-panel-of-normals ${p}.pon.hdf5 \
    --annotated-intervals ${a} \
    --standardized-copy-ratios ${i}/${i}.standardizedCR.tsv \
    --denoised-copy-ratios ${i}/${i}.denoisedCR.tsv

gatk PlotDenoisedCopyRatios \
    --standardized-copy-ratios ${i}/${i}.standardizedCR.tsv \
    --denoised-copy-ratios ${i}/${i}.denoisedCR.tsv \
    --sequence-dictionary ~/rsrch/ref/hg19_k/hg19.dict \
    --minimum-contig-length 46709983 \
    --output ${i}/plots \
    --output-prefix ${i}

# --minimum-contig-length 46709983 \
