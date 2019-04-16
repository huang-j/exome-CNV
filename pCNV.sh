#!/bin/bash
usage() {
    echo "Usage: $0 [-i <string>] [-p <string>] [-a <string>] [-o <string>]"
    echo "        -i input file "
    echo "        -p PoN"
    echo "        -a annotated intervals"
    echo "	  -o output directory"
    exit 1;
}
while getopts ":i:p:a:o:" x; do
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
	o)
	    o=${OPTARG}
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
if [ -z ${o} ]; then
    usage
fi

f=${i##*/}
f=${f%.counts.hdf5}
mkdir ${o}/${f}
mkdir ${o}/${f}/plots


gatk --java-options "-Xmx24G" DenoiseReadCounts \
    -I ${i} \
    --count-panel-of-normals ${p} \
    --annotated-intervals ${a} \
    --standardized-copy-ratios ${o}/${f}/${f}.standardizedCR.tsv \
    --denoised-copy-ratios ${o}/${f}/${f}.denoisedCR.tsv

gatk PlotDenoisedCopyRatios \
    --standardized-copy-ratios ${o}/${f}/${f}.standardizedCR.tsv \
    --denoised-copy-ratios ${o}/${f}/${f}.denoisedCR.tsv \
    --sequence-dictionary ~/rsrch/ref/hg19_k/hg19.dict \
    --minimum-contig-length 46709983 \
    --output ${o}/${f}/plots \
    --output-prefix ${f}
