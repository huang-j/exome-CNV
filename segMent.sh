#!/bin/bash
## Pass in tumor and normals for segmentation
usage() {
    echo "Usage: $0 [-f <string>] [-i <string>]"
    echo "        -f file containing names of samples to be used. tab separated where second column is normals"
    echo "                e.g. MK29-T_S4"
    echo "	  -i interval list"
    exit 1;
}
while getopts ":f:i:" x; do
    case "${x}" in
        f)
            f=${OPTARG}
            ;;
        i)
	    i=${OPTARG}
	    ;;
        *)
            usage
            ;;
    esac
done
if [ -z "${f}" ]; then
    usage
fi
if [ -z ${i} ]; then
    usage
fi
## check for end at the end, for silly truncating while reasons
if [ -z $(cat ${f} | grep 'end') ]; then
  echo $'\n'end >> ${f}
fi

cat ${f} | \
while read CMD; do
  t=$(echo $CMD | awk '{print $1}')
  n=$(echo $CMD | awk '{print $2}')
  ## t : tumor
  ## n : normal
  tbam=$(find ~/rsrch/Agilent -type f | grep ${t}'.sorted.bam$')

  gatk --java-options "-Xmx64g" CollectAllelicCounts \
      -L ${i} \
      -I ${tbam} \
      -R ~/rsrch/ref/hg19_k/hg19.fasta \
      -O ${t}/${t}.allelicCounts.tsv

  if [ ! -z $n ]; then
    nbam=$(find ~/rsrch/Agilent -type f | grep ${n}'.sorted.bam$')
    gatk --java-options "-Xmx64g" CollectAllelicCounts \
        -L ${i} \
        -I ${nbam} \
        -R ~/rsrch/ref/hg19_k/hg19.fasta \
        -O ${t}/${n}.allelicCounts.tsv

    gatk --java-options "-Xmx64g" ModelSegments \
        --denoised-copy-ratios ${t}/${t}.denoisedCR.tsv \
        --allelic-counts ${t}/${t}.allelicCounts.tsv \
        --normal-allelic-counts ${t}/${n}.allelicCounts.tsv \
        --output ${t} \
        --output-prefix ${t}
  else
    gatk --java-options "-Xmx64g" ModelSegments \
        --denoised-copy-ratios ${t}/${t}.denoisedCR.tsv \
        --allelic-counts ${t}/${t}.allelicCounts.tsv \
        --output ${t} \
        --output-prefix ${t}
  fi

  gatk CallCopyRatioSegments \
      --input ${t}/${t}.cr.seg \
      --output ${t}/${t}.called.seg

  gatk PlotModeledSegments \
      --denoised-copy-ratios ${t}/${t}.denoisedCR.tsv \
      --allelic-counts ${t}/${t}.hets.tsv \
      --segments ${t}/${t}.modelFinal.seg \
      --sequence-dictionary ~/rsrch/ref/hg19_k/hg19.dict \
      --minimum-contig-length 46709983 \
      --output ${t}/plots \
      --output-prefix ${t}_clean
done
