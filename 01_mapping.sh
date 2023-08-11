###### RNA-seq 
genome=~/yuyun/genome_data/161010_Chinese_Spring_v1.0_pseudomolecules_parts.fasta
genomesize=~/yuyun/genome_data/161010_Chinese_Spring_v1.0_pseudomolecules_parts.genomesize
gtf=/public/home/zhangyuyun/yuyun/genome_data/IWGSC_v1.1_HC_20170706_parts.gtf
rundir=~/miniconda3/envs/yuyun/bin/
index=/public/home/xieyilin/yilin/genome/iwgsc_cs/161010_Chinese_Spring_v1.0_pseudomolecules_parts

test -d 02_fastp || mkdir 02_fastp
test -d 03_hisat2 || mkdir 03_hisat2

n=0
for i in `cat list`
do
    n=`expr $n + 1`
{
    echo start_${i}
    ${rundir}fastqc -o 02_fastp 01_rawdata/${i}*_1.fq.gz 01_rawdata/${i}*_2.fq.gz &
    ${rundir}fastp -w 8 -l 20 -q 20 -r --cut_right_window_size 4 --cut_right_mean_quality 20 --detect_adapter_for_pe \
        -i 01_rawdata/${i}*_1.fq.gz -I 01_rawdata/${i}*_2.fq.gz -o 02_fastp/${i}_1_val_1.fq.gz -O 02_fastp/${i}_2_val_2.fq.gz
    ${rundir}fastqc -o 02_fastp 02_fastp/${i}_1_val_1.fq.gz 02_fastp/${i}_2_val_2.fq.gz &

    cd 03_hisat2
    ${rundir}hisat2  -x ${index} -1 ../02_fastp/${i}_1_val_1.fq.gz -2 ../02_fastp/${i}_2_val_2.fq.gz -p 10 |\
        ${rundir}samtools sort -@ 10 - > ${i}.sort.bam
    ${rundir}samtools stats ${i}.sort.bam > ${i}.sort.bam.stats &
    ${rundir}samtools view -@ 10 -b -S -q 20 ${i}.sort.bam | ${rundir}bamtools filter -in - -tag NH:1 -out ${i}.sort.q20.uniq.bam
    ${rundir}samtools stats ${i}.sort.q20.uniq.bam > ${i}.sort.q20.uniq.bam.stats

    ${rundir}featureCounts -T 10 -p -a $gtf -o featurecounts.${i}.result ${i}.sort.q20.uniq.bam
    n=`grep 'reads mapped:' ${i}.sort.q20.uniq.bam.stats | cut -f 3`
    tail -n +3 featurecounts.${i}.result | awk -v n=$n '{print $1"\t"($7*2000000000)/(n*$6)}' |\
        sort -k1,1 | cut -f 1,2 > featurecounts.${i}.fpkm.tmp

    scale=`grep 'reads mapped:'  ${i}.sort.q20.uniq.bam.stats | cut -f3 | awk '{print 1000000/$1}'`
    ${rundir}bedtools genomecov -ibam ${i}.sort.q20.uniq.bam  -bg -split -scale ${scale} > ${i}.sort.q20.uniq.rpm.bedgraph
    ${rundir}wigToBigWig -clip ${i}.sort.q20.uniq.rpm.bedgraph ${genomesize} ${i}.sort.q20.uniq.rpm.bw
    rm -rf ${i}.sort.q20.uniq.rpm.bedgraph
    cd ../
    echo finished_${i}
}&
    if [ $n == 3 ]
    then
        echo 'wait'
        n=0
        wait
    fi
done
wait

###### ChIP-seq 
test -d 02_fastqc || mkdir 02_fastqc
test -d 02_fastp || mkdir 02_fastp
test -d 03_bwa || mkdir 03_bwa
test -d 04_macs2 || mkdir 04_macs2
DIR=~/miniconda3/envs/yuyun/bin
function chip(){
    i=$1
    genome=~/yuyun/genome_data/161010_Chinese_Spring_v1.0_pseudomolecules_parts.fasta
    genomesize=~/yuyun/genome_data/161010_Chinese_Spring_v1.0_pseudomolecules_parts.genomesize
    genomelen=14271578887
    index=~/yuyun/genome_data/wheat6_bwaindex/161010_Chinese_Spring_v1.0_pseudomolecules_parts.fasta

    ${DIR}/fastqc -o 02_fastqc  01_rawdata/${i}*_1.fq.gz 01_rawdata/${i}*_2.fq.gz &
    fastp -w 8 -l 20 -q 20 -r --cut_right_window_size 4 --cut_right_mean_quality 20 --detect_adapter_for_pe \
        -i 01_rawdata/${i}*_1.fq.gz -I 01_rawdata/${i}*_2.fq.gz -o 02_fastp/${i}_1_val_1.fq.gz -O 02_fastp/${i}_2_val_2.fq.gz
    ${DIR}/fastqc -o 02_fastp  02_fastp/${i}_1_val_1.fq.gz 02_fastp/${i}_2_val_2.fq.gz &

    cd 03_bwa
    file1=../02_fastp/${i}_1_val_1.fq.gz
    file2=../02_fastp/${i}_2_val_2.fq.gz
    ${DIR}/bwa mem -t 10 ${index} $file1 $file2 | ${DIR}/samtools view -bS -@ 10 - | ${DIR}/samtools sort -@ 10 - > ${i}.sort.bam
    ${DIR}/samtools stats ${i}.sort.bam >  ${i}.sort.bam.stats &
    ${DIR}/samtools view -q 20 -b -@ 10  ${i}.sort.bam | ${DIR}/samtools rmdup -  ${i}.sort.q20.rmdup.bam
    ${DIR}/samtools stats ${i}.sort.q20.rmdup.bam > ${i}.sort.q20.rmdup.bam.stats
    scale=`grep 'reads mapped:'  ${i}.sort.q20.rmdup.bam.stats | cut -f3|awk '{print 1000000/$1}'`
    ${DIR}/bedtools genomecov -ibam  ${i}.sort.q20.rmdup.bam  -bg -split -scale ${scale} > ${i}.sort.q20.rmdup.rpm.bedgraph
    ${DIR}/wigToBigWig -clip ${i}.sort.q20.rmdup.rpm.bedgraph ${genomesize} ${i}.sort.q20.rmdup.rpm.bw &
    ${DIR}/bedtools bamtobed -i ${i}.sort.q20.rmdup.bam | ${DIR}/bedtools sort -i - > ${i}.sort.q20.rmdup.bed &
    cd ../

    cd 04_macs2
    ${DIR}/macs2 callpeak -t ../03_bwa/${i}.sort.q20.rmdup.bam -f BAMPE -g ${genomelen} --nomodel --nolambda -n ${i}_PE
    awk '$8>=10{print $1"\t"$2"\t"$3"\t"$4}' ${i}_PE_peaks.narrowPeak > ${i}_PE_peaks.p10.bed
    cut -f 4 ${i}_PE_peaks.p10.bed | sort | join -1 1 - -2 10 <(grep -v '#' ${i}_PE_peaks.xls  | tail -n +3 | sort -k10,10) |\
        awk '{print $2,$3,$4,$5,$6,$7,$8,$9,$10,$1}' | sort -k9,9nr -k8,8nr > ${i}_PE_peaks.p10.xls #sort by qvalue and then by enrichment
    cut -f 4  ${i}_PE_peaks.p10.bed | sort | join -1 1 - -2 4 <(sort -k4,4 ${i}_PE_summits.bed) |\
        awk '{print $2"\t"$3"\t"$4"\t"$1}' > ${i}_PE_summits.p10.bed
    cd ../

}


j=0
for i in `cat list`
do
    j=`expr $j + 1`
{
    echo "start $i"
    chip ${i} CS
    echo "finished $i"
}&
    if [ $j == 3 ]
    then
        echo 'wait'
        j=0
        wait
    fi
done


