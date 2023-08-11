####### 1. identify the level of pi in each 1k bin 
cut -f 1-4 ../02_pi.vs.nohomo_rmlowcov/wheat6.bin1k.pi.bed > wheat6.bin1k.bed ### remove level0 and level1 coverage bin

dir=/public/home/zhangyuyun/yuyun/dap-seq/wheat_re-seq_SNP/guoweilong
for i in `cat ${dir}/chrlist`
do
    grep $i wheat6.bin1k.bed | bedtools map -a - -b ${dir}/${i}.100bp.pi.part.bed -o sum -c 4 |\
        awk '{if($4==".") print $1,$2,$3,$4,"0";else print $1,$2,$3,$4,$5/10}' | tr ' ' '\t' >> wheat6.bin1k.pi.bed
done

##### quantile pi
cut -f 4,5 wheat6.bin1k.pi.bed | awk '$2>0' > wheat6.bin1k.pi.txt
Rscript work_get_quantile.r wheat6.bin1k.pi.txt 5 up
for j in {1..5}
do
    n1=`awk 'NR=="'$j'"{print $2*1}' wheat6.bin1k.pi.txt.quantile`
    n2=`awk 'NR=="'$j'"{print $3*1}' wheat6.bin1k.pi.txt.quantile`

    awk -v n1=$n1 -v n2=$n2 '$5>=n1&&$5<=n2{print $0"\tlevel'${j}'"}' wheat6.bin1k.pi.bed >> wheat6.bin1k.pi.quantile.bed
done
awk '$5==0{print $0"\tlevel0"}' wheat6.bin1k.pi.bed >> wheat6.bin1k.pi.quantile.bed
bedtools sort -i wheat6.bin1k.pi.quantile.bed > wheat6.bin1k.pi.quantile.bed2 && mv wheat6.bin1k.pi.quantile.bed2 wheat6.bin1k.pi.quantile.bed

##### work_get_quantile.r
info<-commandArgs(T)
data<-read.table(info[1])
level<-as.numeric(info[2])
type<-info[3]
q<-data.frame(quantile(data$V2,c(seq(0,1,length.out=level+1))))
names(q)<-"num"
qtf<-NULL
for (i in 1:level){
    if (type %in% "up"){
        qt<-data.frame(level=paste0("level",i),num1=q[i,1],num2=q[i+1,1])
        qtf<-rbind(qtf,qt)
    }else{
        qt<-data.frame(level=paste0("level",rev(1:level)[i]),num1=q[i,1],num2=q[i+1,1])
        qtf<-rbind(qtf,qt)
    }
}
write.table(qtf,paste0(info[1],".quantile"),row.names=F,col.names=F,quote=F)

####### 2. enrichment analysis: pi vs subg-homo 
genomesize=~/yuyun/genome_data/161010_Chinese_Spring_v1.0_pseudomolecules_parts.genomesize

cat wheat6.bin1k.pi.quantile.bed | \
	awk '{if($6=="level5")print $1,$2,$3,$4,"high-pi";else if($6="level1"||$6=="level0")print $1,$2,$3,$4,"low-pi"}' | \
	tr ' ' '\t' > wheat6.bin1k.pi.quantile.selected.bed
file=wheat6.bin1k.pi.quantile.selected.bed

homo3syn=/public/home/zhangyuyun/yuyun/nucmer/01_CS_ABD/cs_ABDhomo.multiple.RH.homo3-syn.bed
nohomo=/public/home/zhangyuyun/yuyun/nucmer/01_CS_ABD/cs_ABDhomo.multiple.RH.nohomo.bed
ln -s /public/home/zhangyuyun/yuyun/nucmer/01_CS_ABD/cs_ABDhomo.multiple.RH.homo3-syn.bed wheat6.homo.bed
ln -s /public/home/zhangyuyun/yuyun/nucmer/01_CS_ABD/cs_ABDhomo.multiple.RH.nohomo.bed wheat6.nohomo.bed

bedtools intersect -a $file -b $nohomo -wao | awk '$NF>100' | awk '{print $4,$5,$NF,"nohomo"}' > nohomo.temp
bedtools intersect -a $file -b $homo3syn -wao | awk '$NF>100' | awk '{print $4,$5,$NF,"homo"}' > homo.temp

cat *.temp | sort -k1,1 -k3,3nr | sort -k1,1 -u > wheat6.bin1k.pi.quantile.txt

bg=`cat $genomesize | awk '{s=s+$2}END{print s}'`
awk '{a[$5]=a[$5]+$3-$2}END{for(i in a)print i,a[i]}' wheat6.bin1k.pi.quantile.selected.bed > pi-level.length
rm enrichment.plot
for j in nohomo homo
do
	a=`cat wheat6.${j}.bed | awk '{s=s+$3-$2}END{print s}'`
	awk '{a[$2]=a[$2]+$3}END{for(i in a)print i,a[i]}' ${j}.temp | sort -k1,1 | join -1 1 - -2 1 pi-level.length |\
		awk '{print "'$j'",$0,"'$a'","'$bg'",($2/$3)/("'$a'"/"'$bg'")}' >> enrichment.plot.temp 
done
cut -f 2,4 -d ' ' wheat6.bin1k.pi.quantile.txt | sort | uniq -c | sed 's/_div//' | awk '{print $3,$2,$1}' | \
	sort -k1,1 | join -1 1 <(sort -k1,1 enrichment.plot.temp) -2 1 - | awk '$2==$8' | \
	cut -f 8 --complement -d ' ' > enrichment.plot2



