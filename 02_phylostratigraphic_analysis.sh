##### Phylostratigraphic Analysis
grep AT ../Orthogroups.tsv | cut -f 7-9 | sed 's/, /\n/g;s/\t/\n/g' | sed 's/\..*//' | sed '/^$/d' > wheat6.age1.list
grep -v AT ../Orthogroups.tsv | grep -E 'Zm' | cut -f 7-9 | sed 's/, /\n/g;s/\t/\n/g' | sed 's/\..*//' | sed '/^$/d' > wheat6.age2.list
grep -v AT ../Orthogroups.tsv | grep -v Zm | grep LOC | cut -f 7-9 | \
    sed 's/, /\n/g;s/\t/\n/g' | sed 's/\..*//' | sed '/^$/d' > wheat6.age3.list
grep -v AT ../Orthogroups.tsv | grep -v Zm | grep -v LOC | grep Bradi | \
    cut -f 7-9 | sed 's/, /\n/g;s/\t/\n/g' | sed 's/\..*//' | sed '/^$/d' > wheat6.age4.list
grep -v AT ../Orthogroups.tsv | grep -v Zm | grep -v LOC | grep -v Bradi | grep HORVU | \
    cut -f 7-9 | sed 's/, /\n/g;s/\t/\n/g' | sed 's/\..*//'| sed '/^$/d' > wheat6.age5.list
grep -v AT ../Orthogroups.tsv | grep -v Zm | grep -v LOC | grep -v Bradi | grep -v HORVU | grep ScWN | \
    cut -f 7-9 | sed 's/, /\n/g;s/\t/\n/g' | sed 's/\..*//' | sed '/^$/d' > wheat6.age6.list
grep -v AT ../Orthogroups.tsv | grep -v Zm | grep -v LOC | grep -v Bradi | grep -v HORVU | grep -v ScWN | \
    grep -E 'Tu|AET' | cut -f 7-9 | sed 's/, /\n/g;s/\t/\n/g' | sed 's/\..*//' | sed '/^$/d' > wheat6.age7.list
grep -v AT ../Orthogroups.tsv | grep -v Zm | grep -v LOC | grep -v Bradi | grep -v HORVU | grep -v ScWN | grep -v Tu| grep -v AET|\
    grep 'TRID' | cut -f 7-9 | sed 's/, /\n/g;s/\t/\n/g' | sed 's/\..*//' | sed '/^$/d' > wheat6.age8.list
cat wheat6.age[1-8].list <(grep '>' /public/home/zhangyuyun/yuyun/genome_data/IWGSC_v1.1_HC_20170706_pep.fasta.re.filter.longest.fa | sed 's/\..*//;s/>//') |\
    sort | uniq -u | grep -v CSU > wheat6.age9.list

##### enrichment with histone modification targeted genes
dir=~/yuyun/orthofinder/wheat246_Hv_rye_Bd/genome_12species_MSU/gene_age/
bg=104567
rm enrichment.age.txt
for i in {1..9}
do
    n1=`cat ${dir}/wheat6.age${i}.list | wc -l`
    for j in H3K27me3 H3K4me3 H3K9ac
    do
        n2=`cat ../05_histone_triad_nontriad/${j}.targetgene.tssupd3k.list | grep -v CSU | wc -l`
        n3=`cat ${dir}/wheat6.age${i}.list ../05_histone_triad_nontriad/${j}.targetgene.tssupd3k.list | \
            sort | uniq -d | wc -l`

        echo "node"$i $j $n1 $n2 $n3 | awk '{print $0,$5/$3,$5/$4,(($5/$4)/($3/"'$bg'"))}' | tr ' ' '\t' >> enrichment.age.txt
    done
done

###### pvalue
cp /public/home/zhangyuyun/yuyun/dap-seq/02_TF_analysis/28_CS_DAP/06_summary/11_TAD_TE/work_get_chi-squared.r ./
awk '{print $1"%"$2,$5,$4,$3,"'$bg'"}' enrichment.age.txt > chitest.age.txt
Rscript work_get_chi-squared.r chitest.age.txt
sed 's/%/ /' chitest.age.txt.pvalue > chitest.age.pvalue
rm chitest.age.txt.pvalue chitest.age.txt

#### work_get_chi-squared.r
file<-commandArgs(T)
data<-read.table(file)
names(data)<-c("TE","count1","count2","count1","count2")
fres<-NULL
for (i in 1:nrow(data)){
    a=rbind(data[i,c(2,3)],data[i,c(4:5)])
    p=chisq.test(a)$p.value
    res<-data.frame(data[i,1],p)
    fres<-rbind(fres,res)
}
write.table(fres,paste0(file,".pvalue"),quote=F,row.names=F,col.names=F)

