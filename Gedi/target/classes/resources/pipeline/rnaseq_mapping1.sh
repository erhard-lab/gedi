#!/bin/bash

<?JS0

varin("references","Definitions of reference sequences",true);
varin("tmp","Temp directory",true);
varin("mode","SRR/FASTQ",true);
varin("test","Test with the first 10k sequences",false);
varin("nthreads","Number of threads (Default: 8)",false);
varin("minlength","Minimal length of reads to keep (default: 18)",false);
varin("keepUnmapped","Keep the unmapped reads in a fasta file",false);
varin("keepbams","Keep the bam files",false);
varin("keepfastq","Keep the fastq files",false);
varin("sharedMem","Use shared memory for STAR",false);
varin("extract","Extract part of read",false);
varin("nosoft","No softclipping",false);
varin("adapter","Only for single end!",false);
varin("introns","All/Annotated/None",false);
varin("smartseq","*Full length* smart seq library, i.e. trim first 3bp",false);
varin("umi","Reads start with umis (umi-len,spacer-len,umi2-len,spacer2-len; 6,4 for Lexogen UMIs!)",false);
varin("umiregex","UMIs are in Fastq header, specify regex, umi is group1 (only first read for PE)",false);
varin("layout","Specify read layout for umis and linkrs)",false);
varin("keepnoumi","Also keep the cit file without umi deduplication",false);
varin("umiallowmulti","Allow for multimappers when working with UMIs",false);
varin("citparam","e.g. -novar -nosec",false);
varin("samattributes","e.g. nM MD NH",false)
varin("phredfilter","filter by phred, either true or the min score",false)
varin("starindex","folder of a combined STAR index",false);
varin("starparam","additional parameter for STAR",false);
varin("bcid","ids in fastq contain barcodes",false);
varin("cutadaptparam","additional parameter for cutadapt",false);
varin("adapter1","Forward adapter sequence",false);
varin("adapter2","Backward adapter sequence",false);
varin("spike","Spike prefix of chromosomes: Will remove these chromosomes from cit file, but count their stats",false);
varin("SLAMrescue","Indexed pseudogenome for SLAMrescue",false);
varin("trimmomatic","use trimmomatic (for PE)",false);
varin("seqpurge","use seqpurge (for PE)",false);
varin("trim_polyA_R1","Trim poly-A in R1, specify length",false);
varin("trim_polyA_R2","Trim poly-A in R2, specify length",false);
varin("forcecorrect","enforce CorrectCIT also for SE",false);
varin("discardmate","Discard reads 1 or 2 after FastqFilter (default: use both, otherwise specify 1 or 2)",false);


varout("reads","File name containing read mappings");

?>

<?JS
var novar;
var bcid;
var umiregex;
var keepnoumi;
var starparam = starparam?starparam:"";
var cutadaptparam = cutadaptparam?cutadaptparam:"";
var nosec;
var introns;
var keepUnmapped;
var minlength=minlength?minlength:18;
var nthreads = nthreads?nthreads:Math.min(8,Runtime.getRuntime().availableProcessors());
var sharedMem;
var pairedend;
var nosoft;
var smartseq;
var umi;
var umiallowmulti;
var layout;
var spike;
var SLAMrescue;
var trimmomatic;
var discardmate;
var seqpurge;
var trim_polyA_R1=trim_polyA_R1?trim_polyA_R1:0;
var trim_polyA_R2=trim_polyA_R2?trim_polyA_R2:0;

var samattributes = samattributes?samattributes:"nM MD NH"
var phredfilter = phredfilter?phredfilter:false;
if (phredfilter==true) phredfilter=28
var extract;
var extractparam = extract?"-extract "+extract+" ":"";
output.executable=true;

var forcecorrect = forcecorrect?"-f":"";

var intronparam = "";
if (introns && introns.toLowerCase().startsWith("no")) intronparam = "--alignSJDBoverhangMin 9999 --alignSJoverhangMin 9999";
else if (introns && introns.toLowerCase().startsWith("ann")) intronparam = "--alignIntronMax 1";

var citparam = citparam?citparam:"";

var alimode="Local";
if (nosoft) alimode=nosoft;
if (nosoft==true) alimode="EndToEnd";
	


var starindex;
if (!starindex && garray.length==1) starindex=ReadMapper.STAR.getIndex(Genomic.get(garray[0]),null);	
if (!starindex) throw new RuntimeException("Specify a starindex! You can create one by gedi -e GenomicUtils -p -m star -g "+genomes);
	
var fastqname = name+".fastq";

var test;
var keepbams;
var keepfastq;
var fqs;
var minmaq = "";
var checkPElen;
if (new File(tmp+"/"+name).isDirectory() && new File(tmp+"/"+name).listFiles().length>0) log.warning("Non-empty directory "+tmp+"/"+name+" already exists!")
?>


mkdir -p <?JS tmp ?>/<?JS name ?>
cd <?JS tmp ?>/<?JS name ?>

<?JS if (mode=="SRR") { ?>
fasterq-dump --split-files <?JS if(test) print("-X 100000 "); ?><?JS files ?>


<?JS 
	if (pairedend) {
		var f1="_1"
		var f2="_2"
		if (pairedend=="switch") {
			f1="_2"
			f2="_1"
		}
		if (files.contains(" ")) {
			println("cat *_1.fastq > "+name+f1+".fastq");
			println("cat *_2.fastq > "+name+f2+".fastq");
		} else {
			println("mv "+files+"_1.fastq "+name+f1+".fastq");
			println("mv "+files+"_2.fastq "+name+f2+".fastq");
		}
		
		fqs = [name+"_1.fastq",name+"_2.fastq"];
	}
	else {
		if (files.contains(" ")) {
			println("cat *.fastq > "+name+".fastq");
		} else {
			println("mv "+files+".fastq "+name+".fastq");
		}
		fqs = [name+".fastq"];
	}
	fastqname = fqs.join(" ");	
	
} else if (mode=="FASTQ") {  

var ff = [files];
var fqs = [fastqname];

if (pairedend) {
	if (ff.length!=1) throw new RuntimeException("Paired-end reads must be in two files!");
	var fp = FileUtils.findPartnerFile(ff[0],pairedend);
	log.info("Found partner files: "+ff[0]+" "+fp); 
	ff=[ff[0],fp];
	fqs=[name+"_1.fastq",name+"_2.fastq"];
	fastqname = fqs.join(" ");	
}

for (var find=0; find<ff.length; find++) {
	var f1 = ff[find];
	var fq1 = fqs[find];
?>
	<?JS if (!test && f1.endsWith(".gz")) {  ?>
	zcat <?JS f1 ?> > <?JS fq1 ?>
	<?JS } else if (test && f1.endsWith(".gz")) {  ?>
	zcat <?JS f1 ?> | head -n400000 > <? fq1 ?>
	<?JS } else if (!test && f1.endsWith(".bz2")) {  ?>
	bzcat <?JS f1 ?> > <?JS fq1 ?>
	<?JS } else if (test && f1.endsWith(".bz2")) {  ?>
	bzcat <?JS f1 ?> | head -n400000  > <?JS fq1 ?>
	<?JS } else if (!test) {  ?>
	cp <?JS f1 ?> <?JS fq1 ?>
	<?JS } else if (test) {  ?>
	head -n400000 <?JS f1 ?> > <?JS fq1 ?>
	<?JS }
}
} ?>


echo -e "Category\tCount" > <?JS name ?>.reads.tsv
echo -ne "All\t" >> <?JS name ?>.reads.tsv

L=`wc -l <?JS fastqname ?> | head -n1 | awk '{print $1}'`
leftreads=$((L / 4))
echo $leftreads >> <?JS name ?>.reads.tsv

<?JS
var adapter;
if (adapter) { ?>
cutadapt -a <? adapter ?> <? cutadaptparam ?> <?JS fastqname ?> > <?JS fastqname ?>.trimmed
mv <?JS fastqname ?>.trimmed <?JS fastqname ?>

#reaper --nozip --noqc -3p-prefix 1/1/0/0 -swp 1/4/4 -geom no-bc -i <?JS fastqname ?> -basename <?JS name ?>  -3pa <? adapter ?>
#mv <? name ?>.lane.clean <?JS fastqname ?>
#rm <?JS name ?>.lint
<? } ?>

<?JS if (checkPElen) { ?>
EQUAL_LEN=" -checkPElen"
if (( $( x=`head -n2 <? name ?>_1.fastq | tail -n1`; echo ${#x} ) != $( x=`head -n2 <? name ?>_2.fastq | tail -n1`; echo ${#x} ) )); then 
	EQUAL_LEN=""
	echo Not checking for equal read length after trimming! >&2
fi
<?JS } else { ?>
EQUAL_LEN=""
<?JS } ?>

<?JS
var adapter1;
var adapter2;
if (adapter1 && adapter2 && seqpurge) { ?>
SeqPurge -threads <? nthreads ?> -min_len 18 -qcut 0 -a1 <? adapter1 ?> -a2 <? adapter2 ?> -in1 <? name ?>_1.fastq -in2 <? name ?>_2.fastq -out1 tmp1.fastq.gz -out2 tmp2.fastq.gz
gunzip tmp1.fastq.gz
gunzip tmp2.fastq.gz
mv tmp1.fastq <? name ?>_1.fastq 
mv tmp2.fastq <? name ?>_2.fastq
<? } else if (adapter1 && adapter2 && trimmomatic) { ?>
TRIMMER=$( ls -1 `echo $PATH | sed -e 's/:/\/trimmomatic*.jar /'g ` 2>/dev/null | head -n1 )
if [ ! -f $TRIMMER ]; then
	(>&2 echo "Trimmomatic not found, did not trim reads!")
else
	echo -e ">a/1\n<? adapter1 ?>\n>a/2\n<? adapter2 ?>" > adapter.fasta
	java -jar $TRIMMER PE -threads <? nthreads ?> <?JS fastqname ?> -baseout trimmed ILLUMINACLIP:adapter.fasta:2:18:8:1:true
	mv trimmed_1P <? name ?>_1.fastq 
	mv trimmed_2P <? name ?>_2.fastq 
	<? fastqname=name+"_1.fastq "+name+"_2.fastq"; ?>
fi
<? } else if (adapter1 && adapter2) { ?>
cutadapt -a <? adapter1 ?> -A <? adapter2 ?> <? cutadaptparam ?> -o tmp1.fastq -p tmp2.fastq <?JS fastqname ?>
mv tmp1.fastq <? name ?>_1.fastq 
mv tmp2.fastq <? name ?>_2.fastq 
<? fastqname=name+"_1.fastq "+name+"_2.fastq"; ?>
<? }?>

<?JS if (trim_polyA_R1>0) { 
var polyA=StringUtils.repeat("A",trim_polyA_R1)
?>
cutadapt -a <? polyA ?> <? cutadaptparam ?> -o tmp.fastq <? name ?>_1.fastq
mv tmp.fastq <? name ?>_1.fastq
<? }?>

<?JS if (trim_polyA_R2>0) { 
var polyA=StringUtils.repeat("A",trim_polyA_R2)
?>
cutadapt -a <? polyA ?> <? cutadaptparam ?> -o tmp.fastq <? name ?>_2.fastq
mv tmp.fastq <? name ?>_2.fastq
<? }?>


gedi -t . -e FastqFilter -D <?extractparam?>-overwrite<? if (smartseq) print(" -smartseq"); ?><? if (umi) print(" -umi "+umi); ?><? if (umiregex) print(" -umiregex \""+umiregex+"\""); ?><? if (bcid) print(" -keepids "); ?><? if (layout) print(" -layout "+layout); ?> -ld <?JS name ?>.readlengths.tsv $EQUAL_LEN -min <?JS minlength ?> <?JS fastqname ?>


	echo -ne "Trimmed\t" >> <?JS name ?>.reads.tsv
	L=`wc -l <?JS fastqname ?> | head -n1 | awk '{print $1}'`
	leftreads=$((L / 4))
	echo $leftreads >> <?JS name ?>.reads.tsv
	
<?JS if (discardmate) { ?>
rm <? name ?>_<? discardmate ?>.fastq
mv <? name ?>_<? print(new Integer(3-discardmate)) ?>.fastq <? name ?>.fastq
<?JS
	fastqname=name+".fastq"
} ?>


<?JS for each (var index in discardGenomes.keySet()) { ?>


# <? index ?> removal
<?
if (fastqname.contains(" ")) { ?>
bowtie2 -p <? nthreads ?> --un-conc un --local -x <? print(discardGenomes.get(index)) ?> -1 <? name ?>_1.fastq -2 <? name ?>_2.fastq > /dev/null  
	mv un.1 <? name ?>_1.fastq
	mv un.2 <? name ?>_2.fastq
<?
} else { ?>
bowtie2 -p <? nthreads ?> --un un --local -x <? print(discardGenomes.get(index)) ?> -U <? name ?>.fastq > /dev/null  
	mv un <? name ?>.fastq
<? } ?>

	echo -ne "filtered_<? index ?>\t" >> <?JS name ?>.reads.tsv
	L=`wc -l <?JS fastqname ?> | head -n1 | awk '{print $1}'`
	leftreads=$((L / 4))
	echo $leftreads >> <?JS name ?>.reads.tsv
<? } ?>


<? if (spike) { ?>

# mapping
STAR --runMode alignReads --runThreadN <? nthreads ?>  <? intronparam ?> --genomeDir <? starindex ?> <?JS print(sharedMem?"--genomeLoad LoadAndKeep":""); ?> <? starparam ?> --readFilesIn <? fastqname ?> --outSAMmode NoQS --outSAMtype BAM Unsorted --alignEndsType <? alimode ?> --outSAMattributes <? samattributes ?>  <?JS print(keepUnmapped?"--outReadsUnmapped Fastx":""); ?> <?JS print(SLAMrescue?"--outSAMunmapped Within":""); ?>
gedi -e BamFilterSpike Aligned.out.bam <? spike ?> remaining.bam >> <?JS name ?>.reads.tsv
rm Aligned.out.bam
samtools sort -o <? name ?>.bam remaining.bam
rm remaining.bam

<? } else { ?>
<?JS 
var outSAMmode="NoQS"
if (SLAMrescue) outSAMmode="Full"
if (phredfilter) outSAMmode="Full"
?>
# mapping
STAR --runMode alignReads --runThreadN <? nthreads ?>  <? intronparam ?> --genomeDir <? starindex ?> <?JS print(sharedMem?"--genomeLoad LoadAndKeep":""); ?> --limitBAMsortRAM 8000000000 <? starparam ?> --readFilesIn <? fastqname ?> --outSAMmode <? outSAMmode ?> --outSAMtype BAM SortedByCoordinate --alignEndsType <? alimode ?> --outSAMattributes <? samattributes ?> <?JS print(keepUnmapped?"--outReadsUnmapped Fastx":""); ?> <?JS print(SLAMrescue?"--outSAMunmapped Within":""); ?>
mv Aligned.sortedByCoord.out.bam <? name ?>.bam

echo -ne "Unique\t" >> <?JS name ?>.reads.tsv
grep "Uniquely mapped reads number"  Log.final.out | cut -f2 -d'|' | awk '{ print $1}' >> <?JS name ?>.reads.tsv
echo -ne "Multi\t" >> <?JS name ?>.reads.tsv
grep "Number of reads mapped to multiple loci"  Log.final.out | cut -f2 -d'|' | awk '{ print $1}' >> <?JS name ?>.reads.tsv
<? } ?>

<? if (SLAMrescue) { 

?>


# rescue
samtools index <? name ?>.bam
STRAND=$( gedi -e InferStrandness -g <? genomes ?> <? name ?>.bam )

gedi -mem 64G -e ExtractReads -strandness $STRAND -f <? name ?>.bam
samtools view -b -F 4 <? name ?>.bam > <? name ?>.genomemapped.bam
<?JS
var starparam10=starparam.replaceAll("--outFilterMismatchNmax \\d+","--outFilterMismatchNmax 10")
?>
STAR --runMode alignReads --runThreadN <? nthreads ?>  <? intronparam ?> --genomeDir <? SLAMrescueIndex ?> <?JS print(sharedMem?"--genomeLoad LoadAndKeep":""); ?> --limitBAMsortRAM 8000000000 <? starparam10 ?>  --readFilesIn *_unmapped_T2C*fastq --outSAMmode Full --outSAMtype BAM SortedByCoordinate --alignEndsType <? alimode ?> --outSAMattributes <? samattributes ?>
rm  *_unmapped_T2C*fastq  
samtools view -b -F 256 Aligned.sortedByCoord.out.bam > <? name ?>.pseudoMapped.bam
rm Aligned.sortedByCoord.out.bam

gedi -mem 64G -e RescuePseudoReads -genome <? genomes ?> -pseudogenome <? SLAMrescue ?> -strandness $STRAND -origmaps <? name ?>.genomemapped.bam -pseudomaps <? name ?>.pseudoMapped.bam -idMap <? name ?>.idMap
samtools merge -f <? name ?>.bam <? name ?>.genomemapped.bam <? name ?>.pseudoMapped_reverted.bam

echo -ne "Unique rescued\t" >> <?JS name ?>.reads.tsv
grep "Uniquely mapped reads number"  Log.final.out | cut -f2 -d'|' | awk '{ print $1}' >> <?JS name ?>.reads.tsv
echo -ne "Multi rescued\t" >> <?JS name ?>.reads.tsv
grep "Number of reads mapped to multiple loci"  Log.final.out | cut -f2 -d'|' | awk '{ print $1}' >> <?JS name ?>.reads.tsv

<? } ?>

<?JS
if (keepUnmapped) {
	if (fastqname.contains(" ")) {
?>
mv Unmapped.out.mate1 <? name ?>.unmapped_1.fastq
mv Unmapped.out.mate2 <? name ?>.unmapped_2.fastq
<?JS			} else { ?>
mv Unmapped.out.mate1 <? name ?>.unmapped.fastq
<?JS   } 
}
?>

mkdir -p <?JS wd ?>/report

<?JS if (phredfilter) { ?>
gedi -e BamPhredFilter -min <? phredfilter ?> -stats <? name ?>.bam
mv <? name ?>.stats.tsv <?JS wd ?>/report
mv <? name ?>.stats.plots/* <?JS wd ?>/report

<?JS } ?>

samtools index <? name ?>.bam

echo -ne "Mitochondrial\t" >> <?JS name ?>.reads.tsv
gedi -t . -e Bam2CIT<? if (umi || layout || bcid || umiregex) print(" -umi"); ?><? if (umiallowmulti) print(" -umiAllowMulti"); ?> <? minmaq ?> <? citparam ?> <? name ?>.cit <? name ?>.bam | tail -n1 | awk '{ print $2 }'>> <?JS name ?>.reads.tsv
gedi -t . -e CorrectCIT <? forcecorrect ?> <? name ?>.cit
echo -ne "CIT\t" >> <?JS name ?>.reads.tsv
gedi -t . -e ReadCount -g <? genomes ?> -m Weight <? name ?>.cit | tail -n1 | awk '{ print $2 }'>> <?JS name ?>.reads.tsv

mv *.readlengths.* <?JS wd ?>/report

<?JS if (umi || layout || umiregex) { ?>

<?JS if (keepnoumi) { ?>
gedi -t . -e DedupUMI -nocollapse -prefix <? name ?>.nodedup <? name ?>.cit
<?JS } ?>

gedi -t . -e DedupUMI -prefix <? name ?> -mm -plot <? name ?>.cit
echo -ne "Dedup\t" >> <?JS name ?>.reads.tsv
gedi -t . -e ReadCount -g <? genomes ?> -m Weight <? name ?>.cit | tail -n1 | awk '{ print $2 }'>> <?JS name ?>.reads.tsv

mv <?JS name ?>.dedup* <?JS wd ?>/report
<?JS } ?>

mv <?JS name ?>.reads.tsv <?JS wd ?>/report

<?JS if (keepfastq) { ?>
mkdir -p <?JS wd ?>/fastq
mv *.fastq <?JS wd ?>/fastq
gzip <?JS wd ?>/fastq/*.fastq
<?JS } ?>

if [ -f <?JS name ?>.cit ]; then
	if [ -f <?JS name ?>.barcodes.tsv ]; then
		mv <?JS name ?>.barcodes.tsv <?JS wd ?>
	fi
   mv <?JS name ?>.cit* <?JS wd ?>
   <?JS if (keepnoumi) { ?>
   mkdir -p <?JS wd ?>/nodedup
	mv <?JS name ?>.nodedup.cit <?JS wd ?>/nodedup/<?JS name ?>.cit
	mv <?JS name ?>.nodedup.cit.metadata.json <?JS wd ?>/nodedup/<?JS name ?>.cit.metadata.json
	<?JS } ?>

	if [ -d Solo.out ]; then
		mkdir -p <?JS wd ?>/Solo.out
		mv Solo.out <?JS wd ?>/Solo.out/<?JS name ?>
	fi


<?JS if (keepUnmapped) { ?>
   mkdir -p <?JS wd ?>/unmapped
   gzip <?JS name ?>.unmapped*.fastq
   mv <?JS name ?>.unmapped*.fastq.gz <?JS wd ?>/unmapped
<?JS } ?>
<?JS if (keepbams) { ?>
   mkdir -p <?JS wd ?>/bams
   mv <? name ?>.bam* <?JS wd ?>/bams
<?JS } ?>
   cd <?JS wd ?>
   rm -rf <?JS tmp ?>/<?JS name ?>
else
   (>&2 echo "There were some errors, did not delete temp directory!")
   cd <?JS wd ?>
fi



<?JS
var png = "report/"+name+".reads.png";
var table = "report/"+name+".reads.tsv"
processTemplate("plot_mappingstatistics.R",output.file.getParentFile().getParentFile()+"/report/"+name+".reads.R");
?>
Rscript report/<?JS name ?>.reads.R
echo '{"plots":[{"section":"Mapping statistics","id":"mapping<? print(StringUtils.toJavaIdentifier(name)) ?>","title":"<? name ?>","description":"How many reads are removed by adapter trimming and rRNA mapping, how many are mapped to which reference and retained overall.","img":"<? name ?>.reads.png","script":"<? name ?>.reads.R","csv":"<? name ?>.reads.tsv"}]}' > report/<? name ?>.reads.report.json 

