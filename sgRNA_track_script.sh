###sgRNA tracks by pure Sequence analysis
mkdir sgRNA_analysis
cd sgRNA_analysis

##C. elegans genome WS270 / Essentially any version since WB235 will be the same##
wget ftp://ftp.wormbase.org/pub/wormbase/releases/WS270/species/c_elegans/PRJNA13758/c_elegans.PRJNA13758.WS270.genomic.fa.gz

##Generate bed with all NGGs
zcat c_elegans.PRJNA13758.WS270.genomic.fa.gz | awk '{if($0 ~ />/){print "\n"$0;}else{printf $0}} END{print ""}' | tail -n+2 | awk '{if($0 ~ />/){split($0,chr,">");name=chr[2];}else{for(i=1; i<=(length($0)-22); i++){seq=substr($0,i+21,2); if(seq == "GG"){print name"\t"i-1"\t"(i+19)"\t.\t.\t+"}}}}' - > NGG_Cel_WS270_Fwd.bed
##And then in reverse
zcat c_elegans.PRJNA13758.WS270.genomic.fa.gz | awk '{if($0 ~ />/){print "\n"$0;}else{printf $0}} END{print ""}' | tail -n+2 | perl -pe '{if(!(/>/m)){chomp; tr/atcgATCG/tagcTAGC/; $_=reverse; $_.="\n"}}' | awk '{if($0 ~ />/){split($0,chr,">");name=chr[2];}else{for(i=1; i<=(length($0)-22); i++){seq=substr($0,i+21,2); if(seq == "GG"){print name"\t"(length($0)-i-19)"\t"(length($0)-i+1)"\t.\t.\t-"}}}}' - > NGG_Cel_WS270_Rev.bed
##Combine them
cat NGG_Cel_WS270_Fwd.bed NGG_Cel_WS270_Rev.bed | sort -k1,1 -k2,2n > NGG_Cel_WS270.bed

#Uncompress genome and index it for next steps
gzip -d c_elegans.PRJNA13758.WS270.genomic.fa.gz
bwa index c_elegans.PRJNA13758.WS270.genomic.fa

##Get sequences
bedtools getfasta -fi c_elegans.PRJNA13758.WS270.genomic.fa -bed NGG_Cel_WS270.bed -s > NGG_Cel_WS270.fasta

##Check for uniqueness of 20 mers
../QueryInGenomeWithMMinWindow.pl c_elegans.PRJNA13758.WS270.genomic.fa NGG_Cel_WS270.fasta 20 0 1 20 | awk -F"\t" '{if($3==1){print $1"\n"$2}}' > UniqueNGGs-0MM.fasta

##Map and make track
bwa aln -n 0 -o 0 c_elegans.PRJNA13758.WS270.genomic.fa UniqueNGGs-0MM.fasta | bwa samse c_elegans.PRJNA13758.WS270.genomic.fa - UniqueNGGs-0MM.fasta > Cel_gen_0MM_NGG-seqs.sam
grep "XT:A:U" Cel_gen_0MM_NGG-seqs.sam | awk -F"\t" '{OFS="\t";if($2==16){str="-"}else{str="+"}print $3,($4-1),($4+19),$1,$6,str}' | sort -k1,1 -k2,2n > Cel_gen_0MM_NGG-seqs.bed

###1 mismatch by perl strategy onlyjajja
##Convert into 18-mers then add 1 MM
#1 MM
awk '{if($0 ~ />/){name=$0;}else{for(i=1; i<=(length($0)-19); i++){print name"_"i"separatortoto"substr($0,i,20)"\n"substr($0,i+1,18)}}}' UniqueNGGs-0MM.fasta > NGGs-18mers.fasta

##Look for uniqueness in 18 mers with perl script (hash method)
../QueryInGenomeWithMMinWindow2.pl c_elegans.PRJNA13758.WS270.genomic.fa NGGs-18mers.fasta 18 0 1 18 18.tab

##Convert into fasta files sequences that appear only once, place it as 20mer; Ts have been removed before.
awk -F"\t" '{if($3 == 1){print $1}}' 18.tab | awk -F"separatortoto" '{print $1"\n"$2}' - > Unique-18mers.fasta

##Filter for up to one MM as a 20mer
../QueryInGenomeWithMMinWindow.pl c_elegans.PRJNA13758.WS270.genomic.fa Unique-18mers.fasta 20 1 2 19 | awk -F"\t" '{if($3==1){print $1"\n"$2}}' > Unique18mersNGGs-1MM.fasta

##Map and make track
bwa aln -n 0 -o 0 c_elegans.PRJNA13758.WS270.genomic.fa Unique18mersNGGs-1MM.fasta | bwa samse c_elegans.PRJNA13758.WS270.genomic.fa - Unique18mersNGGs-1MM.fasta > Cel_gen_1MM_NGG-seqs.sam
grep "XT:A:U" Cel_gen_1MM_NGG-seqs.sam | awk -F"\t" '{OFS="\t";if($2==16){str="-"}else{str="+"}print $3,($4-1),($4+19),$1,$6,str}' | sort -k1,1 -k2,2n > Cel_gen_1MM_NGG-seqs.bed

###2 mismatches by blast strategy only
##Create DB and filter for alignements larger than 17
#2 MM

##After installing blast, it is required to create a database of the genome to map against
makeblastdb -in c_elegans.PRJNA13758.WS270.genomic.fa -dbtype 'nucl' -title c_elegans.WS270

##Blast filtered unique 18 mer sequences using "blast short algorithm"
blastn -db c_elegans.PRJNA13758.WS270.genomic.fa -query Unique18mersNGGs-1MM.fasta -task "blastn-short" -out Blast18mers-1MM.tab -outfmt 6 -evalue 1 -num_threads 20

##Get names of seqs that are unique in a set that matches at least 18 bp and start the aligment before the third position
awk -F"\t" '{if((($3*$4)/100)>17){if($7 < 3){array[$1]++}}} END{for (key in array){if(array[key]==1){print key}}}' Blast18mers-1MM.tab > Blast_18mers_filt-2MM.txt

##Parse output to create a file with the indexes of the unique sequences
awk '{if($0 ~ />/){name=$0;}else{print name"\t"$0}}' Unique18mersNGGs-1MM.fasta | awk -F">" '{print $2}' | awk -F"\t" '{if(array[$1] != 0){print ">"$1"\n"array[$1]}else{array[$1]=$2}}' - Blast_18mers_filt-2MM.txt > Cel_gen_2MM_NGG-seqs.fasta

##Map unique to obtain position
bwa aln -n 0 -o 0 c_elegans.PRJNA13758.WS270.genomic.fa Cel_gen_2MM_NGG-seqs.fasta | bwa samse c_elegans.PRJNA13758.WS270.genomic.fa - Cel_gen_2MM_NGG-seqs.fasta > Cel_gen_2MM_NGG-seqs.sam
grep "XT:A:U" Cel_gen_2MM_NGG-seqs.sam| awk -F"\t" '{OFS="\t";if($2==16){str="-"}else{str="+"}print $3,($4-1),($4+19),$1,$6,str}' | sort -k1,1 -k2,2n > Cel_gen_2MM_NGG-seqs.bed

###3 mismatches by blast strategy only
## Use previous alignement and filter from it
#3 MM

##Get names of seqs that are unique in a set that matches at least 18 bp and start the aligment before the third position
awk -F"\t" '{if((($3*$4)/100)>16){if($7 < 4){array[$1]++}}} END{for (key in array){if(array[key]==1){print key}}}' Blast18mers-1MM.tab > Blast_18mers_filt-3MM.txt

##Parse output to create a file with the indexes of the unique sequences
awk '{if($0 ~ />/){name=$0;}else{print name"\t"$0}}' Unique18mersNGGs-1MM.fasta | awk -F">" '{print $2}' | awk -F"\t" '{if(array[$1] != 0){print ">"$1"\n"array[$1]}else{array[$1]=$2}}' - Blast_18mers_filt-3MM.txt > Cel_gen_3MM_NGG-seqs.fasta

##Map unique to obtain position
bwa aln -n 0 -o 0 c_elegans.PRJNA13758.WS270.genomic.fa Cel_gen_3MM_NGG-seqs.fasta | bwa samse c_elegans.PRJNA13758.WS270.genomic.fa - Cel_gen_3MM_NGG-seqs.fasta > Cel_gen_3MM_NGG-seqs.sam
grep "XT:A:U" Cel_gen_3MM_NGG-seqs.sam | awk -F"\t" '{OFS="\t";if($2==16){str="-"}else{str="+"}print $3,($4-1),($4+19),$1,$6,str}' | sort -k1,1 -k2,2n > Cel_gen_3MM_NGG-seqs.bed

###4 mismatches by 16-mer and blast strategy
## Use previous unique mer fasta and convert it into 16 mer (2MM-16-2MM)
#4 MM
awk '{if($0 ~ />/){name=$0;}else{for(i=1; i<=(length($0)-19); i++){print name"_"i"separatortoto"substr($0,i,20)"\n"substr($0,i+2,16)}}}' Cel_gen_3MM_NGG-seqs.fasta > 3MM-NGGs-16mers.fasta

##Look for uniqueness in 18 mers with perl script (hash method)
../QueryInGenomeWithMMinWindow2.pl c_elegans.PRJNA13758.WS270.genomic.fa 3MM-NGGs-16mers.fasta 16 0 1 16 16.tab

##Convert into fasta files sequences that appear only once, place it as 20mer; Ts have been removed before.
awk -F"\t" '{if($3 == 1){print $1}}' 16.tab | awk -F"separatortoto" '{print $1"\n"$2}' - > Unique-16mers-3MM.fasta

##Blast filtered unique 16 mer sequences using "blast short algorithm"
blastn -db c_elegans.PRJNA13758.WS270.genomic.fa -query Unique-16mers-3MM.fasta -task "blastn-short" -out Blast16mers-3MM.tab -outfmt 6 -evalue 1 -num_threads 20

##Get names of seqs that are unique in a set that matches at least 16 bp and start the aligment before the third position
awk -F"\t" '{if((($3*$4)/100)>15){if($7 < 5){array[$1]++}}} END{for (key in array){if(array[key]==1){print key}}}' Blast16mers-3MM.tab > Blast_16mers_filt-4MM.txt

##Parse output to create a file with the indexes of the unique sequences
awk '{if($0 ~ />/){name=$0;}else{print name"\t"$0}}' Unique-16mers-3MM.fasta | awk -F">" '{print $2}' | awk -F"\t" '{if(array[$1] != 0){print ">"$1"\n"array[$1]}else{array[$1]=$2}}' - Blast_16mers_filt-4MM.txt > Cel_gen_4MM_NGGs-seqs.fasta

##Map unique to obtain position
bwa aln -n 0 -o 0 c_elegans.PRJNA13758.WS270.genomic.fa Cel_gen_4MM_NGGs-seqs.fasta | bwa samse c_elegans.PRJNA13758.WS270.genomic.fa - Cel_gen_4MM_NGGs-seqs.fasta > Cel_gen_4MM_NGGs-seqs.sam
grep "XT:A:U" Cel_gen_4MM_NGGs-seqs.sam| awk -F"\t" '{OFS="\t";if($2==16){str="-"}else{str="+"}print $3,($4-1),($4+19),$1,$6,str}' | sort -k1,1 -k2,2n > Cel_gen_4MM_NGGs-seqs.bed




###SFU track fashion
###In a first step, for each PAM site in the C. elegans genome we kept the corresponding adjacent 20-base guide only if its GC content was between 20% and 80% and no poly-T tracts of length 5 or longer were present

mkdir replicate_SFU_track
cd replicate_SFU_track

##WGET WS250
wget ftp://ftp.wormbase.org/pub/wormbase/species/c_elegans/PRJNA13758/sequence/genomic/c_elegans.PRJNA13758.WS250.genomic.fa.gz
wget ftp://ftp.wormbase.org/pub/wormbase/species/c_elegans/PRJNA13758/gff/c_elegans.PRJNA13758.WS250.annotations.gff3.gz


##Find all NGGs first in fwd
zcat c_elegans.PRJNA13758.WS250.genomic.fa.gz | awk '{if($0 ~ />/){print "\n"$0;}else{printf toupper($0)}} END{print ""}' | tail -n+2 | awk '{if($0 ~ />/){name=$0;}else{for(i=1; i<=(length($0)-22); i++){seq=substr($0,i+21,2); if(seq == "GG"){print name"_"i"_+_separatortoto"substr($0,i,23)"\n"substr($0,i,20)}}}}' - > NGG_Cel_WS250_Fwd.fasta
##And then in reverse
zcat c_elegans.PRJNA13758.WS250.genomic.fa.gz | awk '{if($0 ~ />/){print "\n"$0;}else{printf toupper($0)}} END{print ""}' | tail -n+2 | perl -pe '{if(!(/>/m)){chomp; tr/atcgATCG/tagcTAGC/; $_=reverse; $_.="\n"}}' | awk '{if($0 ~ />/){name=$0;}else{for(i=1; i<=(length($0)-22); i++){seq=substr($0,i+21,2); if(seq == "GG"){print name"_"(length($0)-i+1)"_-_separatortoto"substr($0,i,23)"\n"substr($0,i,20)}}}}' - > NGG_Cel_WS250_Rev.fasta
##Combine them
cat NGG_Cel_WS250_Fwd.fasta NGG_Cel_WS250_Rev.fasta > NGG_Cel_WS250.fasta

##GC-Filt
awk '{if($0 ~ />/){name=$0;}else{sumA=0;sumT=0;sumC=0;sumG=0;sumN=0;seq=$0;k=length(seq); for (i=1;i<=k;i++) {if (substr(seq,i,1)=="T") sumT+=1; else if (substr(seq,i,1)=="A") sumA+=1; else if (substr(seq,i,1)=="G") sumG+=1; else if (substr(seq,i,1)=="C") sumC+=1; else if (substr(seq,i,1)=="N") sumN+=1}; GCcontent=(sumC+sumG)/k*100; if((GCcontent >= 20) && (GCcontent <= 80)){print name"\n"seq}}}' NGG_Cel_WS250.fasta > NGG_Cel_WS250-GCFilt.fasta

##Poly-T/a tracks filt
awk '{if($0 ~ />/){name=$0;}else{print name"\t"$0}}' NGG_Cel_WS250-GCFilt.fasta | grep -v "TTTTT" | grep -v "AAAAA" | awk -F"\t" '{print $1"\n"$2}' > NGG_Cel_WS250-GC_TFilt.fasta

###Guides for which the seed region (defined as 12 bases at the 3’ end plus PAM) was not unique in the genome were eliminated.
#Uncompress genome and index it for next steps
gzip -d c_elegans.PRJNA13758.WS250.genomic.fa.gz
##convert to uppercase genome
awk '{print toupper($0)}' c_elegans.PRJNA13758.WS250.genomic.fa > tmp
mv tmp c_elegans.PRJNA13758.WS250.genomic.fa

bwa index c_elegans.PRJNA13758.WS250.genomic.fa

grep ">" NGG_Cel_WS250-GC_TFilt.fasta | awk -F"separatortoto" '{print $1"separatortoto"$2"\n"substr($2,9,15)}' - > NGG_Cel_WS250-GC_TFilt-15mers.fasta

../../QueryInGenomeWithMMinWindow.pl c_elegans.PRJNA13758.WS250.genomic.fa NGG_Cel_WS250-GC_TFilt-15mers.fasta 15 1 13 13 | awk -F"\t" '{if($3==1){print $1}}' | awk -F"separatortoto" '{print $1""substr($2,21,3)"\n"substr($2,1,20)}' - > NGG_Cel_WS250-GC_TFilt-Unique15mers.fasta

##The guide + PAM sequences were then aligned to the whole genome with bwa aln (Li and Durbin 2009 Bioinformatics) allowing an edit distance of 3. Guides mapping to multiple locations in the genome were eliminated.
bwa aln -n 3 -k 3 c_elegans.PRJNA13758.WS250.genomic.fa NGG_Cel_WS250-GC_TFilt-Unique15mers.fasta | bwa samse c_elegans.PRJNA13758.WS250.genomic.fa - NGG_Cel_WS250-GC_TFilt-Unique15mers.fasta > Cel_rep_sfu_NGG-seqs.sam
grep "XT:A:U" Cel_rep_sfu_NGG-seqs.sam | awk -F"\t" '{OFS="\t";if($2==16){str="-"}else{str="+"}print $3,($4-1),($4+19),$1,$6,str}' | sort -k1,1 -k2,2n > Cel_rep_sfu_NGG-seqs.bed

##Compare results
wget http://genome.sfu.ca/crispr/NGG_guides_WS250.bed.gz

wc -l Cel_rep_sfu_NGG-seqs.bed
zcat NGG_guides_WS250.bed.gz | wc -l

##Do 4 MM but with blast
#
##After installing blast, it is required to create a database of the genome to map against
makeblastdb -in c_elegans.PRJNA13758.WS250.genomic.fa -dbtype 'nucl' -title c_elegans.WS250

##Blast filtered unique 16 mer sequences using "blast short algorithm"
blastn -db c_elegans.PRJNA13758.WS250.genomic.fa -query NGG_Cel_WS250-GC_TFilt-Unique15mers.fasta -task "blastn-short" -out Blast15mers.tab -outfmt 6 -evalue 1 -num_threads 20

##Get names of seqs that are unique in a set that matches at least 16 bp and start the aligment before the third position
awk -F"\t" '{if((($3*$4)/100)>15){if($7 < 5){array[$1]++}}} END{for (key in array){if(array[key]==1){print key}}}' Blast15mers.tab > Blast_15mers_filt-4MM.txt

##Parse output to create a file with the indexes of the unique sequences
awk '{if($0 ~ />/){name=$0;}else{print name"\t"$0}}' NGG_Cel_WS250-GC_TFilt-Unique15mers.fasta | awk -F">" '{print $2}' | awk -F"\t" '{if(array[$1] != 0){print ">"$1"\n"array[$1]}else{array[$1]=$2}}' - Blast_15mers_filt-4MM.txt > Cel_gen_4MM_NGGs-repSFU.fasta

##Map unique to obtain position
bwa aln -n 0 -o 0 c_elegans.PRJNA13758.WS270.genomic.fa Cel_gen_4MM_NGGs-repSFU.fasta | bwa samse c_elegans.PRJNA13758.WS270.genomic.fa - Cel_gen_4MM_NGGs-repSFU.fasta > Cel_gen_4MM_NGGs-repSFU.sam
grep "XT:A:U" Cel_gen_4MM_NGGs-repSFU.sam | awk -F"\t" '{OFS="\t";if($2==16){str="-"}else{str="+"}print $3,($4-1),($4+19),$1,$6,str}' | sort -k1,1 -k2,2n > Cel_gen_4MM_NGGs-repSFU.bed





