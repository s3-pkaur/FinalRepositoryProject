# FINAL REPOSITORY PROJECT 
### Lab 3 ###

Firstly, uncompress the Proteome file 
``` bash
gunzip proteomes/*.gz
```
Then place all of the proteins in the same file 
``` bash
cat  proteomes/*.faa > allprotein.fas
```
After placing all of the proteins on the same file, create a BLAST database with all of the proteomes. This BLAST database will then be used to find potential homologs of the query sequence 
``` bash
makeblastdb -in allprotein.fas -dbtype prot
```
Create working directory
```bash
mkdir ~/lab03-$MYGIT/mygene
```
Download the query sequence 
```bash
ncbi-acc-download -F fasta -m protein "NP_001121697.2"
```
Once the query sequence has been downloaded, preform a BLAST search 
```bash
blastp -db ../allprotein.fas -query NP_001121697.2.fa -outfmt 0 -max_hsps 1 -out solutecarrier.blastp.typical.out
```
After preforming the initial BLAST search, filter through the ouput file for putative homologs that have high-scoring alignments with an e-value cutoff of less than 1e-30 using awk 
```bash
awk '{if ($6< 1e-30)print $1 }' solutecarrier.blastp.detail.out > solute_carrier.blastp.detail.filtered.out
```
Find the number of paralogs in each species
```bash
grep -o -E "^[A-Z]\.[a-z]+" solute_carrier.blastp.detail.filtered.out  | sort | uniq -c
```
### Lab 4 ###
Create a working directory for solute carrier gene family in lab 4 
```bash
mkdir ~/lab04-$MYGIT/mygene
```
After creating a working directory, using a seqkit command obtain the sequences from the BLAST output file 
```bash
seqkit grep --pattern-file ~/lab03-$MYGIT/mygene/solute_carrier.blastp.detail.filtered.out ~/lab03-$MYGIT/allprotein.fas | seqkit grep -v -p "carpio" > ~/lab04-$MYGIT/mygene/solutecarrier.homologs.fas
```
Preform a multiple sequence alignment using the muscle program via this command 
```bash
-align ~/lab04-$MYGIT/mygene/solutecarrier.homologs.fas -output ~/lab04-$MYGIT/mygene/solutecarrier.homologs.al.fas
```
After preforming a multiple sequence alignment, view the alignment using the alv command. You can also view the alignment under the majority option with the second command 
```bash
alv -kli  ~/lab04-$MYGIT/mygene/solutecarrier.homologs.al.fas | less -RS
```
```bash
alv -kli --majority ~/lab04-$MYGIT/mygene/solutecarrier.homologs.al.fas | less -RS
```
Once the results have been viewed, convert the file into a pdf
```bash
 Rscript --vanilla ~/lab04-$MYGIT/plotMSA.R  ~/lab04-$MYGIT/mygene/solutecarrier.homologs.al.fas
```
Calculate the width/length of the alignment using the alignbuddy program 
```bash
alignbuddy  -al  ~/lab04-$MYGIT/mygene/solutecarrier.homologs.al.fas
```
Calculate the length of the alignment after removing any gaps
```bash
alignbuddy -trm all  ~/lab04-$MYGIT/mygene/solutecarrier.homologs.al.fas | alignbuddy  -al
```
Calculate the length of the alignment after removing any position that is invariant or completely conserved 
```bash 
alignbuddy -dinv 'ambig' ~/lab04-$MYGIT/mygene/solutecarrier.homologs.al.fas | alignbuddy  -al
```
After calculating the width/length of the alignment using the alignbuddy program, calculate the average percent identity with the help of t-coffee 
```bash
t_coffee -other_pg seq_reformat -in ~/lab04-$MYGIT/mygene/solutecarrier.homologs.al.fas -output sim
```
Now, calculate the average percent identity using the alignbuddy program 
```bash
 alignbuddy -pi ~/lab04-$MYGIT/mygene/solutecarrier.homologs.al.fas | awk ' (NR>2)  { for (i=2;i<=NF  ;i++){ sum+=$i;num++} }
     END{ print(100*sum/num) } '
```
### Lab 5 ###
Firstly, create a working directory
```bash 
mkdir ~/lab05-$MYGIT/mygene
```
Navigate to the working directory
```bash
cd ~/lab05-$MYGIT/mygene
```
Using IQ-TREE find the maximum likelihood tree estimate. First, a tree search is preformed where the optimal amino acid substitution model and amino acid frequencies are found. As that is happening, the branch lengths are estimates in addition to ultrafast bootstrap support levels
```bash
iqtree -s ~/lab05-$MYGIT/mygene/solutecarrier.homologsf.al.fas -bb 1000 -nt 2
```
Once this tree is created, we can view it in the nw_display format with the command below 
```bash
~/lab05-$MYGIT/mygene/solutecarrier.homologsf.al.fas.treefile
```
Finally, convert this into a pdf 
```bash
Rscript --vanilla ~/lab05-$MYGIT/plotUnrooted.R ~/lab05-$MYGIT/mygene/solutecarrier.homologsf.al.fas.treefile ~/lab05-$MYGIT/mygene/solutecarrier.homologsf.al.fas.treefile.pdf 0.4 15
```
Using the gotree program, midpoint rotting the tree that was made earlier using IQ-TREE 
```bash
gotree reroot midpoint -i ~/lab05-$MYGIT/mygene/solutecarrier.homologsf.al.fas.treefile -o ~/lab05-$MYGIT/mygene/solutecarrier.homologsf.al.mid.treefile
```
In order to visualize this a bit better, output the midpoint rooted tree as a graphic 
```bash
nw_order -c n ~/lab05-$MYGIT/mygene/solutecarrier.homologsf.al.mid.treefile | nw_display -w 1000 -b 'opacity:0' -s  >  ~/lab05-$MYGIT/mygene/solutecarrier.homologsf.al.mid.treefile.svg -
```
Finally, convert the svg image into a pdf
```bash
convert  ~/lab05-$MYGIT/mygene/solutecarrier.homologsf.al.mid.treefile.svg  ~/lab05-$MYGIT/mygene/solutecarrier.homologsf.al.mid.treefile.pdf
```
We can convert this midpoint rooted tree into a cladogram with the help of the following command
```bash
nw_order -c n ~/lab05-$MYGIT/mygene/solutecarrier.homologsf.al.mid.treefile | nw_topology - | nw_display -s  -w 1000 > ~/lab05-$MYGIT/mygene/solutecarrier.homologsf.al.midCl.treefile.svg -

convert ~/lab05-$MYGIT/mygene/solutecarrier.homologsf.al.midCl.treefile.svg ~/lab05-$MYGIT/mygene/solutecarrier.homologsf.al.midCl.treefile.pdf
```
### Lab 6 ###
First, perform a gene and species tree reconciliation with the help of a software pacakage known as Notung. This command also saves the output as a png 
```bash
java -jar ~/tools/Notung-3.0_24-beta/Notung-3.0_24-beta.jar -s ~/lab05-$MYGIT/species.tre -g ~/lab06-$MYGIT/mygene/solutecarrier.homologsf.al.mid.treefile --reconcile --speciestag prefix --savepng --events --outputdir ~/lab06-$MYGIT/globins/
```
Convert the reconcilied gene tree file from newick format to RecPhyloXML format by making creating a RecPhyloXML object 
```bash
python2.7 ~/tools/recPhyloXML/python/NOTUNGtoRecPhyloXML.py -g ~/lab06-$MYGIT/mygene/solutecarrier.homologsf.al.mid.treefile.rec.ntg --include.species
```
Using thirdkind, preform a gene-reconicliation-within species tree graphic 
```bash
thirdkind -Iie -D 40 -f ~/lab06-$MYGIT/mygene/solutecarrier.homologsf.al.mid.treefile.rec.ntg.xml -o  ~/lab06-$MYGIT/mygene/solutecarrier.homologsf.al.mid.treefile.rec.svg
```
Finally convert svg graphic into a pdf to view the contents a bit easier
```bash
convert  -density 150 ~/lab06-$MYGIT/mygene/solutecarrier.homologsf.al.mid.treefile.rec.svg ~/lab06-$MYGIT/mygene/solutecarrier.homologsf.al.mid.treefile.rec.pdf
```
### Lab 8 ###
Create a working directory and navigate to it 
```bash
mkdir ~/lab08-$MYGIT/mygene && cd ~/lab08-$MYGIT/mygene
```
Make a copy of the unaligned sequence by removing the asterisk (stop codon) with substitutiong any asterisk with nothing 
```bash
sed 's/*//' ~/lab04-$MYGIT/mygene/solutecarrier.homologs.fas > ~/lab08-$MYGIT/mygene/solutecarrier.homologs.fas
```
Since the Pfam database has already been downloaded, run RPS-Blast using the unaligned sequences as query and the Pfam database with an e-value of .0000000001 as the cutoff 
```bash
rpsblast -query ~/lab08-$MYGIT/mygene/solutecarrier.homologs.fas -db ~/data/Pfam/Pfam -out ~/lab08-$MYGIT/mygene/solutecarrier.rps-blast.out -outfmt "6 qseqid qlen qstart qend evalue stitle" -evalue .0000000001
```
Copy the final gene tree from lab 5, this tree will serve where the predicted Pfam domains will be placed 
```bash
cp ~/lab05-$MYGIT/mygene/solutecarrier.homologsf.al.fas.treefile ~/lab08-$MYGIT/mygene
```
Using the script provided by D4. Rest in R, plot the predicted pfam domains that were obtained from RPS-Blast next to their associated protein on the phylogenetic tree that was copied from lab 5. This commnad will also save the file as a pdf for easy viewing purposes 
```bash
Rscript --vanilla ~/lab08-$MYGIT/plotTreeAndDomains.r ~/lab08-$MYGIT/mygene/solutecarrier.homologsf.al.fas.treefile ~/lab08-$MYGIT/mygene/solutecarrier.rps-blast.out ~/lab08-$MYGIT/mygene/solutecarrier.tree.rps.pdf
```
View the tab delimited annotations from the RPS-Blast using the following command
```bash
mlr --inidx --ifs "\t" --opprint  cat ~/lab08-$MYGIT/mygene/solutecarrier.rps-blast.out | tail -n +2 | less -S
```
Examine the distribution of Pfam domains acorss the proteins. With the following command you can see which proteins have no annotations 
```bash
cut -f 1 ~/lab08-$MYGIT/mygene/solutecarrier.rps-blast.out | sort | uniq -c
```
Using the following command, examine which Pfam domain has the most annotations 
```bash
cut -f 6 ~/lab08-$MYGIT/mygene/solutecarrier.rps-blast.out | sort | uniq -c
```
Using the following command, examine which protein has the longest annotated protein domain 
```bash
awk '{a=$4-$3;print $1,'\t',a;}' ~/lab08-$MYGIT/mygene/solutecarrier.rps-blast.out |  sort  -k2nr
```
Finally with the following command examine which protein has a domain annotation with the best e-value
```bash
cut -f 1,5 -d $'\t' ~/lab08-$MYGIT/mygene/solutecarrier.rps-blast.out
```
# THE END 
