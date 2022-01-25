
### MEERKAT MITOGENOME ASSEMBLY - NOVOPplasty3.7. ###

# https://github.com/ndierckx/NOVOPlasty
# DO NOT filter or quality trim the reads!!! Use the raw whole genome dataset (Only adapters should be removed)!
# You can subsample to speed up the process and to reduce the memory requirements
# perl NOVOPlasty3.7 -c config.txt


# Cutadapt v.2.4
# Trim Galore v.0.6.4

# remove adapters
PATH_TO_FILES=/home/derezanin/species_comp/macrogen_sequences/HN00104407

trim_galore --paired --phred33 -j 8 --basename meerkat_adpt_cut -o adapt_cut/ \
$PATH_TO_FILES/30318_1.fastq.gz $PATH_TO_FILES/30318_2.fastq.gz

time NOVOPlasty3.7 -c config1_meerkat.txt
time NOVOPlasty3.7 -c config2_meerkat.txt
time NOVOPlasty3.7 -c config3_meerkat.txt
# time: 270 min (11.10.19)

# without max mem set, uses up to 150GB RAM, max mem set to 20 - subsamples read set but doesn't circulize the genome, outputs 2 contigs
# use whole read set

# downsampling read set to 20%
# changed python version in the downsample.py script to /usr/local/bin/python2.7, due to the error with python2.6.(default on our server)
time /home/derezanin/software/MITObim/misc_scripts/downsample.py -s 20 --interleave \
-r meerkat_adpt_cut_R1_val_1.fq.gz -r meerkat_adpt_cut_R2_val_2.fq.gz | gzip > meerkat_dwnsmp20_intlvd.fq.gz
# time: 288 min (16.10.19)


## Bowtie2 mapping ##

# mapping downsampled read_set back to meerkat circ. mitoassemblies

READS=/data/fg2/derezanin/species_comp/feliformes/herpestidae/meerkat_mitogenome/novoplasty_run/adapt_cut

# bowtie2-build refs
MEERKAT_ASM=/data/fg2/derezanin/species_comp/feliformes/herpestidae/meerkat_mitogenome/novoplasty_run/finished_mitos

time bowtie2 -q --end-to-end -p 6 --no-unal -x $MEERKAT_ASM/smamong_meerkat_asm --interleaved $READS/meerkat_dwnsmp20_intlv.fq.gz \
-S sma_mong_meerkatdwnsmp20.sam

time bowtie2 -q --end-to-end -p 6 --no-unal -x $MEERKAT_ASM/eupmong_meerkat_asm --interleaved $READS/meerkat_dwnsmp20_intlv.fq.gz \
-S eup_mong_meerkatdwnsmp20.sam

# repeated novoplasty meerkat mitoassembly run for the sm.a.mongoose (only as a seed seq. in config, removed as ref) - same output


# renamed Meerkat_sma_mong_seed-only_consensus.fasta to Meerkat.fa in mafft_msa dir

# align 22 mitogenomes with MAFFT (MAFFT v7.407 (2018/Jul/23), installed in conda env "partition finder")

mafft --auto --thread 10 all_mitos.fasta --phylipout > all_mitos.phy
# doesn't work , logfile ->  mafft_auto.log

# Accuracy-oriented methods:
# *G-INS-i (suitable for sequences of similar lengths; recommended for <200 sequences; 
# iterative refinement method incorporating global pairwise alignment information):
# mafft --globalpair --maxiterate 1000 input [> output] 
# ginsi input [> output] 

# --globalpair - All pairwise alignments are computed with the Needleman-Wunsch algorithm. 
# Suitable for a set of globally alignable sequences. Applicable to up to ~200 sequences. 
# A combination with --maxiterate 1000 is recommended (G-INS-i). Default: off (6mer distance is used) 
mafft --large --globalpair --thread 8 --nuc all_mitos.fasta --phylipout > all_mitos.out 
# Unknown option:  all_mitos.fasta
# /home/derezanin/miniconda3/envs/partition_finder/bin/mafft: Cannot open --phylipout.
# run without specifying the output format
mafft --large --globalpair --thread 8 --nuc all_mitos.fasta > all_mitos.out 
# outformat needs to be specified  before input file
mafft --large --globalpair --thread 8 --clustalout --nuc all_mitos.fasta > all_mitos.clust

mafft --large --globalpair --thread 10 --phylipout --nuc all_mitos.fasta > all_mitos.phylip
# no output
mafft --globalpair --thread 10 --phylipout --nuc all_mitos.fasta > all_mitos.phylip

# --maxiterate 1000 recommended for both G-INS-i and E-INS-i methods
mafft --large --globalpair --thread 10 --maxiterate 1000 --phylipout --nuc all_mitos.fasta > all_mitos_i1000.phy
# Iterative refinment is not supported for --memsavetree
# run without --large
mafft --globalpair --thread 10 --maxiterate 1000 --phylipout --nuc all_mitos.fasta > all_mitos_i1000.phy

# The default gap scoring scheme has been changed in version 7.110 (2013 Oct).
# It tends to insert more gaps into gap-rich regions than previous versions.
# To disable this change, add the --leavegappyregion option.
mafft --globalpair --thread 10 --maxiterate 1000 --phylipout --leavegappyregion --nuc all_mitos.fasta > all_mitos_i1000gappy.phy


# *E-INS-i (suitable for sequences containing large unalignable regions; recommended for <200 sequences):
# mafft --ep 0 --genafpair --maxiterate 1000 input [> output]
# einsi input [> output]
# For E-INS-i, the --ep 0 option is recommended to allow large gaps. 
mafft --ep 0 --genafpair --thread 10 --maxiterate 1000 --phylipout --nuc all_mitos.fasta > all_mitos_einsi.phy

mafft --ep 0 --genafpair --thread 10 --maxiterate 1000 --phylipout --leavegappyregion --nuc all_mitos.fasta > all_mitos_einsi_gappy.phy

# rename mito sequences, run mafft again
#ginsi
mafft --globalpair --thread 10 --phylipout --nuc merged_mitos.fa > merged_mitos.phy
mafft --globalpair --thread 10 --maxiterate 1000 --phylipout --nuc merged_mitos.fa > merged_mitos_i1000.phy
mafft --globalpair --thread 10 --maxiterate 1000 --phylipout --leavegappyregion --nuc merged_mitos.fa > merged_mitos_i1000gappy.phy

#einsi
mafft --ep 0 --genafpair --thread 10 --maxiterate 1000 --phylipout --nuc merged_mitos.fa > merged_mitos_einsi.phy
mafft --ep 0 --genafpair --thread 10 --maxiterate 1000 --phylipout --leavegappyregion --nuc merged_mitos.fa > merged_mitos_einsi_gappy.phy

# check msa in Jalview (v.2.11.0-j8 installed locally)
# einsi_gappy msa looks the best (least gaps, better resolved regions around indels)
# try to run partition finder on the msa

PART_FINDER=/home/derezanin/software/partitionfinder-2.1.1/

$PART_FINDER/PartitionFinder.py /home/derezanin/species_comp/feliformes/herpestidae/meerkat_mitogenome/partition_finder -p 8 --raxml
# doesnt run, phylip format without whitespace
# shorten names of species, run again
mafft --ep 0 --genafpair --thread 10 --maxiterate 1000 --phylipout --leavegappyregion --nuc merged_mitos.fa > new_runs/merged_mitos_einsi_gappy.phy
# same partition finder error, phylip out in mafft not adequate, write to mafft dev
# get fasta output from mafft
mafft --ep 0 --genafpair --thread 10 --maxiterate 1000 --leavegappyregion --nuc merged_mitos.fa > new_runs/merged_mitos_einsi_gappy.fa
# convert fasta aln to phylip aln with script: https://github.com/npchar/Phylogenomic/blob/master/fasta2relaxedPhylip.pl

./fasta2relaxedPhylip.pl -f merged_mitos_einsi_gappy.fa -o merged_mitos_einsi_gappy.phy
# works, but changes species order every time
# run PF again
# error hosei separated with space, should be underscore, fix, run again

# ERROR    | 2019-11-06 11:18:09,491 | You cannot estimate a Maximum Likelihood (ML) starting tree (the default behaviour) when you have columns missing from your data block definitions, 
# because the method we use to estimate the ML tree requires all sites in the alignment to be assigned to a data block. We recommend that you either remove the sites you don't want from y
# our alignment or (if possible) include the missing sites in appropriate data blocks. Failing that, you can use the --no-ml-tree command line option. In this case, a NJ (PhyML) or MP(Rax
# ML) starting tree will be estimated for your analysis.
$PART_FINDER/PartitionFinder.py /home/derezanin/species_comp/feliformes/herpestidae/meerkat_mitogenome/partition_finder --no-ml-tree -p 8 --raxml
# error GLIBC.2.14 not found, recompile partition finder from source file from Github - doesn't work
# install glibc in home folder, build from there 
# https://stackoverflow.com/questions/35616650/how-to-upgrade-glibc-from-version-2-12-to-2-14-on-centos
# hard without sudo permissions, get new raxml version
#reinstalled raxml, overwritten executables in PF programs
# ran PF again - everything works

# cat fixed ref mitos and meerkat consensus 
# run mafft again
# convert to phylip format
# input data blocks in PF config file to build ML tree with raxml
mafft --ep 0 --genafpair --thread 10 --maxiterate 1000 --leavegappyregion --nuc merged_mitos.fa > merged_mitos_einsi_gappy_fixed.fa
# sort and cut D-loop in Jalview/Aliview

./fasta2relaxedPhylip.pl -f merged_mitos_einsi_gappy_fixed_srt.fa -o merged_mitos_einsi_gappy_fixed_srt.phy
$PART_FINDER/PartitionFinder.py /home/derezanin/species_comp/feliformes/herpestidae/meerkat_mitogenome/partition_finder --no-ml-tree -p 8 --raxml

# assign all sites, remove \3 from trna/rrna lines, try out ML tree 
$PART_FINDER/PartitionFinder.py /home/derezanin/species_comp/feliformes/herpestidae/meerkat_mitogenome/partition_finder -p 8 --raxml
# check starting tree


# run raxml separately with 1000 bootstrap replicates
/home/derezanin/software/standard-RAxML/raxmlHPC-PTHREADS-SSE3 -f a -T 12 -p 250187 -x 2511987 -N 1000 --bootstop-perms=1000 -m GTRGAMMA -s merged_mitos_einsi_gappy_fixed_srt.phy \
-n raxml_results1000.tree 

# -f a - rapid Bootstrap analysis and search for best­ scoring ML tree in one program run
# -T # of threads
# -p random # seed for parsimony inference
# -x Specify random # seed and turn on rapid bootstrapping CAUTION:   unlike   in   previous   versions   of   \
# RAxML   will   conduct   rapid   BS  replicates under the model of rate heterogeneity you specified via ­m and not by default under CAT
# -b random # seed for bootstraping
# -N # of replicates
# -bootstop­-perms=number specify the number of permutations to be conducted for the bootstopping/bootstrap convergence test, allowed min 100, set to 1000
# -m selection model
# -t specify a user starting tree in Newick format
# -w absolute path to the dir to write output files 

# specify starting tree from PF run
/home/derezanin/software/standard-RAxML/raxmlHPC-PTHREADS-SSE3 -f a -T 12 -p 250187 -x 2511987 -N 1000 -m GTRGAMMA \
-t analysis/start_tree/RAxML_result.BLTREE -s merged_mitos_einsi_gappy_fixed_srt.phy \
-n raxml_results_starttree.tree 
# starting tree ignored by rapid bootstrapping

# compare raxml trees in Figtree - bootstrap values - done
# Herpestes auropunctatus doesn't cluster with H. javanica, check sequences in Geneious, align RefSeq version of H.auropunctatus

# run mafft without H. auropunctatus
mafft --ep 0 --genafpair --thread 8 --maxiterate 1000 --leavegappyregion --nuc merged_mitos_no_ind_mong.fa > merged_mitos_einsi_gappy_no_ind.fa

# sort in Jalview

# convert to phylip
./fasta2relaxedPhylip.pl -f merged_mitos_einsi_gappy_no_ind_srt.fa -o merged_mitos_einsi_gappy_no_ind_srt.phy


# try tree with GTRGAMMAI on all mitos
/home/derezanin/software/standard-RAxML/raxmlHPC-PTHREADS-SSE3 -f a -T 8 -p 250187 -x 2511987 -N 1000 --bootstop-perms=1000 -m GTRGAMMAI -s merged_mitos_einsi_gappy_fixed_srt.phy \
-n raxml_results1000gtrgammai.tree 

# try with GTRCAT
/home/derezanin/software/standard-RAxML/raxmlHPC-PTHREADS-SSE3 -f a -T 8 -p 250187 -x 2511987 -N 1000 --bootstop-perms=1000 -m GTRCAT -s merged_mitos_einsi_gappy_fixed_srt.phy \
-n raxml_results1000gtrcat.tree 
# GTR + Optimization of substitution rates + Optimization of site-specific evolutionary rates 
# which are categorized into numberOfCategories distinct rate categories for greater computational efficiency. 
# Final tree might be evaluated under GTRGAMMA, depending on the tree search option.
time bowtie2 -q --phred33 --end-to-end -p 6 --no-unal -x meerkat_fixed_consensus --interleaved \
/home/derezanin/species_comp/feliformes/herpestidae/meerkat_mitogenome/novoplasty_run/adapt_cut/meerkat_dwnsmp20_intlv.fq.gz \
-S sma_mong_meerkatdwnsmp20.sam

