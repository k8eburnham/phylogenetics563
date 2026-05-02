# 2026-02-09: Data Acquisition and Filtering

## Data description: 

UCE (ultraconserved element) DNA data from Borowiec et al. (2024) focusing on ant phylogeny - see access information below
I know phyluce was used upstream to generate these unaligned reads because the authors noted their use of a standard UCE workflow, and the file structure in their data was indicative of phyluce usage. 

Species included (as named in fasta files): Acropyga_sp_FMNHINS0003471605, Anochetus_mayri_FMNHINS0003165063, Camponotus_floridanus_FMNHINS0003145311, Dolichoderus_bispinosus_FMNHINS0003165075, Formica_occulta_FMNHINS0003471213, Iberoformica_subrufa_CASENT0270631, Myrmica_wheeleri_FMNHINS0003471188, Ooceraea_australis_CASENT0106146, Pheidole_flavens_FMNHINS0003145600, Solenopsis_invicta_FMNHINS0003144951 (represents 10 genuses of ants)

## Steps taken to aquire and filter data:

### Data acquisition: 

Downloaded `2-loci.zip` from Dryad (doi: 10.5061/dryad.547d7wmhb) (url: https://datadryad.org/dataset/doi:10.5061/dryad.547d7wmhb) and extracted 3 unaligned loci: `uce-149`, `uce-741`, and `uce-910`. Loci were selected because their UCE codes were all 3 digits (for consistent argument passing) and they all contained reads for the 10 species listed above

### Data filtering: 

Each file above contained information on far too many taxon, so a custom script (scripts/filter_ants.py) was used to extract information from only the 10 species of interest (list above). Files were verified to contain 10 species using `grep -c ">"` in the terminal, because each species name was preceeded by a ">" in the fasta files. Unfiltered reads are named uce-xxx.unaligned.fasta (where xxx is the 3-digit UCE code for each of the 3 loci), and filtered reads (only containing 10 species of interest) are named filtered-uce-xxx.unaligned.fasta. These fasta files still need to be aligned with clustalw/muscle/mafft and combined into one post-alignment file. The filtered and unfiltered files are both stored in data/unaligned_files.

# 2026-02-16: Sequence alignment using muscle and clustalw; comparison

## Software selection criteria

### Muscle: 

Muscle is a newer MSA software that I want to try because it's widely used in the field and relatively easy to download. This is an advantage because it makes my analysis more comparable to existing papers and more relevant to scientists who are interested in phylogenetics. This alignment software works quickly, and is known in the field of phylogenetics to produce more accurate alignments than older softwares. One downside to this software is that it was difficult to install and run on my operating system.

### ClustalW: 

ClustalW is an alignment software that aligns DNA sequences iteratively, which can take a long time on large datasets. This software generates .dnd files, which might be helpful during downstream distance analysis or tree building. I wanted an older alignment software that I could compare my muscle results with. If there isn't a significant difference between the muscle and clustalw alignment products, I'll look into a third alignment software. This software is advantageous because it was very easy to download and run for me, but its age is a disadvantage.

### Other options:

I considered downloading cactus from the in-class alignathon paper (https://genome.cshlp.org/content/24/12/2077), but it was difficult to download (it required the docker application), and it was optimized for mammal species (I'm working with ants).

## Downloading the alignment software (done in terminal)
```
Muscle: conda install -c bioconda muscle
ClustalW: conda install -c bioconda clustalw
```
## Running the alignments

### Muscle: 
```
muscle -align filtered_uce-xxx.unaligned.fasta -output filtered_uce-xxx.aligned_softwarename.fasta (run from inside phylogenetics563/data/unaligned folder in terminal)
```
### ClustalW: 
```
clustalw -INFILE=filtered_uce-xxx.unaligned.fasta -OUTFILE=filtered_uce-xxx.aligned_softwarename.fasta  -ALIGN -OUTPUT=FASTA (run from inside phylogenetics563/data folder in terminal)
```
#### Output files were named filtered_uce-xxx.aligned_softwarename.fasta (where xxx is the 3-digit UCE identifier for each of the 3 loci) and organized into a new folder, data/aligned_files, after being generated

## Comparing the results of alignments using https://alignmentviewer.org/

### at locus 149: 

clustalW looks better across all species (average gap% is 34 with min 9 and max 56, as comapred to average of 51, min of 33, and max of 67 for muscle)

### at locus 741: 

clustalW looks better across all species (average gap% is 45 with min 29 and max 62, as comapred to average of 49, min of 33, and max of 64 for muscle)

### at locus 910: 

clustalW looks better across all species (average gap% is 47 with min 22 and max 65, as comapred to average of 50, min of 27, and max of 67 for muscle)

#### based on these results, I should use the clustalW-generated aligned fasta files going forward

### Note to self: these aligned files still need to be concatenated before generating distance-based, parsimony-based, ML, or Bayesian (but not coalescent) trees

# 2026-03-03: Generating distance- and parsimony-based trees using R (using clustalw-aligned fasta files; see rationale above)

## Distance-based tree generation using R package ape (run from main phylogenetics 563 folder)

### algorithm description: 

Widely used phylogenetic software that can be used to generate distance-based trees. Uses the classical neighbor-joining algorithm and assumes specific substitution models (eg JC69, K80) sufficiently describe sequence divergence. Limitations include the unreliability of distance-based methods and inability to efficiently execute on extremely large trees

### reproducible script (run by opening R in a terminal window from the main project folder):

1) install packages if necessary
```
install.packages("adegenet", dep=TRUE)
install.packages("phangorn", dep=TRUE)
```
2) load the packages
```
library(ape)
library(adegenet)
library(phangorn)
```
3) read the fasta files and sort them so the species are in the same order in every file
```
f149 <- read.dna("data/aligned_files/filtered_uce-149.aligned_clustalw.fasta", format = "fasta")
f149[sort(rownames(f149)), ]
f741 <- read.dna("data/aligned_files/filtered_uce-741.aligned_clustalw.fasta", format = "fasta")
f741[sort(rownames(f741)), ]
f910 <- read.dna("data/aligned_files/filtered_uce-910.aligned_clustalw.fasta", format = "fasta")
f910[sort(rownames(f910)), ]
```
4) concatenate the loci so each species has one long sequence
```
combined_dna <- cbind(f149, f741, f910, fill.with.gaps = TRUE)
dim(combined_dna) #check the dimensions (should be 10 species x total length)
D <- dist.dna(combined_dna, model="TN93") #calculate the distances (the TN93 model of evolution allows for different rates of transitions and transversions and heterogeneous base frequencies, which are common in UCEs)
tre <- nj(D) #get the NJ tree
tre <- ladderize(tre) #reorganize the internal structure of the tree to get the ladderize effect when plotted
```
5) visualize the tree
```
plot(tre, cex=.6)
```
6) visualize and plot the rooted tree
```
rtre = root(tre, outgroup = "Ooceraea_australis_CASENT0106146", resolve.root = TRUE)
plot(rtre)
nodelabels(rtre$node.label)

output_destination <- file.path("figures", "distance_based_tree_rooted_final.pdf") #or send it out to a pdf in the figures folder
pdf(file = output_destination, width = 8, height = 10)
plot(rtre, main = "Rooted distance-based tree using ape software", cex = 0.8, edge.width = 2, edge.color = "darkblue")
add.scale.bar()
dev.off()
```
## Parsimony-based tree generation using R package phangorn

### algorithm description: 

Software with wider capabilities than ape (including maximum likelihood & simple bayesian methods, but we'll be using it for maximum parsimony here). Typically assumes sites evolve independently from one another and the Q matrix is constant across the entire tree. Slower-running than RAxML and IQ-TREE (other max likelihood methods, will cover later)

### reproducible script:

1) install and load the packages if you haven't yet (see distance-based tree steps 1 and 2 above)

2) create the combined_dna object if you haven't yet (see distance-based tree steps 3 and 4 above)

3) generate an initial tree and compute its parsimony score
```
dna2 <- as.phyDat(combined_dna)
tre.ini <- nj(dist.dna(combined_dna,model="raw"))
parsimony(tre.ini, dna2) #output: 7335
tre.pars <- optim.parsimony(tre.ini, dna2) #search for the tree with maximum parsimony, will print something like: Final p-score 7255 after 2 nni operations 
```
4) visualize the tree
```
plot(tre.pars, cex=0.6)
```
5) root the tree and send it out to a file in the figures folder
```
rtre.pars = root(tre.pars, outgroup = "Ooceraea_australis_CASENT0106146", resolve.root = TRUE)
output_destination <- file.path("figures", "parsimony_based_tree_rooted_final.pdf")
pdf(file = output_destination, width = 8, height = 10)
plot(rtre.pars, main = "Rooted parsimony-based tree using phyluce software", cex = 0.8, edge.width = 2, edge.color = "darkblue")
add.scale.bar()
dev.off()
```
### note to self: both of these trees were the same before rooting, but they aren't after

# 2026-03-19: Maximum likelihod tree generation using RAxML

## algorithm description:

RAxML is a software that searches for a tree with the maximum likelihood. It uses a maximum parsimony tree as its starting point (defaults to using 20 starting trees), which are often fast to compute, relatively close to the maximum likelihood tree, and easy to tweak. The user can specify the model of evolution, and the software assumes all branches and sites are evolving at the same rate. It also assumes that input sequences are aligned, homologous, and lack duplicated regions. It has the capability to flag gappy regions, identical sequences, and/or unexpected characters. It searches the tree space using lazy subtree rearrangement, and it converges when the likelihood score stops increasing statistically significantly

## reproducible script - run from the data/RAxML project subfolder

1) Install RAxML using this link (https://github.com/amkozlov/raxml-ng) and move it to the bin so it can be executed from any folder
```
cd Downloads
ls # check the raxml-ng executable is there
sudo cp raxml-ng /usr/local/bin # I needed to override my mac's permissions with sudo to do this
```
2) Check the alignment of my data (using clustalw-aligned fasta files for the rest of this exercise; see rationale above)
```
raxml-ng --check --msa concatenated_uce.fasta --model partitions.txt 
# using partitions.txt as the --model bc it combines the 3 loci & GTR+G is specified within it 
# using GTR+G model because it has flexibility in base frequency and substitution rate while accommodating a gamma distrubtion (TN93 was preferred for distance and parsimony trees because GTR+G wasn't directly applicable to both)
```
3) Run RAxML
```
raxml-ng --msa concatenated_uce.fasta --model partitions.txt --prefix my_uce_tree --threads 4 #using 4 because I checked how many cores my computer has
```
4) Ensure the correct files were saved in RAxML folder (my_uce_tree.raxml.bestTree, my_uce_tree.raxml.mlTrees, my_uce_tree.raxml.bestModel, my_uce_tree.raxml.log)

5) Plot in R just to see, not saving (will save a better version later, see below)
```
library(ape)
tre = read.tree(file="my_uce_tree.raxml.bestTree")
plot(tre)
```
6) re-run RAxML with non-parametric bootstrapping using the flag --all
```
raxml-ng --all --msa concatenated_uce.fasta --model partitions.txt --bs-trees 100 --prefix final_uce_ml_bootstrap --threads 4 #chose 100 boostraps for higher statistical confidence
```
7) plot the bootstrap tree in R, root it (very important!!), and save pdf to figures
```
library(ape)
tre = read.tree(file="final_uce_ml_bootstrap.raxml.support")
plot(tre)
nodelabels()

rtre = root(tre, outgroup = "Ooceraea_australis_CASENT0106146", resolve.root = TRUE)
plot(rtre)
nodelabels(rtre$node.label)

output_destination <- file.path("../../figures", "RAxML_boostrap_tree_rooted_final.pdf")
pdf(file = output_destination, width = 8, height = 10) #move 2 levels out from RAxML folder, into figures folder
plot(rtre, main = "Rooted Maximum Likelihood Tree using RAxML Boostrapping", cex = 0.8, edge.width = 2, edge.color = "darkblue")
add.scale.bar()
dev.off()
```
# 2026-04-09: MrBayes

## algorithm description, pros and cons

### How it works: 

MrBayes uses markov chain monte carlo to sample the space of possible evolutionary trees and builds a posterior probability distribution of how likely it is that specific clades are related, then reports a tree based on that posterior probability. This is different than maximum likelihood methods, which try to find the single best tree given the data.

### Assumptions: 

MrBayes assumes that the user-inputted model of evolution is correct, that rates of evolution remain constant over the tree and time, and that changes at one DNA site don't affect changes at another site.

### Limitations: 

MrBayes is computationally expensive, so it takes much longer to run than maximum likelihood methods, even on my relatively small dataset. The software originally estimated that it would take 23 hours to run, but it ended up taking less than an hour for me - maybe because I used multiple cores. MRBayes also requires the user to monitor ASDSF and ESS statistics to know when to stop the run, instead of having an internal on/off switch. The software can also be sensitive to priors, especially if the dataset is small.

### use guide I'm using: 

https://evomics.org/learning/phylogenetics/mrbayes/ (I couldn't access the one linked in the homework guidelines)

## reproducible script

1) convert concatenated, aligned .fasta file to .nex file; nexus format is input standard for MrBayes (run from phylogenetics563 folder)
```
pip install biopython
python -c "from Bio import SeqIO; SeqIO.convert('data/RAxML/concatenated_uce.fasta', 'fasta', 'data/MrBayes/concatenated_uce.nex', 'nexus', molecule_type='DNA')"
```
2) create MrBayes block text file with appropriate parameters

``` begin mrbayes;
 set autoclose=yes; 
 prset brlenspr=unconstrained:exp(10.0); #I expect my branch lengths to be a mix of short and long, since I have some closely and some more distantly related species, so the MrBayes default is appropriate
 prset shapepr=exp(1.0); #exponential prior for gamma shape
 prset tratiopr=beta(1.0,1.0); #uniform Beta prior for Ti/Tv ratio
 prset statefreqpr = fixed(empirical); #my data has a high AT skew, so I want my model to fix the state frequencies to those calculated directly from my alignment
 lset nst=6 rates=gamma ngammacat=4; #my data is 22k base pairs, so it warrants using the GTR model, allowing for 6 different rates of substitutions (transitions, transversions, etc)
 mcmcp ngen=1000000 samplefreq=10 printfreq=100 nruns=2 nchains=4 savebrlens=yes; #using 1 million generations to mirror industry standard, using 2 runs to allow analysis of if the software arrived at the same tree twice, using 4 chains to allow for metropolis coupling, where one chain stays cold
 outgroup Ooceraea_australis_CASENT0106146;
 mcmc;
 sumt;
end;
```

3) append MrBayes block to the end of nexus file - I did this manually by copy-pasting the mb block (above) to the end of the nexus file in data/MrBayes because the cat function wasn't working for me

The final nexus file with the MrBayes block for use in running MrBayes is called data/MrBayes/concatenated_uce_with_mbb.nex. The outgroup specified in the MrBayes block above and in this file has been updated to be biologically accurate.

4) Run MrBayes (run from data/MrBayes folder)
```
mpirun -np 8 mb concatenated_uce_with_mbb.nex #using all 8 of my computer's cores to speed it up
```
5) Assess if the chain has converged using the Tracer app and the .p files

#are the file names the same as last time before I changed the outgroup? oops oh well
#this should be redone since I reran mrbayes
My tracer plot has a hairy caterpiller look, an ESS of 1100/1600 (good threshold is 200), and good visual mixing of the two runs. It doesn't really show a signficant burn in period, but instead it levels off very quickly. The marginal density plots look similar, but not identical.

# 2026-04-21: Astral (coalescent)

## algorithm description, pros and cons

### How it works: 

Astral treats gene trees as raw data (instead of nucelotide sequence information) to find a species tree from the tree space that best explains the gene trees

### Assumptions: 

It assumes the gene trees you input are accurate. Any species trees generated from inaccurate gene trees are not guaranteed to be correct. It also assumes there isn't any horizontal gene flow between species, and that the primary species signal onflict comes from incomplete lineage sorting

### Limitations: 

Excessive incomplete lineage sorting, if present in the data, can cause astral to inaccurately predict species relationships. Also, if genes aren't evolving independently, the results can be biased. Finally, astral is very sensitive to errors in the accuracy of the gene trees it receives, since these are its raw data.

### instructions I'm following to install and run: 

https://github.com/smirarab/ASTRAL/blob/master/astral-tutorial.md 

## reproducible script

1) install software in terminal
```
java -jar astral.5.7.8.jar
```
2) locate and prepare input files (gene trees) and decide where to store output files (species trees)

run these from the data folder; these commands generate 3 gene trees for use in Astral using RAxML on the .fasta files created by clustalw alignment (necessary because my previous RAxML and MrBayes work was on concatenated files, which doesn't work with Astral)
using standard substitution model GTRGAMMA and random seed 12345 for now
```
raxml-ng --msa aligned_files/filtered_uce-149.aligned_clustalw.fasta --model GTR+G --prefix astral/project_data/locus149 --seed 12345 --threads 2
raxml-ng --msa aligned_files/filtered_uce-741.aligned_clustalw.fasta --model GTR+G --prefix astral/project_data/locus741 --seed 12345 --threads 2
raxml-ng --msa aligned_files/filtered_uce-910.aligned_clustalw.fasta --model GTR+G --prefix astral/project_data/locus910 --seed 12345 --threads 2
```
3) combine the raxml output files (individual gene trees) into a combined gene tree for use in astral (also run from my data folder in terminal)
```
cat astral/project_data/*.raxml.bestTree > astral/project_data/combined_gene_trees.tre
```
4) run Astral on the combined gene trees - run this in terminal from the data/Astral folder!
```
java -jar astral.5.7.8.jar -i project_data/combined_gene_trees.tre -o project_data/astral_output_species_tree.tre
```
5) Analyze the results

My final normalized quartet score (printed out in terminal by running astral) was 0.8920634920634921, which suggests some minor incomplete lineage sorting, but it's a typical amount for empirical datasets and indicates that the primary species signal is very clear

6) Root and save the output tree (located at data/Astral/project_data/astral_output_species_tree.tre) to the figures folder using R
```
library(ape)
tre = read.tree(file="data/Astral/project_data/astral_output_species_tree.tre")
rtre = root(tre, outgroup = "Ooceraea_australis_CASENT0106146", resolve.root = TRUE)
output_destination <- file.path("figures", "Astral_tree_rooted_final.pdf")
pdf(file = output_destination, width = 8, height = 10) 
plot(rtre, main = "Rooted Coalescent Tree using Astral", cex = 0.8, edge.width = 2, edge.color = "darkblue")
#note to self: R says branch lengths are ignored when doing this plotting
add.scale.bar()
dev.off()
```
