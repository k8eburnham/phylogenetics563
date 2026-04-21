# 2026-02-09: Data Acquisition and Filtering

## Data description: 

### UCE data from Borowiec et al. (2024) focusing on ant phylogeny; I know phyluce was used upstream to generate these unaligned reads because the authors noted their use of a standard UCE workflow, and the file structure in their data was indicative of phyluce usage. 
### Species included (as named in fasta files): Acropyga_sp_FMNHINS0003471605, Anochetus_mayri_FMNHINS0003165063, Camponotus_floridanus_FMNHINS0003145311, Dolichoderus_bispinosus_FMNHINS0003165075, Formica_occulta_FMNHINS0003471213, Iberoformica_subrufa_CASENT0270631, Myrmica_wheeleri_FMNHINS0003471188, Ooceraea_australis_CASENT0106146, Pheidole_flavens_FMNHINS0003145600, Solenopsis_invicta_FMNHINS0003144951

## Steps taken to aquire and filter data:

### Data acquisition: Downloaded `2-loci.zip` from Dryad (doi:10.5061/dryad.547d7wmhb) and extracted 3 unaligned loci: `uce-149`, `uce-741`, and `uce-910`. Loci were selected because their UCE codes were all 3 digits and they all contained reads for the 10 species listed above
### Data filtering: Each file above contained information on far too many taxon, so a custom script (filter_ants.py) was used to extract information from only the 10 species of interest (list above). Files were verified to contain 10 species using `grep -c ">"`. Unfiltered reads are named uce-xxx.unaligned.fasta, and filtered reads (only containing 10 species of interest) are named filtered-uce-xxx.unaligned.fasta. These fasta files still need to be aligned with clustalw/muscle/mafft and combined into one post-alignment file.


# 2026-02-16: Sequence alignment using muscle and clustalw; comparison

## Software selection criteria

### Muscle: This is a newer MSA software that I want to try because it's widely used in the field and relatively easy to download.
### ClustalW: This software generates .dnd files, which might be helpful during downstream distance analysis or tree building. I wanted an older alignment software that I could compare my muscle results with. If there isn't a significant difference between the muscle and clustalw alignment products, I'll look into a third alignment software.
### No softwares from alignathon paper (https://genome.cshlp.org/content/24/12/2077): The cactus software was an attractive option, but it was optimized for mammal species, and I'm working with ants. I had also difficulty trying to download it (it required the docker application).

## Downloading the alignment software

### Muscle: conda install -c bioconda muscle
### ClustalW: conda install -c bioconda clustalw

## Running the alignments

### Muscle: muscle -align input.fasta -output output.fasta (run from inside phylogenetics563/data folder in terminal)
### ClustalW: clustalw -INFILE=input.fasta -OUTFILE=output.fasta  -ALIGN -OUTPUT=FASTA (run from inside phylogenetics563/data folder in terminal)
### Output files were named filtered_uce-***.aligned_softwarename.fasta and organized into the aligned_files folder

## Comparing the results of alignments using https://alignmentviewer.org/

### at locus 149: clustalW looks better across all species (average gap% is 34 with min 9 and max 56, as comapred to average of 51, min of 33, and max of 67 for muscle)
### at locus 741: clustalW looks better across all species (average gap% is 45 with min 29 and max 62, as comapred to average of 49, min of 33, and max of 64 for muscle)
### at locus 910: clustalW looks better across all species (average gap% is 47 with min 22 and max 65, as comapred to average of 50, min of 27, and max of 67 for muscle)

### based on these results, I should use the clustalW-generated aligned fasta files going forward
### other existing mechanisms to compare alignment quality: gap consolidation, average identity in the "Information Content" or "Logo" section, handling of ends

## Note to self: these files still need to be combined into one master file for each locus before I build my phylogenetic tree (concatenation) - is this just for distance/parsimony methods, or ML too?


# 2026-03-03: Generating distance- and parsimony-based trees using R (using clustalw-aligned fasta files; see rationale above)

## Distance-based tree generation using R package ape (run from main phylogenetics 563 folder)

### algorithm description: Widely used phylogenetic software that can be used to generate distance-based trees. Uses the classical neighbor-joining algorithm and assumes specific substitution models (eg JC69, K80) sufficiently describe sequence divergence. Limitations include the unreliability of distance-based methods and inability to efficiently execute on extremely large trees

### #install packages if necessary
### install.packages("adegenet", dep=TRUE)
### install.packages("phangorn", dep=TRUE)

### #load the packages
### library(ape)
### library(adegenet)
### library(phangorn)

### #read the fasta files and sort them so the species are in the same order in every file
### f149 <- read.dna("data/aligned_files/filtered_uce-149.aligned_clustalw.fasta", format = "fasta")
### f149[sort(rownames(f149)), ]
### f741 <- read.dna("data/aligned_files/filtered_uce-741.aligned_clustalw.fasta", format = "fasta")
### f741[sort(rownames(f741)), ]
### f910 <- read.dna("data/aligned_files/filtered_uce-910.aligned_clustalw.fasta", format = "fasta")
### f910[sort(rownames(f910)), ]

### #concatenate the loci so each species has one long sequence
### combined_dna <- cbind(f149, f741, f910, fill.with.gaps = TRUE)
### #check the dimensions (should be 10 species x total length)
### dim(combined_dna)
### #calculate the distances (the TN93 model of evolution allows for different rates of transitions and transversions, heterogeneous base frequencies, and between-site variation of the substitution rate)
### D <- dist.dna(combined_dna, model="TN93")
### #get the NJ tree
### tre <- nj(D)
### #reorganize the internal structure of the tree to get the ladderize effect when plotted
### tre <- ladderize(tre)

### #visualize the tree
### plot(tre, cex=.6)
### #or send it out to a pdf in the figures folder
### output_destination <- file.path("figures", "distance_based_tree.pdf")
### pdf(file = output_destination, width = 8, height = 10)
### plot(tre, main = "distance-based tree using ape software", cex = 0.8, edge.width = 2, edge.color = "darkblue")
### add.scale.bar()
### dev.off()

## Parsimony-based tree generation using R package phangorn

### algorithm description: Software with wider capabilities than ape (including maximum likelihood & simple bayesian methods, but we'll be using it for maximum parsimony here). Typically assumes sites evolve independently from one another and the Q matrix is constant across the entire tree. Slower-running than RAxML and IQ-TREE (other max likelihood methods, will cover later)

### load the packages if you haven't yet (see above)
### dna2 <- as.phyDat(combined_dna) #this is the object created using the ape software - create it if you haven't yet
### #generate an initial tree and compute its parsimony score
### tre.ini <- nj(dist.dna(combined_dna,model="raw"))
### parsimony(tre.ini, dna2) #output: 7335
### #search for the tree with maximum parsimony
### tre.pars <- optim.parsimony(tre.ini, dna2)
### #will print something like: Final p-score 7255 after 2 nni operations 

### #visualize the tree
### plot(tre.pars, cex=0.6)
### #or send it out to a pdf in the figures folder
### output_destination <- file.path("figures", "parsimony_based_tree.pdf")
### pdf(file = output_destination, width = 8, height = 10)
### plot(tre, main = "parsimony-based tree using phyluce software", cex = 0.8, edge.width = 2, edge.color = "darkblue")
### add.scale.bar()
### dev.off()

### note to self: both of these trees are the same, yay!

# 2026-03-19: Maximum likelihod ___ using RAxML on my data

## algorithm description, pros and cons -- **improve this description for final 

### Starting tree: maximum parsimony tree; defaults to using 20 starting trees
### Fast to compute, good enough, maybe not the most accurate so we’ll tweak it

### Model of evolution: chosen based on the data? Or not, sounds like it uses GTR, which is one of the more complex models
### CRL: Done internally, just like IQ-tree
### All ML methods have those 3 assumptions - all branches, sites evolve the same

### Data quality: assumes input sequences are aligned, homologous, no duplicate sequences (orthology step?), can flag gappy regions, ID sequences and/or unexpected characters

### Convergence (end point): searches tree space using lazy subtree rearrangement; stops when likelihood score stops increasing statistically significantly; unclear on methods to ensure reliable results

## reproducible script

### 1. Install RAxML using this link (https://github.com/amkozlov/raxml-ng) and move it to the bin so it can be executed from any folder

cd Downloads
ls # check the raxml-ng executable is there
sudo cp raxml-ng /usr/local/bin # I needed to override my mac's permissions with sudo to do this

### 2. Check the alignment of my data (using clustalw-aligned fasta files for the rest of this exercise; see rationale above)

raxml-ng --check --msa concatenated_uce.fasta --model partitions.txt # using partitions.txt as the --model bc it combines the 3 loci & GTR+G is specified within it #could use some justification for the GTR+G model over LG+G8+F

### 3. Run RAxML

raxml-ng --msa concatenated_uce.fasta --model partitions.txt --prefix my_uce_tree --threads 4 #using 4 because I checked how many cores my computer has

### 4. Ensure the correct files were saved in RAxML folder: my_uce_tree.raxml.bestTree, my_uce_tree.raxml.mlTrees, my_uce_tree.raxml.bestModel, my_uce_tree.raxml.log

### 5. Plot in R just to see, not saving

library(ape)
tre = read.tree(file="my_uce_tree.raxml.bestTree")
plot(tre)

### 6. re-run RAxML with non-parametric bootstrapping using the flag --all

raxml-ng --all --msa concatenated_uce.fasta --model partitions.txt --bs-trees 10 --prefix uce_ml_bootstrap --threads 4 
## note to self: for final, increase bootstrapping number to 100 or use --bs-trees autoMRE to let the program decide

### 7. plot the tree in R and root it (very important!!) - not saving to figures now, but should for final project

library(ape)
tre = read.tree(file="primatesAA-aligned-muscle-raxml-boostrap.raxml.support")
plot(tre)
nodelabels()

rtre = root(tre, outgroup = "Acropyga_sp_FMNHINS0003471605", resolve.root = TRUE) # choosing a species because using a node number was throwing an error - would appreciate help on confirming if this is the right species to use as outgroup for rooting purposes
plot(rtre)
nodelabels(rtre$node.label)


# 2026-04-09: MrBayes

## algorithm description, pros and cons

### How it works: MrBayes uses markov chain monte carlo to sample the space of possible evolutionary trees and builds a posterior probability distribution of how likely it is that specific clades are related, then reports a tree based on that posterior probability. This is different than maximum likelihood methods, which try to find the single best tree given the data.
### Assumptions: MrBayes assumes that the user-inputted model of evolution is correct, that rates of evolution remain constant over the tree and time, and that changes at one DNA site don't affect changes at another site.
### Limitations: MrBayes is computationally expensive, so it takes much longer to run than maximum likelihood methods, even on my relatively small dataset. The software originally estimated that it would take 23 hours to run, but it ended up taking less than an hour for me - maybe because I used multiple cores. MRBayes also requires the user to monitor ASDSF and ESS statistics to know when to stop the run, instead of having an internal on/off switch. The software can also be sensitive to priors, especially if the dataset is small.
### use guide I'm using: https://evomics.org/learning/phylogenetics/mrbayes/ (I couldn't access the one linked in the homework guideliens)

## reproducible script

### 1. convert concatenated, aligned .fasta file to .nex file; nexus format is input standard for MrBayes (run from phylogenetics563 folder)

pip install biopython
python -c "from Bio import SeqIO; SeqIO.convert('data/RAxML/concatenated_uce.fasta', 'fasta', 'data/MrBayes/concatenated_uce.nex', 'nexus', molecule_type='DNA')"

### 2. create MrBayes block text file with appropriate parameters

begin mrbayes;
 set autoclose=yes; 
 prset brlenspr=unconstrained:exp(10.0); #I expect my branch lengths to be a mix of short and long, since I have some closely and some more distantly related species, so the MrBayes default is appropriate
 prset shapepr=exp(1.0); #exponential prior for gamma shape
 prset tratiopr=beta(1.0,1.0); #uniform Beta prior for Ti/Tv ratio
 prset statefreqpr = fixed(empirical); #my data has a high AT skew, so I want my model to fix the state frequencies to those calculated directly from my alignment
 lset nst=6 rates=gamma ngammacat=4; #my data is 22k base pairs, so it warrants using the GTR model, allowing for 6 different rates of substitutions (transitions, transversions, etc)
 mcmcp ngen=1000000 samplefreq=10 printfreq=100 nruns=2 nchains=4 savebrlens=yes; #using 1 million generations to mirror industry standard, using 2 runs to allow analysis of if the software arrived at the same tree twice, using 4 chains to allow for metropolis coupling, where one chain stays cold
 outgroup Anochetus_mayri_FMNHINS0003165063; #belongs to a different subfamily than the rest of my species
 mcmc;
 sumt;
end;

### 3. append MrBayes block to the end of nexus file - I did this manually by copy-pasting the mb block to the end of the nexus file because the cat function wasn't working for me

### 4. Run MrBayes - using all 8 of my computer's cores to speed it up

mpirun -np 8 mb concatenated_uce_with_mbb.nex

### 5. Assess if the chain has converged using Tracer and the .p files

### my tracer plot has a hairy caterpiller look, an ESS of 1100/1600 (good threshold is 200), and good visual mixing of the two runs. It doesn't really show a signficant burn in period, but instead it levels off very quickly. The marginal density plots look similar, but not identical.

# 2026-04-21: Astral (coalescent)

## algorithm description, pros and cons

### How it works: Astral treats gene trees as raw data (instead of nucelotide sequence information) to find a species tree from the tree space that best explains the gene trees
### Assumptions: It assumes the gene trees you input are accurate. Any species trees generated from inaccurate gene trees are not guaranteed to be correct. It also assumes there isn't any horizontal gene flow between species, and that the primary species signal onflict comes from incomplete lineage sorting
### Limitations: Excessive incomplete lineage sorting, if present in the data, can cause astral to inaccurately predict species relationships. Also, if genes aren't evolving independently, the results can be biased. Finally, astral is very sensitive to errors in the accuracy of the gene trees it receives, since these are its raw data.
### instructions I'm following to install and run: https://github.com/smirarab/ASTRAL/blob/master/astral-tutorial.md 

## reproducible script

### 1. software install instructions

#### add steps from instructions linked above?

### 2. locate and prepare input files (gene trees) and decide where to store output files (species trees)

#### run these from the data folder; these commands generate 3 gene trees for use in Astral using RAxML on the .fasta files created by clustalw alignment (necessary because my previous RAxML and MrBayes work was on concatenated files, which doesn't work with Astral)
#### using standard substitution model GTRGAMMA and random seed 12345 for now

raxml-ng --msa aligned_files/filtered_uce-149.aligned_clustalw.fasta --model GTR+G --prefix astral/project_data/locus149 --seed 12345 --threads 2
raxml-ng --msa aligned_files/filtered_uce-741.aligned_clustalw.fasta --model GTR+G --prefix astral/project_data/locus741 --seed 12345 --threads 2
raxml-ng --msa aligned_files/filtered_uce-910.aligned_clustalw.fasta --model GTR+G --prefix astral/project_data/locus910 --seed 12345 --threads 2

#### combine the raxml output files (individual gene trees?) into a combined gene tree for use in astral (also run from my data folder in terminal)

cat astral/project_data/*.raxml.bestTree > astral/project_data/combined_gene_trees.tre

#### run Astral on the combined gene trees - from the data/Astral folder!

java -jar astral.5.7.8.jar -i project_data/combined_gene_trees.tre -o project_data/astral_output_species_tree.tre

### 3. Analyze the results

#### My final normalized quartet score (printed out in terminal by running astral) was 0.8920634920634921, which suggests some minor incomplete lineage sorting, but it's a typical amount for empirical datasets and indicates that the primary species signal is very clear
#### optional: visualize the output tree (data/Astral/project_data/astral_output_species_tree.tre) using R
