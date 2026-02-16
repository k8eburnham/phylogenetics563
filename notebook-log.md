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

## Note to self: these files still need to be combined into one master file for each locus before I build my phylogenetic tree (concatenation)