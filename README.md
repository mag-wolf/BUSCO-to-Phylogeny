# BUSCO-to-Phylogeny
-------------------------------

A pipeline created to automatize all my steps for phylogenetic analyses using the BUSCO toolkit.
Does require perl5, python3 and anaconda environments!

by Magnus Wolf 2021 (magnus.wolf@senckenberg.de)
-------------------------------
M.Sc. Magnus Wolf

PhD-Student

Senckenberg Biodiversity and Climate Research Center (SBiK-F)

Group: Evolutionary Vertebrate Genomics

Georg-Voigt-Straße 14-16, floor 4, room 4.07

60325 Frankfurt am Main, Germany

Installation:
-------------------------------
1.) Download the tarball

2.) copy the tarball and rename it to what you desire

3.) extract the renamed tarball

    tar -xzvf renamed.tar.gz renamed

4.) install dependencies via conda:

    conda create --name MAFFTenv
    conda install -n MAFFTenv -c bioconda mafft 
    conda create --name CLIPKITenv
    conda install -n CLIPKITenv -c jlsteenwyk clipkit
    conda create --name IQTREEenv
    conda install -n IQTREEenv -c bioconda iqtree
    conda create --name BUSCOenv
    conda install -n BUSCOenv -c bioconda -c conda-forge busco=5.1.2  #change version number if a newer version is available.

Now you are ready to start.

Usage:
-------------------------------
1.) Gather whole genome assemblies in fasta format for all species you want to have in your tree. Make sure to include an
outgroup! Copy them in the renamed folder that you extracted previously. For convenience, I would recommend naming the 
fasta files after the scientific species name. E.g.: Balaenoptera_physalus.fasta, Balaenoptera_musculus.fasta etc.... 

2.) Open the script BUSCO-to-Phylogeny.sh with a text editor of your choice (e.g. nano).

3.) Edit general dependencies. Especially the working directory.

4.) Edit BUSCO-to-Phylogeny dependencies and parameters. Especially the OrthoDB database and the augustus 
training species you want to use.

5.) Edit options for BUSCO-to-PHYLOGENY. Call everything "TRUE" that you want to use. The pipeline 
contains 10 subparts that can be run independently if all other subparts are called "FALSE". By
leaving it as it is, everything will run one by one. 

6.) Now simply run:
    
    bash BUSCO-to-Phylogeny.sh

I suggest piping screen outputs to an error log by adding: 

    2>&1 | tee error.log
 

Here a list of all subparts:
###

    runbusco                   #run the BUSCO tool on every assembly for annotation and getting single copy orthologs.
    findsharedscos             #find BUSCO genes that are shared between all of the species you provided
    makealignments             #make alignments of all gene sequences using mafft
    trimmgenealignments        #trimm the alignments using clipkit
    constructgenetrees         #constructing genestrees from trimmed alignments using iqtree
    filtergenetrees            #filter genetrees and alignments based on the maximum likelihood distance (throw out too conserved and too variable genes)
    concatgenealignments       #concatenate gene alignments into one big matrix using FASconCAT
    trimmsupermatrix           #trimm the concatenated alignment using clipkit
    constructsupermtree        #constructing a tree from the concatenated alignment using iqtree
    constructsuperttree        #constructing a consensus tree based on all constructed genes trees using Astral 
    
###



Citation:
-------------------------------
If you use this pipeline for your research I would greatly appreciate if you cite the
study during which I created this pipeline:

Clément Schneider, Christian Woehle, Carola Greve, Cyrille A D'Haese, Magnus Wolf, Michael Hiller, Axel Janke, Miklós Bálint, Bruno Huettel, Two high-quality de novo genomes from single ethanol-preserved specimens of tiny metazoans (Collembola), GigaScience, Volume 10, Issue 5, May 2021, giab035, https://doi.org/10.1093/gigascience/giab035



Good luck!
-------------------------------
