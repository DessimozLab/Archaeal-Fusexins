Fsx1 project Repo
-----------------
The python scripts contained in this repo were used to generate the KMER spectra, structural comparison and HMM comparison results. There is also code used in phylogenetic tree annotation and for generating the annotation of proteomes using HMMER, pfam and ArCOG databases. Each notebook contains a first explanatory markdown cell and comments in the code to help users replicate the analysis.

The Datasource folder also contains files relevant to the manuscript such as movies, tables etc.


HHblits
-------
Ipython notebook used with HHblits. Used for finidng Fsx1 and comparing Fsx1 profiles to known fusexin profiles. ( Figure 1.a , supp Figure 1 ).

Kmer
----
Used for finding IMEs using kmer spectra and extracting transfered regions. ( Figure 4.a , Supp Figure 13 )

Mobilome_ORFs
-------------
Script and notebook to accompany HMMER. These can be used to annotate a proteome using PFAM and arcog domains. ( Source data table 2 )

Structures
---------
Contains the set of monomeric structures used for structure-based phylogenetics figures ( Figure 4.c ).


Protein-comparison-tool
-----------------------
This script runs Fatcat and TMalign as well as FastmME and Clustalo to create trees and sequence alignments based on structural comparison. ( Figure 4.c ).


Hmmer_vs_metaclust
------------------
This notebook was used to parse the results of the phmmer search against the metaclust database. The HMM used to search the DB is also included in this folder. ( Source data table 1 )

Synteny_analysis
----------------
Data, scripts, and notebooks to generate IMEs clustering based on Jaccard indexes (Supplementary Fig. 8) and  synteny plots (Fig. 8; Supplementary Fig. 16 ).

Phylo_Eco_Mapping
-----------------
Data, scripts and structural models to produce the phylogenetic tree with environmental information and selected trimers with electrostatic surfaces to be calculated with APBS (Fig. 6). Also containing sequences, alingment and phylogenetic tree of Fsx1 and HAP2/GCS1 homologs (Suppelemtary Figure 10).

Datasource
-----------
Contains supplementary files produced by bioinformatics analysis, primers, synthesized sequences, movies, sequence identifiers etc.

A note on software and system requirements
----------------------
This software and the scripts used for this project have only been tested on Ubuntu 18.04 and 20.04 environments
At various steps of the project compelementary bioinformatics strategies were used to generate models or search for homology. Furthermore, several external tools are called by these scripts and notebooks. Installation instuctions and dependencies for the software used in this project can be found at the following locations:

HH-Suite:[https://github.com/soedinglab/hh-suite] ( v3.3.0 )
Emboss:[http://emboss.sourceforge.net/download/] ( v6.5.7 )
Blast suite:[https://blast.ncbi.nlm.nih.gov/Blast.cgi] ( v2.12.0 )
Hmmer:[http://hmmer.org/] ( v3.3.2 ) 
OpenMM suite and pdbfixer:[https://openmm.org/] ( v7.6.0 and 1.8.1 respectively )
ClustalO:[http://www.clustal.org/omega/] ( v1.2.2 )
iqtree:[http://www.iqtree.org/] ( v2 )
MAD root:[https://github.com/davidjamesbryant/MADroot] ( v1.0 )
AlphaFold:[https://github.com/deepmind/alphafold] ( v2.1.1 )
TMalign:[https://zhanggroup.org/TM-align/] (Version 20210224)
FATCAT:[https://fatcat.godziklab.org/] ( v2.0 )
I-tasser:[https://zhanggroup.org/I-TASSER/] ( v5.1 )
Modeller:[https://salilab.org/modeller/] ( v10.2 )
TOPCONS:[https://topcons.cbr.su.se/] ( v2.0 )
McScan:[https://github.com/tanghaibao/mcscan] ( v0.8 )
FastME:[http://www.atgc-montpellier.fr/fastme/] ( v2.0 )

The metaclust and uniclust30 databases used with Hmmer and HHblits can be found on the mmseqs web page at [https://metaclust.mmseqs.org/] and [https://uniclust.mmseqs.com/] respectively. The Pfam database is available at Pfam:[http://pfam.xfam.org/]( v33 ). 





