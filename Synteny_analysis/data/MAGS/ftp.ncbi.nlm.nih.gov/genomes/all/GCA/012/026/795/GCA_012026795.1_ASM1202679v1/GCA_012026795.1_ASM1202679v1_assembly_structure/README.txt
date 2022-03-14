################################################################################
README for assembly_structure directories under:
          ftp://ftp.ncbi.nlm.nih.gov/genomes/all/
          ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/
          ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/

Last updated: August 9, 2016
################################################################################

============================
Assembly structure directory
============================
An assembly structure directory is provided to report detailed information 
about the internal structure of the assembly. The directory contains AGP files 
that define how component sequences are organized into scaffolds and/or 
chromosomes. Other files define how scaffolds and chromosomes are organized 
into non-nuclear and other assembly-units, and how any alternate or patch 
scaffolds are placed relative to the chromosomes.

The assembly structure directory is named as:
[assembly accession.version]_[assembly name]_assembly_structure

=============
Organization:
=============
Each assembly_structure directory will contain one or more assembly-unit 
directories. Many assemblies consist of a single assembly-unit, the Primary 
Assembly; other assemblies may be comprised of multiple assembly-units. For 
example, the human GRCh38 assembly contains a Primary Assembly, a series of 
ALT_REF_LOCI assembly-units and a non-nuclear (organelle) assembly-unit.

     Primary_Assembly/
     ALT_REF_LOCI_*/
     non-nuclear/

Each assembly-unit directory contains the following files:
     component_localID2acc
     scaffold_localID2acc
     join_certificate.xml (optional - only present for some assemblies from the 
                           Genome Reference Consortium)

Each assembly-unit directory also contains one or more of the following 
directories (depending on the particular assembly):
     assembled_chromosomes/
     placed_scaffolds/
     unlocalized_scaffolds/
     unplaced_scaffolds/
     alt_scaffolds/ (only in alternate loci and patch assembly-units)
     pseudoautosomal_region/ (only for mammals)

The content of the assembled_chromosomes, placed_scaffolds, 
unlocalized_scaffolds, unplaced_scaffolds, alt_scaffolds and 
pseudoautosomal_region directories is:

assembled_chromosomes/
     chr2acc
     FASTA/
          chr*.fna.gz
     AGP/
          chr*.comp.agp.gz
          chr*.agp.gz

placed_scaffolds/
     FASTA/
          chr*.placed.scaf.fna.gz
     AGP/
          chr*.placed.scaf.agp.gz

unlocalized_scaffolds/
     unlocalized.chr2scaf
     FASTA/
          chr*.unlocalized.scaf.fna.gz
     AGP/
          chr*.unlocalized.scaf.agp.gz

unplaced_scaffolds/
     FASTA/
          unplaced.scaf.fna.gz
     AGP/
          unplaced.scaf.agp.gz

alt_scaffolds/
     FASTA/
          alt.scaf.fna.gz
     AGP/
          alt.scaf.agp.gz
     alt_scaffold_placement.txt
     alignments/
          {scaffold accession.version}_{chromosome accession.version}.asn
          {scaffold accession.version}_{chromosome accession.version}.gff

pseudoautosomal_region/
     par_align.asn
     par_align.gff
     
assembly_structure directories for assemblies with alternate loci or patch 
scaffolds may also contain some additional files:

genomic_regions_definitions.txt 
     This file is provided for assemblies that have defined REGIONS on the 
     Primary Assembly for which alternative loci or patch scaffolds have been 
     provided. The file reports:
     region name | chromosome accession.version | start | stop

all_alt_scaffold_placement.txt
     A file providing the alternate locus scaffold placements for all the 
     alternate assembly-units. See below for the file format.

------
Notes:
------
1. The sequences of the placed scaffolds are redundant with the
sequences of the assembled chromosomes. The placed scaffolds are 
provided for users who prefer to work with scaffolds rather than with
chromosomes.
2. Eukaryote genome assemblies may include an assembly-unit named
"non-nuclear" which contains data from organelle genomes, for example
the mitochondrion or chloroplast.
3. If the assembly is comprised of more than one assembly-unit, the
names for the assembly-units, other than a "non-nuclear"
assembly-unit, are supplied by the submitter.
4. The chromosome-from-scaffold AGP file (chr?.agp.gz), and the
placed_scaffolds directory, may be omitted if the chromosome is
assembled directly from components, or if the chromosome is a complete
sequence with no gaps.
5. The file suffix .agp.gz indicates AGP files. See format specification: 
https://www.ncbi.nlm.nih.gov/genome/assembly/agp/AGP_Specification.shtml

=====================
Description of files:
=====================
1. Files containing genomic sequences in nucleotide fasta format (unmasked)
FILENAME                         CONTENT
chr?.fna.gz                      chromosome sequence
chr?.placed.scaf.fna.gz          placed scaffold sequences
chr?.unlocalized.scaf.fna.gz     unlocalized scaffold sequences
unplaced.scaf.fna.gz             unplaced scaffold sequences
alt.scaf.fna.gz                  alternate loci or patch scaffold sequences

2. AGP files
See format specification: 
https://www.ncbi.nlm.nih.gov/genome/assembly/agp/AGP_Specification.shtml
The AGP files in this directory tree use GenBank or RefSeq accession.versions 
as the identifiers for components, scaffolds, and chromosomes.

FILENAME                         CONTENT
chr?.comp.agp.gz                 chromosome from component AGP
chr?.agp.gz                      chromosome from scaffold AGP
chr?.placed.scaf.agp.gz          placed scaffold from component AGP
chr?.unlocalized.scaf.agp.gz     unlocalized scaffold from component AGP
unplaced.scaf.agp.gz             unplaced scaffold from component AGP
alt.scaf.agp.gz                  alternate loci or patch scaffold from
                                 component AGP

3. component_localID2acc
A two column file associating the submitter component ID with the 
accession.version. 'na' is shown in the ID column if the submitter did not 
provide a name for the component.

4. scaffold_localID2acc
A two column file associating the submitter scaffold ID with the 
accession.version. 'na' is shown in the ID column if the submitter did not 
provide a name for the scaffold.

5. chr2acc
A two column file associating the chromosome, or linkage group name, with the 
accession.version.

6. unlocalized.chr2scaf
A two column file giving the chromosome or linkage group assignment for each 
unlocalized scaffold.

7. join_certificate.xml
This file provides data on joins in the assembly that were curated by the 
Genome Reference Consortium (GRC). This file will not be present for assemblies 
submitted by other groups.

8. genomic_regions_definitions.txt
A file defining the regions on the primary assembly for which alternate loci or 
patch scaffolds are available. May also include pseudoautosomal regions, 
centromere regions or heterochromatin regions if these have been defined for an 
assembly.

The file is tab delimited (including a #header) with the following columns:
col 1: region_name: name for the genomic region
col 2: chromosome: accession.version for the chromosome or 
       unlocalized/unplaced scaffold
col 3: start: the starting position on the chromosome or scaffold
       (in 1 base coordinates)
col 4: stop: the ending position on the chromosome or scaffold
       (in 1 base coordinates)

9. alt_scaffold_placement.txt & all_alt_scaffold_placement.txt
Files associating alternate loci or patch scaffolds with the corresponding 
primary assembly chromosome, providing the location on the chromosome, the 
genomic region name, and the length of any unaligned tails.

The file is tab delimited (including a #header) with the following columns:
col 1: alt_asm_name: name of the assembly-unit that includes the alternate 
       scaffold
col 2: prim_asm_name: name of the primary assembly-unit on which the alternate 
       scaffold is being placed
col 3: alt_scaf_name: name of the alternate scaffold being placed
col 4: alt_scaf_acc: accession.version of the alternate scaffold being placed
col 5: parent_type: type of object on which the alternate scaffold is being 
       placed, either CHROMOSOME or SCAFFOLD
col 6: parent_name: name of the object on which the alternate scaffold is being 
       placed (can be either a chromosome or a scaffold)
col 7: parent_acc: accession.version of the sequence on which the alternate 
       scaffold is being aligned
col 8: region_name: name of the genomic region on the parent within which the 
       alternate scaffold is placed
col 9: ori: orientation of the alignment, '+', '-' or 'b' (mixed)
col10: alt_scaf_start: start of the placement on the alternate scaffold
       (in 1 base coordinates)
col11: alt_scaf_stop: end of the placement on the alternate scaffold
       (in 1 base coordinates)
col12: parent_start: start of the placement on the parent sequence
       (in 1 base coordinates)
col13: parent_stop: end of the placement on the parent sequence 
       (in 1 base coordinates)
col14: alt_start_tail: number of bases at the start of the alternate scaffold 
       not involved in the placement
col15: alt_stop_tail: number of bases at the end of the alternate scaffold not 
       involved in the placement

Note: Every alternate scaffold associated with the assembly-unit will be listed 
in this file. Any alternate scaffold that has no placement will have 'na' in 
columns 5 to 15. Any alternate scaffold that has a chromosome assignment, but 
no alignment, would have the chromosome name in column 6 and 'na' in columns 7 
to 15.

10. alignments/{scaffold accession.version}_{chromosome accession.version}.asn
Files providing alignments of the alternate loci or patch scaffolds to the 
corresponding primary assembly chromosome, in ASN.1 format. These alignments 
indicate how the alternate loci and patch scaffold sequences differ from the 
chromosomes of the primary assembly.

11. alignments/{{scaffold accession.version}_{chromosome accession.version}.gff
Files providing alignments of the alternate loci or patch scaffolds to the 
corresponding primary assembly chromosome, in CIGAR format embedded within a 
GFF format file. These alignments indicate how the alternate loci and patch 
scaffold sequences differ from the chromosomes of the primary assembly.

12. par_align.asn
A file providing alignments between each pseudoautosomal region (PAR) on the X 
chromosome and the corresponding PAR on the Y chromosome, in ASN.1 format. 

13. par_align.gff
A file providing alignments between each pseudoautosomal region (PAR) on the X 
chromosome and the corresponding PAR on the Y chromosome, in CIGAR format 
embedded within a GFF format file.

14. patch_type
A file providing the patch type for each of the scaffolds in a patch assembly-
unit.

The file is tab delimited (including a #header) with the following columns:
col 1: alt_scaf_name: local name for the patch scaffold
col 2: alt_scaf_acc: the accession.version for the patch scaffold
col 3: patch_type: FIX or NOVEL (defined below)

============
Definitions:
============
Assembly:
A set of chromosome assemblies, unlocalized and unplaced sequences and 
alternate loci used to represent an organisms genome. Most current assemblies 
are a haploid representation of an organisms genome, although some loci may be 
represented more than once (see Alternate locus, below). This representation 
may be obtained from a single individual (e.g. chimp or mouse) or multiple 
individuals (e.g. human Genome Reference Consortium assembly). Except in the 
case of organisms that have been bred to homozygosity, the haploid assembly 
does not typically represent a single haplotype, but rather a mixture of 
haplotypes.

Chromosome Assembly: 
A relatively complete pseudo-molecule assembled from smaller sequences
(components) that represent a biological chromosome. Relatively complete 
implies that some gaps may still be present in the assembly, but independent 
measures suggest that most of the sequence is represented by sequenced bases. 
Completeness is submitter defined.

Unlocalized sequence:
A sequence found in an assembly that is associated with a specific chromosome 
but cannot be ordered or oriented on that chromosome. 

Unplaced sequence:
A sequence found in an assembly that is not associated with any chromosome.  

Primary assembly:
An assembly-unit representing the collection of assembled chromosomes,
unlocalized and unplaced sequences that, when combined, should represent a 
non-redundant haploid genome. This excludes any alternate loci.

Alternate locus:
A sequence that provides an alternate representation of a locus found in the 
primary assembly. These sequences do not represent a complete chromosome 
sequence although there is no hard limit on the size of the alternate locus; 
currently these are less than 1 Mb.

Alternate locus group:
An assembly-unit consisting of scaffolds from different loci that are 
considered to be part of the same haplotype (e.g. mouse 129/Sv group).

Genomic region:
A defined span on the primary assembly for which alternate loci or patch 
scaffolds are available. Genomic regions may be named after a gene or gene 
cluster, or may be given arbitrary region numbers.

Major release:
The formal release of a genome assembly, e.g. GRCh38.

Minor release:
A release of a genome assembly including patches that occurs between
major releases.

Genome Patch:
A sequence contig/scaffold that corrects sequence in a major release of the 
genome, or adds sequence to it.

FIX patch:
A patch that corrects sequence or reduces an assembly gap in a given major 
release. FIX patch sequences are meant to be incorporated into the primary or 
existing alt-loci assembly units at the next major release, and their 
accessions will then be deprecated.

NOVEL patch:
A patch that adds sequence to a major release. Typically, NOVEL patch sequences 
are meant to be incorporated into the assembly as new alternate loci at the 
next major release, and their accessions will not be deprecated.

________________________________________________________________________________
National Center for Biotechnology Information (NCBI)
National Library of Medicine
National Institutes of Health
8600 Rockville Pike
Bethesda, MD 20894, USA
tel: (301) 496-2475
fax: (301) 480-9241
e-mail: info@ncbi.nlm.nih.gov
________________________________________________________________________________

