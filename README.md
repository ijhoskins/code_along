# Code Along 11/17/2021

UT Austin Code Along: samtools and SAM/BAM format.

## Installation

If you are on Windows, install a -nix based virtual machine.

Then, to install samtools, execute the following lines in your Terminal application:

```
wget https://github.com/samtools/samtools/releases/download/1.14/samtools-1.14.tar.bz2
tar -xvzf samtools-1.14.tar.bz2
cd samtools-1.14/
./configure --without-curses --disable-bz2 --disable-lzma
make
sudo make install
```

Test that samtools was installed, then create a sandbox (temporary directory):

```
samtools
TEMPDIR = "/tmp"
```

Clone the Code Along repository:
```
git clone https://github.com/ijhoskins/code_along
cd code_along
```

## samtools basics

Take a look at the SAM metadata stored in the header. This normally tells you how the reads were aligned and what kind of processing reads have undergone.

```
samtools view -H examples/CBS_sim.bam
```

Convert between SAM and BAM format and filter reads based on alignment criteria.

Calling samtools view on a BAM produces the SAM output. Pipe the output to head to see the first two records:

```
samtools view examples/CBS_sim.bam | head -n 2
```

Determine what the SAM flags 83 and 163 mean.

```
samtools flags
samtools flags 83
samtools flags 163
```

Generate a quick summary of the reads in the file:

```
samtools flagstat examples/CBS_sim.bam
```

Generate a more detailed summary of the aligments:

```
samtools stats examples/CBS_sim.bam > $TEMPDIR/CBS_sim_stats.txt
```

Extract all R1 sequences:

```
samtools flags READ1
samtools view -f 64 examples/CBS_sim.bam | cut -f 10
```

Determine how the R1s aligned to the reference using the CIGAR field.

The CIGAR field describes the position of insertions and deletions in the alignment. It also reports clipped segments, which can occur in local alignment mode.

```
samtools flags READ1
samtools view -f 64 examples/CBS_sim.bam | cut -f 6
```

## Visualizing alignments and writing reads from SAM/BAM format.

To visualize SAM/BAM alignments in genome browsers such as [IGV](https://software.broadinstitute.org/software/igv/), the alignments must be sorted by position and indexed:

```
samtools sort examples/CBS_sim.bam > $TEMPDIR/CBS_sim_coord_sort.bam
samtools index $TEMPDIR/CBS_sim_coord_sort.bam
ls -l $TEMPDIR/CBS_sim_coord_sort.bam*
samtools view -H $TEMPDIR/CBS_sim_coord_sort.bam
```

It is often useful to sort the reads by query/read name. For example, this allows the user to conveniently extract information from read pairs simultaneously in custom analysis programs. Sorting by read name is also required for writing reads from BAM format to FASTQ:

```
# This is the proper way to write reads from BAM to FASTQ
samtools sort -n examples/CBS_sim.bam > $TEMPDIR/CBS_sim_qname_sort.bam

samtools fastq -n -1 $TEMPDIR/CBS_sim.R1.fq -2 $TEMPDIR/CBS_sim.R2.fq $TEMPDIR/CBS_sim_qname_sort.bam 
```

## Advanced samtools

Often samtools output can be combined with shell pipelines to quickly summarize data.

One case is to count read names in a file. This is useful for determining if any reads have additional alignments (supplementary or secondary alignments) without the SAM flag.

For example, SAM flags may turn invalid if alignments are filtered by the user and not subsequently aligned to re-generate the SAM flags. In this case, "orphan" reads (e.g. a read with an an unmapped mate that was filtered) must be handled by the user.

```
samtools view examples/CBS_sim.bam | cut -f 1 | sort | uniq -c | sort -k1nr | head
```

Alignment tags describe certain properties of alignments such as edit distance (number of edit operations between the read and the reference, i.e. mismatches and InDel stretches), and the position and identity of mismatches and InDels (CIGAR field and MD tag). This information is useful in variant calling applications.

Say a developer wants to simulate variants in alignments and subsequently edits the reads by introducing mismatches or InDels. In this case, reads can written out to FASTQ and realigned to regenerate the NM and MD tags. Alternatively, samtools can be used to recalculate NM and MD alignment tags:

```
samtools calmd examples/CBS_sim.bam examples/CBS.fa > $TEMPDIR/CBS_sim_calmd.bam
samtools view -H $TEMPDIR/CBS_sim_calmd.bam
```

## Python and R interfaces to samtools

If you would rather work strictly with Python or R packages, there are samtools APIs for each language (pysam and Rsamtools).

https://pysam.readthedocs.io/en/latest/

https://bioconductor.org/packages/release/bioc/html/Rsamtools.html

## Additional resources

Read the [samtools docs](http://samtools.sourceforge.net) and the [SAM format specification](https://samtools.github.io/hts-specs/SAMv1.pdf) for more detailed information.

The use of BAM files is encouraged to limit the use of your file space. 

### End-to-end analysis workflows

For RNA-seq analysis workflows, I recommend the following guides on R Bioconductor:

http://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html

http://bioconductor.org/help/course-materials/2017/OSU/B3_RNASeq_Workflow.html

### High-throughput data analysis and statistics

For a nice statistics refresher and more advanced data analysis techniques in R:

http://genomicsclass.github.io/book/

Zelig extends common statistical inference and modeling in R:

http://docs.zeligproject.org/index.html


