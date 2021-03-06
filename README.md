Tuxedo.ScriptCommander
=====================

**Finding Differentially Expressed Genes with Tuxedo protocol**
You can make Tuxedo protocol scripts easily. Use this~:D
**Requirement**
--------------
- python
- **tophat**
- **bowtie**
- **cufflinks** packages 
- **samtools** (please check release date)
- **BWA** is required for calculating fragment size(Calculating tophat parameters)
- **gawk**

**Install**
--------------
```
git clone https://github.com/wjlim/TuxedoProtocolPIPE.git
export PATH=$PWD/TuxedoProtocolPIPE:$PATH (bash shell)
```

**Usage**
--------------
```
**Usage: Tuxedo.ScriptCommander.py -f refence_fasta -m mRNA_fasta -a annotation(gff,gtf) -p total_available_CPUs -i input_info1(comma sep) -i input_info2(comma sep) ...**
Options:
  -h, --help            show this help message and exit
  -f REFERENCE, --reference=REFERENCE
                        reference fasta file
  -m TRANSCRIPTOME, --mRNA=TRANSCRIPTOME
                        mRNA sequences fasta file
  -a ANNOTATION, --annotation=ANNOTATION
                        reference annotation file(gff, gtf)
  -p CPU, --num-threads=CPU
                        Total cpu count you can use[default:8]
  -o OUTFIX, --outfix=OUTFIX
                        User define outfix name [default:PIPE]
  -r, --running         If you want to run Whole script via
                        ./Tuxedo.ScriptCommander.py, add this option
  -i INPUTS, --input=INPUTS
                        Define input files and conditions with multiple -i
                        option    [-i condition1,prefix1,fastq_file1(s) -i
                        condition2,prefix2,fastq_file2(s) ...
```
- example
```
Tuxedo.ScriptCommander.py\
 -f reference/Arabidopsis_thaliana.TAIR10.23.dna.genome.fa\
 -m reference/Arabidopsis_thaliana.TAIR10.23.cds.all.fa\
 -a reference/Arabidopsis_thaliana.TAIR10.23.gtf\
 -p 4\
 -o TEST\
 -i CONTROL,COL_1,raw_data/COL_1.1.fastq.gz,raw_data/COL_1.2.fastq.gz\
 -i CONTROL,COL_2,raw_data/COL_2.1.fastq.gz,raw_data/COL_2.2.fastq.gz\
 -i CONTROL,COL_3,raw_data/COL_3.1.fastq.gz,raw_data/COL_3.2.fastq.gz\
 -i KO,BP12_1,raw_data/BP12_1.1.fastq.gz,raw_data/BP12_1.2.fastq.gz\
 -i KO,BP12_2,raw_data/BP12_2.1.fastq.gz,raw_data/BP12_2.2.fastq.gz\
 -i KO,BP12_3,raw_data/BP12_3.1.fastq.gz,raw_data/BP12_3.2.fastq.gz\
 -r
```
Information
--------------
    Reference fasta file, Transcriptome sequence file and annotation file are required.
    *bowtie2-build Reference.fa Reference* and *bwa index Transcriptome.fa* excute before running this scripts (*ln -s Reference.fa Reference* first)
    Please define OUTFIX as project name.
    If you want to run whole scripts, then add -r option.
    This script parsing with multiple -i options.
    Inputs arguments are delemited with comma.
    If you have replicates data, then put a same condition name.
    You can change detail options in "script" by edit template variables
