#!/usr/bin/env python
#Won-Jun Lim
import sys, os, re, time, subprocess, glob
from optparse import OptionParser
from itertools import *

def myFileIndexingCheck(REF,mRNA):
    if not os.access(mRNA+".bwt", os.F_OK):
        print "Please Indexing Reference file ( bwa index %s )" % mRNA
        sys.exit()

    if not os.access(REF.replace(".fa","")+".rev.2.bt2", os.F_OK):
        print "Please Indexing (bowtie2-build Reference.fa Reference) and make symbolic links (ln -s Reference.fa Reference)" % REF
        sys.exit()

def myMkdir(PATH):
    if not os.access(PATH, os.W_OK):
        os.mkdir(PATH,0751)

def myTime():
    return time.strftime('%Y-%m-%d\t%p %H:%M:%S')

def myHumanReadableSort(LIST):
    conv = lambda x:int(x) if x.isdigit() else x
    KEY = lambda x:[conv(c) for c in re.split('([0-9]+)',x)]
    return sorted(LIST,key=KEY)

def myInputCONDITIONHandler(INPUT_list):
    INPUT_map = {}
    CONDITION_names = set()
    for INPUT in imap(lambda x:x.split(','),INPUT_list):
        # {CONDITION:{PREFIX:(fastq1, fastq2)},...}
        INPUT_map.setdefault(INPUT[0],{}).setdefault(INPUT[1],INPUT[2:])
        CONDITION_names.add(INPUT[0])
    return INPUT_map, myHumanReadableSort(list(CONDITION_names))

def myDistributableCPUcountCalc(CPUs, INPUT_counts):
    if CPUs/INPUT_counts < 4:
        return (False, CPUs)
    else:
        return (True, CPUs/INPUT_counts)

def myJobExecuteAndSuspending(SCRIPTS,Distributable=True):
    def myJob(script):
        print "Script: %s is running start at %s" % (script, myTime())
        process = subprocess.Popen('bash %s' % script, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        stdout, stderr = process.communicate()
        logFile = open(script+'.log','wa')
        errFile = open(script+'.err','wa')
        logFile.write(stdout)
        errFile.write(stderr)
        logFile.close()
        errFile.close()
        return process

    if Distributable:
        processIDs = []
        for script in SCRIPTS:
            process = myJob(script)
            processIDs.append(process)

        for process in processIDs:
            process.wait()

    elif type(SCRIPTS) == str:
        process = myJob(SCRIPTS)
        process.wait()

    else:
        for script in SCRIPTS:
            process = myJob(script)
            process.wait()

def mainScriptMaker(REF,mRNA,ANNOT,CPUs,OUTFIX,OUTPATH,INPUT_list,INPUT_map,CONDITIONS,CPU,Distributable):
    bwaMappingScriptName = '{OUTPATH}/scripts/run.{PREFIX}.bwa.mRNA.sh'
    tophatMappingScriptName = '{OUTPATH}/scripts/run.{PREFIX}.tophat.sh'
    cufflinksScriptName = '{OUTPATH}/scripts/run.{PREFIX}.cufflinks.sh'
    cuffmergeScriptName = '{OUTPATH}/scripts/run.{OUTFIX}.cuffmerge.sh'
    cuffdiffScriptName = '{OUTPATH}/scripts/run.{OUTFIX}.cuffdiff.sh'

    bwaMappingScript = '''#!/bin/bash
mkdir -p {OUTPATH}/mRNA_mapping

bwa mem -t {CPU} {mRNA} {FASTQ} | samtools view -Sb - > {OUTPATH}/mRNA_mapping/{PREFIX}.transcriptome.bwa.bam

tlen=`samtools view {OUTPATH}/mRNA_mapping/{PREFIX}.transcriptome.bwa.bam|awk '$3 !~ /\*/ && $5 == 60 && $7 == "=" && $9 > 0' | head -n 50000 | awk  '{{a+=$9;b+=$9^2;c+=length($10)}} END{{read=c/NR;printf "%.f^%.f",a/NR-(read*2),sqrt(b/NR - (a/NR)^2)}}' 2>/dev/null`
avg=`echo $tlen|cut -d ^ -f 1`
sd=`echo $tlen|cut -d ^ -f 2`

echo -e "sample\\tAverage Insert\\tStandard Dev\\n{PREFIX}\\t$avg\\t$sd" > {OUTPATH}/logs/{PREFIX}.frag.stat
'''
    tophatMappingScript = '''#!/bin/bash
mkdir -p {OUTPATH}/Tophat

tlen=`bwa mem -t {CPU} {mRNA} {FASTQ}|samtools view -S -|awk '$3 !~ /\*/ && $5 == 60 && $7 == "=" && $9 > 0' | head -n 50000 | awk  '{{a+=$9;b+=$9^2;c+=length($10)}} END{{read=c/NR;printf "%.f^%.f",a/NR-(read*2),sqrt(b/NR - (a/NR)^2)}}' 2> /dev/null`
avg=`echo $tlen|cut -d ^ -f 1`
sd=`echo $tlen|cut -d ^ -f 2`

tophat2\
 -o {OUTPATH}/Tophat/{PREFIX}.tophat\
 --num-thread={CPU}\
 -G {ANNOT}\
 --mate-inner-dist=$avg\
 --mate-std-dev=$sd\
 --rg-id={PREFIX}\
 --rg-sample={PREFIX}\
 --rg-library={PREFIX}.lib\
 --library-type=fr-unstranded\
 {REF} {FASTQ}
'''
    cufflinksScript = '''#!/bin/bash
mkdir -p {OUTPATH}/cufflinks

cufflinks\
 -o {OUTPATH}/cufflinks/{PREFIX}.cufflinks\
 --num-threads={CPU}\
 --GTF-guide={ANNOT}\
 --frag-bias-correct {REF}\
 --multi-read-correct\
 --label={PREFIX}\
 --library-type=fr-unstranded\
 {tophatOutput}
'''
    cuffmergeScript = '''#!/bin/bash
mkdir -p {OUTPATH}/cuffmerge

cuffmerge\
 -o {OUTPATH}/cuffmerge\
 --ref-gtf={ANNOT}\
 --ref-sequence={REF}\
 --num-threads={CPUs}\
 {CUFFLINKS_OUTFILEs_list}
'''
    cuffdiffScript = '''#!/bin/bash
mkdir -p {OUTPATH}/cuffdiff

cuffdiff\
 -o {OUTPATH}/cuffdiff\
 -u\
 --frag-bias-correct={REF}\
 --num-threads={CPU}\
 -labels={LABELs}\
 {OUTPATH}/cuffmerge/merged.gtf\
 {inputArgsPackingWithComma}
'''

    def myScriptMaker(template, scriptName ,arguments_map):
        def myScriptWriter(template,scriptName):
            script = open(scriptName,'w')
            script.write(template)
            script.close()

        script = template.format(**arguments_map)
        scriptName = scriptName.format(**arguments_map)

        myScriptWriter(script, scriptName)
        return scriptName

# scripts lists to return execution function
    BWAscripts = []
    TOPHATscripts = []
    CUFFLINKSscripts = []
    CUFFMERGEscripts = []
    CUFFDIFFscripts = []
    CUFFLINKSoutputs = []
    PREFIX_list = []

# Local arguments formatting with templates
    for CONDITION in CONDITIONS:
        for PREFIX in myHumanReadableSort(INPUT_map[CONDITION].keys()):
            FASTQ = ' '.join(map(lambda x:PWD+x,INPUT_map[CONDITION][PREFIX]))
            tophatOutput = '%s/Tophat/%s.tophat/accepted_hits.bam' % (OUTPATH, PREFIX)

            bwaMappingScriptNameFormatted = myScriptMaker(bwaMappingScript,bwaMappingScriptName,locals())
            tophatMappingScriptNameFormatted = myScriptMaker(tophatMappingScript,tophatMappingScriptName,locals())
            cufflinksScriptNameFormatted = myScriptMaker(cufflinksScript, cufflinksScriptName,locals())

            BWAscripts.append(bwaMappingScriptNameFormatted)
            TOPHATscripts.append(tophatMappingScriptNameFormatted)
            CUFFLINKSscripts.append(cufflinksScriptNameFormatted)

            CUFFLINKSoutputs.append('%s/cufflinks/%s.cufflinks/transcripts.gtf'%(OUTPATH,PREFIX))
            PREFIX_list.append(PREFIX)

    CUFFLINKS_OUTFILEs_list = '%s/logs/cufflinks_out.ls' % OUTPATH
    ofile = open(CUFFLINKS_OUTFILEs_list,'w')
    ofile.write('\n'.join(CUFFLINKSoutputs))
    ofile.close()

    cuffmergeScriptNameFormatted = myScriptMaker(cuffmergeScript,cuffmergeScriptName,locals())

    LABELs = ','.join(CONDITIONS)
    inputArgsPacking_list = []
    for CONDITION in CONDITIONS:
        PREFIXs = map(lambda x:'%s/Tophat/%s.tophat/accepted_hits.bam' % (OUTPATH,x) ,myHumanReadableSort(INPUT_map[CONDITION].keys()))
        inputArgsPacking_list.append(','.join(PREFIXs))

    inputArgsPackingWithComma = ' '.join(inputArgsPacking_list)
    cuffdiffScriptNameFormatted = myScriptMaker(cuffdiffScript,cuffdiffScriptName,locals())
    return BWAscripts,TOPHATscripts,CUFFLINKSscripts,cuffmergeScriptNameFormatted, cuffdiffScriptNameFormatted

def myRun(BWAscripts, TOPHATscripts, CUFFLINKSscripts, CUFFMERGEscripts, CUFFDIFFscripts, Distributable):
    startTime = time.time()
    myJobExecuteAndSuspending(BWAscripts,Distributable)
    myJobExecuteAndSuspending(TOPHATscripts,Distributable)
    myJobExecuteAndSuspending(CUFFLINKSscripts,Distributable)
    myJobExecuteAndSuspending(CUFFMERGEscripts,False)
    myJobExecuteAndSuspending(CUFFDIFFscripts,False)
    print "Total Process Time : %.2gs" % (endTime-startTime)

if __name__=='__main__':

    parser = OptionParser("Usage: %prog\
    -f refence_fasta -m mRNA_fasta -a Annotation(gff,gtf)\
    -p total_available_CPUs -i input_info1(comma sep) -i input_info2(comma sep) ...")
    parser.add_option("-f","--reference",dest="reference",help="reference fasta file")
    parser.add_option("-m","--mRNA",dest="transcriptome",help="mRNA sequences fasta file")
    parser.add_option("-a","--annotation",dest="annotation",help="reference annotation file(gff, gtf)")
    parser.add_option("-p","--num-threads",type=int,dest="CPU",default=8,help="Total cpu count you can use[default:8]")
    parser.add_option("-o","--outfix",default="PIPE",dest="OUTFIX",help="User define outfix name [default:PIPE]")
    parser.add_option("-r","--running",action="store_true",default=False,dest="running",help="If you want to run Whole script via %s, add this option" % sys.argv[0])
    # parser.add_option("-s","--stat",action="store_true",default=False,dest="statistics",help="If you want to know simple stat of outputs, add this option")
    parser.add_option("-i","--input",action="append",dest="inputs",help="Define input files and conditions with multiple -i option\
    [-i condition1,prefix1,fastq_file1(s) -i condition2,prefix2,fastq_file2(s) ...")
    (opts, args) = parser.parse_args()

    if len(sys.argv) < 5:
        print "%s -h or --help for help massage" % sys.argv[0]
        sys.exit()

# Arguments
    PWD=os.getcwd()+'/'
    REF = PWD + opts.reference.replace(".fa","")
    mRNA = PWD + opts.transcriptome
    ANNOT = PWD + opts.annotation
    CPUs = opts.CPU
    OUTFIX = opts.OUTFIX
    INPUT_list = opts.inputs
    RUN = opts.running
    OUTPATH = PWD + OUTFIX + '.TuxedoProtocol'

    INPUT_map, CONDITIONS = myInputCONDITIONHandler(INPUT_list)

# Reference indexing check
    myFileIndexingCheck(REF,mRNA)

# Make outpath
    myMkdir(OUTPATH)
    myMkdir(OUTPATH+'/scripts')
    myMkdir(OUTPATH+'/logs')

# CPU capacity test
    Distributable, CPU = myDistributableCPUcountCalc(CPUs, len(INPUT_list))

# Make Running Scripts
    BWAscripts,TOPHATscripts,CUFFLINKSscripts,CUFFMERGEscripts,CUFFDIFFscripts = \
    mainScriptMaker(REF,mRNA,ANNOT,CPUs,OUTFIX,OUTPATH,INPUT_list,INPUT_map,CONDITIONS,CPU,Distributable)

#Job Execution
    if RUN: myRun(BWAscripts, TOPHATscripts, CUFFLINKSscripts, CUFFMERGEscripts, CUFFDIFFscripts, Distributable)
