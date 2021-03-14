#!/usr/bin/env python3
#######################################################################
# File Name: rna-seq.py
# Author: Zhao Chen
# mail: WHULily1996@whu.edu.cn
# Usage: 
# Description: 
# History: 
#######################################################################
import os,sys
# import pandas as pd
# import numpy as np

samtoolsDir = "/data/zhaochen/project/zhangbo/rna-seq/samtools"
cufflinksDir = "/data/zhaochen/project/zhangbo/rna-seq/cufflinks"

def cufflinks(sample,sortbam):
    cmd = "cufflinks -p 21 -o %s/%s %s" % (cufflinksDir,sample,sortbam)
    os.system(cmd)
    return cmd

def assemblies():
    cmd ="ls -R %s/*/transcripts.gtf > %s/assemblies.txt" % (cufflinksDir,cufflinksDir)
    os.system(cmd)
    return cmd

def cuffmerge():
    gtf = "/data/zhaochen/reference/mm10/Annoation/mm10_annotation.gtf"
    fasta = "/data/zhaochen/reference/mm10/Sequence/GRCm38.p6.genome.fa"
    mergetxt = os.path.join(cufflinksDir,"assemblies.txt")
    cmd = "cuffmerge -p 21 -o %s -g %s -s %s %s" % (cufflinksDir,gtf,fasta,mergetxt)
    os.system(cmd)
    return cmd

def cuffdiff(sample_list,bam_list):
    fasta = "/data/zhaochen/reference/mm10/Sequence/GRCm38.p6.genome.fa"
    cmd = "cuffdiff -p 21 -o %s/diff_out -b %s -L %s -u %s/merged.gtf %s" % (cufflinksDir,fasta,sample_list,cufflinksDir,bam_list)
    os.system(cmd)
    return cmd
    # if Error: number of labels must match number of conditions, labels split with ","(comma)

if __name__=="__main__":
    step = 2

    if step < 1:
        for root, dirs, files in os.walk(samtoolsDir):
            for file in files:
                if os.path.splitext(file)[1] == ".bam":
                    sample = file.split("_")[0]
                    file = os.path.join(root, file)
                    print(file)
                    #cufflinks(sample=sample,sortbam=file)


    if step < 2:
        #assemblies()
        cuffmerge()

    if step < 3:
        bam_list = []
        sample_list = []
        for root, dirs, files in os.walk(samtoolsDir):
            for file in files:
                if os.path.splitext(file)[1] == ".bam":
                    sample = file.split("_")[0]
                    sample_list.append(sample)
                    file = os.path.join(root, file)
                    bam_list.append(file)

        bam_list = sorted(bam_list)
        bam_list = " ".join(bam_list)
        sample_list = sorted(sample_list)
        sample_list = ",".join(sample_list)
        print(sample_list)
        print(bam_list)
        cuffdiff(sample_list=sample_list,bam_list=bam_list)


