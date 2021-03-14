#!/usr/bin/env python3
#######################################################################
# File Name: chip-seq.py
# Author: Zhao Chen
# mail: WHULily1996@whu.edu.cn
# Usage: 
# Created Time: 2021年01月18日 星期一 22时13分34秒
# Description: 
# History: 
#######################################################################
import os,sys
import pandas as pd
import numpy as np

rawDir = "/nfs-data2/zhaochen/zhangbo/chip-seq"
cleanDir = "/data/zhaochen/project/zhangbo/chip-seq/clean"
cutadaptDir = "/data/zhaochen/project/zhangbo/chip-seq/cutadapt"
alignDir = "/data/zhaochen/project/zhangbo/chip-seq/align"
samtoolsDir = "/data/zhaochen/project/zhangbo/chip-seq/samtools"
macs2Dir = "/data/zhaochen/project/zhangbo/chip-seq/macs2"
macs2broadDir = "/data/zhaochen/project/zhangbo/chip-seq/macs2broad"
bedDir = "/data/zhaochen/project/zhangbo/chip-seq/samtools/bed"
bwDir = "/data/zhaochen/project/zhangbo/chip-seq/bigwig"

def fastp(sample_R1,sample_R2):
    name_R1 = os.path.basename(sample_R1)
    name_R2 = os.path.basename(sample_R2)
    output_R1 = os.path.join(cleanDir, name_R1.replace("_1.fq.gz", "_R1.fq.gz"))
    output_R2 = os.path.join(cleanDir, name_R2.replace("_2.fq.gz", "_R2.fq.gz"))
    cmd = "fastp -w 21 -i %s -o %s -I %s -O %s" % (sample_R1,output_R1,sample_R2,output_R2)
    os.system(cmd)
    return cmd

def cutadapt_PE(sample_R1, sample_R2):
    """
    PE mains pair end
    SE mains singal end
    """
    adapter_1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
    adapter_2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
    name_R1 = os.path.basename(sample_R1)
    name_R2 = os.path.basename(sample_R2)
    output_R1 = os.path.join(cutadaptDir, name_R1.replace("_R1.fastq.gz", "_R1.fq.gz"))
    output_R2 = os.path.join(cutadaptDir, name_R2.replace("_R2.fastq.gz", "_R2.fq.gz"))
    cmd = "cutadapt -j 28 -a %s -A %s -u 5 -u -35 -U 5 -U -35 -m 30 -o %s -p %s %s %s" % \
          (adapter_1, adapter_2, output_R1, output_R2, sample_R1, sample_R2)
    os.system(cmd)
    return cmd

def bowtie2(sample,clean_R1,clean_R2):
	index = "/data/zhaochen/reference/mm10/bowtie2Index/mm10"
	sample = sample.split("_")[0]
	prefix = os.path.join(alignDir,sample)
	bowtie2_summary = os.path.join("/data/zhaochen/project/zhangbo/chip-seq/bowtie2/bowtie2_summary",sample + "_bowtie2.txt")
	cmd = "bowtie2 -p 21 --end-to-end --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 -x %s -1 %s -2 %s " \
		  "-S %s.sam > %s" % (index,clean_R1,clean_R2,prefix,bowtie2_summary)
	#print(cmd)
	os.system(cmd)
	return cmd

def samtools(sample,samfile):
	view_bam = os.path.join(samtoolsDir,sample + "_view.bam")
	cmd_view = "/usr/local/bin/samtools view -@ 8 -bS -o %s %s" % (view_bam,samfile)
	os.system(cmd_view)
	sort_bam = view_bam.replace("_view.bam","_sort.bam")
	cmd_sort = "/usr/local/bin/samtools sort -@ 8 -O BAM -o %s %s" % (sort_bam,view_bam)
	os.system(cmd_sort)
	rmdup_bam = view_bam.replace("_view.bam","_rmdup.bam")
	cmd_rmdup = "/usr/local/bin/samtools rmdup %s %s" % (sort_bam,rmdup_bam)
	os.system(cmd_rmdup)
	cmd_index = "/usr/local/bin/samtools index -@ 8 -b %s" % (rmdup_bam)
	os.system(cmd_index)

def bamCoverage(rmdup_bam,sample):
	bigwig = os.path.join(bwDir,sample + ".bw")
	cmd = "bamCoverage -b %s -bs 1 -p 30 --effectiveGenomeSize 2652783500 --extendReads --normalizeUsing RPKM -o %s" % (rmdup_bam,bigwig)
	os.system(cmd)
	return cmd

def flagstat(sample,rmdup_bam):
	flagstat = os.path.join(samtoolsDir,sample + "_flagstat.txt")
	cmd = "/usr/local/bin/samtools flagstat -@ 7 %s > %s" % (rmdup_bam,flagstat)
	os.system(cmd)
	return cmd

def macs2(inputfile,outname):
	cmd = "macs2 callpeak --outdir %s -t %s -g mm -p 1e-5 --broad --broad-cutoff 1e-5 --keep-dup all -f BAMPE -n %s" % \
		  (macs2broadDir,inputfile,outname)
	os.system(cmd)
	return cmd

def bedgraph(bamfile,sample):
	bedfile = os.path.join(bedDir,sample + ".bed")
	cmd1 = "bedtools bamtobed -bedpe -i %s > %s" % (bamfile,bedfile)
	cleanbed = os.path.join(bedDir,sample + ".clean.bed")
	cmd2 = "awk '$1==$4 && $6-$2 < 1000 {print $0}' %s > %s" % (bedfile,cleanbed)
	fragmentbed = os.path.join(bedDir,sample + ".fragments.bed")
	cmd3 = "cut -f 1,2,6 %s | sort -k1,1 -k2,2n -k3,3n > %s" % (cleanbed,fragmentbed)
	fragmentbg = os.path.join(bedDir,sample + ".fragments.bedgraph")
	cmd4 = "bedtools genomecov -bg -i %s -g /data/zhaochen/reference/mm10/mm10.chrom.sizes > %s" % (fragmentbed,fragmentbg)
	#os.system(cmd1)
	os.system(cmd2)
	os.system(cmd3)
	os.system(cmd4)

def seacr():
	cmd = "bash /data/zhaochen/software/SEACR-master/SEACR_1.3.sh "
	os.system(cmd)
	return cmd

def annotationPeaks(broadpeak,sample):
	annotation = os.path.join(macs2broadDir,sample + "_peak_annotation.txt")
	cmd = "annotatePeaks.pl %s mm10 > %s" % (broadpeak,annotation)
	os.system(cmd)
	return cmd

if __name__=="__main__":
	step = 8
	filelist = []
	for root, dirs, files in os.walk(rawDir):
		for file in files:
			if "gz" in file:
				sample = file.split("_")[0] + "_" + file.split("_")[1]
				if sample not in filelist:
					filelist.append(sample)

	if step < 1:
		for sample in filelist:
			sample_R1 = os.path.join(rawDir, sample + "_1.fq.gz")
			sample_R2 = os.path.join(rawDir, sample + "_2.fq.gz")
			#print(sample)
			#cutadapt_PE(sample_R1=sample_R1,sample_R2=sample_R2)
			#fastp(sample_R1=sample_R1, sample_R2=sample_R2)

	if step < 2:
		for sample in filelist:
			clean_R1 = os.path.join(cleanDir,sample + "_R1.fq.gz")
			clean_R2 = os.path.join(cleanDir,sample + "_R2.fq.gz")
			bowtie2(sample=sample,clean_R1=clean_R1,clean_R2=clean_R2)
			
	if step < 3:
		for root,dirs,files in os.walk(alignDir):
			for file in files:
				if ".sam" in file:
					sample = file.split(".")[0]
					file = os.path.join(root,file)
					print(sample)
					print(file)
					samtools(samfile=file,sample=sample)

	if step > 4:
		for root,dirs,files in os.walk(samtoolsDir):
			for file in files:
				if os.path.splitext(file)[1] == ".bam":
					sample = file.split("_")[0]
					file = os.path.join(root,file)
					#print(file)
					bamCoverage(rmdup_bam=file,sample=sample)
					#flagstat(sample=sample,rmdup_bam=file)

	if step < 5:
		for root,dirs,files in os.walk(samtoolsDir):
			for file in files:
				if os.path.splitext(file)[1] == ".bam":
					sample = file.split("_")[0]
					file = os.path.join(root,file)
					print(sample)
					print(file)
					bedgraph(bamfile=file,sample=sample)
					#macs2(inputfile=file,outname=sample)
					#seacr()

	if step < 6:
		for root,dirs,files in os.walk(macs2broadDir):
			for file in files:
				if os.path.splitext(file)[1] == ".broadPeak":
					sample = file.split("_")[0]
					file = os.path.join(root,file)
					print(sample)
					print(file)
					annotationPeaks(broadpeak=file,sample=sample)