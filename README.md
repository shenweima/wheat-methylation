# wheat-methylation
Analysis wheat methylation data using [BS-Seeker2](https://www.ncbi.nlm.nih.gov/pubmed/24206606) and [CGmapTools](https://www.ncbi.nlm.nih.gov/pubmed/28968643). These two softwares were developed by Dr. Guo Weilong and Thank him and his student Huang Xiangyi for their help.

The data used in this study came from [SRP133674](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA436361). If you use this data in your study, you could refer this information.  International Wheat Genome Sequencing Consortium (IWGSC). et al., "Shifting the limits in wheat research and breeding using a fully annotated reference genome.", Science, 2018 Aug 16;361(6403)

#### Download data
we can get website of fastq files from [PRJNA436361](https://www.ebi.ac.uk/ena/data/view/PRJNA436361) and saved it in download.txt
```shell
for i in $(cat download.txt); do axel -n 20 -a $i; done
```
#### Index genome
```shell
python2 ./bs_seeker2-build.py -f /data2/Index/161010_Chinese_Spring_v1.0_pseudomolecules_parts.fasta --aligner=bowtie2 -d /data2/Fshare/FastaAndIndex/bs_seeker2_index/
#Considering this step is very slow，we can set 'bowtie2-build' to 'alias bowtie2-build='bowtie2-build --large-index --threads 10'. Thus, we can use mutiple CPUs
```
#### Filter low quality reads
```shell
for i in $(cat input.txt); do fastp -p -w 8 -l 30 -i ${i}_1.fastq.gz -I ${i}_2.fastq.gz -o ${i}_1.filter.fq.gz -O ${i}_2.filter.fq.gz -h ${i}.html; done
```
#### Mapping reads to wheat genome
Here, The paired-end reads were viewed as single end reads.Large fastq file was splitted using seqkit tool and antisense of second reads were obtained using seqkit tool.
```python
#!/usr/bin/python
# -*- coding: utf-8 -*-
__author__ = 'shengwei ma'
__author_email__ = 'shengweima@icloud.com'

import subprocess
from concurrent.futures import ThreadPoolExecutor as Pool


def task(fq):
    proc = subprocess.Popen(['/root/software/BSseeker2/bs_seeker2-align.py', '-i', fq, '-g', '/data2/Fshare/FastaAndIndex/IWGSC_v1.0_bwa/161010_Chinese_Spring_v1.0_pseudomolecules_parts.fasta', '-d', '/data2/Fshare/FastaAndIndex/bs_seeker2_index/', '--temp_dir=/data2/user_data/temp/', '--XSteve', '--aligner=bowtie2', '--bt2-p', '4', '--bt2--end-to-end', '--bt2--very-sensitive', '--bt2--dovetail', '-o', fq.split('fq')[0] + 'bam'], shell=False)
    proc.wait()
    print ('Finish mapping' + fq)

    
if __name__ == '__main__':
    pool = Pool(10)
    lin = [line.strip() for line in open('input_fastq.txt', 'r')]
    # this input file contains fastq file name 
    for li in lin:
        dir0 = li.split('.')[0]
        print (li)
        print (dir0)
        if '_1' in li:
            proc = subprocess.Popen(['seqkit', 'split', '-p', '10', '-j', '10', '-O', dir0, 'fastq/' + li], shell=False)
            proc.wait()
        if '_2' in li:
            anti = li.split('.')[0] + '.filter.fq.gz'
            print (anti)       
            proc = subprocess.Popen(['seqkit','seq','-p','-r','-j', '16', '-o', anti, 'fastq/' + li], shell=False)
            proc.wait()
            print ('Finsh antisense')
            proc = subprocess.Popen(['seqkit', 'split', '-p', '10', '-j', '10', '-O', dir0, anti], shell=False)
            proc.wait()
        infiles = ['001','002','003','004','005','006','007','008','009','010']
        infiles = [dir0 + '/' + dir0 + '.filter.part_' + n + '.fq.gz' for n in infiles]

        r = pool.map(task, infiles)
        proc = subprocess.Popen('samtools merge -@ 10 ' + dir0 + '.bam ' + dir0 + '/*.bam', shell=True)
        proc.wait()
```
#### Merge bams
Use `samtools merge` to merge bam files
```shell
samtools merge -@ 6 out.bam in1.bam in2.bam in3.bam
```
#### Sort bams
```shell
samtools sort -@ 10 -o out.bam in.bam
```
#### Split bam file by chromosome
```shell
for i in `samtools view -H cs_leaf2_sorted.bam | awk -F"\t" '/@SQ/{print $2}' |  cut -d":" -f2`; do samtools view -h -F 0x4 cs_leaf2_sorted.bam $i | samtools view -hbS - > cs_leaf2_sorted.$i.bam ; done
```
#### Call methylation
Use the python script to call methylation in parallel
this step will be generated three files(bw, CGmap,ATCGmap)
```python
#!/usr/bin/python
# -*- coding: utf-8 -*-
__author__ = 'shengwei ma'
__author_email__ = 'shengweima@icloud.com'

import subprocess
from concurrent.futures import ThreadPoolExecutor as Pool


def task(bam):
    print (bam)
    proc = subprocess.Popen(['/root/software/BSseeker2/bs_seeker2-call_methylation.py', '-i', bam, '-o', bam.rstrip(".bam"), '--sorted', '-d', '/data2/Fshare/FastaAndIndex/bs_seeker2_index/161010_Chinese_Spring_v1.0_pseudomolecules_parts.fasta_bowtie2/', '-x'], shell=False)
    proc.wait()
    print ('Finsh mapping ' + bam)

    
if __name__ == '__main__':
    pool = Pool(10)
    bams = []
    cs_bam = ['cs_leaf2_sorted']
    infiles = ["chr1A_part1", "chr1A_part2", "chr1B_part1", "chr1B_part2", "chr1D_part1", "chr1D_part2", "chr2A_part1", "chr2A_part2", "chr2B_part1", "chr2B_part2", "chr2D_part1", "chr2D_part2", "chr3A_part1", "chr3A_part2", "chr3B_part1", "chr3B_part2", "chr3D_part1", "chr3D_part2", "chr4A_part1", "chr4A_part2", "chr4B_part1", "chr4B_part2", "chr4D_part1", "chr4D_part2", "chr5A_part1", "chr5A_part2", "chr5B_part1", "chr5B_part2", "chr5D_part1", "chr5D_part2", "chr6A_part1", "chr6A_part2", "chr6B_part1", "chr6B_part2", "chr6D_part1", "chr6D_part2", "chr7A_part1", "chr7A_part2", "chr7B_part1", "chr7B_part2", "chr7D_part1", "chr7D_part2", "chrUn"]
    for cs in cs_bam:
        for i in infiles:
            bams.append(cs + '.' + i + '.bam')
    r = pool.map(task, bams)
```
#### Combine output files
```shell
cat cs_leaf2_sorted.chr*.CGmap.gz > cs_leaf2_sorted.CGmap.gz
cat cs_leaf2_sorted.chr*.bw.gz > cs_leaf2_sorted.bw.gz
cat cs_leaf2_sorted.chr*.ATCGmap.gz > cs_leaf2_sorted.ATCGmap.gz
```
#### Combine two parts of the chromosome into a whole chromosome
Here, we use a python script to combine it.
```python
#!/usr/bin/python
# -*- coding: utf-8 -*-
__author__ = 'shengwei ma'
__author_email__ = 'shengweima@icloud.com'


chr = [['chr1A', 471304005], ['chr1B', 438720154], ['chr1D', 452179604], ['chr2A', 462376173], ['chr2B', 453218924],
       ['chr2D', 462216879], ['chr3A', 454103970], ['chr3B', 448155269], ['chr3D', 476235359], ['chr4A', 452555092],
       ['chr4B', 451014251], ['chr4D', 451004620], ['chr5A', 453230519], ['chr5B', 451372872], ['chr5D', 451901030],
       ['chr6A', 452440856], ['chr6B', 452077197], ['chr6D', 450509124], ['chr7A', 450046986], ['chr7B', 453822637],
       ['chr7D', 453812268]]

with open('cs_leaf2.CGmap', 'r') as f:
    for line in f:
        line = line.replace('_part1', '')
        line = line.strip().split('\t')
        if line[0].endswith('part2'):
            for i in chr:
                if line[0].split('_')[0] == i[0]:
                    line[2] = int(line[2]) + int(i[1])
                    line[0] = line[0].split('_')[0]
        for m in line[:-1]:
            print (str(m) + '\t'),
        print (line[-1] + '\n'),
```
Now we can do downstream analysis using [CGmapTools](https://www.ncbi.nlm.nih.gov/pubmed/28968643).
