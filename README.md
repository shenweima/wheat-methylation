# wheat-methylation
Analysis wheat methylation data using [BS-Seeker2](https://www.ncbi.nlm.nih.gov/pubmed/24206606) and [CGmapTools](https://www.ncbi.nlm.nih.gov/pubmed/28968643). These two softwares were developed by Dr. Guo Weilong and Thank him and his student Huang Xiangyi for their help.

The data used in this study came from [SRP133674](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA436361). If you use this data in your study, you could refer this information.  International Wheat Genome Sequencing Consortium (IWGSC). et al., "Shifting the limits in wheat research and breeding using a fully annotated reference genome.", Science, 2018 Aug 16;361(6403)

#### Download data
we can get download fastq files from [PRJNA436361](https://www.ebi.ac.uk/ena/data/view/PRJNA436361) and saved it in download.txt
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
Here, The paired-end reads were viewed as single end reads. Large fastq file was splited using seqkit tool and antisense of second reads were obtained using seqkit tool.
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