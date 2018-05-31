# HiCDB – user guide
Authors: Fengling Chen，Guipeng Li
## Overview
&emsp;&emsp;HiCDB is an open-source MATLAB/GNU OCTAVE code that detects the contact domain boundaries (CDBs) from Hi-C contact matrix. HiCDB.m takes raw or normalized contact matrix and outputs conservation annotated CDBs (option: singlemap) or differential CDBs (option: comparemap). visHiCDB.m takes raw or normalized contact matrix and HiCDB results and outputs visualization of CDBs on single Hi-C map (option: singlemap) or differential CDBs on two Hi-C maps (option: comparemap). You could run HiCDB and visHiCDB in MATLAB enviroment or in command line.  
&emsp;&emsp;Here is the general features of HiCDB software.  
<div align=center><img width="420" height="300" src="https://github.com/ChenFengling/HiCDB/blob/master/images/pipeline.png"/></div>

&emsp;&emsp;Here is the general steps of how we detect CDBs.  
<div align=center><img width="580" height="380" src="https://github.com/ChenFengling/HiCDB/blob/master/images/pipe.png"/></div>

### Requirements and install
install MATLAB or GNU OCTAVE and simply download HiCDB and run the scripts.
```shell
git clone https://github.com/ChenFengling/HiCDB.git
```
### Quick start
Unzip the testdata.tar.gz, you will  find the dense format Hi-C data of hESC (Doxin et al.). 
```shell
tar -zxvf testdata.tar.gz
```
To run HiCDB or visHiCDB in command line with nohup, simply add the HiCDB software path and use the code as in MATLAB environment.  
##### Run HiCDB to get the CDBs.
```shell
nohup matlab -r "addpath(genpath('/home/fchen/HiCDB/'));HiCDB({'/home/fchen/HiCDB/h1_rep1/'},40000,'hg19','ref','hg19');exit;" > mylog.txt < /dev/null &
```
This will take the  intra-chromosome matrix ('chr1.matrix',...,'chr23.matrix') in '/home/fchen/HiCDB/h1_rep1/' as input and set the resolution as 40000,chrsizes as 'hg19', the CTCF motif ref as 'hg19' and output the contact domain boundaries. "</dev/null &" shoud be added to avoid MATLAB entering an infinite loop.  
##### Run visHiCDB to display the region chr17:67100000-71100000.
```shell
 nohup matlab -r "addpath(genpath('/home/fchen/HiCDB/'));fig=visHiCDB({'/home/fchen/HiCDB/h1_rep1/chr17.matrix'},{'/home/fchen/HiCDB/h1_rep1/CDB.txt'},40000,17,67100000,71100000);exit;" > mylog2.txt < /dev/null &
```
You will get this output. The dot is CDB detected(light blue:not conserved CDB;dark blue:conserved CDB)  
<div align=center><img width="500" height="500" src="https://github.com/ChenFengling/HiCDB/blob/master/images/17_67100000_71100000_HiCmap.png"/></div>

#####  Get .bed file
As "chrX" is named as "chr23" and as "23" in the output CDB.txt file. You could use the following shell code to change CDB.txt into .bed file.    
```shell
awk -v OFS="\t" '{ print "chr"$1,$2,$3,$4,$5}' CDB.txt >CDB.bed
sed -i  's/chr23/chrX/g' CDB.bed
```  
## 1. Run HiCDB 
### Input
**hicfile:** The directory of all intra-chromosome matrix of a sample. The intra-chromosome matrix must be named as "chr+number.matrix" according to the chromosome order like 'chr1.matrix','chr2.matrix',...,'chr23.matrix'. As HiCDB matches "chr\*.matrix" to recognize the Hi-C matrix, avoid to use the "chr\*.matrix" as the name of other files. The intra-chromosome matrix could be in a dense (a NxN matrix) or sparse (a Kx3 table,Rao et al.) format. hicfile should be set as {'SAMPLE_DIR'} when option is "singlemap" or {'SAMPLE_DIR1','SAMPLE_DIR2'} when option is ‘comparemap’. This is required.  
Dense format contains	the contact	frequencies	of the Hi-C NxN	matrix.  
Sparse format (Rao et al.) has three fields: i, j, and M_i,j. (i and j are written as the left edge of the bin at a given resolution; for example, at 10 kb resolution, the entry corresponding to the first row and tenth column of the matrix would correspond to M_i,j, where i=0, j=90000). As the Hi-C matrix is symmetric, only the upper triangle of the matrix is saved in sparse format. An example is as following:     

| | | |
|-|-|-|
|  50000   |  50000    |  1.0 |
|  60000   |  60000    |  1.0 |
|  540000  |  560000   |  1.0 |
|  560000  | 560000    | 59.0 |
|  560000  |  570000   |  1.0 |
|  560000  |  600000   |  1.0 |
|  560000  |  700000   |  1.0 |
|  690000  |  710000   |  1.0 |
|  700000  |  710000   |  1.0 |
|  710000  | 710000    | 66.0 |

**resolution:** resolution of Hi-C matrix. This is required.  
**chrsizes:** Ordered chromosome sizes of the genome. Optional setting is ‘hg19’, ‘hg38’, ‘mm9’, ‘mm10’ or any other chromosome size files which can be generated following the instructions in annotation/README.md. This is required.  
**ref:** ref should be set when you want to get a cutoff using a CTCF motif or the option is 'comparemap'. Optional ref is ‘hg19’, ‘hg38’, ‘mm9’, ‘mm10’ or any other custom motif locus files which can be generated from instructions in annotation/README.md. Only ‘hg19’ and ‘hg38’ can be annotated with conservation.  
### Run
&emsp;&emsp;To run HiCDB, type **“HiCDB(hicfile,resolution,chrsizes,options & parameter values)”** in MATLAB. A detailed description of how to use HiCDB of all the options and of how to modify parameter values can be found by typing **“help HiCDB”** in MATLAB or **matlab -r "addpath(genpath('HiCDB_DIR/'));help HiCDB;exit;"** in command line. See also some examples below.
### Examples
#### 1. Output all the local maximum peaks and let customers to decide the cutoff.  
```matlab
HiCDB({'/home/data/sample'},10000,'hg19');
HiCDB({'/home/data/sample'},10000,'hg19','outdir','/home/data/sample/outputs/');       
```
#### 2. Use GSEA-like methods to decide the cutoff . 
```matlab
HiCDB({'/home/data/sample'},10000,'hg38','ref','hg38');
HiCDB({'/home/data/sample'},10000,'hg19','ref','custom_motiflocs.txt');
```
#### 3. To detect differential CDBs
```matlab
HiCDB({'/home/data/sample1','/home/data/sample2'},10000,'hg19','option','comparemap','ref','hg19');
```  
#### 4. To detect differential CDBs with replicates
```matlab
HiCDB({{'sample1_rep1/','sample1_rep2/'},{'sample2_rep1/','sample2_rep2/'}},'hg19',10000,'option','comparemap','ref','hg19');
```
### Output(s)
**1.CDB.txt:**

| chr | start|  end|    LRI|    avgRI   | conserve_or_not |     consistent_or_differential |
|:------:|:----:|:------:|:----:|:------:|:----:|:----:|
| 19    | 53100000      | 53140000      | 0.394707211   | 0.647392804 | 0 |     1 |
| 16    | 5060000         | 5100000       | 0.342727704 | 0.663101081 | 1 |     1 |
| 19    | 19620000  | 19660000  |       0.329837698 | 0.609237673 |     1 |     0 |

**2. localmax.txt:** all the local maximum peaks detected before cutoff decision. User can decide custom CDB cutoff upon this file.  
**3. EScurve.png:** CTCF motif enrichment on ranked local maximum peaks.  
These output files can be found in custom output directory or default directory namely the directory of the first sample.      
## 2. Run visHiCDB
### Input
**hicfile:** the same as in function HiCDB  
**CDBfile:** CDBfile sould be a cell array storing the CDB location. The CDB files should be formatted as the output files of function HiCDB.  
**resolution:** resolution of Hi-C map.  
**chr,startloc,endloc:** observation locus on Hi-C map.  
### Run
&emsp;&emsp;To run visHiCDB, type **“fig=visHiCDB(hicfile,CDBfile,resolution,chr,startloc,endloc,options & parameter values)”** in the command line with input matrix, CDB file, options and parameter values. A detailed description of how to use visHiCDB of all the options and of how to modify parameter values can be found by typing **“help visHiCDB”** in MATLAB or **matlab -r "addpath(genpath('HiCDB_DIR/'));help visHiCDB;exit;"** in command line. See also some examples below.  
### Examples
#### 1.Show CDB on single Hi-C map
```matlab
fig=visHiCDB({'/home/data/sample/chr18.matrix'},{'CDB1.txt'},10000,18,25000000,31150000);
```
#### 2. Show differential CDBs on Hi-C maps 
```matlab
fig=visHiCDB({'/home/data/sample1/chr18.matrix','/home/data/sample2/chr18.matrix'},{'/home/data/sample1/CDB1.txt','/home/data/sample2/CDB2.txt'},10000,18,25000000,31150000,'option','comparemap');
```
### Output(s)
**HiCmap.pdf：** a pdf containing figure showing CDBs on single Hi-C map or different kinds of CDBs between two Hi-C maps.  
These output files can be found in custom output directory or default directory namely the directory of the first sample. 


