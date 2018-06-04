function []=HiCDB(hicfile,resolution,chrsizes,varargin)
%   Detect CDBs and differential CDBs on Hi-C matrix
%   Implemented by
%   Fengling Chen
%   May 2017
%   Tsinghua University
%   https://github.com/ChenFengling/HiCDB
%   Any suggestions and remarks might be addressed to Fengling Chen:cfl15@mails.tsinghua.edu.cn
%
%HiCDB using Hi-C raw contact matrix to detect Hi-C contact domain boundaries(CDBs).
%      it outputs annotated CDBs, differential CDBs  on the chosen options
%
%    *Possible outputs*:
%    1.CDB.txt: 
%    #chr   start      end        LRI          avgRI  conserve_or_not  consistent_or_differential 
%    19	53100000	53140000	0.394707211	0.647392804     0           1
%    16	5060000  	5100000	    0.342727704	0.663101081     1           1
%    19	19620000	19660000	0.329837698	0.609237673     1           1
%    7	71660000  	71700000	0.316790592	0.518786767     0           1
%    19	3420000	    3460000	    0.314223401	0.599579273     0           0
%    2.localmax.txt: all the local maximum peaks detected before cutoff
%    decision. User can decide custum CDB cutoff upon this file.
%    3.EScurve.png: CTCF motif enrichment on ranked local maximum peaks.
%    4.aRI.txt: average RI score for each  genomic bin.  
%    5.LRI.txt: LRI score for each genomic bin.
%    
%    HiCDB(hicfile,resolution,chrsizes,'ref',refvalue),where hicfile,resolution and chrsizes are required
%    for calculate the CDBs. hicfile shound be a cell array storing the directory of all intra-chromosome Hi-C matrixes with 
%    sparse or dense format. The intra-chromosome matrix must be named as "chr+number.matrix" 
%    according to the chromosome order like 'chr1.matrix','chr2.matrix',..., 'chr23.matrix'. 
%    As HiCDB matches "chr*.matrix" to recognize the Hi-C matrix, avoid to use the "chr*.matrix" as 
%    the name of other files. The intra-chromosome matrix could be in a dense (a NxN matrix) or sparse (a Kx3 table,Rao et al.) format. 
%    If ref is not set, this function will
%    output all the local maximum peaks. If ref is set, this function will
%    output local maximum peaks and final CDBs.If ref is 'hg38' or 'hg19',
%    CDBs will also be annotated as conserved or not conserved.
% 
%    HiCDB(hicfile,resolution,chrsizes,'option','comparemap','ref',refvalue), where
%    hicfile shoule be a cell arrary containing the directories of comparing sample files.
%    If you don't have replicate, set hicfile as {{'SAMPLE1_DIR'},{'SAMPLE2_DIR'}}.This function will
%    first perform CDB detection on each sample and then compare the difference between their final CDBs by intersection.
%    If  replicates is provided, set hicfile as {{'SAMPLE1_REP1_DIR','SAMPLE1_REP2_DIR'},{'SAMPLE2_REP1_DIR','SAMPLE2_REP2_DIR'}}.
%    The function will find CDBs on each sample with  merged Hi-C matrix, calculate aRI score on each replicates, 
%    then decide a CDB as differetial or not by statistical test on  aRI scores of each CDB. 
%    ref must be set when option is 'comparemap'. If ref is 'hg38' or 'hg19', CDBs will also be annotated
%    as conserved or not conserved.
% 
%    [ ... ] = 	HiCDB(..., 'PARAM1',val1, 'PARAM2',val2, ...) specifies
%     optional parameter name/value pairs to detect the CDBs or to find
%     differential CDBs and/or to save the data desired directory. 
%     Parameters are : 
%    
%
%    'option'     -  Choices are :
%         'singlemap'    - detect CDBs on one sample. Default value is 'singlemap'
%         'comparemap'   - compare CDBs of two samples. This option will
%         use GSEA-like method to decide cut off by default.
% 
% 
%    'ref'   - reference CTCF motif locs on the genome. If it is set, the output will use
%     the GSEA-like methods to decide the cutoff. Default is 'no'. Choices are :
%        'no'
%        'hg19'
%        'hg38'
%        'mm9'
%        'mm10'
%         or other customfile  for example 'genome.txt' made following the
%         instruction in annotation/README.md.
%          Example for 'genome.txt':
%              #chr   motifcenterlocus
%              10	15100928
%              10	15188593
%     
%     'outdir'   -The output directory. Default will be the directory of
%     the first sample.
%
%     'mind'    - Minimum separation between local maximum peaks (measured by bin), specified 
%                 as a positive integer scalar. Use this argument to have findpeaks 
%                 ignore small peaks that occur in the neighborhood of a
%                 larger peak. Recommend to set as 40000/resolution.
% 
%     'wd'      - The smallest window sizes. 
%     
%     'wdsize'  - The number of different window size. The whole window size 
%                  scale will be wd:(wd+wdsize-1).
%     
%      default value for 'mind','wd','wdsize' on different resolution
%     resolution mind  wd wdsize
%     10k         4    3   6
%     40k         2    1   6
%      5k         8    6   8
%
%           
% 
%    Examples:
%
%       1. Output all the local maximum peaks and let customers to decide the cutoff.
%       HiCDB({'sample1/'},10000,'hg19');
%       HiCDB({'sample1/'},10000,'hg19','outdir','sample1/outputs/');       
%       2. Use GSEA-like methods to decide the cutoff
%       HiCDB({'sample1/'},10000,'hg38','ref','hg38');
%       HiCDB({'sample1/'},10000,'custom_chrsize.txt','ref','custom_motiflocs.txt');
%       3. To detect differential CDBs
%       HiCDB({'sample1/','sample2/'},'hg19',10000,'option','comparemap','ref','hg19');
%       4. To detect differential CDBs with replicates 
%       HiCDB({{'sample1_rep1/','sample1_rep2/'},{'sample2_rep1/','sample2_rep2/'}},'hg19',10000,'option','comparemap','ref','hg19');
%
%   See also visHiCDB
%

%% Inputs
if length(hicfile)==1||length(hicfile)==2
    if length(hicfile)==1
    pathstr = hicfile{1};
    else
    pathstr = hicfile{1}{1};   
    end
else
   error('invalid cell number for hicfile') 
end

p = inputParser;
p.FunctionName = 'HiCDB';

defaultOption = 'singlemap';
validOption = {'singlemap','comparemap'};
checkOption = @(x) any(validatestring(x,validOption));

defaultOutdir = pathstr;
defaultRef='no'; 

switch resolution
    case 10000
        defaultMind=4;
        defaultWd=3;
        defaultWdsize=6;
    case 40000
        defaultMind=2;
        defaultWd=1;
        defaultWdsize=6;
    case 5000
        defaultMind=8;
        defaultWd=6;
        defaultWdsize=8;
    otherwise 
        error(' The default value is not set on this Hi-C resolution, please use help HiCDB to see how to set parameters including mind,wd and wdsize.'); 
end

addRequired(p,'hicfile');
addRequired(p,'resolution',@isnumeric);
addRequired(p,'chrsizes')
addParamValue(p,'option',defaultOption,checkOption)
addParamValue(p,'outdir',defaultOutdir)
addParamValue(p,'ref',defaultRef)
addParamValue(p,'mind',defaultMind)
addParamValue(p,'wd',defaultWd)
addParamValue(p,'wdsize',defaultWdsize)


p.parse(hicfile,resolution,chrsizes,varargin{:})

option = p.Results.option
outdir = p.Results.outdir
ref = p.Results.ref
mind = p.Results.mind;
wd = p.Results.wd;
wdsize = p.Results.wdsize;


validateattributes(option,{'char'},{},'HiCDB','option')
validateattributes(outdir,{'char'},{},'HiCDB','outdir')
validateattributes(ref,{'char'},{},'HiCDB','ref')
validateattributes(mind,{'numeric'},{'positive'},'HiCDB','mind')
validateattributes(wd,{'numeric'},{'positive'},'HiCDB','wd')
validateattributes(wdsize,{'numeric'},{'positive'},'HiCDB','wdsize')



if strcmp(option,'singlemap')==1
    if ~(length(hicfile)==1)
       error('The length of hicfile should be 1 when option equals ''singlemap'' ');
    end      
elseif strcmp(option,'comparemap')==1
     if ~(length(hicfile)==2)
       error('The length of hicfile should be 2 when option equals ''comparemap'' ');
     end
     if strcmp(ref,'no')==1
       error('the ref parameter should be set when option equals ''comparemap'' ');
     end
else
     error('Invalid option setting,option should be setted as ''singlemap'' or ''comparemap'' ');
end

if ~strcmp(ref,'no')==1
           switch ref
               case 'hg38'
                   motiffile='annotation/CTCF_hg38_site.txt';
               case 'hg19'
                   motiffile='annotation/CTCF_hg19_site.txt';
               case 'mm9'
                   motiffile='annotation/CTCF_mm9_site.txt';
               case 'mm10'
                   motiffile='annotation/CTCF_mm10_site.txt';
               otherwise
                   motiffile=ref;
           end        
end

switch chrsizes
	case 'hg38'
	   chrsizes=dlmread('annotation/hg38.chrom.sizes.txt');
        case 'hg19' 
           chrsizes=dlmread('annotation/hg19.chrom.sizes.txt');
	case 'mm9'
	   chrsizes=dlmread('annotation/mm9.chrom.sizes.txt');
        case  'mm10'
	   chrsizes=dlmread('annotation/mm10.chrom.sizes.txt');
   otherwise
	   chrsizes=dlmread(chrsizes);
end
%
if strcmp(option,'singlemap')==1
        hicfile=hicfile{1};
        tmp=dir([hicfile,'/chr*.matrix']);
        totalchr=length(tmp);
        allpeaks=[];
	allLRI=[];
	allaRI=[];
        disp('* Detect all local maximum peaks on each chromosome.')
        for i=1:totalchr
          disp(['Processing chr' num2str(i)]);
          hicmap=[hicfile,'/chr',num2str(i),'.matrix'];  %
          [LRI,aRI,localmax]= chrpeaks(hicmap,resolution,chrsizes(i),'KR',mind,wd,wdsize,i);
          allpeaks=[allpeaks;localmax];
	      allLRI=[allLRI;LRI];
	      allaRI=[allaRI;aRI];
        end
        [~,peakidx]=sort(allpeaks(:,4),'descend');
        allpeaks=allpeaks(peakidx,:);
        disp(' ')
        disp('* Write outputs: localmax.txt.')	
        dlmwrite([outdir '/localmax.txt'],allpeaks,'delimiter','\t','precision', 9);
        dlmwrite([outdir '/aRI.txt'],allaRI,'delimiter','\t','precision', 9)
	dlmwrite([outdir '/LRI.txt'],allLRI,'delimiter','\t','precision', 9);


        if ~strcmp(ref,'no')==1
        disp(' ')
        disp('* Decide cut off using enrichmentscore. Write outputs:EScurve.png.')
        [CDB,fig]=EScutoff(allpeaks,motiffile,outdir);
        print(fig,'-dpng',[outdir '/EScurve.png']);
        disp(' ')
          if strcmp(ref,'hg19')==1||strcmp(ref,'hg38')==1
          disp('* Annotate CDBs with CDB conservation. Write outputs:CDB.txt.')
          disp('CDB.txt format:#chr start end LRI avgRI conserve_or_not')
              if strcmp(ref,'hg19')==1
                 conserve=load('annotation/conserved_CDB_hg19.bed');
              else
                 conserve=load('annotation/conserved_CDB_hg38.bed');
              end
              conserve=[conserve(:,1),(conserve(:,2)+conserve(:,3))/2];
              CDB=[CDB,zeros(size(CDB,1),1)];
              for conchr=1:23
              contmp=conserve(conserve(:,1)==conchr,2);
              chrpos=find(CDB(:,1)==conchr);
              CDBcenter=(CDB(chrpos,2)+CDB(chrpos,3))/2;
              pos=min(pdist2(contmp,CDBcenter)) < 40000+resolution/2;
              CDB(chrpos,6)=pos; 
              end 
          else
          disp('* CDB annotation step is only for hg19 and hg38! Write outputs:CDB.txt.')  
          disp('CDB.txt format:#chr start end LRI avgRI')
          end
        dlmwrite([outdir '/CDB.txt'],CDB,'delimiter','\t','precision', 9);
        end

        
elseif strcmp(option,'comparemap')==1
     if length(hicfile{1})==1
            hicfile1=hicfile{1}{1};
            hicfile2=hicfile{2}{1};
            % sample 1
            disp(['A. Processing the first sample ',hicfile1]);
            tmp=dir([hicfile1,'/chr*.matrix']);
            totalchr=length(tmp);
            allpeaks=[];
            allaRI=[];
            allLRI=[];
            disp(' ')
            disp('* Detect all local maximum peaks on each chromosome.')
            for i=1:totalchr
              disp(['Processing chr' num2str(i)]);
              hicmap=[hicfile1,'/chr',num2str(i),'.matrix'];  %
              [LRI,aRI,localmax]= chrpeaks(hicmap,resolution,chrsizes(i),'KR',mind,wd,wdsize,i);
              allpeaks=[allpeaks;localmax];
              allaRI=[allaRI;aRI];
              allLRI=[allLRI;LRI]; 
            end
            [~,peakidx]=sort(allpeaks(:,4),'descend');
            allpeaks=allpeaks(peakidx,:);
            disp(' ')
            disp('* Write outputs: localmax1.txt ')
            dlmwrite([outdir '/localmax1.txt'],allpeaks,'delimiter','\t','precision', 9);
           dlmwrite([outdir '/aRI1.txt'],allaRI,'delimiter','\t','precision', 9);
           dlmwrite([outdir '/LRI1.txt'],allLRI,'delimiter','\t','precision', 9);
            disp(' ')
            disp('* Decide cut off using enrichmentscore. Write outputs: EScurve1.png')
            [CDB1,fig]=EScutoff(allpeaks,motiffile,outdir);
            print(fig,'-dpng',[outdir '/EScurve1.png']);

           % sample 2
            disp(' ')
            disp(['B. Processing the second sample ' hicfile2]);
            tmp=dir([hicfile2,'/chr*.matrix']);
            totalchr=length(tmp);
            allpeaks=[];
            allaRI=[];
            allLRI=[];
            disp(' ')
            disp('* Detect all local maximum peaks on each chromosome.')
            for i=1:totalchr
              disp(['Processing chr' num2str(i)]);
              hicmap=[hicfile2,'/chr',num2str(i),'.matrix'];  %
              [LRI,aRI,localmax]= chrpeaks(hicmap,resolution,chrsizes(i),'KR',mind,wd,wdsize,i);
              allpeaks=[allpeaks;localmax];
             allaRI=[allaRI;aRI];
             allLRI=[allLRI;LRI];
            end
            [~,peakidx]=sort(allpeaks(:,4),'descend');
            allpeaks=allpeaks(peakidx,:);
            disp(' ')
            disp('* Write outputs: localmax2.txt ')
            dlmwrite([outdir '/localmax2.txt'],allpeaks,'delimiter','\t','precision', 9);
            dlmwrite([outdir '/aRI2.txt'],allaRI,'delimiter','\t','precision', 9);
            dlmwrite([outdir '/LRI2.txt'],allLRI,'delimiter','\t','precision', 9);
        disp(' ')
            disp('* Decide cut off using enrichmentscore. Write outputs: EScurve2.png')
            [CDB2,fig]=EScutoff(allpeaks,motiffile,outdir);
            print(fig,'-dpng',[outdir '/EScurve2.png']);

            % compare two samples
            disp(' ')
            disp('C.Find differential CDBs.')
              if strcmp(ref,'hg19')==1||strcmp(ref,'hg38')==1
              disp('* Annotate differential CDBs with CDB conservation. Write outputs:CDB1.txt,CDB2.txt.')
              disp('CDB.txt format:#chr start end LRI avgRI conserve_or_not consistent_or_differential')
                  if strcmp(ref,'hg19')==1
                     conserve=load('annotation/conserved_CDB_hg19.bed');
                  else
                     conserve=load('annotation/conserved_CDB_hg38.bed');
                  end
                  conserve=[conserve(:,1),(conserve(:,2)+conserve(:,3))/2];
                  CDB1=[CDB1,zeros(size(CDB1,1),2)];
                  CDB2=[CDB2,zeros(size(CDB2,1),2)];
                  for conchr=1:23
                  contmp=conserve(conserve(:,1)==conchr,2);
                  %conserve  or not  
                  chrpos1=find(CDB1(:,1)==conchr);
                  CDBcenter1=(CDB1(chrpos1,2)+CDB1(chrpos1,3))/2;
                  pos1=min(pdist2(contmp,CDBcenter1)) < 40000+resolution/2;
                  CDB1(chrpos1,6)=pos1; 

                  chrpos2=find(CDB2(:,1)==conchr);
                  CDBcenter2=(CDB2(chrpos2,2)+CDB2(chrpos2,3))/2;
                  pos2=min(pdist2(contmp,CDBcenter2))< 40000+resolution/2 ;
                  CDB2(chrpos2,6)=pos2; 

                  %differential 
                  pos=min(pdist2(CDBcenter2(pos2),CDBcenter1(pos1)))<= 40000;  
                  CDB1(chrpos1(pos1),7)=pos;
                  pos=min(pdist2(CDBcenter1(pos1),CDBcenter2(pos2)))<= 40000; 
                  CDB2(chrpos2(pos2),7)=pos;

                  pos=min(pdist2(CDBcenter2(~pos2),CDBcenter1(~pos1)))<= resolution+resolution/2;  
                  CDB1(chrpos1(~pos1),7)=pos;
                  pos=min(pdist2(CDBcenter1(~pos1),CDBcenter2(~pos2)))<= resolution+resolution/2; 
                  CDB2(chrpos2(~pos2),7)=pos;              
                  end  
              else 
              disp('CDB annotation step is only for hg19 and hg38!')
              disp('Write outputs without conservation annotation: CDB1.txt, CDB2.txt')      
              disp('CDB.txt format:#chr start end LRI avgRI consistent_or_differential')
                  CDB1=[CDB1,zeros(size(CDB1,1),1)];
                  CDB2=[CDB2,zeros(size(CDB2,1),1)];
                  for conchr=1:23
                  %differential 
                  chrpos1=find(CDB1(:,1)==conchr);
                  CDBcenter1=(CDB1(chrpos1,2)+CDB1(chrpos1,3))/2;
                  chrpos2=find(CDB2(:,1)==conchr);
                  CDBcenter2=(CDB2(chrpos2,2)+CDB2(chrpos2,3))/2;
                  pos=min(pdist2(CDBcenter2,CDBcenter1))<= resolution+resolution/2;  
                  CDB1(chrpos1,6)=pos;
                  pos=min(pdist2(CDBcenter1,CDBcenter2))<= resolution+resolution/2; 
                  CDB2(chrpos2,6)=pos;              
                  end  
             end
             dlmwrite([outdir '/CDB1.txt'],CDB1,'delimiter','\t','precision', 9);
             dlmwrite([outdir '/CDB2.txt'],CDB2,'delimiter','\t','precision', 9);
     else
             n1=length(hicfile{1});
             n2=length(hicfile{2});
             tmp=dir([hicfile{1}{1},'/chr*.matrix']);
             totalchr=length(tmp);
             allpeaks={[],[]};
             aRI=[];
	     aRI1=[];
	     aRI2=[];
	     LRI1=[];
	     LRI2=[];
             disp('A. Calculate the aRI score of each Hi-C replicate.')
             for  i=1:totalchr
                 disp(['Processing chr',num2str(i)]);
                 [aRItmp,aRI1tmp,aRI2tmp,LRI1tmp,LRI2tmp,resultA,resultB]=chrpeaks_diff(hicfile,resolution,chrsizes(i),mind,wd,wdsize,i);
                 allpeaks{1}=[allpeaks{1};resultA];
                 allpeaks{2}=[allpeaks{2}; resultB];
                 %size_merge=rbind(size_merge,cbind(out$resultA$aRI(,4),out$resultB$aRI(,4)))
                 N=ceil(chrsizes(i)/resolution)-1;
                 loci=[i*ones(1,N);(1:N)*resolution-resolution/2;[(1:N-1)*resolution+resolution/2,chrsizes(i)]]';
                 aRI=[aRI;[loci,aRItmp(1:N,:)]];
                 aRI1=[aRI1;aRI1tmp];
		 LRI1=[LRI1;LRI1tmp];
                 aRI2=[aRI2;aRI2tmp];
		 LRI2=[LRI2;LRI2tmp];
             end
             % first sample
             disp('B. Get CDBs in the first condition.');
             [~,peakidx]=sort(allpeaks{1}(:,4),'descend');
             allpeaks{1}=allpeaks{1}(peakidx,:);
             disp(' ')
             disp('* Write outputs: localmax.txt.')
             dlmwrite([outdir '/localmax1.txt'],allpeaks{1},'delimiter','\t','precision', 9);
             dlmwrite([outdir '/aRI1.txt'],aRI1,'delimiter','\t','precision', 9);
             dlmwrite([outdir '/LRI1.txt'],LRI1,'delimiter','\t','precision', 9);
	     disp(' ')
             disp('* Decide cut off using enrichmentscore. Write outputs:EScurve1.png.')
             [CDB1,fig]=EScutoff(allpeaks{1},motiffile,outdir);
             print(fig,'-dpng',[outdir '/EScurve1.png']);
             %
             disp('C. Get CDBs in the second condition.');
             [~,peakidx]=sort(allpeaks{2}(:,4),'descend');
             allpeaks{2}=allpeaks{2}(peakidx,:);
             disp(' ')
             disp('* Write outputs: localmax.txt.')
             dlmwrite([outdir '/localmax2.txt'],allpeaks{2},'delimiter','\t','precision', 9);
             dlmwrite([outdir '/aRI2.txt'],aRI2,'delimiter','\t','precision', 9);
	     dlmwrite([outdir '/LRI2.txt'],LRI2,'delimiter','\t','precision', 9);
	     disp(' ')
             disp('* Decide cut off using enrichmentscore. Write outputs:EScurve2.png.')
             [CDB2,fig]=EScutoff(allpeaks{2},motiffile,outdir);
             print(fig,'-dpng',[outdir '/EScurve2.png']);
             %%
             disp('D. Differential CDBs testing.')
             %merged peaks
             disp('* Write outputs: allaRI.txt. genomic locus + normalized aRI score in each replicates' )
	     dlmwrite([outdir '/allaRI.txt'],aRI,'delimiter','\t','precision', 9);
	     aRIsize=nanmean(aRI(:,4:(3+n1+n2)));
             norm=repmat(aRIsize,[size(aRI,1),1]);
             aRI(:,4:(3+n1+n2))=aRI(:,4:(3+n1+n2))./norm*nanmean(aRI(:,4));              
             CDB=unique([CDB1(:,1:3);CDB2(:,1:3)],'rows','stable');
             ind1=ismember(CDB,CDB1(:,1:3),'rows');
             ind2=ismember(CDB,CDB2(:,1:3),'rows');
             CDB=[CDB,NaN(size(CDB,1),4+n1+n2)];
             CDB(ind1,4:5)=CDB1(:,4:5);
             CDB(ind2,6:7)=CDB2(:,4:5);
             CDB=sortrows(CDB,[1 2]);
             [C,ia,ib]=intersect(CDB(:,1:3),aRI(:,1:3),'rows');
             CDB(ia,8:(7+n1+n2))=aRI(ib,4:(3+n1+n2));
             n=size(CDB,1);
             idx=zeros(n,3);
             idx(:,1)=[0;CDB(2:n,2)==CDB(1:(n-1),3)];
             idx(:,2)=isnan(CDB(:,6));
             idx(:,3)=isnan(CDB(:,4));
             pos=find(idx(:,1));
            CDB_other=CDB;
            CDB_other([pos;pos-1],:)=[];
            CDB_tmp=[CDB(pos-1,4:7),CDB(pos,4:7)];
            CDB_tmp(isnan(CDB_tmp(:,1)),1:2)=CDB_tmp(isnan(CDB_tmp(:,1)),5:6);
            CDB_tmp(isnan(CDB_tmp(:,3)),3:4)=CDB_tmp(isnan(CDB_tmp(:,3)),7:8);
            CDB_overlap=NaN(length(pos),size(CDB,2));
            CDB_overlap(:,1)=CDB(pos-1,1);
            CDB_overlap(:,2:3)=(CDB(pos-1,2:3)+CDB(pos,2:3))/2;
            CDB_overlap(:,4:7)=CDB_tmp(:,1:4);
            CDB_overlap(:,8:(7+n1+n2))=(CDB(pos-1,8:(7+n1+n2))+CDB(pos,8:(7+n1+n2)))/2;
            CDB=[CDB_other;CDB_overlap];
            % CDB normalize ?
            %CDBsize=nanmean(CDB(:,8:(7+n1+n2)));
	    %norm=repmat(CDBsize,[size(CDB,1),1]);
	    %CDB(:,8:(7+n1+n2))=CDB(:,8:(7+n1+n2))./norm*nanmean(CDB(:,7));
	    %
	    diff=sum(CDB(:,(8+n1):(7+n1+n2)),2)-sum(CDB(:,8:(7+n1)),2);
            CDB=[CDB,diff];
            pvalue=[];
            for i=1:size(CDB,1)
                    [~,p]=ttest2(CDB(i,(8+n1):(7+n1+n2)),CDB(i,8:(7+n1)));
                    pvalue=[pvalue;p];
            end
            CDB=[CDB,pvalue];
            disp('* Write outputs: CDB_merged.txt. genomic locus + score in merged data + normalized aRI score in each replicates +difference +pvalue')
            dlmwrite([outdir '/CDB_merged.txt'],CDB,'delimiter','\t','precision', 9);
            
            CDB_A=CDB(((pvalue<0.05 & diff< -quantile(abs(diff),0.5))|diff<quantile(diff,0.1)) & isnan(CDB(:,6)),:);
            CDB_B=CDB(((pvalue<0.05 & diff> quantile(abs(diff),0.5))|diff>quantile(diff,0.9)) & isnan(CDB(:,4)),:);
            
            if strcmp(ref,'hg19')==1||strcmp(ref,'hg38')==1
              disp('* Annotate differential CDBs with CDB conservation. Write outputs:CDB1.txt,CDB2.txt.')
              disp('CDB.txt format:#chr start end LRI avgRI conserve_or_not consistent_or_differential')
                  if strcmp(ref,'hg19')==1
                     conserve=load('annotation/conserved_CDB_hg19.bed');
                  else
                     conserve=load('annotation/conserved_CDB_hg38.bed');
                  end
                  conserve=[conserve(:,1),(conserve(:,2)+conserve(:,3))/2];
                  CDB1=[CDB1,zeros(size(CDB1,1),2)];
                  CDB2=[CDB2,zeros(size(CDB2,1),2)];
                  for conchr=1:23
                  contmp=conserve(conserve(:,1)==conchr,2);
                  %conserve  or not  
                  chrpos1=find(CDB1(:,1)==conchr);
                  CDBcenter1=(CDB1(chrpos1,2)+CDB1(chrpos1,3))/2;
                  pos1=min(pdist2(contmp,CDBcenter1)) < 40000+resolution/2;
                  CDB1(chrpos1,6)=pos1; 

                  chrpos2=find(CDB2(:,1)==conchr);
                  CDBcenter2=(CDB2(chrpos2,2)+CDB2(chrpos2,3))/2;
                  pos2=min(pdist2(contmp,CDBcenter2))< 40000+resolution/2 ;
                  CDB2(chrpos2,6)=pos2; 
                  
                  
                  %differential 
                  chrpos3=find(CDB_A(:,1)==conchr);
                  CDBcenter3=(CDB_A(chrpos3,2)+CDB_A(chrpos3,3))/2;
                  chrpos4=find(CDB_B(:,1)==conchr);
                  CDBcenter4=(CDB_B(chrpos4,2)+CDB_B(chrpos4,3))/2;
                  
                  pos=min(pdist2(CDBcenter3,CDBcenter1))<= resolution/2+1;  
                  CDB1(chrpos1,7)=1-pos;
                  pos=min(pdist2(CDBcenter4,CDBcenter2))<= resolution/2+1; 
                  CDB2(chrpos2,7)=1-pos;         
             end  
              else 
              disp('CDB annotation step is only for hg19 and hg38!')
              disp('Write outputs without conservation annotation: CDB1.txt, CDB2.txt')      
              disp('CDB.txt format:#chr start end LRI avgRI consistent_or_differential')
                  CDB1=[CDB1,zeros(size(CDB1,1),1)];
                  CDB2=[CDB2,zeros(size(CDB2,1),1)];
                  for conchr=1:23
                  %differential 
                  chrpos1=find(CDB1(:,1)==conchr);
                  CDBcenter1=(CDB1(chrpos1,2)+CDB1(chrpos1,3))/2;
                  chrpos2=find(CDB2(:,1)==conchr);
                  CDBcenter2=(CDB2(chrpos2,2)+CDB2(chrpos2,3))/2;
                  chrpos3=find(CDB_A(:,1)==conchr);
                  CDBcenter3=(CDB_A(chrpos3,2)+CDB_A(chrpos3,3))/2;
                  chrpos4=find(CDB_B(:,1)==conchr);
                  CDBcenter4=(CDB_B(chrpos4,2)+CDB_B(chrpos4,3))/2;
                  
                  pos=min(pdist2(CDBcenter3,CDBcenter1))<= resolution/2+1;  
                  CDB1(chrpos1,6)=1-pos;
                  pos=min(pdist2(CDBcenter4,CDBcenter2))<= resolution/2+1; 
                  CDB2(chrpos2,6)=1-pos;              
                  end  
             end         
            dlmwrite([outdir '/CDB1.txt'],CDB1,'delimiter','\t','precision', 9);
            dlmwrite([outdir '/CDB2.txt'],CDB2,'delimiter','\t','precision', 9);
     end
end

end

function [allLRI,allscore,boundary] = chrpeaks(hicfile,resolution,N,normalize,mind,wd,wdsize,chr)
% Detect local maximum peaks on a single chromosome  
% Implemented by
% Fengling Chen,Guipeng Li
% Tsinghua University
% https://github.com/ChenFengling

% wdsize=6 is the best in most exmaples we test
%wdsize=6; 
im=read2dense(hicfile,N,resolution);
[imnew,gapidx]=KRnorm(im);
disp('find local maximum');
oridx=1:size(im,1);
oridx(gapidx)=[];
imnew(gapidx,:)=[];
imnew(:,gapidx)=[];

% calculate RIP 
tmpx = log(1+imnew);     
tmpn = size(tmpx,1);
tads = zeros(wdsize,tmpn);

for w = wd:(wd+wdsize-1);
    xx = (w+1):(tmpn-w-1);
   
    sumup = 0;
    for i=-w:-1
        for j=(i+1):0;
            sumup = sumup+tmpx(xx+i+ tmpn*( xx+j-1));
        end
    end
    
    sumdown = 0;
    for i=1:w
        for j=(i+1):(w+1);
            sumdown = sumdown+tmpx(xx+i+ tmpn*( xx+j-1));
        end
    end

    sumrect = 0;
    for i=-w+1:0
        for j=1:w+i;
            sumrect = sumrect+tmpx(xx+i+ tmpn*( xx+j-1));
        end
    end

    tadscores = (sumup+sumdown-sumrect)./(sumup+sumdown+sumrect);
    tadscores = [repmat(tadscores(1),w,1); tadscores'; repmat(tadscores(end),w+1,1)];
    tads(w-wd+1,:) = tadscores';
end

% normalize
%meanylower=zeros(wdsize,1);
%for i=1:wdsize
%[~,ylower] = envelope(tads(i,:));
%[~,ylower2] = envelope(ylower);
%meanylower(i)=nanmean(ylower2);
%end
%tads=tads./repmat(meanylower,1,tmpn);

% mean
tadscores = mean(tads);

% calculate envelope
[~,ylower] = envelope(tadscores);
[~,ylower2] = envelope(ylower);
peakscore=tadscores-ylower2;
LRI=peakscore;
% find peaks 
[~,TAD_boundaries] = findpeaks(tadscores,'MINPEAKDISTANCE',mind);

% rescue  neighbouring peak with  high score
[~,TAD_boundaries2] = findpeaks(tadscores);
highcut=sort(peakscore(TAD_boundaries),'descend');
highcut=highcut(round(tmpn/(500000/resolution)));
logical1=peakscore>highcut;
logical2=zeros(1,tmpn);
logical2(TAD_boundaries2)=1;
logical3=zeros(1,tmpn);
logical3(TAD_boundaries)=1;
rescue=find(logical1==1&logical2==1&logical3==0);
TAD_boundaries=[TAD_boundaries,rescue];
peakscore=peakscore(TAD_boundaries);
peakavgRI=tadscores(TAD_boundaries);
boundary=[TAD_boundaries;peakscore;peakavgRI];

% change to original  boundary
boundary(1,:)=oridx(boundary(1,:));

% remove repeat region near 
rRa = min(pdist2(gapidx',boundary(1,:)'));
rpflt = rRa>wd;
boundary=boundary(:,rpflt);
[~,peakidx]=sort(boundary(2,:),'descend');
boundary=boundary(:,peakidx);
boundary=boundary(:,~isnan(boundary(2,:)));
boundary=[repmat(chr,1,size(boundary,2));boundary(1,:)*resolution-resolution/2;boundary(1,:)*resolution+resolution/2;boundary(2:end,:)];
boundary=boundary';

% aRI
tmp=NaN(size(im,1),1);
tmp(oridx)=tadscores;
tadscores=tmp;
Ns=size(im,1)-1;
allscore=[repmat(chr,1,Ns);(1:Ns)*resolution-resolution/2;[(1:Ns-1)*resolution+resolution/2,N];tadscores(1:Ns,1)'];
allscore=allscore';

% LRI
tmp=NaN(size(im,1),1);
tmp(oridx)=LRI;
LRI=tmp;
allLRI=[repmat(chr,1,Ns);(1:Ns)*resolution-resolution/2;[(1:Ns-1)*resolution+resolution/2,N];LRI(1:Ns,1)'];
allLRI=allLRI';
end

function [aRI,aRI1,aRI2,LRI1,LRI2,resultA,resultB] = chrpeaks_diff(hicfile,resolution,N,mind,wd,wdsize,chr)
% step.1 calucate aRI score
    n1=length(hicfile{1});
    n2=length(hicfile{2});
	sampleA={};
	sampleB={};
	gapidxA={};
	gapidxB={};
	sumKR=[];
	aRI=[];
	mergedA=zeros(ceil(N/resolution),ceil(N/resolution),'single');
	mergedB=zeros(ceil(N/resolution),ceil(N/resolution),'single');
	% KR norm + depthnorm
	% sample A
	disp('step.1 calucate aRI score on each replicate')
	for i=1:n1
        out=read2dense([hicfile{1}{i},'/chr',num2str(chr),'.matrix'],N,resolution);
	    mergedA=mergedA+out;
	    [out,gapidx]=KRnorm(out);
	    sumKR(i)=sum(sum(out));
	    sampleA{i}=out/sumKR(i)*sumKR(1);
	    gapidxA{i}=gapidx;
    end 
	%sampleA=getNorm(sampleA,plt=1,fout ="sampleA_MA_before.png")
	for i=1:n1
		  aRI(:,i)=calcu_aRI(sampleA{i},gapidxA{i},wd,wdsize);
        end 

	% sample B
        for i=1:n2
	    out=read2dense([hicfile{2}{i},'/chr',num2str(chr),'.matrix'],N,resolution);
	    mergedB=mergedB+out;
	    [out,gapidx]=KRnorm(out);
	    sumKR(n1+i)=sum(sum(out));
	    sampleB{i}=out/sumKR(n1+i)*sumKR(1);
	    gapidxB{i}=gapidx;
        end 
	%sampleA=getNorm(sampleA,plt=1,fout ="sampleA_MA_before.png")
        for i=1:n2
		  aRI(:,n1+i)=calcu_aRI(sampleB{i},gapidxB{i},wd,wdsize);
        end 
        disp('step2. calculate local peaks in merged samples');
	% step2.calculate CDBs in merged samples
	[LRI1,aRI1,resultA]=chrpeaks(mergedA,resolution,N,'KR',mind,wd,wdsize,chr);
	[LRI2,aRI2,resultB]=chrpeaks(mergedB,resolution,N,'KR',mind,wd,wdsize,chr);
end


function [tadscores]=calcu_aRI(imnew,gapidx,wd,wdsize)
 N=size(imnew,1);
 oridx=1:N;
 oridx(gapidx)=[];
 imnew(gapidx,:)=[];
 imnew(:,gapidx)=[];

% calculate RIP 
tmpx = log(1+imnew);     
tmpn = size(tmpx,1);
tads = zeros(wdsize,tmpn);

for w = wd:(wd+wdsize-1);
    xx = (w+1):(tmpn-w-1);
   
    sumup = 0;
    for i=-w:-1
        for j=(i+1):0;
            sumup = sumup+tmpx(xx+i+ tmpn*( xx+j-1));
        end
    end
    
    sumdown = 0;
    for i=1:w
        for j=(i+1):(w+1);
            sumdown = sumdown+tmpx(xx+i+ tmpn*( xx+j-1));
        end
    end

    sumrect = 0;
    for i=-w+1:0
        for j=1:w+i;
            sumrect = sumrect+tmpx(xx+i+ tmpn*( xx+j-1));
        end
    end

    tadscores = (sumup+sumdown-sumrect)./(sumup+sumdown+sumrect);
    tadscores = [repmat(tadscores(1),w,1); tadscores'; repmat(tadscores(end),w+1,1)];
    tads(w-wd+1,:) = tadscores';
end

tadscores = mean(tads);

tmp=NaN(N,1);
tmp(oridx)=tadscores;
tadscores=tmp;
end


function [boundary,fig] = EScutoff(maxall,motiffile,outdir)
% Decide a boundary cut off by  CTCF notif enrichment
%
% boundary: boundary with  LRI score
% motiffile: CTCF sites
% resolution: resolution of HiC matrix
%
% Implemented by  Fengling Chen
% Tsinghua University
% https://github.com/ChenFengling/
motiffile=dlmread(motiffile);
[~,peakidx]=sort(maxall(:,4),'descend');
maxall=maxall(peakidx,:);
chr=unique(maxall(:,1),'stable');
allpos=zeros(size(maxall,1),1);
for i=1:length(chr)
    motif=motiffile(motiffile(:,1)==chr(i),:);
    chrpos=find(maxall(:,1)==chr(i));
    localmax=maxall(chrpos,:);
    pos=min(pdist2(motif(:,2),(localmax(:,2)+localmax(:,3))/2)) < 20000;
    allpos(chrpos)=pos; 
end

allpos=double(allpos);
pcoe=maxall(allpos==1,4)/sum(maxall(allpos==1,4));
ncoe=1/sum(allpos==0);
allpos(allpos==0)= -ncoe;
allpos(allpos==1)= pcoe;
score= cumsum(allpos);
[~,I]=max(score);
fig=figure('visible','off');
plot(score);
boundary=maxall(1:I,:);
end

function [im] = read2dense(hicfile,N,resolution)
  % read hicmatrix
   if size(hicfile,1)==1
   triple=dlmread(hicfile);
   else
	   triple=hicfile;
   end
   
   if size(triple,1)==size(triple,2)
    im=triple;
    N=size(im,1);
   else
    N=ceil(N/resolution);
    im = zeros(N,N,'single');
    im(triple(:,1)/resolution+1+N*(triple(:,2)/resolution))=triple(:,3);
    im = im+im';
    im = double(im - diag(diag(im)));
  end

end

function [imnew,gapidx]=KRnorm(im)     
% sumim=sum(im>0);
%cutoff = quantile(sumim(sumim>0),0.05);
%gapidx=find(sum(im)<cutoff);
  % 80% nearest region should have reads
   pos=find(sum(im>0)==0);
   pos2=find(sum(im>0)~=0);
   im2=im;
   im2(pos,:)=[];
   im2(:,pos)=[];
   n=size(im2,1);
   A=ones(41,n);
   for i=1:20
   A(21+i,1:(n-i))=diag(im2,-i);
   A(21-i,(i+1):n)=diag(im2,i);
   end
   gapidx=find(sum(A>0)<35);
   gapidx=[pos,pos2(gapidx)];



  if(sum(sum(round(im)~=im))==0)
  disp('Your input matrix is raw matrix, perform KR normalization');
  imnew = bnewt2(im);
  else 
  disp('Your input matrix is normalized matrix, skip KR normalization');
  imnew=im;
  end
end


function [upperenv lowerenv] = envelope(sig, method)
% Find upper and lower envelopes of a given signal
% The idea is from Envelope1.1 by Lei Wang, but here it works well when the signal contains
% successive equal samples and also includes first and last samples of the signal in the envelopes.
% inputs:
%   sig: vector of input signal
%   method: method of interpolation (defined as in interp1)
% outputs:
%   upperenv: upper envelope of the input signal
%   lowerenv: lower envelope of the input signal
if nargin == 1 
    method = 'linear';
end
upperind = find(diff(sign(diff(sig))) < 0) + 1;
lowerind = find(diff(sign(diff(sig))) > 0) + 1;
f = 1;
l = length(sig);
try
    upperind = [f upperind l];
    lowerind = [f lowerind l];
catch 
    upperind = [f; upperind; l];
    lowerind = [f; lowerind; l];
end
xi = f : l;
upperenv = interp1(upperind, sig(upperind), xi, method, 'extrap');
lowerenv = interp1(lowerind, sig(lowerind), xi, method, 'extrap');
end


function [Anew,KRnorm,res] = bnewt2(A0,tol,x0,delta,Delta,fl)
n0 = size(A0,1);
KR0 = sum(A0)>100;
A = A0(KR0,KR0);

% BNEWT A balancing algorithm for symmetric matrices
%
% X = BNEWT(A) attempts to find a vector X such that
% diag(X)*A*diag(X) is close to doubly stochastic. A must
% be symmetric and nonnegative.
%
% X0: initial guess. TOL: error tolerance.
% delta/Delta: how close/far balancing vectors can get
% to/from the edge of the positive cone.
% We use a relative measure on the size of elements.
% FL: intermediate convergence statistics on/off.
% RES: residual error, measured by norm(diag(x)*A*x - e).
% Initialise
n = size(A,1); e = ones(n,1); res=[];
if nargin < 6, fl = 0; end
if nargin < 5, Delta = 2; end
if nargin < 4, delta = 0.05; end
if nargin < 3, x0 = e; end
if nargin < 2, tol = 1e-5; end
% Inner stopping criterion parameters.
g=0.9; etamax = 0.05;
eta = etamax; stop_tol = tol*.5;
x = x0; rt = tol^2; v = x.*(A*x); rk = 1 - v;

rho_km1 = rk'*rk; rout = rho_km1; rold = rout;
MVP = 0; % We’ll count matrix vector products.
i = 0; % Outer iteration count.
if fl == 1, fprintf('it in. it res\n'), end
while rout > rt % Outer iteration
    i = i + 1; k = 0; y = e;
    innertol = max([eta^2*rout,rt]);
    while rho_km1 > innertol %Inner iteration by CG
        k = k + 1;
        if k == 1
            Z = rk./v; p=Z; rho_km1 = rk'*Z;
        else
            beta=rho_km1/rho_km2;
            p=Z + beta*p;
        end
        % Update search direction efficiently.
        w = x.*(A*(x.*p)) + v.*p;
        alpha = rho_km1/(p'*w);
        ap = alpha*p;
        % Test distance to boundary of cone.
        ynew = y + ap;
        if min(ynew) <= delta
            if delta == 0, break, end
            ind = find(ap < 0);
            gamma = min((delta - y(ind))./ap(ind));
            y = y + gamma*ap;
            break
        end
        if max(ynew) >= Delta
            ind = find(ynew > Delta);
            gamma = min((Delta-y(ind))./ap(ind));
            y = y + gamma*ap;
            break
        end
        y = ynew;
        rk = rk - alpha*w; rho_km2 = rho_km1;
        Z = rk./v; rho_km1 = rk'*Z;
    end
    x = x.*y; v = x.*(A*x);
    rk = 1 - v; rho_km1 = rk'*rk; rout = rho_km1;
    MVP = MVP + k + 1;

    % Update inner iteration stopping criterion.
    rat = rout/rold; rold = rout; res_norm = sqrt(rout);
    eta_o = eta; eta = g*rat;
    if g*eta_o^2 > 0.1
        eta = max([eta,g*eta_o^2]);
    end
    eta = max([min([eta,etamax]),stop_tol/res_norm]);
    if fl == 1
        fprintf('%3d %6d %.3e \n',i,k, r_norm);
        res=[res; r_norm];
    end
end


KRnorm = ones(n0,1);
KRnorm(KR0) = x;
s = sum(sum(A0))/n0;
KRnorm = KRnorm*sqrt(s);
Anew = A0.*(KRnorm*KRnorm');
%fprintf('Matrix-vector products = %6d\n', MVP)
end


