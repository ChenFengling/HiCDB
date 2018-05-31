function [fig] = visHiCDB(hicfile,CDBfile,resolution,chr,startloc,endloc,varargin)
%   Visualization of CDBs or differential CDBs on Hi-C maps
%   Implemented by
%   Fengling Chen
%   May 2017
%   Tsinghua University
%   https://github.com/ChenFengling/HiCDB
%   Any suggestions and remarks might be addressed to Fengling Chen:cfl15@mails.tsinghua.edu.cn
%
%visHiCDB uses Hi-C contact matrix and CDBs as input and outputs figures
%    of CDBs or differential CDBs on Hi-C maps.
%
%    fig=HiCDB(hicfile,CDBfile,resolution,chr,startloc,endloc) will use singlemap as a default option.
%      hicfile is the file name of  the intra-chromosome matrix.
%      The intra-chromosome matrix could be in a dense (a NxN matrix) or sparse (a Kx3 table,Rao et al.) format.
%      Show CDBs on one sample,set hicfile as 'SAMPLE_File'. CDBfile is the file name of the CDB files.
%     The CDB file should be formated as the output file of HiCDB. Resolution is the resolution of Hi-C matrix.
%     chr,startloc and endloc represent the observation locus. They are all
%     numeric. This function outputs a pdf figure showing CDBs on a Hi-C
%     map on desired locus. KR normalization is performed if the matrix is raw. The
%     running time will be saved a lot if input Hi-C matrix are already normalized.
%    
%    fig=visHiCDB(hicfile,CDBdir,resolution,chr,startloc,endloc,'option','comparemap'), where
%    hicfile shoule be a cell arrary contaioning two Hi-C matrix like {'SAMPLE1_FILE','SAMPLE2_FILE'}.
%    CDBfile should be a cell arrary storing file name of the two CDB files of two samples like {'SAMPLE1_CDB','SAMPLE2_CDB'}.
%    This function outputs a pdf figure showing different kinds of CDBs
%    between these two samples.
%
%    [ ... ] = 	visHiCDB(..., 'PARAM1',val1, 'PARAM2',val2, ...) specifies
%     optional parameter name/value pairs to show the CDBs or show
%     differential CDBs on Hi-C maps and/or to save the figure in desired directory. 
%     Parameters are : 
%    
%
%    'option'     -  Choices are :
%         'singlemap'    - CDBs on one Hi-C map. Default value is 'singlemap'
%         'comparemap'   - compare CDBs on two Hi-C maps. 
%    
%
%     'outdir'   -The output direction. Default will be the directory of
%     the first sample.        
%     
%    Examples:
%
%       1. Show CDB on single Hi-C map
%       fig=visHiCDB({'sample1/chr18.matrix'},{'CDB1.txt'},40000,18,25000000,31150000);
%       2. Show differential CDBs
%       fig=visHiCDB({'sample1/chr18.matrix','sample2/chr18.matrix'},{'CDB1.txt','CDB2.txt'},40000,18,25000000,31150000,'option','comparemap');
% 
%   See also HiCDB
%

%  input
if length(hicfile)==1||length(hicfile)==2
    pathstr = fileparts(hicfile{1});
else
   error('invalid cell number for hicfile.') 
end

if length(hicfile)~=length(CDBfile)
   error('unequal hicfile and CDB file length.')
end

p = inputParser;
p.FunctionName = 'visHiCDB';

defaultOption = 'singlemap';
validOption = {'singlemap','comparemap'};
checkOption = @(x) any(validatestring(x,validOption));

defaultOutdir = pathstr;

defaultNormalized='raw';
validNormalized = {'raw','KR'};
checkNormalized = @(x) any(validatestring(x,validNormalized));


addRequired(p,'hicfile');
addRequired(p,'CDBfile');
addRequired(p,'resolution',@isnumeric);
addRequired(p,'chr',@isnumeric);
addRequired(p,'startloc',@isnumeric);
addRequired(p,'endloc',@isnumeric);
addParamValue(p,'option',defaultOption,checkOption)
addParamValue(p,'outdir',defaultOutdir)
addParamValue(p,'normalized',defaultNormalized,checkNormalized)

p.parse(hicfile,CDBfile,resolution,chr,startloc,endloc,varargin{:})

option = p.Results.option
outdir = p.Results.outdir
normalized = p.Results.normalized

validateattributes(option,{'char'},{},'visHiCDB','option')
validateattributes(outdir,{'char'},{},'visHiCDB','outdir')
validateattributes(normalized,{'char'},{},'visHiCDB','normalized')

tX=round(startloc/resolution+0.5) :round(endloc/resolution);
%%
if strcmp(option,'singlemap')==1
    triple=dlmread(hicfile{1});
    if size(triple,1)==size(triple,2)
        im=triple;
        N=size(im,1);
    else
        N = max(triple(:,1));
        N = max([triple(:,2);N])/resolution+1;
        im = zeros(N,N,'single');
        im(triple(:,1)/resolution+1+N*(triple(:,2)/resolution))=triple(:,3);
        im = im+im';
        im = double(im - diag(diag(im)));
    end
    if strcmp(normalized,'KR')~=1
        imnew = bnewt2(im);
    else 
        imnew=im;
    end

    CDB=dlmread(CDBfile{1}); 
    CDB=CDB(CDB(:,1)==chr,:); 

    fig=figure('visible','off');
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 7 6.8])
    imagesc(log(1+imnew(tX,tX)));
  
    if size(CDB,2)==5
    boundary=(CDB(:,2)+CDB(:,3))/resolution/2;
    hold on
    flt = boundary>=tX(1) & boundary <=tX(end);
    tTAD = boundary-tX(1)+1;
    scatter(tTAD(flt)+0.5,tTAD(flt)+0.5,600,'w.');
    else
        if(length(CDB(:,6)==1)~=0)
        boundary=CDB(CDB(:,6)==1,:);
        boundary=(boundary(:,2)+boundary(:,3))/resolution/2;
        hold on
        flt = boundary>=tX(1) & boundary <=tX(end);
        tTAD = boundary-tX(1)+1;
        scatter(tTAD(flt)+0.5,tTAD(flt)+0.5,600,[0,32,96]/255,'.');   
        end
        if(length(CDB(:,6)==0)~=0)
        boundary=CDB(CDB(:,6)==0,:);
        boundary=(boundary(:,2)+boundary(:,3))/resolution/2;
        hold on
        flt = boundary>=tX(1) & boundary <=tX(end);
        tTAD = boundary-tX(1)+1;
        scatter(tTAD(flt)+0.5,tTAD(flt)+0.5,600,[91,155,213]/255,'.');   
        end     
    
    end
    saveas(fig,[outdir  '/' num2str(chr) '_' num2str(startloc) '_' num2str(endloc) '_HiCmap.fig']);
    print(fig,'-dpdf',[outdir  '/' num2str(chr) '_' num2str(startloc) '_' num2str(endloc) '_HiCmap.pdf']);
else
%%
   disp('*Processing sample 1*')
    triple=dlmread(hicfile{1});
    if size(triple,1)==size(triple,2)
        im=triple;
        N=size(im,1);
    else
        N = max(triple(:,1));
        N = max([triple(:,2);N])/resolution+1;
        im = zeros(N,N,'single');
        im(triple(:,1)/resolution+1+N*(triple(:,2)/resolution))=triple(:,3);
        im = im+im';
        im = double(im - diag(diag(im)));
    end
    if strcmp(normalized,'KR')~=1
       imnew = bnewt2(im);
    else 
        imnew=im;
    end

     depth=sum(sum(imnew)); 
    CDB=dlmread(CDBfile{1}); 
    CDB=CDB(CDB(:,1)==chr,:); 
    
    fig=figure('visible','off');
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 9 3.9])
    subplot(1,2,1)
    imagesc(log(1+imnew(tX,tX)));
    
    if size(CDB,2)==5
    boundary=(CDB(:,2)+CDB(:,3))/resolution/2;
    hold on
    flt = boundary>=tX(1) & boundary <=tX(end);
    tTAD = boundary-tX(1)+1;
    scatter(tTAD(flt)+0.5,tTAD(flt)+0.5,200,'w.');
    elseif size(CDB,2)==6
        if(length(CDB(:,6)==1)~=0)
        boundary=CDB(CDB(:,6)==1,:);
        boundary=(boundary(:,2)+boundary(:,3))/resolution/2;
        hold on
        flt = boundary>=tX(1) & boundary <=tX(end);
        tTAD = boundary-tX(1)+1;
        scatter(tTAD(flt)+0.5,tTAD(flt)+0.5,200,[91,155,213]/255,'.');   
        end
        if(length(CDB(:,6)==0)~=0)
        boundary=CDB(CDB(:,6)==0,:);
        boundary=(boundary(:,2)+boundary(:,3))/resolution/2;
        hold on
        flt = boundary>=tX(1) & boundary <=tX(end);
        tTAD = boundary-tX(1)+1;
        scatter(tTAD(flt)+0.5,tTAD(flt)+0.5,200,[255,255,104]/255,'.');   
        end       
    elseif size(CDB,2)==7
        if(length(CDB(:,6)==1 & CDB(:,7)==1)~=0)
        boundary=CDB(CDB(:,6)==1 & CDB(:,7)==1,:);
        boundary=(boundary(:,2)+boundary(:,3))/resolution/2;
        hold on
        flt = boundary>=tX(1) & boundary <=tX(end);
        tTAD = boundary-tX(1)+1;
        scatter(tTAD(flt)+0.5,tTAD(flt)+0.5,200,[0,32,96]/255,'.');   
        end
        if(length(CDB(:,6)==0 & CDB(:,7)==1)~=0)
        boundary=CDB(CDB(:,6)==0 & CDB(:,7)==1,:);
        boundary=(boundary(:,2)+boundary(:,3))/resolution/2;
        hold on
        flt = boundary>=tX(1) & boundary <=tX(end);
        tTAD = boundary-tX(1)+1;
        scatter(tTAD(flt)+0.5,tTAD(flt)+0.5,200,[91,155,213]/255,'.');   
        end  
        if(length(CDB(:,6)==1 & CDB(:,7)==0)~=0)
        boundary=CDB(CDB(:,6)==1 & CDB(:,7)==0,:);
        boundary=(boundary(:,2)+boundary(:,3))/resolution/2;
        hold on
        flt = boundary>=tX(1) & boundary <=tX(end);
        tTAD = boundary-tX(1)+1;
        scatter(tTAD(flt)+0.5,tTAD(flt)+0.5,200,[237,125,49]/255,'.');   
        end
        if(length(CDB(:,6)==0 & CDB(:,7)==0)~=0)
        boundary=CDB(CDB(:,6)==0 & CDB(:,7)==0,:);
        boundary=(boundary(:,2)+boundary(:,3))/resolution/2;
        hold on
        flt = boundary>=tX(1) & boundary <=tX(end);
        tTAD = boundary-tX(1)+1;
        scatter(tTAD(flt)+0.5,tTAD(flt)+0.5,200,[255,255,104]/255,'.');   
        end  
    end
    
    disp('*Processing sample 2*')
    triple=dlmread(hicfile{2});
    if size(triple,1)==size(triple,2)
        im=triple;
        N=size(im,1);
    else
        N = max(triple(:,1));
        N = max([triple(:,2);N])/resolution+1;
        im = zeros(N,N,'single');
        im(triple(:,1)/resolution+1+N*(triple(:,2)/resolution))=triple(:,3);
        im = im+im';
        im = double(im - diag(diag(im)));
    end

    if strcmp(normalized,'KR')~=1
       imnew = bnewt2(im);
    else 
        imnew=im;
    end
   
    imnew=imnew/sum(sum(imnew))*depth;
    CDB=dlmread(CDBfile{2}); 
    CDB=CDB(CDB(:,1)==chr,:); 

    subplot(1,2,2)
    imagesc(log(1+imnew(tX,tX)));
    
    if size(CDB,2)==5
    boundary=(CDB(:,2)+CDB(:,3))/resolution/2;
    hold on
    flt = boundary>=tX(1) & boundary <=tX(end);
    tTAD = boundary-tX(1)+1;
    scatter(tTAD(flt)+0.5,tTAD(flt)+0.5,200,'w.');
    elseif size(CDB,2)==6
        if(length(CDB(:,6)==1)~=0)
        boundary=CDB(CDB(:,6)==1,:);
        boundary=(boundary(:,2)+boundary(:,3))/resolution/2;
        hold on
        flt = boundary>=tX(1) & boundary <=tX(end);
        tTAD = boundary-tX(1)+1;
        scatter(tTAD(flt)+0.5,tTAD(flt)+0.5,200,[91,155,213]/255,'.');   
        end
        if(length(CDB(:,6)==0)~=0)
        boundary=CDB(CDB(:,6)==0,:);
        boundary=(boundary(:,2)+boundary(:,3))/resolution/2;
        hold on
        flt = boundary>=tX(1) & boundary <=tX(end);
        tTAD = boundary-tX(1)+1;
        scatter(tTAD(flt)+0.5,tTAD(flt)+0.5,200,[153,255,153]/255,'.');   
        end       
    elseif size(CDB,2)==7
        if(length(CDB(:,6)==1 & CDB(:,7)==1)~=0)
        boundary=CDB(CDB(:,6)==1 & CDB(:,7)==1,:);
        boundary=(boundary(:,2)+boundary(:,3))/resolution/2;
        hold on
        flt = boundary>=tX(1) & boundary <=tX(end);
        tTAD = boundary-tX(1)+1;
        scatter(tTAD(flt)+0.5,tTAD(flt)+0.5,200,[0,32,96]/255,'.');   
        end
        if(length(CDB(:,6)==0 & CDB(:,7)==1)~=0)
        boundary=CDB(CDB(:,6)==0 & CDB(:,7)==1,:);
        boundary=(boundary(:,2)+boundary(:,3))/resolution/2;
        hold on
        flt = boundary>=tX(1) & boundary <=tX(end);
        tTAD = boundary-tX(1)+1;
        scatter(tTAD(flt)+0.5,tTAD(flt)+0.5,200,[91,155,213]/255,'.');   
        end  
        if(length(CDB(:,6)==1 & CDB(:,7)==0)~=0)
        boundary=CDB(CDB(:,6)==1 & CDB(:,7)==0,:);
        boundary=(boundary(:,2)+boundary(:,3))/resolution/2;
        hold on
        flt = boundary>=tX(1) & boundary <=tX(end);
        tTAD = boundary-tX(1)+1;
        scatter(tTAD(flt)+0.5,tTAD(flt)+0.5,200,[112,48,160]/255,'.');   
        end
        if(length(CDB(:,6)==0 & CDB(:,7)==0)~=0)
        boundary=CDB(CDB(:,6)==0 & CDB(:,7)==0,:);
        boundary=(boundary(:,2)+boundary(:,3))/resolution/2;
        hold on
        flt = boundary>=tX(1) & boundary <=tX(end);
        tTAD = boundary-tX(1)+1;
        scatter(tTAD(flt)+0.5,tTAD(flt)+0.5,200,[153,255,153]/255,'.');   
        end  
    end

    saveas(fig,[outdir '/' num2str(chr) '_' num2str(startloc) '_' num2str(endloc) '_HiCmap_compare.fig']);
    print(fig,'-dpdf',[outdir '/' num2str(chr) '_' num2str(startloc) '_' num2str(endloc) '_HiCmap_compare.pdf']);
end

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

