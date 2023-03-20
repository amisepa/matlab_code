function flag = bivariate_outliers(X)

% routine that find the bivariate outliers using orthogonal projection and
% box plot rule 

% find the centre of the data cloud using mid-covariance determinant
n   = size(X,1);
result = mcdcov(X,'cor',1,'plots',0,'h',floor((n+size(X,2)*2+1)/2));
center = result.center;

% orthogonal projection to the lines joining the center
% followed by outlier detection using box plot rule

gval = sqrt(chi2inv(0.975,2)); % in fact depends on size(X,2) but here always = 2
for i=1:n % for each row
    dis = NaN(n,1);
    B   = (X(i,:)-center)';
    BB  = B.^2;
    bot = sum(BB);
    if bot~=0
        for j=1:n
            A = (X(j,:)-center)';
            dis(j)= norm(A'*B/bot.*B);
        end
        % IQR rule
        [ql,qu]=idealf(dis);
        record{i} = (dis > median(dis)+gval.*(qu-ql)) ; % + (dis < median(dis)-gval.*(qu-ql));
    end
end

try
    flag = nan(n,1);
    flag = sum(cell2mat(record),2); % if any point is flagged
    
catch ME  % this can happen to have an empty cell so loop
    flag = nan(n,size(record,2));
    index = 1;
    for s=1:size(record,2)
        if ~isempty(record{s})
            flag(:,index) = record{s};
            index = index+1;
        end
    end
    flag(:,index:end) = [];
    flag = sum(flag,2);
end

%% Subfunctions

%CLASSSVD performs the singular value decomposition of a matrix with more
% rows than columns (uses svd.m)
%
% Required input: 
%          x : data matrix of size n by p where n>p
%
% This function is part of LIBRA: the Matlab Library for Robust Analysis,
% available at: 
%              http://wis.kuleuven.be/stat/robust.html
%
% Written by Sabine Verboven, Mia Hubert
% Last Update: 06/07/2004

function [P,T,L,r,centerX,cX]=classSVD(x)
[n,p] = size(x); 

if n==1
    error('The sample size is 1. No SVD can be performed.')
end
cX=mean(x);
centerX=x-repmat(cX,n,1); 
[U,S,loadings]=svd(centerX./sqrt(n-1),0); 
eigenvalues=diag(S).^2;
tol = max([n p])*eigenvalues(1)*eps;
r= sum(eigenvalues>tol);
L=eigenvalues(1:r);
P=loadings(:,1:r);
T=centerX*P; 

% Compute the ideal fourths for data in x
% The estimate of the interquartile range is:
% IQR=qu-ql;
% Adapted from Rand Wilcox's idealf R function, described in
% Rand Wilcox, Introduction to Robust Estimation & Hypothesis Testing, 3rd
% edition, Academic Press, Elsevier, 2012

% Cyril Pernet & Guillaume Rousselet, v1 - September 2012
% ---------------------------------------------------
%  Copyright (C) Corr_toolbox 2012

function [ql,qu]=idealf(x)

j=floor(length(x)/4 + 5/12);
y=sort(x);
g=(length(x)/4)-j+(5/12);
ql=(1-g).*y(j)+g.*y(j+1);
k=length(x)-j+1;
qu=(1-g).*y(k)+g.*y(k-1);