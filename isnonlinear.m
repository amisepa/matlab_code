%% Nonlinearity test for univariate time series. 
% Includes the following nonliearity tests:
%   Keenan,D. M.(1985). A Tukey nonadditivity-type test for time series
%       nonlinearity, Biometrika 72, 39-44.
%   Mohammadi S.(2019). Neural network for univariate and multivariate
%       nonlinearity tests. Stat Anal Data Min: The ASA DataSci Journal.
%       13:50-70. https://doi.org/10.1002/sam.11441
%   Ramsey,J. B.(1969). Tests for specification errors in classical linear
%       least squares regression analysis, J. R. Stat. Soc B 31, 350-371.
%   Terasvirta,T., C. Lin, and C. W. J. Granger(1993), Power of the neural
%       network linearity test, J. Time Ser. Anal. 14, 209-220.
%   Tsay, R. S.(1986). Nonlinearity tests for times series, Biometrika 73,
%       4, 61-466.
%
% Input:
%   - y is a univariate time series in the column vector form.
% 
% Output:
%   - Probability value of the tests. Reject linearity(H0) if pval is
%           less than 5% (or any predetermined level of significance).
% 
% This Matlab code has been used in S. Mohammadi(2019).McLeod-Li and Arch test for
% nonlinearity can be done by MATLAB built-in functions, namely lbqtest and
% archtest. Please pay attention to that in the nonlinearity test by lbtest
% you should use squares of residual series of the time series.
% 
% 
% Copyright, Shapour Mohammadi 2020.shmohmad@ut.ac.ir

function [PvalRamsey,PvalKeenan,PvalTras,PvalTsay] = isnonlinear(y)

n = length(y);
nlags = 5;      %number of lags

% aiccriter = zeros(5,1);

% Find optimal lag
% for lag = 1:nlags
%     xlagmat = lagmatrix(y,1:lag);
%     xlagmat = [ones(n-lag,1) xlagmat(lag+1:end,:)];
%     [mohstats] = moholsforann(y(lag+1:end,1),xlagmat);
%     resmdl1 = mohstats{1,1};
%     aiccriter(lag,1) = log((resmdl1'*resmdl1)/(n-lag))+2*lag/(n-lag);
% end
% optlag = find(aiccriter==min(aiccriter));
% X = lagmatrix(y,1:optlag);

% Try bypass because lagmatrix not available
% for i = 1:nlags
%    X(:,i) = y(nlags-i+1:n-i);
% end

index0 = isnan(X(:,end));
index = find(index0==1);
maxindex = max(index)+1;
X = X(maxindex:end,:);
y = y(maxindex:end,:);
[n, K] = size(X);

[mohstats] = moholsforann(y,X);
resid = mohstats{1,1};
yfit = mohstats{2,1};

%Nonlinearity Tests:

% Ramsey(1969) Test
SSR0 = sum(resid.^2);
[mohstatsRamsey] = moholsforann(resid,[X yfit.^2 yfit.^3]);
residRamsey=mohstatsRamsey{1,1};
SSR1=sum(residRamsey.^2);
FRamsey=((SSR0-SSR1)/(2+K+1))/(SSR1/(n-2*K-2-1));
PvalRamsey=1-fcdf(FRamsey,2+K+1,n-2*K-2-1);


%Keenan Test
[mohstatsKeenan]=moholsforann(yfit.^2,X);
residKeenan=mohstatsKeenan{1,1};
[mohstatsKeenan2b]=moholsforann(resid,residKeenan);
betasKeenan=mohstatsKeenan2b{3,1};
etahato=betasKeenan(2,1);
etahat=etahato*(sum(residKeenan.^2))^0.5;
FKeenan=(etahat^2)*(n-K-2)/(sum(resid.^2)-etahat^2);
PvalKeenan=1-fcdf(FKeenan,1,n-K-2);



% Tsay Test(1986)
CRterms=[];
for i=1:K
    for j=i:K
        CRterm=X(:,i).*X(:,j);
        CRterms=[CRterms CRterm];
    end
end
[mohstatsTsay1]=moholsforann(CRterms,X);
residTsay1=mohstatsTsay1{1,1};

[~,m]=size(CRterms);
[mohstatsTsay2]=moholsforann(resid,residTsay1);
residTsay2=mohstatsTsay2{1,1};

FTsay=((resid'*residTsay1)*(residTsay1'*residTsay1)^(-1)*(resid'*residTsay1)'...
    /m)/((residTsay2'*residTsay2)/(n-m-K-1));
PvalTsay=1-fcdf(FTsay,m,n-m-K-1);




% Trasvirta et. al. (1993) test V23(Volterra 2nd and 3rd degree series)
CRtermsTras2=[];
for i=1:K
    for j=i:K
        CRtermTras2=X(:,i).*X(:,j);
        CRtermsTras2=[CRtermsTras2 CRtermTras2];
    end
end
CRtermsTras3=[];
for i=1:K
    for j=i:K
        for l=j:K
            CRtermTras3=X(:,i).*X(:,j).*X(:,l);
            CRtermsTras3=[CRtermsTras3 CRtermTras3];
        end
    end
end
SSR0Tras=sum(resid.^2);
[~,mTras]=size([CRtermsTras2 CRtermsTras3]);
[mohstatsTras]=moholsforann(resid,[X CRtermsTras2 CRtermsTras3]);
residTras=mohstatsTras{1,1};
SSR1Tras=sum(residTras.^2);

FTras=((SSR0Tras-SSR1Tras)/(mTras))/(SSR1Tras/(n-mTras-K-1));
PvalTras=1-fcdf(FTras,mTras,n-mTras-K-1);



% fast function for ordinary least squares estimation.
function [mohstats]=moholsforann(y,x)
x=[ones(length(y),1), x];

betaols=x\y;
yfit=x*betaols;
e=y-yfit;

mohstats{1,1}=e;
mohstats{2,1}=yfit;
mohstats{3,1}=betaols;



