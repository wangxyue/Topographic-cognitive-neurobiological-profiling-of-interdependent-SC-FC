function [ICC] = Liao_icc(x)

%% Calculate the intra-class correlation with the formula:
% ICC=(F-1)/(F+(k-1)), where F=MSb/MSw, k=n
% Note: each unit in the matrix is a single measurement of each rater
%       formula corresponds to ICC(1,1)
%
% Input: x, dim: n-by-k    - includes data of different targets (subjects) within repeat
%        measurements, each row means each target assessed by k raters
%        e.g. (a1,a2,a3,...,ak) % target a
%             (b1,b2,b3,...,bk) % target b
%             (..,..,...     )
%             (n1,n2,n3,...,nk) % target n
%        n    - number of targets sampled from a large population
%        k    - number of raters sampled from a population
%        numS - number of all samples in total
%
% Output: ICC - intra-class correltation ICC(1,1)

%
% Xuhong Liao 2012/03/26
% Modified by Xuhong Liao, 2012/04/05
% Ref: Shrout1979, IPN_icc.m of Zuo X. N.

[n,k]=size(x);                       % n, number of targets; k, number of raters
numS = n*k;                          % total number of measurements

Vb = n-1;                            % freedom degree of between targets;
Vw = numS-n;                         % freedom degree of within targets =n*(k-1);

ave = sum(x(:))/numS;                % average of all groups

aveT = sum(x,2)/k;                     % average of each target across k raters;

SSb = sum((aveT-ave).^2)*k;            % variation between different targets;

%SSw = sum(sum((Indata-ave).^2,1),2)-SSb;      % variation within each target;
temp = (x-repmat(aveT,1,k)).^2;
SSw = sum(temp(:));                        % variations within each target;

MSb = SSb/Vb;
MSw = SSw/Vw;
F = MSb/MSw;

ICC = (F-1) / (F+(k-1));
end

