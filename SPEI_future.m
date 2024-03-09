%Calculate SPEI from 1-D precipitation and evapotranspiration arrays
%Diogo S. A. Araújo
%Florida Institute of Technology
%Inputs:
%   1) Pr: Vector with precipitation data series.
%   2) Evap: Vector with potential evapotranspiration data series.
%   3) Scale: scalar which represents the frequency, e.g. 1, 3, 6, 9 or 12
%   months
%   4) Distribution: select a distribution to fit Pr-Evap data -
%   'GeneralizedExtremeValue', 'LogLogistic'
%   5) parm: Distribution parameters from fitted distribution, e.g. output
%   from SPEI function.
%   Note: both vectors must have same size, units and refer to the same period.
%Outputs:
%   1) spei: Vector with SPEI values from desired time scale.

function [spei] = SPEI_future(pr,evap,scale,dist,parm)

%Difference of precipitation and evapotranspiration. Check the sign for
%evap
Data = pr - evap;

% check if data lenth is multiple of 12
if(rem(length(Data),12)~=0)
    error('data record is not multiple of 12 months....')
end

% Data setting to scaled dataset
A1=[];
for is=1:scale, A1=[A1,Data(is:length(Data)-scale+is)];end
XS=sum(A1,2);

%Pre-locating space for variable U
%spei = zeros(size(pr,1),size(pr,2),size(pr,3)-scale+1);

%Gamma parameter estimation and tranform    
if strcmp(dist,'GeneralizedExtremeValue')
    for stp =1:12
        tind=stp:12:length(XS); %Creates a time index
        Xn=XS(tind); %Gets value from Xn
        parm0 = parm; %Get the parameters for GEV historic distribution;
        Dist_xs= cdf(dist,Xn,parm0(stp,1),parm0(stp,2),parm0(stp,3)); %calculates accumulated probabilty from gamma distribution;
        spei(tind)= norminv(Dist_xs); %returns the SPI value from the inverse of accumulated probability;
        %spei(spei < -5) = -5; %SPEI limit values
        %spei(spei > 5) = 5; %SPEI limit values
    end

%Gamma parameter estimation and tranform    
elseif strcmp(dist,'LogLogistic')
    
    for stp =1:12
        tind=stp:12:length(XS); %Creates a time index
        Xn=XS(tind); %Gets value from Xn
        parm0 = squeeze(parm(stp,:)); %Gets the initial value for the parameters
        valid_answer = find((Xn >= parm0(3) & parm0(2)<0)|(Xn <= parm0(3) & parm0(2)>0)); %Search for valid answers interval
        %Calculation of probability of excedance 
        for t = 1:length(Xn)
           if ismember(t,valid_answer)
               Dist_xs(t) = 1-(1+((Xn(t)-parm0(3))/parm0(1))^(-parm0(2)))^(-1);
           else
               Dist_xs(t) = NaN;
           end
        end      
        spei(tind)=norminv(Dist_xs); %returns the SPEI value from the inverse of accumulated probability;
        clear Dist_xs; %Avoid missuse of data from previous loops
    end
    
end
end

