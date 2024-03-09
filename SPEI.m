%Calculate SPEI from 1-D precipitation and potential evapotranspiration arrays
%Diogo S. A. Araújo
%Florida Institute of Technology
%Inputs:
%   1) Pr: Vector with precipitation data series.
%   2) Evap: Vector with potential evapotranspiration data series.
%   3) Scale: scalar which represents the frequency, e.g. 1, 3, 6, 9 or 12
%   months
%   4) Dist: select a distribution to fit Pr-Evap data -
%   'GeneralizedExtremeValue', 'LogLogistic' or 'NonParametric'
%   Note: both vectors must have same size, units and refer to the same period.
%Outputs:
%   1) spei: Vector with SPEI values from desired time scale.
%   2) dist_parm: Parameters of fitted distribution.

function [spei, dist_parm] = SPEI(pr,evap,scale,dist)

%Difference of precipitation and potential evapotranspiration. Check the sign for
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
        %parm = fitdist(Xn,dist); %calculating the parameters for a GEV distribution;
        parm = gevfit_lmom(Xn); %calculating the parameters for a GEV distribution;
        dist_parm(stp,1) = parm(1); %saving the parameters in dist_parm array;
        dist_parm(stp,2) = parm(2); %saving the parameters in dist_parm array;
        dist_parm(stp,3) = parm(3); %saving the parameters in dist_parm array;
        Dist_xs= cdf(dist,Xn,parm(1),parm(2),parm(3)); %calculates accumulated probabilty from gamma distribution;
        spei(tind)= norminv(Dist_xs); %returns the SPI value from the inverse of accumulated probability;
    end
    
    %Gamma parameter estimation and tranform
elseif strcmp(dist,'LogLogistic')
    
    LogLogistic3_pdf = @(x,alpha,beta,gamma) (beta/alpha)*(((x-gamma)/alpha).^(beta-1)).*((1+((x-gamma)/alpha).^beta).^-2);
    LogLogistic3_cdf = @(x,alpha,beta,gamma) (1+(alpha/(x-gamma)).^beta).^-1;
    
    o = statset('mlecustom');
    o.FunValCheck = 'off';
    
    for stp =1:12
        tind=stp:12:length(XS); %Creates a time index
        Xn=XS(tind); %Gets value from Xn
        parm0 = par_3_LogLogistic(Xn); %Calculate the initial value for the parameters
        %parm = mle(Xn,'pdf',LogLogistic3_pdf,'cdf',LogLogistic3_cdf,'start',parm0,'Optimfun','fminsearch','options',o); %calculating the parameters for a log-logistic distribution;
        dist_parm(stp,1) = parm0(1); %saving the parameters in dist_parm array;
        dist_parm(stp,2) = parm0(2); %saving the parameters in dist_parm array;
        dist_parm(stp,3) = parm0(3); %saving the parameters in dist_parm array;
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

    %Gamma parameter estimation and tranform
elseif strcmp(dist,'NonParametric')
    
    for i=1:length(Data)-scale+1
        data_sc(i)=sum(Data(i:i+scale-1));
    end
    
    for m =1:12   
        clear d data_rank p
        
        d=data_sc(m:12:length(data_sc));%data collected over the years for period sc
        
        %rank data
        for i=1:length(d)
            
            data_rank(i,1)=sum(d<=d(i));
        end
        
        %calculate the empirical probability
        p=(data_rank(:,1)-0.44)/(length(data_rank(:,1))+0.12);
        
        %calculate standardized index
        spei(m:12:length(data_sc))=norminv(p);     
    end
end

%Limiting to [-5,5] range.
spei(spei < -5) = -5; % Limiting results to -5 SPEI values
spei(spei > 5) = 5; % Limiting results to 5 SPEI values

end

