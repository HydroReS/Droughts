function [SI_drought,t_start,t_end]= SI_to_SIdroughts(SI,SI_thres)
%
% 
% identifies start/end of drought events and creates monthly times series
% equal to SPI length but considering SPI values only for the months that
% are part of a drought event 
%
% F. Marra, Aug 2020, modified from  
% - SI_to_droughts  
% - drought_charact_monthly_and_accum.m
% by Diogo and Thymios 
% 
%Input: 1) A vector with Standardized Index time series (e.g SPI)
%       2) A threshold value for SI to be used for drought identification
%       (usually SI_thres = -1)
%
%Output: A vector for each of the following variables:
%
%       1) SPI_drought: SPI time series with SPI values only when within a drought event
%       2) Start time index
%       3) End time index
%

%% identifies start and end of droughts
k1=find(SI<SI_thres | SI == SI_thres); %find index of SI that satify conditions of drought identification (i.e. < SI_thres)
k2=find(SI>0);

if isempty(k1)

	t_start = [];
	t_end = [];
	SI_drought = zeros(size(SI));
else

	for i=1:length(k1)

    	idx_start=find(k2<k1(i),1,'last');%find positive SI before the "drought incidence...SI<SI_thres"
    	if isempty(idx_start) %start of drought event at the beginning of the SI series
        	t_start(i)=1;
    	else
        	t_start(i)=k2(idx_start)+1;
    	end

    	idx_end=find(k2>k1(i),1,'first');
    	if isempty(idx_end) %end of drought event at the end of the SI series
        	t_end(i)=length(SI);
    	else
        	t_end(i)=k2(idx_end)-1;
    	end

	end

	%% removes duplicates
	t_start = unique(t_start,'stable');
	t_end = unique(t_end,'stable');

	%% gets SPI for drought events only
	SI_drought=zeros(size(SI));
	for t = 1:length(t_start)
	    if isnan(t_start(t)); continue; end
	    SI_drought(t_start(t):t_end(t))=SI(t_start(t):t_end(t));
	end
end

end
    
    
    
    
