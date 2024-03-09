%function that identifies start/end of drought events and calculates drought characteristics (mean/max intensity, duration and severity).
%identification of drought events is performed based on a threshold value of the
%standardized index.
%
%Input: 1) A vector with Standardized Index time series (e.g SPI)
%       2) A threshold value for SI to be used for drought identification
%       (usually SI_thres = -1)
%       3)Option "fig_gen" for generating a figure or not (1=generate figure,
%       0= no figure)
%
%Output: A vector for each of the following variables:
%
%       1) Start time index
%       2) End time index
%       3) Duration (in # of months)
%       4) Mean intensity
%       5) Max intensity
%       6) Severity
%
%
%

function [t_start,t_end,D,Im,Imax,S]= SI_to_droughts(SI,SI_thres,fig_gen)

%% Identify start/end index of drought events

k1=find(SI<SI_thres); %find index of SI that satify conditions of drought identification (i.e. < SI_thres)
k2=find(SI>=0);

if isempty(k1)
	t_start = [];
	t_end = [];
	D = [];
	Im = [];
	Imax = [];
	S = [];
	return
end


for i=1:length(k1)

idx_start=find(k2<k1(i),1,'last');%find positive SI before the "drought incidence...SI<SI_thres"
if(isempty(idx_start)==1)%start of drought event at the beginning of the SI series
    t_start(i)=1;
else
    t_start(i)=k2(idx_start)+1;
end

idx_end=find(k2>k1(i),1,'first');
if(isempty(idx_end)==1)%end of drought event at the end of the SI series
    t_end(i)=length(SI);
else
    t_end(i)=k2(idx_end)-1;
end

end

%remove duplicates
t_start=unique(t_start,'stable');
t_end=unique(t_end,'stable');


%% Calculate drought characteristics
for i=1:length(t_start)
D(i) = 1+ (t_end(i) - t_start(i)); %duration of drought events

Im(i)=mean(SI(t_start(i):t_end(i))); %Mean intensity (mean SI value)

Imax(i)=min(SI(t_start(i):t_end(i))); %Max (negative) intensity 

S(i) = D(i)*Im(i); %Severity of drought event
end

%% Generate figure

if fig_gen==1
    
    figure
    
    plot(SI,'b-'); hold on
    rf=refline(0,0);set(rf,'color','k','linewidth',2)
    for i=1:length(t_start)
          area((t_start(i):t_end(i)),(SI(t_start(i):t_end(i))),'FaceColor','r')
    end
    
    rf_1=refline(0,-1);set(rf_1,'color','y','linewidth',2,'linestyle','--')

    rf_2=refline(0,-2);set(rf_2,'color','r','linewidth',2,'linestyle','--')

legend([rf_1,rf_2],{'moderate','severe'})

ylabel('Standardized Index','fontsize',20)
xlabel('Time (months)','fontsize',20)
end
    
    
    
    
