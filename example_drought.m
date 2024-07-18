%% Example Script of how to use the functions in Drought Repository

%Creating Synthetic Precipitation data 
%Random creation of precipitation for 480 months ranging from 0 to 150
% mm
pr_hist = (150-0).*rand(480,1) + 0;
%Random creation of precipitation for 960 months ranging from 0 to 160
% mm
pr_fut = (160-0).*rand(960,1) + 0;

%Creating Synthetic PET data
%Random creation of PET for 480 months ranging from 0 to 200
% mm
pet_hist = (200-0).*rand(480,1) + 0;
%Random creation of PET for 960 months ranging from 0 to 210
% mm
pet_fut = (210-0).*rand(960,1) + 0;

%Time accumulation of the drought index in months(e.g. 3, 6, 12, etc)
scale = 12; %12-month timescale
% CDF for characterization of Extremes in SPEI calculation
dist = 'GeneralizedExtremeValue'; %GEV distribution

% Calculate SPI based on McKee et. al.(1993) and Stagge et. al.(2015) for
% historic period and saving parameters and zero probabilities
[spi_hist,gamma_par,zero_prob] = SPI_modified(pr_hist,scale);
% future period using parameters and zero probabilities from historic
% period
spi_fut = SPI_future(pr_fut,scale,gamma_par,zero_prob);

% Calculate SPEI based on Vicente-Serrano et. al.(2010) and Stagge et. al.(2015)
% for historic period and saving parameters and zero probabilities
[spei_hist,gev_par] = SPEI(pr_hist,pet_hist,scale,dist);
% future period using parameters and zero probabilities from historic
% period
spei_fut = SPEI_future(pr_fut,pet_fut,scale,dist,gev_par);

% Calculating drought characteristics
% t_start: starting index of each drought event
% t_end: ending index of each drought event
% D: durations of each drought event
% Im: mean intensities of each drought event
% Imax: max intensities of each drought event
% S: severities of each drought event
gen_figure = 1; %0 - don't create figure, %1-create figure
drought_thresh = -1; %Index threshold for drought identification (eg. 0, -1, -2)
%Historic SPEI
[t_start_hist_spei,t_end_hist_spei,D_hist_spei,Im_hist_spei,Imax_hist_spei,S_hist_spei] = ...
    SI_to_droughts(spei_hist,drought_thresh,gen_figure); 
%Historic SPI
[t_start_hist_spi,t_end_hist_spi,D_hist_spi,Im_hist_spi,Imax_hist_spi,S_hist_spi] = ...
    SI_to_droughts(spi_hist,drought_thresh,gen_figure);
%Future SPEI
[t_start_fut_spei,t_end_fut_spei,D_fut_spei,Im_fut_spei,Imax_fut_spei,S_fut_spei] = ...
    SI_to_droughts(spei_fut,drought_thresh,gen_figure);
%Future SPI
[t_start_fut_spi,t_end_fut_spi,D_fut_spi,Im_fut_spi,Imax_fut_spi,S_fut_spi] = ...
    SI_to_droughts(spi_fut,drought_thresh,gen_figure);


