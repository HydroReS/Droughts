% Diogo S A Araujo
% Florida Institute of Techonology
% February 2020
% Modified code from SPI Programmed by Taesam Lee,  Dec.03,2009
%% Calculate Standardized Precipitation Index (SPI) of a data set of precipitation
%
%Input: 
% 1)Data: a column vector of monthly time-series of precipitation;
% 2)scale: an integer defining the length of time window i.e. in number of
% concecutive time steps.

%Output: 
% 1) Standardized Precipitation Index (SPI): vector with SPI time series;
% 2) Distribuition Gamma parameters alpha and beta: matrix with the
% parameters calculated for alpha and beta.

%% NOTE: SPI calculation require multi-year (e.g. 30 or more) of data to be reasonable
% and the code assumes complete years in record (i.e. record = 12month * number of years
%% function

function [U, gamma_parm, zero_prob] = SPI_modified(Data,scale)
% check if data lenth is multiple of 12
if(rem(length(Data),12)~=0)
    error('data record is not multiple of 12 months....')
end

% Data setting to scaled dataset
A1=[];
for is=1:scale, A1=[A1,Data(is:length(Data)-scale+is)];end
XS=sum(A1,2);

%Pre-locating space for variable U and gamma_parm
%U = zeros(size(Data,1),size(Data,2),size(Data,3)-scale+1);


%Gamma parameter estimation and tranform
for stp =1:12
    tind=stp:12:length(XS); %Creates a time index
    Xn=XS(tind); %Gets value from Xn 
    [zeroa]=find(Xn==0);
    Xn_nozero=Xn;Xn_nozero(zeroa)=[];
    q=length(zeroa)/length(Xn);
    parm = gamfit(Xn_nozero); %calculating the parameters for a gamma distribution;
    gamma_parm(stp,1:2) = parm; %saving the parameters in gamma_parm array;
    zero_prob(stp) = q; % saving zero probability
    Gam_xs=q+(1-q)*gamcdf(Xn,parm(1),parm(2)); %calculates accumulated probabilty from gamma distribution;
    U(tind)=norminv(Gam_xs); %returns the SPI value from the inverse of accumulated probability;
end
