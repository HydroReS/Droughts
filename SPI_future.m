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
% 3) parm: a cell array containing the probility distribution
% coefficients of the historic period for each pixel
% 4) zero_prob: array containing the probability of zero for the historic
% period

%Output:
% 1) Standardized Precipitation Index (SPI): vector with SPI time series;

%% NOTE: SPI calculation require multi-year (e.g. 30 or more) of data to be reasonable
% and the code assumes complete years in record (i.e. record = 12month * number of years
%% function

function [U] = SPI_future(Data,scale,parm,zero_prob)

% check if data lenth is multiple of 12
if(rem(length(Data),12)~=0)
    error('data record is not multiple of 12 months....')
end

% Data setting to scaled dataset
A1=[];
for is=1:scale, A1=[A1,Data(is:length(Data)-scale+is)];end
XS=sum(A1,2);

%Pre-locating space for variable U
%U = zeros(size(Data,1)-scale+1,size(Data,2));

%Gamma parameter estimation and tranform
for stp =1:12
    tind=stp:12:length(XS); %Creates a time index
    Xn=XS(tind); %Gets value from Xn 
    q = zero_prob(stp); % saving zero probability
    % In case there are no zeros in the historic period, but there are in
    % the future
    if q ==0
       Xn(Xn==0) = 0.0000000001;
    end
    Gam_xs=q+(1-q)*gamcdf(Xn,parm(stp,1),parm(stp,2)); %calculates accumulated probabilty from gamma distribution;
    U(tind)=norminv(Gam_xs); %returns the SPI value from the inverse of accumulated probability;
end
