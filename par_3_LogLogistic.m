%Calculate  three parameter Log-Logistic Distribution parameteres from
%Probability Weighted Moments (PWMs) using the unbiased estimator 
%described in Hosking (1986)
%Diogo S. A. Araújo
%Florida Institute of Technology
%Inputs:
%   1) Vector with data
%Outputs:
%   1) Vector with scale, shape and location parameters. [alpha, beta, gamma]

function parm = par_3_LogLogistic(data)

%Sort the vector into descending order
data = sort(data,'descend');

%Initializing PWMs
w0 = 0;
w1 = 0;
w2 = 0;

%Calculating PWMs for 0 order
for stp = 1:length(data)
    w0 = (nchoosek(length(data)-stp,0)*data(stp))*(1/length(data))/nchoosek(length(data)-1,0) + w0;
end

%Calculating PWMs for 1 order
for stp = 1:length(data)-1
    w1 = (nchoosek(length(data)-stp,1)*data(stp))*(1/length(data))/nchoosek(length(data)-1,1) + w1;
end

%Calculating PWMs for 2 order
for stp = 1:length(data)-2
    w2 = (nchoosek(length(data)-stp,2)*data(stp))*(1/length(data))/nchoosek(length(data)-1,2) + w2;
end

beta = (2*w1 - w0)/(6*w1-w0-6*w2);

alpha = (w0-2*w1)*beta/(gamma(1+1/beta)*gamma(1-1/beta));

gamma_parm = w0 - alpha*gamma(1+1/beta)*gamma(1-1/beta);

parm = [alpha, beta, gamma_parm];

end