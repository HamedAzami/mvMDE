function Out_mvMDE=mvMDE(X,m,c,stau,Scale)
% This function calculates multivariate multiscale dispersion entropy (mvMDE) using normal cumulative distribution function (NCDF).
%
% mvMDE deals with the shortcomings of mvMDE_I, mvMDE_II, and mvMDE_III although each of these four methods has its own pros and cons.
%
% Inputs:
% X: multivariate signal - a matrix of size nvar (the number of channels) x nsamp (the number of sample points for each channel)
% m: scalar embedding value
% c: number of classes (we usullay set c=5, 6, or 7)
% stau: scalar time lag  value (it is usually equal to 1)
% Scale: maximum number of scale factors
%
% Output:
% Out_mvMDE: a vector of size 1 * Scale - the mvMDE of X
%
% Ref:
% [1] H. Azami, A. Fernandez, and J. Escudero, "Multivariate Multiscale Dispersion Entropy of Biomedical Times Series", Entropy, 21(9), 913, 2019.
% [2] H. Azami, M. Rostaghi, D. Abasolo, and J. Escudero,"Refined composite multiscale dispersion entropy and its
% application to biomedical signals", IEEE Transactions on Biomedical Engineering, 64(12), 2872-2879, 2017.
%
% Hamed Azami and Javier Escudero
% hmd.azami@gmail.com and javier.escudero@ed.ac.uk
%
%  10-Sep-19
%%
Out_mvMDE=NaN*ones(1,Scale);
% multivariate dispersion entropy at scale 1.
Out_mvMDE(1)=mvDE_NCDF(X,m,c,stau);


% multivariate dispersion entropy at temporal scales 2 to maximum scale factor.
sigma=std(X,0,2);
mu=mean(X,2);

for j=2:Scale
    Xs = Multi(X,j);
    Out_mvMDE(j)=mvDE_NCDF_ms(Xs,m,c,mu,sigma,stau);
end


function M_Data = Multi(Data,S)

%  generate the consecutive coarse-grained time series
%  Input:   Data: time series;
%           S: scale factor

% Output:
%           M_Data: coarse-grained time series at scale factor S

L = size(Data,2);
J = fix(L/S);

for j=1:size(Data,1)
    for i=1:J
        M_Data(j,i) = mean(Data(j,(i-1)*S+1:i*S));
    end
end
