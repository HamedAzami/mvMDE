function [Out_mvDE, npdf]=mvDE_NCDF_ms(X,m,c,mu,sigma,stau)
% This function calculates multivariate dispersion entropy (mvMDE) using normal cumulative distribution function (NCDF) with defined mean (mu) and standard deviation (sigma) values.
%
% mvDE deals with the shortcomings of mvDE_I, mvDE_II, and mvDE_III although each has its own pros and cons. 
%
% Inputs:
% X: multivariate signal - a matrix of size nvar (the number of channels) x nsamp (the number of sample points for each channel)
% m: scalar embedding value
% c: number of classes (we usullay set c=5, 6, or 7)
% stau: scalar time lag  value (it is usually equal to 1)
%
% Output:
% Out_mvMDE: a scaler value - the mvDE of X
%
% Ref:
% [1] H. Azami, A. Fernandez, and J. Escudero, "Multivariate Multiscale Dispersion Entropy of Biomedical Times Series", Entropy, 2019.
% [2] H. Azami, M. Rostaghi, D. Abasolo, and J. Escudero,"Refined composite multiscale dispersion entropy and its
% application to biomedical signals", IEEE Transactions on Biomedical Engineering, 64(12), 2872-2879, 2017.
%
% Hamed Azami and Javier Escudero
% hmd.azami@gmail.com and javier.escudero@ed.ac.uk
%
%  10-Sep-19
%%
CH=size(X,1); %number of channels
N=size(X,2); % Length of each signal

for i_CH=1:CH
    x=X(i_CH,:);
    %time series is normalized between (0,1)
    % % %Normal cumulative distribution function (NCDF) to map x into y from 0 to 1 %% Normalize with normcdf
    
    x=normcdf(x,mu(i_CH),sigma(i_CH));
    
   x(x==1)=1-1e-10;
    x(x==0)=1e-10;
    
    X(i_CH,:)=round((x*c)+0.5);
end

tau=stau*ones(1,CH);
M=m*ones(1,CH);
S_M=sum(M);

y=[];
for j_embd=1:CH
    for i_embd=1:N-m+1
        temp1(i_embd,:)=X(j_embd,i_embd:tau(j_embd):i_embd+M(j_embd)-1);
    end
    y=horzcat(y,temp1);
    temp1=[];
end

%time series is normalized between (0,1)
% % %Normal cumulative distribution function (NCDF) to map x into y from 0 to 1 %% Normalize with normcdf


%  generate all possible dispersion pattern
all_pattern=[1:c]';
for f=2:m
    temp=all_pattern;
    all_pattern=[];
    j=1;
    for w=1:c
        [a,b]=size(temp);
        all_pattern(j:j+a-1,:)=[temp,w*ones(a,1)];
        j=j+a;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pattern_num=c^m;% obtain the number of all possible dispersion pattern
key=zeros(1,pattern_num);
for i_p=1:pattern_num
    for ii_p=1:m
        key(i_p)=key(i_p)*10+all_pattern(i_p,ii_p);
    end
end

pdf=zeros(1,pattern_num); %initialize frequency array
Comb_M_m=nchoosek(S_M,m);

pv=zeros(Comb_M_m,N-m+1);
aa=combnk(1:S_M,m);

for i_pv=1:Comb_M_m
    for iii_p=m:-1:1
        pv(i_pv,:)=pv(i_pv,:)+y(:,aa(i_pv,iii_p))'*10^(iii_p-1);
    end
end


for i=1:pattern_num
    
    [a ,ignore]=find(pv==key(i));
    pdf(i)=length(a); %
    
end

npdf=pdf/((N-m+1)*Comb_M_m); % normalize the frequency array to obtain probability density function. 

p=npdf(npdf~=0);

Out_mvDE=-sum(p .* log(p));
end