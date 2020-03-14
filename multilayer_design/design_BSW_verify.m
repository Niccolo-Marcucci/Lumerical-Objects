close all;
clear;
addpath('functions');

n_Ti2O2 = 2.53+1i*1e-4;           	% high refractive index
n_AlO   = 1.65+1i*1e-4;             % low refractive index
n_SiO2  = 1.46+1i*1e-4;             % low refractive index
n_pmma  = 1.48+1i*1e-4;
n_Ta2O5 = 2.08+1i*1e-4;            	% high refractive inde

theta = linspace(40,56,1e4);                    % angle vector
lambda=570e-9;

% the stack is fixed to be S-(B-A)xN-B-AlO-last-tail-air
load('best_param4');%num2str(kkk)
% N=8;  
% best_thicknesses;
% best_indeces;

pol = 'p';
scale_period = 1;

tail = best_thicknesses(1);
dlast = best_thicknesses(2);
detch = best_thicknesses(3);
dsecondlast_A = best_thicknesses(4)*scale_period;
dsecondlast_B = best_thicknesses(5)*scale_period;
dB = best_thicknesses(6)*scale_period;
dA = best_thicknesses(7)*scale_period;

n_in = best_indeces(1);        
n_out = best_indeces(2);       
nB = n_AlO;best_indeces(3);          
nA = best_indeces(4);          
nlast = best_indeces(5);      
ntail = best_indeces(6);       
netch = best_indeces(7);    

N=14;
n_Ti2O2 = 2.53+1i*1e-4;           	% high refractive index
n_AlO   = 1.65+1i*1e-4;             % low refractive index
n_SiO2  = 1.46+1i*1e-4;             % low refractive index
n_pmma  = 1.48+1i*1e-4;
n_Ta2O5 = 2.08+1i*1e-4;            	% high refractive inde
dlast = 60e-9;% parameters(2)*0.85;
detch = 15e-9;%parameters(3);
dsecondlast_B = 133e-9;%parameters(4)*scale_period;
dsecondlast_A = 57e-9;%parameters(5)*scale_period;
dA = 57e-9;%parameters(6)*scale_period;
dB = 133e-9;
% the stack is fixed to be S-(A-B)xN-A-etch-last-tail-air
N_layers=2*N+6;    

n1 = zeros(N_layers,1);             % with all layer including pmma
d1 = zeros(N_layers,1);

n_eff = Inf*ones(3,1);

n1(1) = n_in;
n1(2:2:2*N) = nA;
n1(3:2:2*N+1) = nB;
n1(end-4) = nA;
n1(end-3) = netch;
n1(end-2) = nlast;
n1(end-1) = ntail;
n1(end) = n_out;

d1(1) = 1e-6;
d1(2:2:2*N) = dA;
d1(3:2:2*N+1) = dB;
d1(end-5) = dsecondlast_B;
d1(end-4) = dsecondlast_A;
d1(end-3) = detch;
d1(end-2) = dlast;
d1(end-1) = tail;
d1(end) = 3e-6;

for k = 2:3
    n = [n1(1:end-k) ; n1(end)];
    d = [d1(1:end-k) ; d1(end)];
    
    [dr,nr,~,~] = prepare_multilayer(d,n);
    
    [R,r,t] = reflectivity(lambda,theta,dr,nr,pol);
    [pks,idxs] = findpeaks(1-R);
    [~,pk_ix] = max(pks);
    idx = idxs(pk_ix);
    n_eff(k)= sin(theta(idx)*pi/180)*n_in;

    figure(1)
    hold on
    plot(theta,R)
%     R(idx)
%     t(idx)
    [z1, ~, P ] = field_distribution(lambda,theta(idx),d,n,...
                                                r(idx),t(idx),pol);
    figure(2)
    hold on
    plot(z1,2*P)
        
end  
n_eff(1)=n_eff(2);

n_eff(1)/n_eff(3)

[d1,n1] = prepare_multilayer(d1,n1);

[z1, nz] = field_distribution(lambda,theta(idx),d1,n1);
figure(2);
plot(z1,real(nz)*400);
legend('Field with last layer', 'Field without last layer', 'Refreactive index x 400');
nicePlot
figure(1);
legend('With last layer', 'Without last layer');
nicePlot

% savedata('designs/best_design_TM5',idx_layers,d_layers,n_eff1,n_eff2,n_eff3);
