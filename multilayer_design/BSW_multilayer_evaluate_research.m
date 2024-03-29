% Copyright 2020 Niccolò Marcucci <niccolo.marcucci@polito.it>
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
%
%     http://www.apache.org/licenses/LICENSE-2.0
%
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clear;
addpath('functions');

n_Ti2O2 = 2.53+1i*1e-4;           	% high refractive index
n_AlO   = 1.65+1i*1e-4;             % low refractive index
n_SiO2  = 1.46+1i*1e-4;             % low refractive index
n_pmma  = 1.48+1i*1e-4;
n_Ta2O5 = 2.08+1i*1e-4;            	% high refractive inde

theta = linspace(40,56,1e4);                    % angle vector
lambda=580e-9;

% the stack is fixed to be S-(B-A)xN-B-AlO-last-tail-air
load('best_param6');
% N=9;  
% best_thicknesses;
% best_indeces;

pol = 'p';  
scale_period = 1.0;

tail = best_thicknesses(1);
dlast = best_thicknesses(2);
detch = best_thicknesses(3);
dsecondlast_A = best_thicknesses(4)*scale_period;
dsecondlast_B = best_thicknesses(5)*scale_period;
dB = best_thicknesses(6)*scale_period;
dA = best_thicknesses(7)*scale_period;

n_in = best_indeces(1);        
n_out = best_indeces(2);       
nB = best_indeces(3);          
nA = best_indeces(4);          
nlast = best_indeces(5);      
ntail = best_indeces(6);       
netch = best_indeces(7);    

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

for k = 1:2
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
    [z1, ~, P ] = field_distribution(lambda,theta(idx),d,n,pol, 1e3);
    figure(2)
    hold on
    plot(z1,P)
        
end  
% n_eff(1)=n_eff(2);

n_eff(1)/n_eff(2)

[d1,n1] = prepare_multilayer(d1,n1);

[z1, nz] = field_distribution(lambda,theta(idx),d1,n1);
figure(2);
plot(z1,real(nz-1)*400);
legend('Field with last layer', 'Field without last layer',...
                                                '(n_z- 1) x 400');
nicePlot
figure(1);
legend('With last layer', 'Without last layer');
nicePlot

idx_layers=n1;
d_layers=d1;
n_eff1=n_eff(1);
n_eff2=n_eff(2);
n_eff3=n_eff(3);
name='designs/design_TM8_SiO2';
% stopBeforeSaving(name)
% save(name,'idx_layers','d_layers','n_eff1','n_eff2','n_eff3')