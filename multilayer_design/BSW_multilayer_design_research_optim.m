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

clear, close all
addpath('functions');

n_TiO2  = 2.4444+1i*1e-4;           	% high refractive index
n_AlO   = 1.6282+1i*1e-4;             % low refractive index
n_SiO2  = 1.4615+1i*1e-4;             % low refractive index
n_pmma  = 1.48+1i*1e-4;
n_Ta2O5 = 2.08+1i*1e-4;            	% high refractive inde
n_high  = 3.00+1i*1e-4;
n_low   = 2.00+1i*1e-4;


theta = linspace(40,56,1e3);        % angles vector
lambda = 600e-9;
c0 = physconst('lightspeed');

n_in =  real(n_SiO2);
n_out = 1;
nA = n_TiO2;
nB = n_SiO2;
nlast = n_SiO2;
ntail = 1;
netch = n_AlO;

pol  = 'p';

N=9;

theta_lim = asin(n_out/n_in) + 1.5/180*pi;
theta = linspace(theta_lim/pi*180-3,56,1e3);        % angles vector

kB = 2*pi/lambda*real(nB)*cos(40/180*pi);
kA = 2*pi/lambda*real(nA)*cos(40/180*pi);
FF = 0.5;
period = 1*pi/(kB*FF+kA*(1-FF));
dB_0=period*FF;
dA_0=period*(1-FF);

param_file='best_param8';

best_thicknesses = zeros(7,1);
best_indeces = zeros(7,1);
%{
best_thicknesses(1) = 0;            % second layer above etch stop
best_thicknesses(2) = dA_0;         % first layer above etch stop
best_thicknesses(3) = 15e-9;        % etch stop layer
best_thicknesses(4) = dA_0;         % layer below etch stop
best_thicknesses(5) = dB_0;         % second layer below etch stop
best_thicknesses(6) = dB_0;         % dielectric B
best_thicknesses(7) = dA_0;         % dielectric A
%}

% copied from another deign
best_thicknesses(1) = 0;            % second layer above etch stop
best_thicknesses(2) = 69e-9;         % first layer above etch stop
best_thicknesses(3) = 10e-9;        % etch stop layer
best_thicknesses(4) = 69e-9;         % layer below etch stop
best_thicknesses(5) = 155e-9;         % second layer below etch stop
best_thicknesses(6) = 155e-9;         % dielectric B
best_thicknesses(7) = 69e-9;         % dielectric A

best_indeces(1) = n_in;             % substrate refractive index
best_indeces(2) = n_out;            % air refractive index
best_indeces(3) = nB;               % second material
best_indeces(4) = nA;               % first material
best_indeces(5) = nlast;            % first layer above etch stop 
best_indeces(6) = ntail;            % second layer above etch stop
best_indeces(7) = netch;            % etch stop


%%%
% load('best_param6');
n_in = best_indeces(1);        
n_out = best_indeces(2);       
nB = best_indeces(3);          
nA = best_indeces(4);          
nlast = best_indeces(5);      
ntail = best_indeces(6);       
netch = best_indeces(7);    

% the stack is fixed to be S-(A-B)xN-A-AlO-last-tail-air
N_layers = 2*N+6;

n1 = zeros(N_layers,1);             % with all layer including pmma
d1 = zeros(N_layers,1);

n_eff = Inf*ones(3,1);
R = zeros(3,length(theta));
Xi = ones(10,1);                    % partition function contributions
Partition = 0;                      % total partition
i=1;
rng(mod(now*1e10,2^32));
tol=1;

tic;
f = waitbar(0,'Please wait...');
N_max = 1e4;

q = 1;
RRR = zeros(N_max,2);
PPP = zeros(N_max,3);

n1(1) = n_in;
n1(2:2:2*N) = nA;
n1(3:2:2*N+1) = nB;
n1(end-4) = nA;
n1(end-3) = netch;
n1(end-2) = nlast;
n1(end-1) = ntail;
n1(end) = n_out;
    
for j = 1: N_max

    value = best_thicknesses*1e9;
    rangee = value * 0.5;
    if j == 1
        parameters = value*1e-9;
    else
        % parameters = round(value-rangee.*(rand(7,1)-0.5))*1e-9;
        parameters = (value-rangee.*(rand(7,1)-0.5))*1e-9;
    end
    
    % constrains
    % parameters(4) = parameters(7);
    parameters(5) = parameters(6);
    parameters(3)=10e-9;
    parameters(1) = 0;
    parameters(2) = 60e-9;


    tail = parameters(1);
    dlast= parameters(2);      
    detch = parameters(3);
    dsecondlast_A = parameters(4);
    dsecondlast_B = parameters(5);
    dB = parameters(6);
    dA = parameters(7);
    
    d1(1) = 0;
    d1(2:2:2*N) = dA;
    d1(3:2:2*N+1) = dB;
    d1(end-5) = dsecondlast_B;
    d1(end-4) = dsecondlast_A;
    d1(end-3) = detch;
    d1(end-2) = dlast;
    d1(end-1) = tail;
    d1(end) = 0;
    
    min_r = 0.450;
    max_r = 0.60;
    min_H = 1000;
    max_H = 2000;
    % for k = 2:3
    for k = 2:3
        n = [n1(1:end-k) ; n1(end)];
        d = [d1(1:end-k) ; d1(end)];
        [Rk,r,t] = reflectivity(lambda,theta,d,n,pol);
        R(k,:) = Rk;
        [pks,idxs] = findpeaks(1-Rk);
        [~,pk_ix] = max(pks);
        idx = idxs(pk_ix);
        
        if isempty(idx)
            Xi(k) = 0;
            % Xi(3+k) = 0;
        else
            n_eff(k) = sin(theta(idx)*pi/180)*n_in;

            P_end = abs(t(idx))^2;

            % contributions to the partition function

            Xi(k) = threshold(1-Rk(idx),min_r,max_r);
            % % % Xi(3+k) = threshold(1-n_eff(k), 0,-0.05);
            % Xi(3+k) = threshold(P_end,min_H,max_H);
            if n_eff(k) > 1.008
                RRR(q,k-1) = 1-Rk(idx);
                PPP(q,:) = [dA, dB, dsecondlast_A];
                if k == 3 
                    q=q+1;
                end
            end
        end
    end
    % contrast1a=n_eff(1)-n_eff(2);
    % contrast1b=n_eff(1)/n_eff(2);
    contrast2a=-(n_eff(3)-n_eff(2));
    % contrast2b=n_eff(2)/n_eff(3);
    
    % contribPr =utions to the partition function
    
    min_c=0.080;
    max_c=0.1;
    % Xi(7) = threshold(contrast1a,min_c,max_c);
    % Xi(8) = threshold(contrast2a,min_c,max_c);
    % Xi(9) = threshold(contrast1b,min_c+1,max_c+1);
    % Xi(10) = threshold(contrast2b,1+min_c,1+max_c);

    Xi_new = prod(Xi);
    if (Xi_new > Partition) && all(n_eff > 1.008)
        Partition = Xi_new;
        best_thicknesses = parameters;
        j

        plot(theta,R(2,:),theta,R(3,:));
        nicePlot
        drawnow
    end
    if mod(j,10)==1
        t = toc;
        waitbar(j/N_max,f,['iteration ',num2str(j), ': missing ',...
                                num2str(round((1 - j/N_max)*(t/(j/N_max)))), 's']);
    end 
end
delete(f)

save(param_file,'best_thicknesses','best_indeces','N')

%%
palette = hot(2^9);
figure
subplot(1,2,1)
c =  squeeze(getRGB(RRR(:,1),palette));
scatter3(PPP(:,1)*1e9,PPP(:,2)*1e9,PPP(:,3)*1e9, 10, c,'filled')
hold on
plot3(best_thicknesses(7)*1e9,best_thicknesses(6)*1e9,best_thicknesses(4)*1e9,'+w','linewidth', 5,'MarkerSize',20)
set(gca, 'Color', [0 0 0])
subplot(1,2,2)
c =  squeeze(getRGB(RRR(:,2),palette));
scatter3(PPP(:,1)*1e9,PPP(:,2)*1e9,PPP(:,3)*1e9, 10, c,'filled')
hold on
plot3(best_thicknesses(7)*1e9,best_thicknesses(6)*1e9,best_thicknesses(4)*1e9,'+w','linewidth', 5,'MarkerSize',20)
set(gca, 'Color', [0 0 0])


function rgbImage = getRGB(data,map,range)
    if nargin < 3
        range = data;
    end
    dmin = min(data(:));
    dmax = max(data(:));
    if nargin > 2 
        data =   round( (data - min(range)) / (max(range)-min(range)) *( length(map)-1 ) ) +1 ;  % normalize again inside range and scale to map
    else
        data = round((data - dmin)/(dmax - dmin)*( length(map)-1 ) ) +1 ;                           % normalize between 0 and 1 and scale to map
    end
    rgbImage = ind2rgb(data,map);
end

% arbitrary partition functions
function y =threshold(x, minimum, maximum)
    width = maximum-minimum;
    average = (maximum+minimum)/2;
    KbT = width/5;
    y = 1-1/(exp((x-average)/KbT)+1);
end
function y = range_f(x, minimum, maximum)
    fwhm = maximum-minimum;
    average = (maximum+minimum)/2;
    sigma = fwhm/4;
    y = 1-exp(-(x-average)^2/(2*sigma^2));
end
% random tossing
function y = toss(value, ranges)
    y = round(value-ranges*(rand-0.5));
end

function f_optim(x,n1)
    % constrains
    % x(4) = parameters(7);
    parameters(5) = parameters(6);
    parameters(3)=10e-9;
    parameters(1) = 0;
    parameters(2) = 60e-9;


    tail = parameters(1);
    dlast= parameters(2);      
    detch = parameters(3);
    dsecondlast_B = parameters(5);
    dsecondlast_A = x(3);
    dB = x(2);
    dA = x(1);
    
    d1(1) = 0;
    d1(2:2:2*N) = dA;
    d1(3:2:2*N+1) = dB;
    d1(end-5) = dsecondlast_B;
    d1(end-4) = dsecondlast_A;
    d1(end-3) = detch;
    d1(end-2) = dlast;
    d1(end-1) = tail;
    d1(end) = 0;
    
    min_r = 0.450;
    max_r = 0.60;
    Xi = 1;
    for k = 2:3
        n = [n1(1:end-k) ; n1(end)];
        d = [d1(1:end-k) ; d1(end)];
        [Rk,r,t] = reflectivity(lambda,theta,d,n,pol);
        [pks,idxs] = findpeaks(1-Rk);
        [~,pk_ix] = max(pks);
        idx = idxs(pk_ix);
        
        if isempty(idx)
            Xi = 0;
        else
            n_eff = sin(theta(idx)*pi/180)*n_in;

            P_end = abs(t(idx))^2;

            % contributions to the partition function

            Xi(k) = threshold(1-Rk(idx),min_r,max_r);
        end
    end
    return prod(Xi)
end
