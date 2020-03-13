clear, close all
addpath('functions');


n_Ti2O2 = 2.53+1i*1e-4;           	% high refractive index
n_AlO = 1.65+1i*1e-4;               % low refractive index
n_SiO2 = 1.46;                    	% low refractive index
n_pmma = 1.48+1i*1e-4;
n_Ta2O5 = 2.08+1i*1e-4;            	% high refractive inde

theta = linspace(40,56,1e4);        % angles vector
lambda=570e-9;
c0 = physconst('lightspeed');

n_in = n_SiO2;
n_out = 1;
nA = n_Ti2O2;
nB = n_AlO;
nlast = n_SiO2+1i*1e-4;
ntail = n_pmma;
netch = n_AlO;

pol  = 'p';
N=11;

theta_lim = asin(n_out/n_in) + 0.7/180*pi;
kB=2*pi/lambda*real(nB)*cos(40/180*pi);
kA=2*pi/lambda*real(nA)*cos(40/180*pi);
FF=0.7;
period=pi*0.95/(kB*FF+kA*(1-FF));
dB_0=period*FF;
dA_0=period*(1-FF);
    

best_thicknesses = zeros(7,1);
best_indeces = zeros(7,1);

best_thicknesses(1) = 0;            % second layer above etch stop
best_thicknesses(2) = dB_0;         % first layer above etch stop
best_thicknesses(3) = 15e-9;        % etch stop layer
best_thicknesses(4) = dA_0;         % layer below etch stop
best_thicknesses(5) = dB_0;         % second layer below etch stop
best_thicknesses(6) = dB_0;         % dielectric 1
best_thicknesses(7) = dA_0;         % dielectric 2

best_indeces(1) = n_in;             % substrate refractive index
best_indeces(2) = n_out;            % air refractive index
best_indeces(3) = nB;               % second material
best_indeces(4) = nA;               % first material
best_indeces(5) = nlast;            % first layer above etch stop 
best_indeces(6) = ntail;            % second layer above etch stop
best_indeces(7) = netch;            % etch stop


%%%
% loaddata('best_param2');
% n_in = best_indeces(1);        
% n_out = best_indeces(2);       
% nB = best_indeces(3);          
% nA = best_indeces(4);          
% nlast = best_indeces(5);      
% ntail = best_indeces(6);       
% netch = best_indeces(7);

% the stack is fixed to be S-(A-B)xN-A-AlO-last-tail-air
N_layers=2*N+6;

n1 = zeros(N_layers,1);             % with all layer including pmma
d1 = zeros(N_layers,1);

n_eff = Inf*ones(3,1);
R = zeros(3,length(theta));
Xi = ones(10,1);                    % Partition function contribution
Xi1 = ones(3,1); 
X12 = ones(3,1); 
Partition = 0;
i=1;
rng(now);
tol=1;
% Pm=zeros(11,400);

for j = 1:3
    
    % randomize the parameter. One at a time. 
    i=i+1;
    if i==8                 
        % reset shuffling index
        i=2;
    end
    if i==3                 
        % don't change etch stop thickness (comment if needed)
        i=4;
    end
    if i==4 || i == 5
        % don't change below etch stop (comment if needed)
        i=6;
    end
    
    parameters = best_thicknesses ;

    value = best_thicknesses(i)*1e9;
    switch i
        case 1              
            % pmma thickness
            rangee = value/2;
        case 2              
            % last layer thickness
            rangee = value/2;
        case 3              
            % AlO layer thicknes (etch stop layer)
            rangee = value/2;
        case 4              
            % first below etch stop
            rangee = value*2;
        case 5              
            % second below etch stop
            rangee = value/2;
        case 6
            % B
            rangee = value*0.4;
        case 7
            % A
            rangee = value*0.4;
    end
    parameters(i)=toss(value,rangee)*1e-9;
    
    parameters(4) = parameters(7);
    parameters(5) = parameters(6);
    
    tail = parameters(1);
    dlast= parameters(2);      
    detch = parameters(3);
    dsecondlast_A = parameters(4);
    dsecondlast_B = parameters(5);
    dB = parameters(6);
    dA = parameters(7);
    
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
    d1(end) = 0;
    
    min_r = 0.40;
    max_r = 0.60;
    min_H = 1000;
    max_H = 2000;
    for k = 2:3
        n = [n1(1:end-k) ; n1(end)];
        d = [d1(1:end-k) ; d1(end)];
        [Rk,r,t] = reflectivity(lambda,theta,d,n,pol);
        R(k,:) = Rk;
        [pks,idxs] = findpeaks(1-Rk);
        [~,pk_ix] = max(pks);
        idx = idxs(pk_ix);
        
        if isempty(idx)
            Xi1(k) = 0;
            Xi2(k) = 0;
        else
            n_eff(k) = sin(theta(idx)*pi/180)*n_in;

            P_end = abs(t(idx))^2;

            % contributions to the partition function

            Xi1(k) = threshold(1-Rk(idx),min_r,max_r);

            Xi2(k) = threshold(P_end,min_H,max_H);
        end
    end
%     contrast1a=n_eff(1)-n_eff(2);
%     contrast1b=n_eff(1)/n_eff(2);
    contrast2a=-(n_eff(3)-n_eff(2));
    contrast2b=n_eff(2)/n_eff(3);
    
    % contributions to the partition function
    
    Xi(1:3) = Xi1;
    Xi(4:6) = Xi2;
    
    min_c=0.110;
    max_c=0.130;
%     Xi(7) = threshold(contrast1a,min_c,max_c);
    Xi(8) = threshold(contrast2a,min_c,max_c);
%     Xi(9) = threshold(contrast1b,min_c+1,max_c+1);
    Xi(10) = threshold(contrast2b,min_c+1,max_c+1);
     
    Xi_new = prod(Xi);
%     Pm(1:end-1,j)=Xi;
    if (Xi_new > Partition) && all(n_eff > sin(theta_lim)*n_in)
        Partition = Xi_new;
        best_thicknesses = parameters;
        j

        plot(theta,R(2,:),theta,R(3,:));
        drawnow
    end
%     Pm(end,j)=Partition;
end
% Pm=log10(Pm);
% plot(1:400,Pm(2,:),Pm(4,:),Pm(5,:),Pm(7,:),Pm(8,:),Pm(10,:),Pm(11,:));
% plot(1:400,10^Pm(11,:));

save('best_param4','best_thicknesses','best_indeces','N')
% savedata('best_param3',best_parameters,best_index);

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