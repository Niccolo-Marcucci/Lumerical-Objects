clear;
close all;
addpath('functions');

n_Ti2O2 = 2.53+1i*1e-4;           	% high refractive index
n_AlO   = 1.65+1i*1e-4;             % low refractive index
n_SiO2  = 1.46+1i*1e-4;             % low refractive index
n_pmma  = 1.48+1i*1e-4;
n_Ta2O5 = 2.08+1i*1e-4;            	% high refractive inde

% new design general properties
design_name = "TM9";
design_file = strcat("designs/design_",design_name,".mat");
pol='p';                            % polarisation: 'p' or 's'
design_type='buried';               % either 'buried' or empty

lambda=570e-9;
theta=linspace(30,70,1e5);

% define layers parameters
if strcmp(design_type,'buried')
    n_in=2;
else
    n_in=real(n_SiO2);
end
n_out=1;
d1=zeros(13,1);
n1=ones(13,1);

%% Specify layer properties in group
% n1(1)=n_in;n_Ti2O2;
% n1(2:2:end)=n_Ti2O2;
% n1(3:2:end)=n_AlO;
% n1(end-4)=n_Ti2O2;
% n1(end-3)=n_AlO;%n_SiO2;%n_Ti2O2;%
% n1(end-2)=n_SiO2;%n_Ti2O2;%
% n1(end-1)=n_Ti2O2;%n_SiO2;%
% n1(end)=n_out;

% d1(1)=0.4e-6;
% d1(2:2:end)=105e-9;
% d1(3:2:end)=169.722826087e-9;
% d1(end-4)=100e-9;
% d1(end-3)=430e-9;
% d1(end-2)=20e-9;
% d1(end-1)=85e-9;
% d1(end)=3e-6;

%% Specify layer properties, layer by layer
n1(1)=n_in;
n1(2)=n_Ti2O2;
n1(3)=n_AlO;
n1(4)=n_Ti2O2;
n1(5)=n_AlO;
n1(6)=n_Ti2O2;
n1(7)=n_AlO;
n1(8)=n_Ti2O2;
n1(9)=n_AlO;
n1(10)=n_SiO2;
n1(11)=n_Ti2O2;
n1(12)=1;

d1(1)=0.4e-6;
d1(2)=105e-9;
d1(3)=169.722826087e-9;
d1(4)=99.836956522e-9;
d1(5)=169.7228260866e-9;
d1(6)=99.8369565218e-9;
d1(7)=169.7228260869e-9;
d1(8)=99.8369565218e-9;
d1(9)=429.29891304346e-9;
d1(10)=19.967391304350e-9;
d1(11)=79.869565217390e-9;
d1(12)=249.184782608699e-9;

[d1,n1] = prepare_multilayer(d1,n1);    % removes 0 thickness layers
                                        % if any
for k =1:2
    if strcmp(design_type,'buried')
        n = n1;
        d = d1;
        n(end-2)=n1(end-k);
    else
        n = [n1(1:end-k) ; n1(end)];
        d = [d1(1:end-k) ; d1(end)];
    end
    [dr,nr,~,~] = prepare_multilayer(d,n);  
    
    [R,r,t] = reflectivity(lambda,theta,dr,nr,pol);
    [pks,idxs] = findpeaks(1-R);
    [~,pk_ix] = max(pks);
    idx = idxs(pk_ix);
    n_eff(k)= sin(theta(idx)*pi/180)*n_in;

    figure(1)
    hold on
    plot(theta,R)
    
    [z, ~, P] = field_distribution(lambda,theta(idx),d,n,...
                                                r(idx),t(idx),pol);
    figure(2)
    hold on
    plot(z,P)
        
end  
n_eff(3)=n_eff(2);

n_eff(1)/n_eff(3)

[z, nz] = field_distribution(lambda,theta(idx),d1,n1);
figure(2);
plot(z,real(nz-1)*400);
if strcmp(design_type,'buried')
    legend('Field without buried layer','Field with buried layer',...
                                                '(n_z -1) x 400');
else
    legend('Field with last layer','Field without last layer',...
                                                '(n_z -1) x 400');
end
nicePlot

figure(1);
if strcmp(design_type,'buried')
    legend('Without buried layer','With buried layer');
else
    legend('With last layer','Without last layer');
end
nicePlot
% set(gca,'yscale','log')

idx_layers=n1;
d_layers=d1;
n_eff1=n_eff(1);
n_eff2=n_eff(2);
n_eff3=n_eff(3);

stopBeforeSaving(design_file)
save(design_file,'idx_layers','d_layers','n_eff1','n_eff2','n_eff3')