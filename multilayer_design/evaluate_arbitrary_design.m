close all;
clear;
addpath('functions');

n_Ti2O2 = 2.53+1i*1e-4;           	% high refractive index
n_AlO   = 1.65+1i*1e-4;             % low refractive index
n_SiO2  = 1.46+1i*1e-4;             % low refractive index
n_pmma  = 1.48+1i*1e-4;
n_Ta2O5 = 2.08+1i*1e-4;            	% high refractive inde

n_in=2;
n_out=1;
d1=zeros(21,1);
n1=ones(21,1);

n1(1)=n_in;
n1(2:2:end)=n_Ti2O2;
n1(3:2:end)=n_AlO;
n1(end-3)=n_Ti2O2;
n1(end-2)=n_AlO;
n1(end-1)=n_SiO2;
n1(end)=n_out;

d1(1)=1e-6;
d1(2:2:end)=100e-9;
d1(3:2:end)=160e-9;
d1(end-3)=60e-9;
d1(end-2)=15e-9;
d1(end-1)=60e-9;
d1(end)=3e-6;

pol='p';
lambda=570e-9;

theta=linspace(40,70,1e5);

for k =1:2
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
    plot(z1,P)
        
end  
n_eff(3)=n_eff(2);

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

idx_layers=n1;
d_layers=d1;
n_eff1=n_eff(1);
n_eff2=n_eff(2);
n_eff3=n_eff(3);
name='designs/design_TM_Desc2';
% stopBeforeSaving(name)
% save(name,'idx_layers','d_layers','n_eff1','n_eff2','n_eff3')