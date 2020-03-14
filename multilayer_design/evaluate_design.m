clear
close all
addpath('functions');

load best_design_TE.mat
load field_TE

d1=d_layers;
n1=idx_layers;

lambda=570e-9;
theta = linspace(40,56,1e4);     
pol='s';

for k = 1:2
    n = [n1(1:end-k) ; n1(end)];
    d = [d1(1:end-k) ; d1(end)];
    
    [dr,nr,~,~] = prepare_multilayer(d,n);
    
    [R,r,t] = reflectivity(lambda,theta,dr,nr,pol);
    [pks,idxs] = findpeaks(1-R);
    [~,pk_ix] = max(pks);
    idx = idxs(pk_ix);
    figure(1)
    hold on
    plot(theta,R)
%     R(idx)
%      t(idx)
    [z1, ~, P ] = field_distribution(lambda,theta(idx),d,n,...
    r(idx),t(idx),pol);
    figure(2)
    hold on
    plot(z1,2*P)
end


[d1,n1] = prepare_multilayer(d1,n1);

[z1, nz] = field_distribution(lambda,theta(idx),d1,n1);
figure(2);
plot(z1,real(nz)*400);
legend('Field with last layer', 'Field without last layer', 'Refreactive index x 400');
figure(1);
legend('With last layer', 'Without last layer');
nicePlot
figure(2)
plot(z-1e-6,H1,z-1e-6,H2)
nicePlot
