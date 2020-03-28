clear
close all
addpath('functions');

design_name = "TM_gd3_buriedDBR";
design_file = strcat("designs/design_",design_name,".mat");
load(design_file);
pol='p';

d1=d_layers;
n1=idx_layers;
d1(end)=1e-6;


lambda=570e-9;
theta = linspace(40,70,1e4);  

for k = 1:2
    n = n1;
    d = d1;
    n(end-2)=n1(end-k);
    
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
    plot(z1,P)

end
    
n_eff1/n_eff3

[d1,n1] = prepare_multilayer(d1,n1);

[z1, nz] = field_distribution(lambda,theta(idx),d1,n1);
figure(2);
plot(z1,real(nz-1)*400  );
xlabel('z [m]')
ylabel('Power density')
legend('Field without buried layer', 'Field with buried layer',...
                                        '(n_z -1) x 400');
nicePlot

figure(1);
% ylim([0.2,1])
xlabel('Incidence angle [degrees]')
ylabel('Reflectivity')
legend('Without buried layer', 'With buried layer');
nicePlot
% set(gca,'yscale','log')
% error

folder="~/OneDrive/BSW/Pictures/";
name=strcat(folder,design_name,"_Field_distribution");
stopBeforeSaving(name)
saveas(figure(2),name,'png')
name=strcat(folder,design_name,"_BWS_lines");
saveas(figure(1),name,'png')