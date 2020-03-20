clear, close all
tic
theta_v = linspace(0,90,5e3);                    % angle vector
lambda_v=linspace(400,800,5e3)*1e-9;
addpath('functions');

load designs/design_TM_guided.mat
pol='p';

[d1,n1,~,~] = prepare_multilayer(d_layers,idx_layers);

R=zeros(length(lambda_v),length(theta_v));
parfor i=1:length(lambda_v)
    lambda=lambda_v(i);
    [R(i,:),~]=reflectivity(lambda,theta_v,d1,n1,pol);
end

%%
[theta,lambda]=meshgrid(theta_v,lambda_v);

beta=sin(theta/180*pi)*n1(1)*2*pi./lambda*(d1(3)-d1(2));
c=physconst('lightspeed');
figure('units','normalized','outerposition',[0 0 1 1])
s=surf(beta,2*pi*c./lambda,1-R);
ylabel('angular frequency')
xlabel('beta')
s.EdgeColor='none';
colormap('hot')
colorbar
view(2)
whitebg(figure(1),'black'); 

set(gcf, 'InvertHardCopy', 'off');
folder="~/OneDrive/BSW/Pictures/";
name=strcat(folder,"TM_guided_dispersion_beta");
saveas(figure(1),name,'png')

figure('units','normalized','outerposition',[0 0 1 1])
s=surf(theta,lambda,1-R);
ylabel('wavelength')
xlabel('Angle of incidence')
s.EdgeColor='none';
colormap('hot')
colorbar
view(2)
whitebg(figure(2),'black'); 

set(gcf, 'InvertHardCopy', 'off');
folder="~/OneDrive/BSW/Pictures/";
name=strcat(folder,"TM_guided_dispersion_theta");
saveas(figure(2),name,'png')
%%
toc