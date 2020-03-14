clear, close all
tic
theta_v = linspace(0,90,5e3);                    % angle vector
lambda_v=linspace(460,630,1e3)*1e-9;
addpath('functions');

load best_design_TE.mat
[d1,n1,~,~] = prepare_multilayer(d_layers,idx_layers);

R=zeros(length(lambda_v),length(theta_v));
parfor i=1:length(lambda_v)
    lambda=lambda_v(i);
    [R(i,:),~]=reflectivity(lambda,theta_v,d1,n1,'s');
end

[theta,lambda]=meshgrid(theta_v,lambda_v);
% R_p=squeeze(RT1.Rp);
beta=sin(theta/180*pi)*n1(1)*2*pi./lambda*(d1(3)-d1(2));
c=physconst('lightspeed');
figure
s=surf(beta,2*pi*c./lambda,1-R);
% s=surf(theta,lambda,1-R);
s.EdgeColor='none';
colormap('hot')
colorbar
view(2)
toc