clear, 
% close all

theta_v = [linspace(0,90,2e3)];% , asin(linspace(1,2,2e3))/pi*180
lambda_v=linspace(400,800,2e3)*1e-9;
addpath('functions');

design_name = "TM_guidedDesc";
design_file = strcat("designs/design_",design_name,".mat");
load(design_file)
pol='p';

[d,n,~,~] = prepare_multilayer(d_layers,...
                             idx_layers);%(end:-1:1)

R=zeros(length(lambda_v),length(theta_v));
parfor i=1:length(lambda_v)
    lambda=lambda_v(i);
    [R(i,:),~]=reflectivity(lambda,theta_v,d,n,pol);
end
R(R>1)=1;
%%
[theta,lambda]=meshgrid(theta_v,lambda_v);

beta=real(sin(theta/180*pi))*2*pi*n(1)./lambda;
c=physconst('lightspeed');
f=figure('units','normalized','outerposition',[0 0 1 1]);
omega=2*pi*c./lambda;
s=surf(beta,omega,(1-R));
ylabel('angular frequency')
xlabel('beta')
s.EdgeColor='none';
colormap('hot')
colorbar;
view(2)
whitebg(f,'black');


% the air line is beta*c/(sin(theta_crit)*n_sub), but the denominator
% is equal to sin(pi/2)*n_air by def.

n_eq = real((n(2)*n(3))/sqrt(n(2)^2 + n(3)^2));
theta_crit = real(asin(n(end)/n(1)));
theta_brew =      asin( n_eq /n(1));
% theta_brew*180/pi 
% theta_crit*180/pi

beta = linspace(0,beta(1,end),2);
omega_crit=c*beta/(sin(theta_crit)*n(1));
% omega_brew=c*beta/(sin(theta_brew)*n(1));
omega_brew=c*beta/n_eq;

hold on
plot3(beta,omega_crit,[1 1],'--w','linewidth',3)
if pol == 'p'
    plot3(beta,omega_brew,[1 1],'--g','linewidth',3)
end
ylim([omega(end,end) omega(1,1)])
xlim([0 2.5e7])
%%
set(gcf, 'InvertHardCopy', 'off');
folder="~/OneDrive/BSW/Pictures/";
name=strcat(folder,design_name,"_dispersion_beta");
stopBeforeSaving(name)
saveas(figure(1),name,'png')

theta2=theta_v;
theta2(imag(theta2)~=0)=90+abs(imag(theta2(imag(theta2)~=0)));
[theta,lambda]=meshgrid(real(theta2),lambda_v);

figure('units','normalized','outerposition',[0 0 1 1])
s=surf(theta,lambda,1-R);
ylabel('wavelength')
xlabel('Angle of incidence')
s.EdgeColor='none';
colormap('hot')
colorbar
view(2)
whitebg(figure(2),'black'); 

hold on
plot3(theta_crit *[1 1]*180/pi,[lambda_v(1),lambda_v(end)],...
                                        [1 1],'--w','linewidth',3)
if pol == 'p'
   plot3(real(theta_brew)*[1 1]*180/pi,[lambda_v(1),lambda_v(end)],...
                                        [1 1],'--g','linewidth',3)
end
set(gcf, 'InvertHardCopy', 'off');
folder="~/OneDrive/BSW/Pictures/";
name=strcat(folder,design_name,"_dispersion_theta");
saveas(figure(2),name,'png')
%%
