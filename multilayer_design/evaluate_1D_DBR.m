clear
close all
addpath functions
design_file = 'designs/design_TM_guided.mat';
load(design_file);

theta_v=10;
lambda_v=linspace(400,700,1e5)*1e-9;
d=zeros(201,1);
n=zeros(201,1);

kB = 2*pi/570e-9*n_eff3;
kA = 2*pi/570e-9*n_eff1;
FF = 0.5;
period = 187e-9;%2* pi/(kB*FF+kA*(1-FF));
da=period*FF;
db=period*(1-FF);

d(1)=1e-6;
d(2:2:end)=da;
d(3:2:end)=db;
d(end/2+1/2)=175e-9;
d(end)=1e-6;
n(1)=n_eff3;
n(2:2:end)=n_eff1;
n(3:2:end)=n_eff3;
% n(end)=n_eff3;


[dr,nr,~,~] = prepare_multilayer(d,n);

R=zeros(length(lambda_v),length(theta_v));
for i=1:length(lambda_v)
    lambda=lambda_v(i);
    R(i,:)=reflectivity(lambda,theta_v,dr,nr,'s');
end
plot(lambda_v*1e9,R)
nicePlot
xlabel('wavelength [nm]')
ylabel('reflectance')
xlim([400 700])

[z,nz]=field_distribution(570,0,d,n);
figure
plot(z,nz)
nicePlot

folder="~/OneDrive/BSW/Pictures/";
name=strcat(folder,"TM_guided_DBR_reflectance");
% stopBeforeSaving(name)
% saveas(figure(1),name,'png')