clear
close all
addpath functions
design_name = "TM_guided3";
design_file = strcat("designs/design_",design_name,".mat");
load(design_file);

theta_v=0;
lambda_v=linspace(500,650,1e5)*1e-9;

N_periods=20;
d=zeros(N_periods*4+1,1);
n=zeros(N_periods*4+1,1);

lambda = 570e-9;
nA = n_eff1;
nB = n_eff3; % spacer dielectric
kA = 2*pi/lambda*nA;
kB = 2*pi/lambda*nB;
% FF = 0.5;
% period = pi/(kB*FF+kA*(1-FF));
%        = 2*pi/(kB+kA);
period = round(lambda/(nA+nB)*1e9)*1e-9;
db=period/2;
da=period/2;

SBwidth = lambda*4/pi*...  %Stop Band width
           asin(abs(n_eff1-n_eff3)/(n_eff1+n_eff3));

d(1)=1e-6;
d(2:2:end)=da;
d(3:2:end)=db;
d(end/2+1/2)=round(lambda/nB*1e9)*1e-9;
d(end)=1e-6;
n(1)=nB;
n(2:2:end)=nA;
n(3:2:end)=nB;
% n(end)=n_eff3;


[dr,nr,~,~] = prepare_multilayer(d,n);

R=zeros(length(lambda_v),length(theta_v));
tic
parfor i=1:length(lambda_v)
    R(i,:)=reflectivity(lambda_v(i),theta_v,dr,nr,'s');
end
toc

%% resonance
SBrange = (lambda_v > lambda-SBwidth/2)&(lambda_v < lambda+SBwidth/2);
R_tmp=R;
R_tmp(~SBrange,1)=1;
spectrum = 1-R_tmp;
[pks,idxs] = findpeaks(spectrum);
[peak,pk_ix] = max(pks);
idx = idxs(pk_ix);

min_s =     min(abs(spectrum(1:idx)- peak/2));
min_idx = find((abs(spectrum(1:idx)- peak/2) == min_s));
max_s =     min(abs(spectrum(idx:end)- peak/2));
max_idx = find((abs(spectrum(idx:end)- peak/2) == max_s)) + idx;

Q = 0;%lambda_v(idx)/abs(lambda_v(min_idx) -lambda_v(max_idx));

%% plot
plot(lambda_v*1e9,R)
title(strcat("Period ",string(period*1e9),"nm, Spacer ",...
                       string(d(end/2+1/2)*1e9),...
                       "nm, N^o of rings ", string(N_periods)))
legend(strcat("Resonance at ",...
              string(round(lambda_v(idx)*1e9)),...
              "nm, Q = ",string(round(Q))));
nicePlot
xlabel('wavelength [nm]')
ylabel('reflectance')
xlim([500 650])

% [z,nz]=field_distribution(570,0,d,n,'',5e2);
% figure
% plot(z,nz)
% nicePlot

folder="~/OneDrive/BSW/Pictures/";
name=strcat(folder,design_name,"_DBR_reflectance");
stopBeforeSaving(name)
saveas(figure(1),name,'png')