% Copyright 2020 Niccolò Marcucci <niccolo.marcucci@polito.it>
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
%
%     http://www.apache.org/licenses/LICENSE-2.0
%
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all
addpath('functions');


% design general properties
% design_name = "TM_SiO2TiO2_532_N9";
design_name = "TM_TiO2SiO2_202603";
design_file = strcat("designs/design_",design_name,".mat");
load(design_file)
pol='p';                            % polarisation: 'p' or 's'

% idx_layers(end-1) = 1.57+1i*1e-2;
d_layers(end-1) =0e-9; 


theta_v = [linspace(43,60,2e3)];% , asin(linspace(1,2,2e3))/pi*180];%
lambda_v=linspace(400,850,2e3)*1e-9;

lambdaDBR=265e-9;                   % periodicity of the 2D DBR 
                                    % inscribed on the last layer
lambda=600e-9;                      % wavelength of interest

[d,n,~,~] = prepare_multilayer(d_layers,...
                             idx_layers);

n_=disperse_indices(n,lambda_v);

R=zeros(length(lambda_v),length(theta_v));
parfor i=1:length(lambda_v)
    [R(i,:),~]=reflectivity(lambda_v(i),theta_v,d,n_(:,i),pol);
end
R(R>1)=1;                           % necessary when sampling the 
                                    % dispersion using imaginary 
                                    % angles
%% dispersion in beta
[theta_m,lambda_m]=meshgrid(theta_v,lambda_v);

beta=real(sin(theta_m/180*pi))*2*pi*n(1,:)./lambda_m;
c=physconst('lightspeed');
% {
f=figure('units','normalized','outerposition',[0 0 1 1]);
omega=2*pi*c./lambda_m;
s=surf(beta,omega,(1-R));
ylabel('Angular frequency \omega [rad/s]','fontsize',18)
xlabel('\beta [rad/m]','fontsize',18)
s.EdgeColor='none';
colormap('parula')
colorbar;
view(2)
whitebg(f,'black');
set(gcf, 'InvertHardCopy', 'off');          % for maintaing dark 
                                            % background after saving
%}
%%
% {
n_eff=sin(theta_m/180*pi)*n(1,:);
f=figure('units','normalized','outerposition',[0 0 1 1]);
s=surf(n_eff,lambda_m,(R));
ylabel('Wavelength - nm','fontsize',18)
xlabel('Effective index','fontsize',18)
s.EdgeColor='none';
colormap('parula')
colorbar;
view(2)
yline(lambda,'--m','linewidth',3)
% whitebg(f,'black');
% set(gcf, 'InvertHardCopy', 'off');          % for maintaing dark 
                                            % background after saving
%}
beta = linspace(0,2.5e7,2);

n_eq = real((n(2)*n(3))/sqrt(n(2)^2 + n(3)^2));
theta_brew =      asin( n_eq /n(1));
theta_crit = real(asin(n(end)/n(1)));

omega_crit=c*beta/(sin(theta_crit)*n(1));   % air line
omega_brew=c*beta/n_eq;                     % Brewster line
omega_silica=c*beta/1.46;                   % Silica substrate line
%{
hold on
plot3(beta,omega_crit,[1 1],'--w','linewidth',3)
if pol == 'p'
    plot3(beta,omega_brew,[1 1],'--g','linewidth',3)
end
plot3(beta,2*pi*c/lambda*[1 1],[1 1],'--m','linewidth',2)
plot3(pi/lambdaDBR*[1 1],omega_brew,[1,1],'--y','linewidth',2)
plot3(beta,omega_silica,[1,1],'--r','linewidth',3)

ylim([omega(end,end) omega(1,1)])
xlim([0 2.5e7])

folder="~/OneDrive/BSW/Pictures/";
name=strcat(folder,design_name,"_dispersion_beta");
% stopBeforeSaving(name)
% saveas(figure(1),name,'png')
%}
%% dispersion in theta
theta2=theta_v;
% plot imaginary angles at values > 90
theta2(imag(theta2)~=0)=90+abs(imag(theta2(imag(theta2)~=0)));
[theta_m,lambda_m]=meshgrid(real(theta2),lambda_v);
n_eff = sin(theta_m*pi/180)*n(1,:);

f = figure('units','normalized','outerposition',[0 0 1 1]);
imagesc(theta_v,lambda_v,(R));
set(gca,'ydir','norm')
% s=surf(theta_m,lambda_m,log(1-R));
% s.EdgeColor='none';
ylabel('Wavelength')
xlabel('Angle of incidence')
colormap('parula')
colorbar
% view(2)
whitebg(f,'black'   ); 
set(gcf, 'InvertHardCopy', 'off');

hold on
plot3(theta_crit *[1 1]*180/pi,[lambda_v(1),lambda_v(end)],...
                                        [1 1],'--w','linewidth',3)
if pol == 'p'
   plot3(real(theta_brew)*[1 1]*180/pi,[lambda_v(1),lambda_v(end)],...
                                        [1 1],'--g','linewidth',3)
end
% plot3([theta_v(1) theta_v(end)],lambda*[1 1],[1 1],'--m')
yline(lambda,'--m','linewidth',3)

folder="~/OneDrive/BSW/Pictures/";
name=strcat(folder,design_name,"_dispersion_theta");
% saveas(figure(2),name,'png')
%%

function str = idx2str(value)
    x = real(value);
    if (1.44 <= x) &&  (x <= 1.47)
        str = "SiO_2";
    elseif (1.6 <= x) &&  (x <= 1.66)
        str = "Al_2O_3";
    elseif (2.3 <= x) &&  (x <= 2.6)
        str = "TiO_2";
    elseif x == 1.48
        str = "PMMA";
    elseif (1.95 <= x) &&  (x <= 2.2)
        str = "Ta_2O_5";
    elseif x == 1
        str = "Air";
    else
        str = strcat("fake",string(real(value)),"");
    end
end

function n_out=disperse_indices(n,lambda)
    n_out = ones(length(n),length(lambda));
    n_unique = unique(n);
    for i = 1:length(n_unique)
        switch idx2str(n_unique(i))
            case "SiO_2"
                data = readtable("./Dispersioni/SiO2 150deg PEALD 20260312.txt");
                n_SiO2 = spline(data.nm, data.n, lambda*1e9);
                if imag(n_unique(i)) == 0
                    n_SiO2 = real(n_SiO2);
                end
            case "Al_2O_3"
                data = readtable("./Dispersioni/Al2O3 150deg plasma.txt");
                n_Al2O3 = spline(data.nm, data.n, lambda*1e9);
            case "TiO_2"
                data = readtable("./Dispersioni/TiO2 150deg PEALD 20260312.txt");
                n_TiO2 = spline(data.nm, data.n, lambda*1e9);
            case "PMMA"
                n_PMMA = real(n_unique(i))*ones(1,length(lambda));
            case "Ta_2O_5"
                n_Ta_2O_5 = real(n_unique(i))*ones(1,length(lambda));
            case "Air"
                n_Air = real(n_unique(i))*ones(1,length(lambda));
        end
    end
    for i = 1:length(n)
        switch idx2str(n(i))
            case "SiO_2"
                n_out(i,:) = n_SiO2 + 1i*imag(n(i));
            case "Al_2O_3"
                n_out(i,:) = n_Al2O3 + 1i*imag(n(i));
            case "TiO_2"
                n_out(i,:) = n_TiO2 + 1i*imag(n(i));
            case "PMMA"
                n_out(i,:) = n_PMMA + 1i*imag(n(i));
            case "Ta_2O_5"
                n_out(i,:) = n_Ta_2O_5 + 1i*imag(n(i));
            case "Air"
                n_out(i,:) = n_Air + 1i*imag(n(i));
        end
    end
end