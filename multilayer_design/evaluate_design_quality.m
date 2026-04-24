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
load(design_file);
pol='p';                            % polarisation: 'p' or 's'
design_type='';               % either 'buried' or empty

lambda=600e-9;
theta = linspace(43,60,1e5);

% idx_layers(end-1) = 1.60+1i*1e-4;
%d_layers(end-1) = 30e-9; 
% d_layers(2:9) = 0e-9; 

d1=d_layers;
n1=idx_layers;

[d1,n1] = prepare_multilayer(d1,n1);
n1 = disperse_indices(n1,lambda);

for k = 1:2
    if strcmp(design_type,'buried')
        n = n1;
        d = d1;
        n(end-2)=n1(end-k);
    else
        n = [n1(1:end-k) ; n1(end)];
        d = [d1(1:end-k) ; d1(end)];
    end
    
    [dr,nr,~,~] = prepare_multilayer(d,n);
    
    [R,r] = reflectivity(lambda,theta,dr,nr,pol);
    [pks,idxs] = findpeaks(1-R);
    [~,pk_ix] = max(pks);
    idx = idxs(pk_ix);
    
    n_eff= sin(theta(idx)*pi/180)*n(1);
    
    figure(1)
    hold on
    plot(theta,R)
    text(theta(idx),0.5+R(idx)/2,...
            strcat(" n_{eff}=",string(n_eff)),'fontsize',14);
        
    [z1, ~,P ] = field_distribution(lambda,theta(idx),d,n,pol);
    figure(2)
    hold on
    plot(z1,P)
end

[d1,n1] = prepare_multilayer(d1,n1);

[z1, nz] = field_distribution(lambda,theta(idx),d1,n1,'',1e3);
figure(2);
plot(z1,real(nz-1)*400  );
xlabel('z [m]')
ylabel('Power density')
if strcmp(design_type,'buried')
    legend('Field without buried layer','Field with buried layer',...
                                                '(n_z -1) x 400');
else
    legend('Field with last layer','Field without last layer',...
                                                '(n_z -1) x 400');
end
title(create_design_title(design_file))
nicePlot

figure(1);
ylim([0.2,1])
xlabel('Incidence angle [degrees]')
ylabel('Reflectivity')
if strcmp(design_type,'buried')
    legend('Without buried layer','With buried layer');
else
    legend('With last layer','Without last layer');
end
title(create_design_title(design_file))
nicePlot
% set(gca,'yscale','log')

folder="~/OneDrive/BSW/Pictures/";
name=strcat(folder,design_name,"_Field_distribution");
% stopBeforeSaving(name)
% saveas(figure(2),name,'png')
% name=strcat(folder,design_name,"_BWS_lines");
% saveas(figure(1),name,'png')

%%

function str = idx2str(value)
    x = real(value);
    if (1.441 <= x) &&  (x <= 1.47)
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