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

design_name = "TM_TiO2SiO2_202603";
design_file = strcat("designs/design_",design_name,".mat");
load(design_file);
pol='p';                            % polarisation: 'p' or 's'
design_type='';               % either 'buried' or empty
lambda_DBR = inf;                % determines the beta at which the
                                    % reflectivity is computed (see
                                    % usage).
d1=d_layers;
n1=idx_layers;

beta = pi/lambda_DBR;
lambda = linspace(400,1000,1e4)*1e-9;
K = 2*pi./lambda*n1(1);
theta = asin(beta./K)/pi*180;

for k = 2
    if strcmp(design_type,'buried')
        n = n1;
        d = d1;
        n(end-2)=n1(end-k);
    else
        n = [n1(1:end-k) ; n1(end)];
        d = [d1(1:end-k) ; d1(end)];
    end
    
    R = zeros(1,length(lambda));
    r = zeros(1,length(lambda));
     
    [dr,nr,~,~] = prepare_multilayer(d,n);

    nr_=disperse_indices(nr,lambda);
    for i=1:length(lambda)
        [R(i),r(i)] = reflectivity(lambda(i),theta(i),dr,nr_(:,i),pol);
    end
    [pks,idxs] = findpeaks(1-R);
    [~,pk_ix] = max(pks);
    idx = idxs(pk_ix);
%     n_eff(k)= sin(theta(idx)*pi/180)*nr(1);

    
    figure(1)
    hold on
    plot(lambda*1e9,R)

    [z1, ~, P ] = field_distribution(lambda(idx),theta(idx),d,n,pol);
    figure(2)
    hold on
    plot(z1,P)

end

[d1,n1] = prepare_multilayer(d1,n1);

[z1, nz] = field_distribution(lambda(idx),theta(idx),d1,n1);
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
nicePlot

figure(1);
% plot(570*[1 1],[1e-3 1],'--k')
xlabel('wavelength [nm]')
ylabel('Reflectivity')
if strcmp(design_type,'buried')
    legend('Without buried layer','With buried layer');
else
    legend('With last layer','Without last layer');
end
nicePlot
% set(gca,'yscale','log')

folder="~/OneDrive/BSW/Pictures/";
name=strcat(folder,design_name,"_Field_distribution");
stopBeforeSaving(name)
saveas(figure(2),name,'png')
name=strcat(folder,design_name,"_BWS_lines");
saveas(figure(1),name,'png')


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