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
design_name = "TM_gd3_buriedDBR";
design_file = strcat("designs/design_",design_name,".mat");
load(design_file);
pol='p';                            % polarisation: 'p' or 's'
design_type='buried';               % either 'buried' or empty

lambda=570e-9;
theta = linspace(40,70,1e4);

d1=d_layers;
n1=idx_layers;

[d1,n1] = prepare_multilayer(d1,n1);

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
stopBeforeSaving(name)
saveas(figure(2),name,'png')
name=strcat(folder,design_name,"_BWS_lines");
saveas(figure(1),name,'png')