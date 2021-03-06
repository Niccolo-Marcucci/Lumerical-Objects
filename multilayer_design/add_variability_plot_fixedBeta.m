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

% to be only run after running 'evaluate_design_quality.m'

load(strcat("stability_analysis_onBeta_design",design_name,".mat"))
peak_lambdas=peak_lambdas*1e9;

figure(1), hold on
for k=1:2
    media=mean(peak_lambdas(:,k));
    delta=std(peak_lambdas(:,k));
    idx=(peak_lambdas(:,k)<media+2.5*delta) & ...
                            (peak_lambdas(:,k)>media-2.5*delta);
    d1=fitdist(peak_lambdas(idx,k),'normal');
    d2=fitdist(peak_lambdas(~idx,k),'normal');

    correction1=3;%sum(idx)/length(idx);
    correction2=1;%sum(~idx)/length(idx);
    p1=plot(lambda*1e9,pdf(d1,lambda*1e9)*correction1,...
                'displayname','Variability of the peak');

    text(d1.mu,pdf(d1,d1.mu)*correction1/2,...
        strcat(" \sigma=",string(d1.sigma),"nm"),'fontsize',14);
%     plot(theta,pdf(d2,theta)*correction2,'--',...
%                 'colors',p1.Color,'displayname','Secondary peak');
end
nicePlot


folder="~/OneDrive/BSW/Pictures/";
name=strcat(folder,design_name,"_variabilityOfPeaks_onLambda");
stopBeforeSaving(name)
saveas(figure(1),name,'png')