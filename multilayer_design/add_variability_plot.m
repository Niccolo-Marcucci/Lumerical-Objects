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

load(strcat("stability_analysis_design_",design_name,".mat"))

figure(1), hold on
for k=1:2
    media=mean(peak_thetas(:,k));
    delta=std(peak_thetas(:,k));
    idx=(peak_thetas(:,k)<media+2.5*delta) & ...
                (peak_thetas(:,k)>media-2.5*delta);
    d1=fitdist(peak_thetas(idx,k),'normal');
    d2=fitdist(peak_thetas(~idx,k),'normal');

    correction1=sum(idx)/length(idx);
    correction2=sum(~idx)/length(idx);
    p1=plot(theta,pdf(d1,theta)*correction1,...
                'displayname','Variability of the peak');

    text(d1.mu,pdf(d1,d1.mu)/2,...
                strcat(" \sigma=",string(d1.sigma)),'fontsize',14);
%     plot(theta,pdf(d2,theta)*correction2,'--','colors',...
%                         p1.Color,'displayname','Secondary peak');
end
nicePlot


folder="~/OneDrive/BSW/Pictures/";
name=strcat(folder,design_name,"_variabilityOfPeaks_onTheta");
stopBeforeSaving(name)
saveas(figure(1),name,'png')