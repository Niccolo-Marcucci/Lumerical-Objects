clear
close all
addpath('functions');

design_name = "TM_gd3_buriedDBR";
design_file = strcat("designs/design_",design_name,".mat");
load(design_file);
pol='p';                            % polarisation: 'p' or 's'
design_type='buried';               % either 'buried' or empty
lambda_DBR = 165e-9;                % determines the beta at which the
                                    % reflectivity is computed (see
                                    % usage).
d1=d_layers;
n1=idx_layers;

beta = pi/lambda_DBR;
lambda = linspace(500,640,1e4)*1e-9;
K = 2*pi./lambda*n1(1);
theta = asin(beta./K)/pi*180;

for k = 1:2
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
    for i=1:length(lambda)
        [R(i),r(i)] = reflectivity(lambda(i),theta(i),dr,nr,pol);
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
plot(570*[1 1],[1e-3 1],'--k')
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