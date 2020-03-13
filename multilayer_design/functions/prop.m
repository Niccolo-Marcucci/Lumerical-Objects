% 
% 
function [Pr,Pl] = prop(k,d,costheta)
    kz=k*costheta; 
    Pr = exp(+1i*kz*d);  
    Pl = exp(-1i*kz*d);
end