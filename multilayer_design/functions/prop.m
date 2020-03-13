% 
% 
function P = prop(k,d,costheta)
    kz=k*costheta; 
    P = [exp(+1i*kz*d), 0   ;  
          0 , exp(-1i*kz*d) ];
end