% 
% 
function T = Tij(n_i,n_j,beta,pol)
    costheta_i = sqrt(n_i.^2-beta^2)./n_i;  
    costheta_j = sqrt(n_j.^2-beta^2)./n_j;  
    if pol == 's'
        rij = (n_i*costheta_i-n_j*costheta_j)./...
              (n_i*costheta_i+n_j*costheta_j);
        rji = -rij;
        tji =  rji + 1;
    elseif pol == 'p'
        rij = (n_j*costheta_i-n_i*costheta_j)./...
              (n_j*costheta_i+n_i*costheta_j);
        rji = -rij;
        tji = (rji + 1)*n_j/n_i;
    else 
        error("Invalid Polarization. Valid options are 's' or 'p'")
    end

    T = 1/tji* [ 1 , rji ;
                rji,  1  ];  
    
end