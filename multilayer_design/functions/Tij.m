% 
% 
function T = Tij(n_i,n_j,theta_i,theta_j,pol)
    if pol == 's'
        rij = (n_i*cos(theta_i)-n_j*cos(theta_j))./...
              (n_i*cos(theta_i)+n_j*cos(theta_j));
        rji = -rij;
        tji =  rji + 1;
    elseif pol == 'p'
        rij = (n_j*cos(theta_i)-n_i*cos(theta_j))./...
              (n_j*cos(theta_i)+n_i*cos(theta_j));
        rji = -rij;
        tji = (rji + 1)*n_j/n_i;
    else 
        error("Invalid Polarization. Valid options are 's' or 'p'")
    end

    T = 1/tji* [ 1 , rji ;
                rji,  1  ];  
    
end