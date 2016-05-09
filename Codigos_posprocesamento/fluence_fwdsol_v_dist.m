function phi = fluence_fwdsol_v_dist(mua, mus, g, r, n)

    %%%% Analytical Solution for TRS Semi-Infinite Geometry
    
    c= 2.99792458e10;

    % Speed of light in medium
    v = c/n; 
    
    musp = (1-g)*mus;
    D = v/(3*(musp+mua));
    
    % Place source 1/musp into the medium.
    z_tr = 1/(musp+mua);
    zb = 2/(3*(musp+mua));

    % Fluence as a function of position:
    phi = v ./ (4*pi*D) .* ...
          (exp(-sqrt(3*musp*mua)*sqrt(z_tr^2 + r.^2))   ./    sqrt(z_tr^2 + r.^2) - ...
           exp(-sqrt(3*musp*mua)*sqrt((z_tr+2*zb)^2 + r.^2)) ./ sqrt((z_tr+2*zb)^2 + r.^2));
    

