
    % Pull out common variables from the ephemeris matrix
    %======================================================================
    %toe = sat_ephem(jj,17);           % time of ephemeris
    dt  = dt_input(tt) - dt_ephem(jj); % seconds difference from epoch
    
    a   = 5153.65^2;           % semimajor axis, sqrt(a) = gps_ephem_all(:,5) (meters)
    ecc = 0;             % eccentricity
    n0  = sqrt(muE/a^3);               % nominal mean motion (rad/s)
    n   = n0 + sat_ephem(jj,3);        % corrected mean motion, delta_n = gps_ephem_all(:,3)
    
    
    M = 0.2e1; % OR 2?


    % Compute perigee, true and eccentric anomaly...
    %======================================================================

    % Load argument of perigee to a local variable and add perigee rate, rad
    perigee  = sat_ephem(jj,8); % + perigee_rate * dt;  

    % Compute Eccentric Anomaly, rad
    E    = mean2eccentric(M,ecc);
    cosE = cos(E);  
    sinE = sin(E);

    % Compute true anomaly, rad
    nu    = atan2( sqrt(1 - ecc*ecc).*sinE,  cosE-ecc ); 

    % Compute the argument of latitude, rad 
    u = nu + perigee;  % true anomaly + argument of perigee


    % Compute radius and inclination
    %======================================================================

    r   = a * (1 - ecc*cosE) ;                        % corrected radius  
    inc = 0.96 ;   %  inclination 
                                                               % i_dot = sat_ephem(jj,10)

    cosu = cos(u);    
    sinu = sin(u);  

    % Compute satellite position in orbital plane (Eq. 13)
    %======================================================================
    xo = r * cosu;    % satellite x-position in orbital plane
    yo = r * sinu;    % satellite y-position in orbital plane

    % Corrected longitude of ascending node for node rate and Earth rotation
    %======================================================================
    % Ascending node = ephem_all(jj,6)
    r_rate = 0;
    node = sat_ephem(jj,6) + (r_rate - wE)*dt -  (wE * sat_ephem(jj,17)); % Toe = gps_ephem_all(jj,17)

    % Calculate GPS Satellite Position in ECEF (m)
    %======================================================================
    cosi = cos(inc);    sini = sin(inc);
    coso = cos(node);   sino = sin(node);


    % Satellite position in ECEF (m)
    x(tt,1) = xo*coso - yo*cosi*sino  %x-position  

    x(tt,2) = xo*sino + yo*cosi*coso  %y-position 

    x(tt,3) = yo*sini                 %z-position
    

    % Keep track of health of each satellite
    %======================================================================      
    health(tt,1) = sat_ephem(jj,25); % satellite health (0.00 is useable)




