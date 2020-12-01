% fname = "brdc2450.20n";
fname = "nist2450.20n";

% answer_x = [-15236.226182   -738.304955  21546.270275  ];

seconds_in_day = 60*60*24;
tow = seconds_in_day * 2; % Because we are on Tuesday

[ephem_all,ionoparams] = read_clean_GPSbroadcast(fname,true);

t = [0:1:24*60]*60+172800;  %Time of week in seconds with 1 min spacing
WN = 73*ones(size(t));
t_input = [WN;t].';            %[nx2]

how = mod(t, 172800)/3600; % hours in week

PRN = ephem_all(:,1);

% for pp = 1:length(PRN)
pp = 1;

% This is our prn that we are researching
prn = PRN(pp);

% Now we need to grab position by following the slides

muE = 3.986005e14;     % WGS-84 value, m^3/s^2
wE  = 7.2921151467e-5; % WGS-84 value, rad/s
c   = 2.99792458e8;    % GPS acceptd speed of light, m/s

sz         = size(t_input,1);
x          = ones(sz,3) * NaN;
health     = ones(sz,1) * NaN;

% Pull out ephemerides for PRN in question
kk  = find(ephem_all(:,1) == prn);  % kk is vector containing row numbers of ephem_all that are for sat.no. 'index'
sat_ephem = ephem_all(kk,:);        % sat_ephem is matrix of all ephem data for each entry of sat.no. 'index'

if isempty(kk),return,end

% Compute elapsed times of each ephemeris epoch wrt first entry, seconds
dt_ephem = (sat_ephem(:,19) - sat_ephem(1,19))*604800 + (sat_ephem(:,17) - sat_ephem(1,17)); % This is every two hours in seconds


% Compute elapsed times of each input time wrt first ephemeris entry, seconds
dt_input = (t_input(:,1) - sat_ephem(1,19))*604800 + (t_input(:,2) - sat_ephem(1,17));


% Start main calc loop
% loop through all input times (we'll start with first one)
for tt = 1:sz % loop through all input times
    
    % Pull out most recent ephemeris values
    % jj = max( find(dt_input(tt) >= dt_ephem) ); % sat_ephem(:,17) = toe (sec into GPS week) of each entry
    % jj = row of specific sat. ephem. data with epoch closest to input time
    
    % Pull out nearest ephemeris values
    [mn,jj] = min(abs( dt_input(tt) - dt_ephem ));
    
    % if isempty(jj),continue,end  % no matching ephemeris time found. continue to next input time
    
    
    % Pull out common variables from the ephemeris matrix
    %======================================================================
    toe = sat_ephem(jj,17);           % time of ephemeris
    dt  = t(tt) - toe; % seconds difference from epoch % tow in sec - toe
    
    a   = sat_ephem(jj,5)^2;           % semimajor axis, sqrt(a) = gps_ephem_all(:,5) (meters)
    ecc = sat_ephem(jj,4);             % eccentricity
    n0  = sqrt(muE/a^3);               % nominal mean motion (rad/s)
    n   = n0 + sat_ephem(jj,3);        % corrected mean motion, delta_n = gps_ephem_all(:,3)
    M   = sat_ephem(jj,2) + n*dt;      % mean anomaly, M0 = gps_ephem_all(:,2)
    
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
    
    r   = a * (1 - ecc*cosE) ; % radius
    inc = sat_ephem(jj,7) ;   %  inclination
    i_dot = sat_ephem(jj,10); % for correcting inclination
    
    cosu = cos(u);
    sinu = sin(u);
    
    % Add Corrections
    %======================================================================
    Cuc = ephem_all(jj,11);
    Cus = ephem_all(jj,12);
    Crc = ephem_all(jj,13);
    Crs = ephem_all(jj,14);
    Cic = ephem_all(jj,15);
    Cis = ephem_all(jj,16);
    
    % Corrections
    delta_uj = Cus*sin(2*u)+Cuc*cos(2*u); % Argument of latitude correction
    delta_rj = Crs*sin(2*u)+Crc*cos(2*u); % Argument of latitude correction
    delta_ij = Cis*sin(2*u)+Cic*cos(2*u); % Argument of latitude correction
    
    u = u + delta_uj;
    r = r + delta_rj;
    inc = inc + delta_ij+i_dot*dt;
    
    % Compute satellite position in orbital plane (Eq. 13)
    %======================================================================
    xo = r * cosu;    % satellite x-position in orbital plane
    yo = r * sinu;    % satellite y-position in orbital plane
    
    % Corrected longitude of ascending node for node rate and Earth rotation
    %======================================================================
    % Ascending node = ephem_all(jj,6)
    node = sat_ephem(jj,6) + (sat_ephem(jj,9) - wE)*dt -  (wE * sat_ephem(jj,17)); % Toe = gps_ephem_all(jj,17)
    
    % Calculate GPS Satellite Position in ECEF (m)
    %======================================================================
    cosi = cos(inc);    sini = sin(inc);
    coso = cos(node);   sino = sin(node);
    
    % Satellite position in ECEF (m)
    x(tt,1) = xo*coso - yo*cosi*sino;  %x-position
    
    x(tt,2) = xo*sino + yo*cosi*coso;  %y-position
    
    x(tt,3) = yo*sini;                 %z-position
    
    sprintf('Output: %0.5g %0.5g %0.5g',[x(tt,1)/1000, x(tt,2)/1000, x(tt,3)/1000])
%     sprintf('Answer: %0.5g %0.5g %0.5g',[answer_x(tt,1), answer_x(tt,2), answer_x(tt,3)])
    
    % Keep track of health of each satellite
    %======================================================================
    health(tt,1) = sat_ephem(jj,25); % satellite health (0.00 is useable)
    
end
