% fname = "brdc2450.20n";
fname = "nist2450.20n";
% fname = "BRDX00DLR_S_20200570000_01D_MN.rnx.gz";

seconds_in_day = 60*60*24;
tow = seconds_in_day * 2; % Because we are on Tuesday

[gps_ephem,ionoparams] = read_clean_GPSbroadcast(fname,true);

t = [0:1:24*60]*60+172800;  %Time of week in seconds with 1 min spacing 
WN = 73*ones(size(t));
T_in = [WN;t].';            %[nx2] 

PRN = gps_ephem(:,1);
num_sat = length(PRN);
for kk = 1:num_sat
    [health,pos] = broadcast2pos(gps_ephem,T_in,PRN(kk));
%     POS_ecef(:,:,kk) = pos;     %n x 3 x num_sat
%     [x,y,z] = ECEF2llh(pos(1,:))
    pos(1,:)
    [lat,lon,alt] = ECEF2llh(pos(1,:))
end

%Problem 1 functions
function [lat,lon,h] = ECEF2llh(xyz)
R = 6378.137e3;             %SemiMajor axis of ellipsoid (meters)
f = 1/298.257223563;% flattening parameter of ellipsoid
e2 = 2*f-f^2;       %square of eccentricity of ellipsoid. 

lambda = atan2d(xyz(2),xyz(1));     %Longitude (deg)
p = sqrt(xyz(1)^2+xyz(2)^2);        
r =sqrt(sum(xyz.^2));

phi_gd = asind(xyz(3)/r);   %initial guess for geodetic latitude (deg)
tolerance = 1e-8;    
dif = 1;            
while dif > tolerance %iterate until delta is small
    
    C = R/sqrt(1-e2*sind(phi_gd)^2); %Radius of curvature in the meridian
    tan_phi_gd = (xyz(3)+C*e2*sind(phi_gd))/p;
    phi_gd_new = atand(tan_phi_gd);   %Geodetic latitude (deg)
    dif = abs(phi_gd_new - phi_gd);      %difference between guess and new
    phi_gd = phi_gd_new;            %set old to new
end

h = p/cosd(phi_gd)-C;       %ellipsoidal height (meters) above surface
lat = phi_gd;               %latitude (deg)
lon = lambda;               %longitude (deg)

end