
% fname = "DATA\brdc2450.20n";
% 
% [gps_ephem,ionoparams] = read_clean_GPSbroadcast(fname,true);


fname = "DATA\nist2450.20o";
rinex_file = read_rinex_obs8(fname, 1, 1)