function [ wn,tow ] = cal2gps(ymd )
% function [ wn,tow ] = cal2gps(ymd )
% INPUT  
%   [ymd] = array containing the date as year, month, and day
% OUTPUT 
%   wn = Full GPS week number (not modulo 1024)
%   tow =  time of week in seconds
% Notes: intended for use in integer dates not guaranteed to preserve
% sub-second values.
%
jd = cal2jd(ymd(:,1),ymd(:,2),ymd(:,3));
[wn,tow]=jd2gps(jd);
end

