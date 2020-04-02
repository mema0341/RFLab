%% Problem 2
% ccc % Destroy the evidence
theta = 0:0.001:360; theta = deg2rad(theta); % Define Theta

%% Problem 2a Half Power Beam Width
y = cos(theta) .* cos(2*theta);

hpbw = 0.5;
% hpbw = sqrt(2)/2;

tmp = abs(y-hpbw);
tmp1 = (tmp<0.0001);

% Find the first value where theta = 0.5
tmpx = theta(tmp1);
tmpy = y(tmp1);

th=deg2rad(1:36:480);
specline = ones(size(th))*sqrt(2)/2;

figure
% Plot y, HPBW line, and intersection of the two
plot(theta,y,'LineWidth',2); hold on
plot(th, specline, '--r', 'LineWidth',2)
plot(tmpx, tmpy, '*b', 'LineWidth',4)


answer2a_hpbw_rad = tmpx(1)
answer2a_hpbw_deg = rad2deg(tmpx(1))

% Problem 2a First Null Beam Width
npbw = 0;

% Find the first value where theta = 0.5
tmp = abs(y-npbw);
tmp1 = (tmp<0.0001);

tmpx = theta(tmp1);
tmpy = y(tmp1);

plot(tmpx, tmpy, '*g', 'LineWidth',4)
axis([0, 2*pi, 0, 1])

answer2a_npbw_rad = tmpx(1)
answer2a_npbw_deg = rad2deg(tmpx(1))

% Make Graph Pretty
title({['Probleam 2a: cos(\theta)cos(2\theta)'];['HPBW = ',num2str(round(answer2a_hpbw_rad*2,2)),' [rad], ', num2str(round(answer2a_hpbw_deg*2,2)), ' [deg]']; ['NPBW = ',num2str(round(answer2a_npbw_rad,2)),' [rad], ', num2str(round(answer2a_npbw_deg,2)), ' [deg]']})
xlabel('Theta [radians]')
ylabel('Normalized Radiation Intensity')
axis([0, 2*pi, -1, 1])

%% Problem 2c Half Power Beam Width
y = cos(theta) .* cos(3*theta);

hpbw = 0.5;

tmp = abs(y-hpbw);
tmp1 = (tmp<0.0001);

% Find the first value where theta = 0.5
tmpx = theta(tmp1);
tmpy = y(tmp1);

th=deg2rad(1:36:480);
specline = ones(size(th))*0.5;

figure
% Plot y, HPBW line, and intersection of the two
plot(theta,y,'LineWidth',2); hold on
plot(th, specline, '--r', 'LineWidth',2)
plot(tmpx, tmpy, '*b', 'LineWidth',4)

answer2c_hpbw_rad = tmpx(1)
answer2c_hpbw_deg = rad2deg(tmpx(1))

% Problem 2a First Null Beam Width
npbw = 0;

% Find the first value where theta = 0.5
tmp = abs(y-npbw);
tmp1 = (tmp<0.0001);

tmpx = theta(tmp1);
tmpy = y(tmp1);

plot(tmpx, tmpy, '*g', 'LineWidth',4)
axis([0, 2*pi, 0, 1])

answer2c_npbw_rad = tmpx(1)
answer2c_npbw_deg = rad2deg(tmpx(1))

% Make Graph Pretty
title({['Problem 2c: cos(\theta)cos(3\theta)'];['HPBW = ',num2str(round(answer2c_hpbw_rad,2)),' [rad], ', num2str(round(answer2c_hpbw_deg,2)), ' [deg]']; ['NPBW = ',num2str(round(answer2c_npbw_rad,2)),' [rad], ', num2str(round(answer2c_npbw_deg,2)), ' [deg]']})
axis([0, 2*pi, -1, 1])
xlabel('Theta [radians]')
ylabel('Normalized Radiation Intensity')

%% Problem 2d Half Power Beam Width
y = (cos(theta).^2) .* (cos(3*theta).^2);

hpbw = 0.5;

tmp = abs(y-hpbw);
tmp1 = (tmp<0.0001);

% Find the first value where theta = 0.5
tmpx = theta(tmp1);
tmpy = y(tmp1);

th=deg2rad(1:36:480);
specline = ones(size(th))*0.5;

figure
% Plot y, HPBW line, and intersection of the two
plot(theta,y,'LineWidth',2); hold on
plot(th, specline, '--r', 'LineWidth',2)
plot(tmpx, tmpy, '*b', 'LineWidth',4)

answer2d_hpbw_rad = tmpx(1)
answer2d_hpbw_deg = rad2deg(tmpx(1))

% Problem 2a First Null Beam Width
npbw = 0;

% Find the first value where theta = 0.5
tmp = abs(y-npbw);
tmp1 = (tmp<0.0001);

tmpx = theta(tmp1);
tmpy = y(tmp1);

plot(tmpx, tmpy, '*g', 'LineWidth',4)
axis([0, 2*pi, 0, 1])

answer2d_npbw_rad = tmpx(1)
answer2d_npbw_deg = rad2deg(tmpx(1))

% Make Graph Pretty
title({['Problem 2d: cos(\theta)cos(3\theta)'];['HPBW = ',num2str(round(answer2d_hpbw_rad,2)),' [rad], ', num2str(round(answer2d_hpbw_deg,2)), ' [deg]']; ['NPBW = ',num2str(round(answer2d_npbw_rad,2)),' [rad], ', num2str(round(answer2d_npbw_deg,2)), ' [deg]']})
axis([0, 2*pi, 0, 1])
xlabel('Theta [radians]')
ylabel('Normalized Radiation Intensity')
