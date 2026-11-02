function [x2, z2] = rotateNoseDown(x, z, theta_deg)
% [x2, z2] = rotateNoseDown(x, z, theta_deg) rotates AP (x) / DV (z) coordinates about
% lambda (0,0) by theta_deg nose down: anterior points move ventral. Use a
% negative angle for nose up. Going from skull-groove angle A1 to A2 (e.g. 65 ->
% 75 deg) is rotateNoseDown(x, z, A2 - A1).
%
% Hannah Payne, 2026

c = cosd(theta_deg); s = sind(theta_deg);
x2 =  x*c + z*s;
z2 = -x*s + z*c;
end
