function [a, b] = pixelsToTemplate(X, sec, direction)
% [ml, dv] = pixelsToTemplate(X, sec) converts full-resolution pixel
% positions X = [x y] (x right, y down) in one coronal section image to ML/DV (mm)
% in the frame of a template built by makeAnatomyTemplates.
% Xpix = pixelsToTemplate([ml dv], sec, 'inverse') goes back to pixels (one
% drawn midline only).
%   ML = sec.px * distance from the midline (sec.midline, 2x2 [x y]), positive
%        toward image right, negated if sec.flip
%   DV = sec.px * height along the midline above the center of the outline's
%        dorsal-ventral extent (sec.edge_px), + sec.shift
% sec: one row of the template's sections table (px = mm per pixel, from the
% section's TIFF header). A 4x2 midline (two separated hemispheres in one image)
% is fused first with fuseHemispheres.
%
% Hannah Payne, 2026

if nargin < 3, direction = 'forward'; end
m = sec.midline;   if iscell(m), m = m{1}; end
E = sec.edge_px;   if iscell(E), E = E{1}; end
if size(m, 1) == 4                    % two hemispheres mounted apart, one line each
    assert(strcmp(direction, 'forward'), 'pixelsToTemplate: inverse not available for fused sections');
    E_fused = fuseHemispheres(m, E, E);
    [X, m] = fuseHemispheres(m, E, X);
    E = E_fused;
end
[~, i_top] = min(m(:, 2));            % dorsal end (smaller image y)
top = m(i_top, :);
u = (m(3 - i_top, :) - top) / norm(m(3 - i_top, :) - top);   % along the midline, pointing ventral
v = [u(2) -u(1)];                     % perpendicular, pointing to image right
s_E = (E - top) * u';
c = (min(s_E) + max(s_E)) / 2;        % center of the outline's extent along the midline

if strcmp(direction, 'forward')
    a = sec.px * ((X - top) * v');                 % ml
    if sec.flip, a = -a; end
    b = -sec.px * ((X - top) * u' - c) + sec.shift;   % dv
else
    ml = X(:, 1); dv = X(:, 2);
    if sec.flip, ml = -ml; end
    a = top + (ml / sec.px) .* v + (c - (dv - sec.shift) / sec.px) .* u;   % [x y]
    b = [];
end
end
