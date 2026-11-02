function [Xf, mid_f, img_f] = fuseHemispheres(midline, E, X, img, img_scale)
% [Xf, mid_f, img_f] = fuseHemispheres(midline, E, X, [img, img_scale]) joins a
% section whose two hemispheres were mounted apart, each with its own drawn
% midline (midline = 4x2 [x y], points 1-2 one line, 3-4 the other, full-res
% pixels). E: the section outline [x y] (full-res), around both hemispheres.
%
% Each pixel belongs to the hemisphere on its side of the halfway line between
% the two drawn lines. Each hemisphere is rotated so its own line is parallel to
% the halfway line and moved onto it; along the line, the hemispheres are shifted
% so the centers of their outlines' dorsal-ventral extents meet (the ends of the
% drawn lines are too noisy to align on).
%   X       n x 2 full-res pixel positions (cells, outline); [] to skip
%   img     optional image at 1/img_scale of full resolution
%   Xf      fused positions;  mid_f  fused midline, 2x2 [dorsal; ventral]
%   img_f   fused image (same size as img; gap pixels 0)
%
% Hannah Payne, 2026

lines = {midline(1:2, :), midline(3:4, :)};
top = zeros(2, 2); u = zeros(2, 2); len = zeros(2, 1);
for h = 1:2
    L = lines{h};
    [~, i_top] = min(L(:, 2));                    % dorsal end = smaller image y
    top(h, :) = L(i_top, :);
    u(h, :) = (L(3 - i_top, :) - top(h, :)) / norm(L(3 - i_top, :) - top(h, :));   % pointing ventral
    len(h) = norm(L(2, :) - L(1, :));
end
u_mid = mean(u, 1) / norm(mean(u, 1));
top_mid = mean(top, 1);
side = @(P) (P - top_mid) * [u_mid(2) -u_mid(1)]';   % > 0: image-right of the halfway line
s_line = arrayfun(@(h) side(top(h, :) + 0.5*len(h)*u(h, :)), 1:2);
right_hemi = 1 + (s_line(2) > s_line(1));     % index of the hemisphere on the image right
onSide = @(P, h) (side(P) > 0) == (h == right_hemi);

% Rigid transform per hemisphere: rotate about its line's dorsal end onto u_mid,
% move that end to top_mid, then slide along u_mid so the hemispheres' outline
% extents are centered on each other
R = cell(2, 1); center = zeros(2, 1);
for h = 1:2
    a = atan2(u_mid(2), u_mid(1)) - atan2(u(h, 2), u(h, 1));
    R{h} = [cos(a) -sin(a); sin(a) cos(a)];
    Eh = (E(onSide(E, h), :) - top(h, :)) * R{h}' + top_mid;
    t = (Eh - top_mid) * u_mid';
    center(h) = (min(t) + max(t)) / 2;
end
slide = mean(center) - center;
move = @(P, h) (P - top(h, :)) * R{h}' + top_mid + slide(h)*u_mid;
mid_f = top_mid + [0; mean(len)] .* u_mid;

Xf = [];
if ~isempty(X)
    Xf = zeros(size(X));
    for h = 1:2
        in = onSide(X, h);
        Xf(in, :) = move(X(in, :), h);
    end
end

img_f = [];
if nargin >= 4 && ~isempty(img)
    % Inverse mapping per output pixel: an output pixel on side s came from
    % hemisphere s (pull it back through that hemisphere's inverse move)
    [xo, yo] = meshgrid(1:size(img, 2), 1:size(img, 1));
    Po = ([xo(:), yo(:)] - 0.5)*img_scale + 0.5;       % output pixel centers in full-res coords
    img_f = zeros(size(img), 'like', img);
    for h = 1:2
        Psrc = (Po - top_mid - slide(h)*u_mid) * R{h} + top(h, :);   % R' inverse = R^T
        take = onSide(Po, h) & onSide(Psrc, h);
        vals = interp2(double(img), (Psrc(:, 1) - 0.5)/img_scale + 0.5, (Psrc(:, 2) - 0.5)/img_scale + 0.5, 'linear', 0);
        img_f(take) = vals(take);
    end
end
end
