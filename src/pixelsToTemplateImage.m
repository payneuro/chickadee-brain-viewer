function vals = pixelsToTemplateImage(sec, img, ds, ML, DV)
% vals = pixelsToTemplateImage(sec, img, ds, ML, DV) samples one section image,
% read at 1/ds of full resolution, at template-frame points ML, DV (mm, same
% size; the inverse of pixelsToTemplate). Points off the image are NaN.
% A 4x2 midline (two hemispheres mounted apart) is fused first with
% fuseHemispheres. For an A/B image (one hemisphere), points on the far side of
% the midline from its outline are NaN, which drops the partner hemisphere's ghost.
% sec: one row of the template's sections table.
%
% Hannah Payne, 2026

img = double(img);
m = sec.midline;   if iscell(m), m = m{1}; end
if size(m, 1) == 4
    E = sec.edge_px;   if iscell(E), E = E{1}; end
    [Ef, mf, img] = fuseHemispheres(m, E, E, img, ds);
    sec.midline = {mf}; sec.edge_px = {Ef};
end
Xp = pixelsToTemplate([ML(:) DV(:)], sec, 'inverse');   % full-res pixel of each point
vals = interp2(img, (Xp(:, 1) - 0.5)/ds + 0.5, (Xp(:, 2) - 0.5)/ds + 0.5, 'linear', NaN);
if sec.single_hemi
    o = sec.outline;   if iscell(o), o = o{1}; end
    vals((ML(:) >= 0) ~= (median(o(:, 1), 'omitnan') > 0)) = NaN;   % outlines may hold NaN separators
end
vals = reshape(vals, size(ML));
end
