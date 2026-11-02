function cells = loadDryadCells(data_dir, P)
% cells = loadDryadCells(data_dir, [P]) loads every labeled cell from the Dryad
% data of Applegate et al. 2023 (dryad_applegate2023\anatomy.mat: one table per
% bird), as one table with a bird column: bird, ROI, Long, Transverse,
% Hemisphere, AP, ML, DV, wavelength.
% With P, DV is moved from Marissa's MC029_orig onto MC029_fixed: Dryad DV was set
% from MC029_orig's dorsal surface, so each cell gets z_fixed - z_orig (the
% dorsal surface heights, surfaceHeight), read where the cell lands on MC029 under the
% fitted placement in P (a template's P from makeAnatomyTemplates: section_step,
% ap_mid, ap_shift, default_angle, ref_angle). Cells where either surface is
% missing keep their Dryad DV. AP and ML are unchanged. Without P: Dryad's DV.
%
% Hannah Payne, 2026

anatomy = load(fullfile(data_dir, 'dryad_applegate2023', 'anatomy.mat')).anatomy;
per_bird = cell(numel(anatomy), 1);
for ib = 1:numel(anatomy)
    t = anatomy(ib).data;
    per_bird{ib} = addvars(t, repmat({anatomy(ib).bird}, height(t), 1), 'Before', 1, 'NewVariableNames', 'bird');
end
cells = vertcat(per_bird{:});
if nargin < 2, return; end

D_orig = loadBrain3d(fullfile(data_dir, 'MC029_orig.mat'));
D_fixed = loadBrain3d(fullfile(data_dir, 'MC029_fixed.mat'));
scale = P.section_step / 0.1;
ap_frame = dryadToTemplate(cells.AP, cells.ML, cells.DV, ...   % AP only: scale and shift, no rotation
    struct('ap_mid', P.ap_mid, 'ap_scale', scale, 'ap_shift', P.ap_shift, 'theta', 0, 'c', 0));
ap_ref = rotateNoseDown(ap_frame, cells.DV, -(P.default_angle - P.ref_angle));   % into MC029's own frame
ml_ref = cells.ML / scale;   % MC029's ML is stretched by scale
dv = cells.DV + surfaceHeight(D_fixed, ap_ref, ml_ref) - surfaceHeight(D_orig, ap_ref, ml_ref);
ok = ~isnan(dv);
cells.DV(ok) = dv(ok);
end


function z = surfaceHeight(D, ap, ml)
% Dorsal-most Ibrain voxel at the nearest grid column to each (ap,ml), NaN if
% none
inside = D.Ibrain > 0.5;
[~, top] = max(flip(inside, 3), [], 3);         % first inside voxel from the top
z_map = D.zz(numel(D.zz) + 1 - top);
z_map(~any(inside, 3)) = NaN;
ix = nearestIndex(D.xx, ap);
iy = nearestIndex(D.yy, ml);
z = reshape(z_map(sub2ind(size(z_map), iy, ix)), size(ap));
end


function ind = nearestIndex(grid, vals)
% Same as [~,ind] = min(abs(grid - v)) per value (ties go to the lower index),
% evaluated once per unique value in chunks
[u, ~, j] = unique(vals(:));
ind_u = zeros(size(u));
chunk = 1e4;
for k = 1:chunk:numel(u)
    kk = k:min(k+chunk-1, numel(u));
    [~, ind_u(kk)] = min(abs(grid(:)' - u(kk)), [], 2);
end
ind = ind_u(j);
end
