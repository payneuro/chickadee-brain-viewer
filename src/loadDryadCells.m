function cells = loadDryadCells(data_dir, P)
% cells = loadDryadCells(data_dir, [P]) loads every labeled cell from the Dryad
% data of Applegate et al. 2023 (dryad_applegate2023\anatomy.mat: one table per
% bird), as one table with a bird column: bird, ROI, Long, Transverse,
% Hemisphere, AP, ML, DV, wavelength. Without P: Dryad's DV.
% With P (a template's P from makeAnatomyTemplates), DV becomes depth below the
% top of the cell's own section, undoing how the Dryad sections were placed: on
% MC029_orig rotated P.dryad_tilt deg nose down, in its coronal plane at Dryad
% AP + P.dryad_ap_offset, with the section's top at that slice's dorsal-most
% point + P.dryad_dv_offset. dryadToTemplate then places the cells in the
% template. AP and ML are unchanged.
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

D = loadBrain3d(fullfile(data_dir, 'MC029_orig.mat'));
[ap_u, ~, j] = unique(cells.AP);   % sections: one top each
top = arrayfun(@(a) sliceTop(D, a + P.dryad_ap_offset, P.dryad_tilt), ap_u);
cells.DV = cells.DV - (top(j) + P.dryad_dv_offset);
end


function z = sliceTop(D, ap, tilt)
% Dorsal-most point of D's brain in the coronal plane at ap of D rotated tilt deg
% nose down (nearest voxel), NaN if the plane misses it
[Z, Y] = ndgrid(D.zz, D.yy);
[X, Zr] = rotateNoseDown(ap*ones(size(Z)), Z, -tilt);   % back into D's own frame
dI = D.xx(2) - D.xx(1);
ix = round((X - D.xx(1))/dI) + 1; iy = round((Y - D.yy(1))/dI) + 1; iz = round((Zr - D.zz(1))/dI) + 1;
ok = ix >= 1 & ix <= numel(D.xx) & iy >= 1 & iy <= numel(D.yy) & iz >= 1 & iz <= numel(D.zz);
inside = false(size(Z));
inside(ok) = D.Ibrain(sub2ind(size(D.Ibrain), iy(ok), ix(ok), iz(ok))) > 0.5;
z = max(Z(inside));
if isempty(z), z = NaN; end
end
