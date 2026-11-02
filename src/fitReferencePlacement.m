%% Fit where an anatomy bird's sections sit on the reference brain (volumes only)
% Free parameters:
%   scale     vibratome thickness / nominal (same blade for both brains): MC022's
%             coronal spacing (AP, about the series middle) and MC029's sagittal
%             spacing (ML, about the midline) are both multiplied by it
%   ap_shift  added to the section AP (as opts.ap_shift in makeAnatomyTemplates)
%   tilt      deg nose down, reference -> bird frame (bird_angle - ref_angle)
%   dv        global DV offset added after each section's dorsal-surface alignment
% For every candidate, each section (A/B halves together) is re-aligned in DV to
% the dorsal surface of the transformed reference slice at its AP, as in
% buildBirdTemplate. Score: Dice of section masks vs reference slices, summed
% over sections. Coarse grid search (dv = 0) on a 0.2 mm grid, then fminsearch
% on a 0.1 mm grid (~2 min in all).
% Input: data\<bird>_75.mat (outline_smooth, shift, AP per section), MC029_fixed.mat
% Output: results\<bird>_ref_fit.png, printed parameters
%
% Hannah Payne, 2026

bird = 'MC022';
code_dir = fileparts(fileparts(mfilename('fullpath')));   % repository root
results_dir = fullfile(code_dir, 'results');
if ~isfolder(results_dir), mkdir(results_dir); end
T = load(fullfile(dataDir(), [bird '_75.mat']));
ref = load(fullfile(dataDir(), 'MC029_fixed.mat'));
if ~exist('grid_coarse', 'var')     % set it in the workspace beforehand to search a different grid
    grid_coarse = struct('scale', 0.85:0.05:1.25, 'ap_shift', 0.5:0.25:3.5, 'tilt', -15:3:45);
end

% Sections at their unaligned DV (A/B halves share a group, as in buildBirdTemplate)
sec = T.sections;
scale_built = T.P.section_step / 0.1;           % section spacing the template was built with / 0.1
ap_mid = T.P.ap_mid;                            % scaling keeps the series middle fixed
ap0 = ap_mid + (sec.AP - T.P.ap_shift - ap_mid) / scale_built;   % Dryad section AP (undo the build's placement)
aligned = ~isnan(sec.shift);
group = sec.k; group(sec.single_hemi) = -sec.trueOrder(sec.single_hemi);
gks = unique(group(aligned));
G = struct('ap0', {}, 'outline', {});
for ig = 1:numel(gks)
    k = find(aligned & group == gks(ig), 1);
    o = sec.outline_smooth{k}; o(:, 2) = o(:, 2) - sec.shift(k);
    G(ig).ap0 = ap0(k); G(ig).outline = o;
end

% Coarse grid search
D = prepareGrid(G, ref, 0.2);
[S, A, R] = ndgrid(grid_coarse.scale, grid_coarse.ap_shift, grid_coarse.tilt);
dice_grid = zeros(size(S));
tic
for i = 1:numel(S)
    dice_grid(i) = sectionDice([S(i) A(i) R(i) 0], D, ap_mid);
end
fprintf('grid search: %d candidates, %.0f s\n', numel(S), toc);
[~, i_best] = max(dice_grid(:));
p0 = [S(i_best) A(i_best) R(i_best) 0];

% Refine on the fine grid
D = prepareGrid(G, ref, 0.1);
p_now = [scale_built T.P.ap_shift T.P.default_angle - T.P.ref_angle 0];
dice_now = sectionDice(p_now, D, ap_mid);
p_fit = fminsearch(@(p) -sectionDice(p, D, ap_mid), p0, optimset('TolX', 1e-3, 'TolFun', 1e-5));
dice_fit = sectionDice(p_fit, D, ap_mid);
fprintf('current: scale %.3f, ap_shift %.2f mm, tilt %.1f deg, dv %.2f mm -> Dice %.4f\n', p_now, dice_now);
fprintf('grid:    scale %.3f, ap_shift %.2f mm, tilt %.1f deg, dv %.2f mm\n', p0);
fprintf('fit:     scale %.3f, ap_shift %.2f mm, tilt %.1f deg, dv %.2f mm -> Dice %.4f\n', p_fit, dice_fit);

% Figure: Dice along each parameter through the fit; sagittal overlays
names = {'scale', 'AP shift (mm)', 'tilt (deg nose down)', 'DV offset (mm)'};
sweeps = {0.8:0.02:1.3, -1:0.1:4, -10:1:45, -0.5:0.05:0.5};
fig = figure('Color', 'w');
fig.OuterPosition = [-1910 50 1900 1020];
for j = 1:4
    ax = subplot(2, 4, j); hold(ax, 'on')
    d = zeros(size(sweeps{j}));
    for i = 1:numel(d)
        p = p_fit; p(j) = sweeps{j}(i);
        d(i) = sectionDice(p, D, ap_mid);
    end
    plot(ax, sweeps{j}, d, 'k-')
    xline(ax, p_fit(j), 'r-'); xline(ax, p_now(j), 'b--');
    xlabel(ax, names{j}); ylabel(ax, 'Dice');
    title(ax, sprintf('fit %.3g (red), current %.3g (blue)', p_fit(j), p_now(j)))
    drawnow
end
mls = [0.5 2 3.5 5];
for j = 1:numel(mls)
    ax = subplot(2, 4, 4 + j); hold(ax, 'on')
    plotSagittal(ax, p_fit, D, ap_mid, mls(j));
    axis(ax, 'equal'); ax.XDir = 'reverse'; grid(ax, 'on')
    xlabel(ax, 'Anterior (mm)'); ylabel(ax, 'Dorsal (mm)');
    title(ax, sprintf('ML %.1f, fit: MC029 grey, %s red', mls(j), bird))
    drawnow
end
exportgraphics(fig, fullfile(results_dir, [bird '_ref_fit.png']), 'Resolution', 80);


function D = prepareGrid(G, ref, dI)
% Section masks (unaligned DV) and dorsal surfaces on a grid of spacing dI;
% reference volume and grid for sampling its slices
D.dI = dI;
D.yy = (floor(min(ref.yy)/dI) - 40 : ceil(max(ref.yy)/dI) + 40) * dI;   % padded: MC022 is wider
D.zz = (floor(min(ref.zz)/dI) - 20 : ceil(max(ref.zz)/dI) + 40) * dI;
D.ap0 = [G.ap0];
D.mask = false(numel(D.yy), numel(D.zz), numel(G));
for ig = 1:numel(G)
    o = G(ig).outline; brk = [0; find(isnan(o(:, 1))); size(o, 1) + 1];
    for i = 1:numel(brk) - 1
        p = o(brk(i)+1:brk(i+1)-1, :);
        if size(p, 1) < 3, continue; end
        D.mask(:, :, ig) = D.mask(:, :, ig) | poly2mask((p(:, 2) - D.zz(1))/dI + 1, (p(:, 1) - D.yy(1))/dI + 1, numel(D.yy), numel(D.zz));
    end
end
D.top = columnTop(D.mask);
D.ref = ref;
end


function [M, top] = refSlices(p, D, ap)
% Reference slices (yy x zz x section) at frame APs ap: ML stretched by the
% scale, rotated tilt deg nose down; nearest-voxel sampling
ref = D.ref; dI = ref.xx(2) - ref.xx(1);
[Y, Z, X] = ndgrid(D.yy, D.zz, ap);
[Xn, Zn] = rotateNoseDown(X, Z, -p(3));
iy = round((Y/p(1) - ref.yy(1))/dI) + 1; ix = round((Xn - ref.xx(1))/dI) + 1; iz = round((Zn - ref.zz(1))/dI) + 1;
ok = iy >= 1 & iy <= numel(ref.yy) & ix >= 1 & ix <= numel(ref.xx) & iz >= 1 & iz <= numel(ref.zz);
M = false(size(Y));
M(ok) = ref.Ibrain(sub2ind(size(ref.Ibrain), iy(ok), ix(ok), iz(ok)));
top = columnTop(M);
end


function top = columnTop(M)
% Index of the highest true voxel in each ML column (NaN if none), yy x section
[~, i_top] = max(flip(M, 2), [], 2);
top = size(M, 2) + 1 - squeeze(i_top);
top(squeeze(~any(M, 2))) = NaN;
end


function [d, B] = sectionDice(p, D, ap_mid)
% Dice over all sections after per-section dorsal alignment plus the global dv.
% B: aligned section masks
ap = ap_mid + (D.ap0 - ap_mid) * p(1) + p(2);
[M, rtop] = refSlices(p, D, ap);
B = false(size(D.mask));
for ig = 1:numel(ap)
    diff_top = rtop(:, ig) - D.top(:, ig);
    if nnz(~isnan(diff_top)) < 5, continue; end        % no reference here: section counts as unmatched
    s = round(median(diff_top, 'omitnan') + p(4)/D.dI);
    B(:, :, ig) = circshift(D.mask(:, :, ig), s, 2);
end
d = 2*nnz(B & M) / (nnz(D.mask) + nnz(M));
end


function plotSagittal(ax, p, D, ap_mid, ml)
% Transformed reference (grey) and aligned sections (red outline) in the
% sagittal plane at ml
[~, B] = sectionDice(p, D, ap_mid);
[~, iy] = min(abs(D.yy - ml));
ap = ap_mid + (D.ap0 - ap_mid) * p(1) + p(2);
[ap, o] = sort(ap);
D1 = D; D1.yy = D.yy(iy);                       % reference along the full AP axis, this ML only
R = squeeze(refSlices(p, D1, D.ref.xx));        % zz x AP
imagesc(ax, D.ref.xx, D.zz, double(R), 'AlphaData', 0.4*R);
colormap(ax, gray); clim(ax, [0 1.5]); ax.YDir = 'normal';
contour(ax, ap, D.zz, double(squeeze(B(iy, :, o))), [0.5 0.5], 'r');
end
