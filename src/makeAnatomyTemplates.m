%% Build a template (and Dryad -> template transform) for each anatomy bird
% Each bird's folder "<imaging_root>\<bird> sortedSlices" has the layout read by
% loadAnatomyBird. Add birds to the list as their section data become available.
%
% Output per bird: data\<bird>_75.mat with xx, yy, zz, Ibrain (brain outline only),
% P, sections (per-section AP, DV shift, midline, pixel size, image file, ...) and
% tf (Dryad transform, used by brainViewer through dryadToTemplate). sections.cells:
% each section's labelled cells (all wavelengths) in the template frame.
% injections: the injection site of each wavelength's tracer (Dryad Transverse,
% Long).
% Section images are read at 1/8 through readSectionImage (cached in data\<bird>_ds8).
% Tissue outlines and joins of detached pieces: data\<bird>_outlines.mat (run
% makeSectionOutlines first).
% QC figures (plot_qc = true): results\<bird>_75_qc.png, _midline_qc_1/2.png,
% _auto_qc.png, _smoothing_qc.png
%
% Steps are the local functions below: loadAnatomyBird (section data + pixel size
% and AP of each section), buildBirdTemplate (dorsal-surface alignment and
% outline volume), fitDryadTransform (Dryad DV -> template DV).
% Needs: pixelsToTemplate, dryadToTemplate, rotateNoseDown (shared with brainViewer),
% joinPieces.
%
% Hannah Payne, 2026

birds = {'MC022'};
imaging_root = 'Z:\Hannah\imaging';
draw_corrections = false;   % true: first open sections with a missing/broken midline or outline to draw by hand
                            % (saved to <bird folder>\midlineCorrections.mat, used on every later run)
plot_qc = false;            % true: also draw the QC figures and save them to results\
% ref_angle, ap_shift and section_step: fitReferencePlacement (MC022 section
% volume vs MC029, Dice 0.946; DV offset fixed at 0). section_step is the same
% vibratome for all birds: it scales each bird's AP spacing (about the series
% middle) and MC029's ML (sagittal sections, about the midline).
opts = struct('ref_angle', 55.1, ...  % MC029_fixed's effective skull-groove angle: 75 - 19.9 deg fitted tilt
    'bird_angle', 75, ...             % anatomy birds: 75 deg (capillary on beak groove)
    'ap_shift', 0.5, ...              % added to Dryad section AP (fit 0.50 mm)
    'section_step', 0.1062, ...       % vibratome section thickness (mm): nominal 0.1 x fitted 1.062
    'bin', 0.05, ...                  % ML bin for dorsal-surface matching (mm)
    'close_um', 75, ...               % closing radius for section masks: fills cracks/infoldings < 2x this
    'open_um', 240, ...               % opening radius after closing: removes protrusions < 2x this
    'median_sections', 3);            % template: median over this many neighboring sections along AP
qc_aps = [-1.0 0.0 1.0 2.0 3.0 3.8 4.5 5.5];   % sections shown in the QC figure
% Region of each labelledCells label number (columns: label 1, 2, ...), per
% channel. Not stored with the data: for MC022, inferred by matching per-section
% counts to Dryad. "" = no matching region (cells left out). Replace with the key
% from Marissa when available.
label_key.MC022 = table(["488"; "555"; "647"], ...
    ["DL/CDL" "HD" "NIL_M" "NIL_L" ""; "DL/CDL" "HD" "NIL_M" "NFL" ""; "DL/CDL" "HD" "NIL_M" "NFL" ""], ...
    'VariableNames', {'wavelength', 'roi'});
% For parameter sweeps: a struct opts_override in the workspace replaces opts fields
if exist('opts_override', 'var')
    for f = fieldnames(opts_override)', opts.(f{1}) = opts_override.(f{1}); end
end

code_dir = fileparts(fileparts(mfilename('fullpath')));   % repository root
data_dir = fullfile(code_dir, 'data');
results_dir = fullfile(code_dir, 'results');
if ~isfolder(results_dir), mkdir(results_dir); end
ref = load(fullfile(data_dir, 'MC029_fixed.mat'));
ref.yy = ref.yy * opts.section_step / 0.1;   % MC029's sagittal sections were cut with the same vibratome
cells = loadDryadCells(data_dir);   % Dryad's own DV

for ib = 1:numel(birds)
    dry = cells(string(cells.bird) == birds{ib}, :);   % raw Dryad AP/ML/DV columns
    bird_dir = fullfile(imaging_root, [birds{ib} ' sortedSlices']);
    bird = loadAnatomyBird(bird_dir, dry, data_dir);
    % to draw: still unusable, or fixed automatically (confirm or adjust)
    to_draw = find(arrayfun(@(k) sectionProblem(bird, k), 1:numel(bird.meta)) ~= "" | bird.auto ~= "");
    if draw_corrections && ~isempty(to_draw)
        drawCorrections(bird, to_draw);
        bird = loadAnatomyBird(bird_dir, dry, data_dir);    % reload with the corrections
    end
    T = buildBirdTemplate(bird, ref, opts);
    dry_fixed = loadDryadCells(data_dir, T.P);   % Dryad DV moved onto MC029_fixed, at the placement just used
    dry.DV = dry_fixed.DV(string(dry_fixed.bird) == birds{ib});
    T.tf = fitDryadTransform(T, bird, dry);
    fprintf('%s: %d/%d sections aligned; dorsal fit MAD median %.3f mm\n', bird.name, ...
        nnz(~isnan(T.sections.shift)), height(T.sections), median(T.sections.fit_mad, 'omitnan'));
    fprintf('  Dryad -> template: theta %.2f deg nose down, c %.3f mm, RMS %.3f mm over %d section-channels\n', ...
        T.tf.theta, T.tf.c, T.tf.rms, height(T.tf.per_section));
    T.P.overlap_ref = overlapWithRef(T, ref);
    fprintf('  overlap with %s over the template''s AP range: Dice %.3f\n', 'MC029', T.P.overlap_ref);
    S = T.tf.per_section;
    bad_ml = abs(S.ml_pixels - S.ml_dryad) > 0.5;
    fprintf('  ML check: %d section-channels differ from Dryad by > 0.5 mm (k = %s)\n', nnz(bad_ml), mat2str(unique(S.k(bad_ml))'));
    T.sections.cells = sectionCells(T, bird, label_key.(bird.name));
    T.injections = unique(table(string(dry.wavelength), string(dry.Transverse), dry.Long, ...   % tracer per wavelength (Dryad)
        'VariableNames', {'wavelength', 'Transverse', 'Long'}));
    save(fullfile(data_dir, [bird.name '_75.mat']), '-struct', 'T', '-v7');
    if ~plot_qc, continue; end
    plotTemplateQC(T, bird, ref, opts, dry, qc_aps, fullfile(results_dir, [bird.name '_75_qc.png']));
    plotMidlineQC(T, bird, fullfile(results_dir, [bird.name '_midline_qc']));
    plotAutoQC(bird, fullfile(results_dir, [bird.name '_auto_qc.png']));
    plotSmoothingQC(T, opts, fullfile(results_dir, [bird.name '_smoothing_qc.png']));
end


function plotMidlineQC(T, bird, png_stem)
% Sections whose midline was not a single drawn line: two hemispheres in one
% image (raw image + both lines | fused image + fused midline + cells), and A/B
% image pairs fused in mm (plotPairQC).
ds = 8;   % image read at 1/8 scale
[~, w] = max(double([bird.labelled.wavelength]));
L = bird.labelled(w).cells;
plotPairQC(T, sprintf('%s_2.png', png_stem));
ks = find(contains(T.sections.midline_source, "two hemispheres fused"))';
if isempty(ks), return; end
fig = figure('Position', [30 30 1040*min(2, numel(ks)) 420*ceil(numel(ks)/2)]);
tl = tiledlayout(fig, ceil(numel(ks)/2), 2*min(2, numel(ks)), 'TileSpacing', 'compact');
for k = ks
    sec = T.sections(k, :);
    [img, full_size] = readSectionImage(sec.image_file);
    img = double(img);
    clim_img = prctile(img(:), [1 99.5]);
    E = bird.edge_scale * bird.meta(k).edge(:, [2 1]);
    P = L.pos{k}(L.label{k} > 0, :);
    m = sec.midline{1};
    ext = [1 full_size(2) 1 full_size(1)];
    nexttile(tl); imagesc(ext(1:2), ext(3:4), img, clim_img); colormap(gray); axis image; hold on
    plot(E(:, 1), E(:, 2), 'c-')
    plot(m(1:2, 1), m(1:2, 2), 'y-o', m(3:4, 1), m(3:4, 2), 'm-o', 'LineWidth', 1.5)
    plot(P(:, 1), P(:, 2), 'r.', 'MarkerSize', 3)
    title(sprintf('k=%d raw: two drawn lines', k))
    [Ef, mf, img_f] = fuseHemispheres(m, E, E, img, ds);
    Pf = fuseHemispheres(m, E, P);
    nexttile(tl); imagesc(ext(1:2), ext(3:4), img_f, clim_img); axis image; hold on
    plot(Ef(:, 1), Ef(:, 2), 'c-'); plot(mf(:, 1), mf(:, 2), 'y-o', 'LineWidth', 1.5)
    plot(Pf(:, 1), Pf(:, 2), 'r.', 'MarkerSize', 3)
    title(sprintf('k=%d fused', k))
end
drawnow
exportgraphics(fig, sprintf('%s_1.png', png_stem), 'Resolution', 70);
end


function cells = sectionCells(T, bird, key)
% Each section's labelled cells (all wavelengths) in the template frame: one
% table per section (ml, dv, wavelength, roi), empty if not aligned. key: region
% of each label number per wavelength (label_key above); labels without a
% region are left out.
n = height(T.sections);
cells = repmat({table(zeros(0, 1), zeros(0, 1), strings(0, 1), strings(0, 1), ...
    'VariableNames', {'ml', 'dv', 'wavelength', 'roi'})}, n, 1);
for k = find(~isnan(T.sections.shift))'
    for w = 1:numel(bird.labelled)
        L = bird.labelled(w).cells;
        roi_of = key.roi(key.wavelength == bird.labelled(w).wavelength, :);
        lab = L.label{k};
        keep = lab > 0 & lab <= numel(roi_of);
        keep(keep) = roi_of(lab(keep)) ~= "";
        if ~any(keep), continue; end
        [ml, dv] = pixelsToTemplate(L.pos{k}(keep, :), T.sections(k, :));
        cells{k} = [cells{k}; table(ml, dv, repmat(bird.labelled(w).wavelength, numel(ml), 1), ...
            roi_of(lab(keep))', 'VariableNames', cells{k}.Properties.VariableNames)];
    end
end
end


function plotSmoothingQC(T, opts, png_file)
% Sections whose masks changed most from their outline to the template: outline
% (grey), closed + opened (blue), after the median over neighboring sections
% (red), in the template frame
n_show = 10;
[changed, order] = sort(T.sections.median_change_mm2, 'descend', 'MissingPlacement', 'last');
order = order(~isnan(changed));
key = T.sections.k; key(T.sections.single_hemi) = -T.sections.trueOrder(T.sections.single_hemi);   % A/B halves: one mask
[~, first] = unique(key(order), 'stable'); order = order(first);
ks = order(1:min(n_show, numel(order)))';
fig = figure('Position', [30 30 1700 760]);
tl = tiledlayout(fig, 2, ceil(numel(ks)/2), 'TileSpacing', 'compact');
for k = ks
    nexttile(tl); hold on
    o = cellfun(@(x) [x; NaN NaN], T.sections.outline(key == key(k)), 'UniformOutput', false);
    o = vertcat(o{:});
    plot(o(:, 1), o(:, 2), '-', 'Color', .5*[1 1 1], 'LineWidth', 0.75)
    for c = {T.sections.outline_smooth{k}, T.sections.outline_median{k}; [0 0.4 1], [1 0 0]}
        plot(c{1}(:, 1), c{1}(:, 2), '-', 'Color', c{2}, 'LineWidth', 0.75)
    end
    [~, f] = fileparts(T.sections.image_file(k));
    axis equal; grid on
    title(sprintf('%s (k=%d), median changed %.2f mm^2', f, k, T.sections.median_change_mm2(k)), 'Interpreter', 'none', 'FontSize', 9)
end
legend(nexttile(tl, 1), {'outline', sprintf('closed %g um, opened %g um', opts.close_um, opts.open_um), ...
    sprintf('median of %d sections', opts.median_sections)}, 'Location', 'southwest')
xlabel(tl, 'ML (mm)'); ylabel(tl, 'DV (mm)')
drawnow
exportgraphics(fig, png_file, 'Resolution', 75);
end


function plotAutoQC(bird, png_file)
% Sections fixed automatically: image, outline (detected or stored), the flat
% stretch used for the midline, and the midline
ds = bird.edge_scale;
ks = find(bird.auto ~= "");
if isempty(ks), return; end
fig = figure('Position', [30 30 420*numel(ks) 470]);
tl = tiledlayout(fig, 1, numel(ks), 'TileSpacing', 'compact');
for k = ks
    [img, full_size] = readSectionImage(bird.image_files{k});
    img = double(img);
    nexttile(tl); imagesc([1 full_size(2)], [1 full_size(1)], img, prctile(img(:), [1 99.5])); colormap(gray); axis image; hold on
    E = ds * bird.meta(k).edge(:, [2 1]);
    plot(E([1:end 1], 1), E([1:end 1], 2), 'c-')
    m = bird.meta(k).midLineCorr;
    if contains(bird.auto(k), "flat")
        side = "left"; if bird.ab(k) == "B", side = "right"; end
        [~, fi] = flatSideMidline(E, side, bird.px(k));
        plot(fi.points(:, 1), fi.points(:, 2), 'm.', 'MarkerSize', 5)
    end
    if size(m, 1) == 2, plot(m(:, 1), m(:, 2), 'y-o', 'LineWidth', 1.5); end
    [~, f] = fileparts(bird.image_files{k});
    title(sprintf('%s (k=%d)\n%s', f, k, bird.auto(k)), 'Interpreter', 'none', 'FontSize', 7)
end
legend(nexttile(tl, 1), {'outline', 'flat stretch', 'midline'}, 'TextColor', 'w', 'Color', 'k')
drawnow
exportgraphics(fig, png_file, 'Resolution', 70);
end


function plotPairQC(T, png_file)
% A/B pairs (two images of one section, one hemisphere each) as fused into the
% template frame by renderSection (after the pair's flip), with the pair's
% smoothed outline and the longest-wavelength cells
S = T.sections;
orders = unique(S.trueOrder(S.single_hemi))';
[ML, DV] = meshgrid(-7:0.02:7, -10:0.02:2);
w = string(max(double(vertcat(S.cells{:}).wavelength)));
fig = figure('Position', [30 30 360*numel(orders) 420]);
tl = tiledlayout(fig, 1, numel(orders), 'TileSpacing', 'compact');
for o = orders
    nexttile(tl); hold on
    [canvas, ks] = renderSection(S, o, ML, DV);
    imagesc(ML(1, :), DV(:, 1), canvas, [0 1]); colormap(gray); set(gca, 'YDir', 'normal')
    if ~isempty(ks), oc = S.outline_smooth{ks(1)}; plot(oc(:, 1), oc(:, 2), 'c-'); end   % pair's joint smoothed outline
    for k = ks
        c = S.cells{k}(S.cells{k}.wavelength == w, :);
        plot(c.ml, c.dv, 'r.', 'MarkerSize', 3)
    end
    xline(0, 'y:'); axis image; xlim(ML(1, [1 end])); ylim(DV([1 end], 1))
    ka = find(S.trueOrder == o)';
    [~, names] = fileparts(S.image_file(ka));
    title(sprintf('trueOrder %d: %s\n%s', o, strjoin(names, ' + '), ...
        strjoin(S.midline_source(ka), ' | ')), 'Interpreter', 'none', 'FontSize', 7)
    xlabel('ML (mm)'); ylabel('DV (mm)')
end
drawnow
exportgraphics(fig, png_file, 'Resolution', 70);
end


function plotTemplateQC(T, bird, ref, opts, dry, aps_show, png_file)
% Example sections: rotated reference slice, bird outline, cells from pixels and
% Dryad cells after the fitted transform (longest-wavelength channel)
[~, w] = max(double([bird.labelled.wavelength]));
L = bird.labelled(w).cells;
dw = dry(string(dry.wavelength) == bird.labelled(w).wavelength, :);
fig = figure('Position', [50 50 1600 760]);
for ii = 1:numel(aps_show)
    k = find(abs(bird.section_ap(:) - aps_show(ii)) < 1e-6 & cellfun(@(l) any(l > 0), L.label(:)), 1);
    if isempty(k) || isnan(T.sections.shift(k)), continue; end
    sec = T.sections(k, :);
    subplot(2, 4, ii); hold on
    [Yq, Zq] = meshgrid(ref.yy, ref.zz);
    [Xn, Zn] = rotateNoseDown(sec.AP*ones(size(Zq)), Zq, -(opts.bird_angle - opts.ref_angle));
    Sref = interp3(ref.xx, ref.yy, ref.zz, single(ref.Ibrain), Xn, Yq, Zn, 'nearest', 0);
    contour(Yq, Zq, Sref, [0.5 0.5], 'k-', 'LineWidth', 0.75)
    o = sec.outline{1}; plot(o(:, 1), o(:, 2), '-', 'Color', .6*[1 1 1], 'LineWidth', 0.5)
    [ml, dv] = pixelsToTemplate(L.pos{k}(L.label{k} > 0, :), sec);
    q = abs(dw.AP - bird.section_ap(k)) < 1e-6;
    [~, mlq, dvq] = dryadToTemplate(dw.AP(q), dw.ML(q), dw.DV(q), T.tf);
    plot(ml, dv, 'r.', 'MarkerSize', 5); plot(mlq, dvq, 'b.', 'MarkerSize', 5)
    axis equal; grid on; xline(0, 'k:')
    title(sprintf('AP %.1f (k=%d), shift %.2f, MAD %.2f', aps_show(ii), k, sec.shift, sec.fit_mad))
    xlabel('ML (mm)'); ylabel('DV (mm)')
end
legend(subplot(2, 4, 1), {sprintf('reference rotated %g\\circ', opts.bird_angle - opts.ref_angle), ...
    [bird.name ' outline'], 'cells from pixels', 'Dryad cells, transformed'}, 'Location', 'southwest')
drawnow
exportgraphics(fig, png_file, 'Resolution', 80);
end


function bird = loadAnatomyBird(bird_dir, dry, data_dir)
% bird = loadAnatomyBird(bird_dir, dry, data_dir) collects one anatomy bird's
% section data for buildBirdTemplate / fitDryadTransform, assuming the folder
% layout of "<bird> sortedSlices":
%   metaDataUpdated3.mat              metaData, one entry per section image
%   images\slice_*_640.tif            section images (pixel size from the header)
%   <bird>_<wavelength>_labelledCells.mat   detected cells (pos, label) per section
% dry: this bird's Dryad rows (AP, ML, DV, wavelength), used to find section APs.
%
%   midlineCorrections.mat (optional) hand-drawn midlines/outlines, made with
%                                     draw_corrections = true; replace metaData
% and from data_dir, <bird>_outlines.mat (makeSectionOutlines): tissue outlines,
% which replace metaData.edge (hand-drawn ones take precedence), and the joins
% of detached pieces. In a joined image, metaData midlines and the cells are
% moved with their pieces (joinPoints), so everything is in joined-image pixels;
% hand-drawn midlines are drawn on the joined image already.
% bird: name, meta (edge: [row col] at 1/edge_scale, several polygons separated
% by NaN rows), corrected, edge_corrected, midline_corrected, mirror (per
% section, from midlineCorrections), ab ("A", "B" or ""), single_hemi (A/B
% images: one hemisphere each), auto (automatic midline fixes, "" if none),
% image_files, px (mm/pixel per section), flip (metaData.reverse: ipsilateral
% hemisphere on the image left, so ML is negated to make ipsilateral positive as
% in Dryad; checked against the cells in k = 36 and the A/B images), edge_scale,
% join (per section: lab, steps, peaks; steps empty if not joined), labelled
% (per wavelength), section_ap (AP relative to lambda per section), dir.
%
% Section AP is not stored in metaData. Sections are 0.1 mm apart in trueOrder,
% so AP = ap0 - 0.1*trueOrder; ap0 is the value at which per-section labelled-cell
% counts match the Dryad per-AP counts best.
%
% Hannah Payne, 2026

section_step = 0.1;   % mm between consecutive trueOrder values
[~, folder] = fileparts(bird_dir);
bird.name = strtrim(erase(folder, 'sortedSlices'));
bird.dir = bird_dir;
bird.meta = load(fullfile(bird_dir, 'metaDataUpdated3.mat')).metaData;
n = numel(bird.meta);
% Hand corrections per section (k): midLineCorr, edge (either may be empty =
% keep), mirror (A/B images: mirror the fused section left-right)
bird.edge_corrected = false(1, n);
bird.midline_corrected = false(1, n);
bird.mirror = false(1, n);
corr_file = fullfile(bird_dir, 'midlineCorrections.mat');
corrections = struct('k', {}, 'midLineCorr', {}, 'edge', {}, 'mirror', {});
if isfile(corr_file), corrections = load(corr_file).corrections; end
for c = corrections
    if ~isempty(c.midLineCorr), bird.meta(c.k).midLineCorr = c.midLineCorr; bird.midline_corrected(c.k) = true; end
    if ~isempty(c.edge), bird.meta(c.k).edge = c.edge; bird.edge_corrected(c.k) = true; end
    if isfield(c, 'mirror') && ~isempty(c.mirror), bird.mirror(c.k) = c.mirror; end
end
bird.corrected = bird.edge_corrected | bird.midline_corrected;
bird.edge_scale = 8;   % edges are [row col] in the 1/8-scale image (blob-file pyramid level 4; readSectionImage)

files = dir(fullfile(bird_dir, 'images', 'slice_*_640.tif'));
assert(numel(files) == n, '%s: %d images for %d metaData entries', bird.name, numel(files), n);
bird.image_files = fullfile({files.folder}, {files.name})';   % sorted by name = metaData order
% A/B images (e.g. slice_02A_640.tif): one hemisphere each. The letter only says
% where to start looking for the medial side (A: image left, B: image right);
% once a midline exists, the hemisphere's side comes from its outline.
tok = regexp({files.name}, '\d([AB])_', 'tokens', 'once');
bird.ab = strings(1, n);
bird.ab(~cellfun(@isempty, tok)) = string(cellfun(@(t) t{1}, tok(~cellfun(@isempty, tok)), 'UniformOutput', false));
bird.single_hemi = bird.ab ~= "";

% Tissue outlines and joins (makeSectionOutlines)
O = load(fullfile(data_dir, [bird.name '_outlines.mat'])).S;
assert(numel(O) == n, '%s: %d outlines for %d metaData entries', bird.name, numel(O), n);
bird.join = struct('lab', {O.lab}, 'steps', {O.steps}, 'peaks', {O.peaks});
for k = 1:n
    if ~bird.edge_corrected(k) && ~isempty(O(k).outline)
        polys = O(k).outline(:);
        if bird.single_hemi(k)   % one hemisphere per image: drop the partner hemisphere's ghost (the smaller blob)
            [~, i] = max(cellfun(@(b) polyarea(b(:, 2), b(:, 1)), polys));
            polys = polys(i);
        end
        bird.meta(k).edge = cell2mat(cellfun(@(b) [b; NaN NaN], polys, 'UniformOutput', false));
        bird.meta(k).edge(end, :) = [];
    end
    if O(k).joined && ~bird.midline_corrected(k) && ~isempty(bird.meta(k).midLineCorr)
        bird.meta(k).midLineCorr = joinPoints(bird.meta(k).midLineCorr, bird.join(k), bird.edge_scale);
    end
end

bird.px = cellfun(@tiffPixelSize, bird.image_files);
bird.flip = [bird.meta.reverse];

% Automatic fix where the midline is missing or implausible (hand corrections,
% loaded above, take precedence): A/B images (one hemisphere each, medial side on
% the image left for A, right for B) get a midline along their flat medial side.
% bird.auto says what was done.
bird.auto = strings(1, n);
for k = find(~bird.midline_corrected)
    if sectionProblem(bird, k) ~= "" && bird.single_hemi(k) && ~isempty(bird.meta(k).edge)
        side = "left"; if bird.ab(k) == "B", side = "right"; end
        [bird.meta(k).midLineCorr, info] = flatSideMidline(bird.edge_scale * bird.meta(k).edge(:, [2 1]), side, bird.px(k));
        bird.auto(k) = sprintf("midline on flat %s side (%.2f mm within %.0f um, RMS %.0f um)", ...
            side, info.length_mm, info.tol_um, info.rms_um);
        fprintf('  section %d: %s\n', k, bird.auto(k));
    end
end

lf = dir(fullfile(bird_dir, [bird.name '_*_labelledCells.mat']));
for w = 1:numel(lf)
    bird.labelled(w).wavelength = string(regexp(lf(w).name, '_(\d+)_labelledCells', 'tokens', 'once'));
    bird.labelled(w).cells = load(fullfile(lf(w).folder, lf(w).name)).cells;
    for k = find(~cellfun(@isempty, {bird.join.steps}))     % cells move with their pieces
        bird.labelled(w).cells.pos{k} = joinPoints(bird.labelled(w).cells.pos{k}, bird.join(k), bird.edge_scale);
    end
end

% Section AP: try every ap0 that puts the most-labelled section on a Dryad AP
order = [bird.meta.trueOrder];
n_lab = zeros(numel(bird.labelled), n);
for w = 1:numel(bird.labelled)
    n_lab(w, :) = cellfun(@(l) nnz(l > 0), bird.labelled(w).cells.label);
end
[~, k_max] = max(sum(n_lab, 1));
best = -1;
for ap = unique(dry.AP)'
    ap0 = ap + section_step*order(k_max);
    sec_ap = round((ap0 - section_step*order)/section_step)*section_step;
    n_match = 0;
    for w = 1:numel(bird.labelled)
        dw = dry(string(dry.wavelength) == bird.labelled(w).wavelength, :);
        n_dry = arrayfun(@(a) nnz(abs(dw.AP - a) < 1e-6), sec_ap);
        n_match = n_match + nnz(n_lab(w, :) > 0 & n_lab(w, :) == n_dry);
    end
    if n_match > best, best = n_match; bird.section_ap = sec_ap; ap0_best = ap0; end
end
fprintf('%s: section AP = %.1f - %.1f*trueOrder (labelled-cell counts match Dryad in %d of %d section-channels)\n', ...
    bird.name, ap0_best, section_step, best, nnz(n_lab > 0));
end


function px = tiffPixelSize(file)
% px = tiffPixelSize(file) returns the pixel size (mm) of a TIFF from its header:
% XResolution (pixels per unit), with the unit read from the ImageJ description
% ("unit=micron") or from ResolutionUnit.
%
% Hannah Payne, 2026

info = imfinfo(file);
info = info(1);
unit = '';
if isfield(info, 'ImageDescription')
    tok = regexp(info.ImageDescription, 'unit=(\S+)', 'tokens', 'once');
    if ~isempty(tok), unit = tok{1}; end
end
if isempty(unit), unit = info.ResolutionUnit; end
mm_per_unit = struct('micron', 1e-3, 'um', 1e-3, 'mm', 1, 'Centimeter', 10, 'Inch', 25.4);
unit_key = strrep(unit, char(181), 'u');   % ImageJ may write a micro sign
if ~isfield(mm_per_unit, unit_key)
    error('tiffPixelSize: unknown resolution unit "%s" in %s', unit, file);
end
px = mm_per_unit.(unit_key) / info.XResolution;
end


function T = buildBirdTemplate(bird, ref, opts)
% T = buildBirdTemplate(bird, ref, opts) builds a brain-outline template from one
% anatomy bird's coronal section outlines, aligned dorso-ventrally to a reference
% template (as in Applegate et al. 2023: each section's dorsal surface was aligned
% to the reference brain).
%
% Steps
%   1. Rotate the reference nose down by bird_angle - ref_angle (75 - 55.1 =
%      19.9 deg, fitted by fitReferencePlacement) about lambda (0,0,0), so it is
%      at the bird's skull-groove angle.
%   2. For each section: outline in mm, ML from the hand-drawn midline, height
%      along the midline. Its dorsal surface (highest outline point in each ML bin)
%      is compared with the dorsal surface of the rotated reference slice at the
%      section's AP; the DV shift is the median difference over all shared bins.
%   3. Fill each shifted outline into the planes within half a section of its AP.
%
% bird (struct, from loadAnatomyBird)
%   name        e.g. 'MC022'
%   meta        metaData struct array, one per section image, with fields
%               midLineCorr (2x2 [x y], full-res pixels), edge ([row col] at
%               1/edge_scale), trueOrder
%   section_ap  AP of each section relative to lambda (mm), same order as meta
%   px          mm per full-res pixel, per section (from the TIFF headers)
%   edge_scale  edge is stored at 1/edge_scale of full resolution
%   flip        logical per section: mirror ML (image flipped left-right)
% ref (struct)   reference template (load of e.g. MC029_fixed.mat): xx, yy, zz, Ibrain
%               (yy may be stretched, so its spacing can differ from xx and zz)
% opts (struct)
%   ref_angle   skull-groove angle of the reference (deg)
%   bird_angle  skull-groove angle the bird was sectioned at (deg)
%   ap_shift    added to section_ap to place sections on the reference (mm)
%   section_step  physical section spacing (mm); section_ap (0.1 mm apart) is
%               scaled to it about the middle of the series
%   bin         ML bin for comparing dorsal surfaces (mm)
%   close_um, open_um  closing and opening radii for each section's mask (um)
%   median_sections  sections per AP median of the masks
%
% T: xx, yy, zz (reference grid, padded in ML/DV), Ibrain (logical), P, sections
% (table: k, trueOrder, AP (frame), shift, flip, px, midline, fit_mad, n_bins,
% outline, image_file, join_lab, join_steps, peaks, outline_smooth,
% outline_median, median_change_mm2; sectionCells adds cells). P also holds
% section_step and ap_mid. Use pixelsToTemplate(X, T.sections(k,:))
% to map pixels (of the joined image, if joined) into T.
%
% Hannah Payne, 2026

rot = opts.bird_angle - opts.ref_angle;      % deg nose down, reference -> bird frame
dI = ref.xx(2) - ref.xx(1);
n = numel(bird.meta);
% Frame AP: Dryad AP (0.1 mm per section) + ap_shift, with the spacing scaled to
% section_step about the middle of the series
ap_mid = (min(bird.section_ap) + max(bird.section_ap)) / 2;
frame_ap = ap_mid + (bird.section_ap(:) - ap_mid) * opts.section_step / 0.1 + opts.ap_shift;
sections = table((1:n)', [bird.meta.trueOrder]', frame_ap, nan(n, 1), ...
    logical(bird.flip(:)), logical(bird.single_hemi(:)), bird.px(:), cell(n, 1), cell(n, 1), strings(n, 1), nan(n, 1), zeros(n, 1), ...
    cell(n, 1), string(bird.image_files(:)), ...
    'VariableNames', {'k', 'trueOrder', 'AP', 'shift', 'flip', 'single_hemi', 'px', 'midline', 'edge_px', 'midline_source', ...
    'fit_mad', 'n_bins', 'outline', 'image_file'});
sections.join_lab = {bird.join.lab}';          % joins of detached pieces (renderSection: joinPieces)
sections.join_steps = {bird.join.steps}';
sections.peaks = vertcat(bird.join.peaks);     % log intensity of background and tissue, for display contrast

% A/B pairs: reverse is applied to the fused pair, not to each image (the two
% images' flags differ). The pair is mirrored if the image holding the cells has
% reverse = 1 (ipsilateral to ML > 0, as for whole sections); a hand-set mirror
% inverts that.
n_cells = zeros(1, n);
for w = 1:numel(bird.labelled)
    n_cells = n_cells + reshape(cellfun(@(l) nnz(l > 0), bird.labelled(w).cells.label), 1, []);
end
pair_flip = false(1, n);
order = [bird.meta.trueOrder];
for o = unique(order(bird.single_hemi))
    ks = find(bird.single_hemi & order == o);
    [n_max, i_max] = max(n_cells(ks));
    pair_flip(ks) = xor(n_max > 0 && bird.flip(ks(i_max)), any(bird.mirror(ks)));
end

for k = 1:n
    % Midline: drawn (2x2), or two separated hemispheres (4x2, fused by
    % pixelsToTemplate). Missing or implausible ones are skipped until drawn.
    problem = sectionProblem(bird, k);
    if problem ~= ""
        sections.midline_source(k) = "needs drawing: " + problem;
        fprintf('  section %d: %s; skipped (set draw_corrections = true to draw it)\n', k, problem);
        continue
    end
    m = bird.meta(k).midLineCorr;
    src = "drawn";
    if size(m, 1) == 4, src = src + ", two hemispheres fused"; end
    if bird.single_hemi(k), src = src + ", one hemisphere (A/B image, fused by AP)"; end
    if bird.corrected(k), src = src + ", hand-corrected"; end
    if bird.auto(k) ~= "", src = "auto: " + bird.auto(k); end
    sections.midline_source(k) = src;
    sections.midline{k} = m;
    sections.edge_px{k} = bird.edge_scale * bird.meta(k).edge(:, [2 1]);   % outline, full-res [x y]
    sec0 = sections(k, :); sec0.shift = 0;
    % A/B images: each hemisphere lies on its own side of its midline (ML > 0 or
    % < 0), so the pair forms a whole section; then the pair's flip (above)
    if sec0.single_hemi
        sections.flip(k) = pair_flip(k);
        sec0.flip = pair_flip(k);
    end
    E = densify(sections.edge_px{k}, 0.01 / sec0.px);           % ~10 um spacing
    [eML, eH] = pixelsToTemplate(E, sec0);                   % height above the outline's DV-extent center
    sections.outline{k} = [eML, eH];                             % shift added below
end

% DV shift: the dorsal surface (highest outline point per ML bin) against the
% rotated reference slice's dorsal surface at the same AP. An A/B pair is one
% section: its two hemispheres keep the centers of their outlines' DV extents
% together (as in fuseHemispheres) and get one shift from their combined outline.
usable = ~cellfun(@isempty, sections.outline);
group = sections.k;
pair = usable & sections.single_hemi;
group(pair) = -sections.trueOrder(pair);                     % A/B halves share a group
for gk = unique(group(usable))'
    ks = find(usable & group == gk);
    o = vertcat(sections.outline{ks});
    [bml, hTop] = upperEnvelope(o(:, 1), o(:, 2), opts.bin);
    tAt = interp1(ref.yy, refSliceTop(ref, sections.AP(ks(1)), rot), bml);
    use = ~isnan(tAt) & ~isnan(hTop);
    sections.n_bins(ks) = nnz(use);
    if nnz(use) < 5
        for k = ks', sections.outline{k} = []; end
        continue
    end
    d = tAt(use) - hTop(use);
    sections.shift(ks) = median(d);
    sections.fit_mad(ks) = median(abs(d - median(d)));
    for k = ks', sections.outline{k}(:, 2) = sections.outline{k}(:, 2) + median(d); end
end

% Fill outlines into the reference grid, padded by grid_pad in ML and DV (the
% bird's brain can be wider than the reference). One mask per section (A/B
% halves together, so they join across the midline): filled outline, closed
% with a disk of radius close_um (fills cracks and infoldings narrower than its
% diameter), holes filled, opened with radius open_um (removes protrusions
% narrower than its diameter) -> outline_smooth. Then, so a single bad section
% is outvoted by its neighbors: the median, over the median_sections sections
% nearest in AP (only sections that exist), of each mask's signed distance to
% its edge; inside where < 0 -> outline_median, filled into the planes within
% half a section of the AP. median_change_mm2: area that differs between the
% two. Outlines are [ML DV], several polygons separated by NaN rows.
grid_pad = 1.5;   % mm
ml_max = ceil((max(abs(ref.yy)) + grid_pad)/dI)*dI; yy = -ml_max:dI:ml_max;
zz = flip(ref.zz(1):-dI:ref.zz(1) - grid_pad); zz = [zz(1:end-1), ref.zz, ref.zz(end) + (dI:dI:grid_pad)];
Ibrain = false(numel(yy), numel(ref.xx), numel(zz));
se_close = strel('disk', round(opts.close_um/1000/dI));
se_open = strel('disk', round(opts.open_um/1000/dI));
aligned = ~isnan(sections.shift);
gks = unique(group(aligned))';
masks = cell(size(gks)); sd = cell(size(gks)); ap = zeros(size(gks));
for ig = 1:numel(gks)
    ks = find(aligned & group == gks(ig))';
    mask = false(numel(yy), numel(zz));
    for k = ks, mask = mask | fillPolygons(sections.outline{k}, yy, zz); end
    masks{ig} = imopen(imfill(imclose(mask, se_close), 'holes'), se_open);
    sd{ig} = single(bwdist(masks{ig}) - bwdist(~masks{ig}));   % signed distance, < 0 inside
    ap(ig) = sections.AP(ks(1));
end
sections.outline_smooth = cell(n, 1);
sections.outline_median = cell(n, 1);
sections.median_change_mm2 = nan(n, 1);
for ig = 1:numel(gks)
    ks = find(aligned & group == gks(ig))';
    [~, near] = mink(abs(ap - ap(ig)), opts.median_sections);
    mask = median(cat(3, sd{near}), 3) < 0;
    sections.outline_smooth(ks) = {maskOutline(masks{ig}, yy, zz)};
    sections.outline_median(ks) = {maskOutline(mask, yy, zz)};
    sections.median_change_mm2(ks) = nnz(xor(mask, masks{ig})) * dI^2;
    for p = find(abs(ref.xx - ap(ig)) <= opts.section_step/2 + 1e-9)   % a plane halfway between gets both
        Ibrain(:, p, :) = squeeze(Ibrain(:, p, :)) | mask;
    end
end

T.xx = ref.xx; T.yy = yy; T.zz = zz;
T.Ibrain = Ibrain;
T.P = struct('bird', bird.name, 'dI', dI, 'default_angle', opts.bird_angle, ...
    'ref_angle', opts.ref_angle, 'ap_shift', opts.ap_shift, 'section_step', opts.section_step, 'ap_mid', ap_mid, 'bin', opts.bin, ...
    'built', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm')));
T.sections = sections;
end


function problem = sectionProblem(bird, k)
% Why section k can't be used as is ("" if it can): outline missing or narrower
% than min_width; no midline; or a drawn midline lying off the tissue (median
% distance of its points from the outline, counting points inside as 0). Whole
% sections allow max_gap; A/B images (one hemisphere, midline = medial wall)
% allow max_gap_hemi.
min_width = 1;         % mm
max_gap = 0.2;         % mm
max_gap_hemi = 0.05;   % mm
problem = "";
e = bird.meta(k).edge;
if isempty(e) || bird.px(k) * bird.edge_scale * (max(e(:, 2)) - min(e(:, 2))) < min_width   % max/min skip the NaN separators
    problem = "outline missing or broken"; return
end
m = bird.meta(k).midLineCorr;
if ~ismember(size(m, 1), [2 4]) || size(m, 2) ~= 2
    problem = "no midline"; return
end
if size(m, 1) == 2
    E = densify(bird.edge_scale * e(:, [2 1]), 0.01 / bird.px(k));   % coarse outlines: measure to the polygon, not its vertices
    seg = m(1, :) + linspace(0, 1, 50)' .* (m(2, :) - m(1, :));
    d = min(pdist2(seg, E), [], 2) * bird.px(k);
    d(inpolygon(seg(:, 1), seg(:, 2), E(:, 1), E(:, 2))) = 0;
    gap = median(d);
    if gap > max_gap || (bird.single_hemi(k) && gap > max_gap_hemi)
        problem = sprintf("midline %.2f mm from the outline", gap);
    end
end
end


function X = joinPoints(X, jn, ds)
% Full-res pixel positions X = [x y] in a section image -> the same positions
% after its detached pieces were joined (jn: lab, steps from findJoin, at 1/ds
% scale, [row col]): each point follows the piece it lies on through every step
% that moved that piece. Points off all pieces take the nearest piece.
if isempty(X) || isempty(jn.steps), return; end
P = (X(:, [2 1]) - 0.5)/ds + 0.5;                % [row col] at the scale of lab
[~, idx] = bwdist(jn.lab > 0);
ri = min(max(round(P(:, 1)), 1), size(jn.lab, 1));
ci = min(max(round(P(:, 2)), 1), size(jn.lab, 2));
piece = jn.lab(idx(sub2ind(size(jn.lab), ri, ci)));
for s = jn.steps
    on = ismember(piece, s.moved);
    R = [cosd(s.theta) -sind(s.theta); sind(s.theta) cosd(s.theta)];
    P(on, :) = (P(on, :) - s.ctr) * R' + s.ctr + s.shift;
end
X = (P(:, [2 1]) - 0.5)*ds + 0.5;
end


function [mid, info] = flatSideMidline(E, side, px)
% Midline of a single hemisphere along its flat medial side. For each angle within
% +-max_angle of vertical: the line touching the outline on that side, and the
% longest unbroken stretch of outline within tol of it. The angle with the longest
% stretch wins (ties: smaller RMS); the line is then refitted to that stretch.
% Length = extent of the stretch along the line. tol starts at tol_um and grows by
% tol_step (up to max_tol_um) until that length is at least min_frac of the
% outline's dorsal-ventral height (a straight medial wall, not a short flat patch).
% E: outline [x y] full-res pixels (ordered); side "left"/"right"; px mm/pixel.
% mid: 2x2 [dorsal; ventral]; info: length_mm, rms_um, tol_um, points.
tol_um = 50; tol_step = 25; max_tol_um = 200; max_angle = 30; min_frac = 0.6;
E = densify(E, 0.01 / px);                               % ~10 um spacing
height = (max(E(:, 2)) - min(E(:, 2))) * px;
best = longestStretch(E, side, px, tol_um, max_angle);
while best.len < min_frac * height && tol_um + tol_step <= max_tol_um
    tol_um = tol_um + tol_step;
    best = longestStretch(E, side, px, tol_um, max_angle);
end
Q = E(best.idx, :);
c = mean(Q, 1);
[V, ~] = eig(cov(Q)); d = V(:, 2)';                       % refit line to the stretch
t = (Q - c) * d';
ends = c + [min(t); max(t)] .* d;
[~, order] = sort(ends(:, 2));                            % dorsal (smaller y) first
mid = ends(order, :);
info = struct('length_mm', best.len, 'rms_um', best.rms, 'tol_um', tol_um, 'points', Q);
end


function best = longestStretch(E, side, px, tol_um, max_angle)
% Over angles: longest unbroken outline stretch within tol_um of the line touching
% the outline on the given side, length = its extent along the line (see
% flatSideMidline)
best = struct('len', -1, 'rms', inf);
for th = -max_angle:1:max_angle
    nrm = [cosd(th) -sind(th)];                          % normal, pointing image right
    along = [sind(th) cosd(th)];                         % line direction
    proj = E * nrm';
    if side == "left", dist = proj - min(proj); else, dist = max(proj) - proj; end
    near = dist * px * 1000 < tol_um;
    lab = bwlabel([near; near]);                          % doubled to handle wrap-around
    runs = unique(lab(lab > 0))';
    for r = runs
        idx = find(lab == r);
        if numel(idx) > numel(near), idx = idx(1:numel(near)); end
        idx = mod(idx - 1, numel(near)) + 1;
        len = range(E(idx, :) * along') * px;
        rms = sqrt(mean(dist(idx).^2)) * px * 1000;
        if len > best.len + 1e-9 || (abs(len - best.len) < 1e-9 && rms < best.rms)
            best = struct('len', len, 'rms', rms, 'idx', idx);
        end
    end
end
end


function drawCorrections(bird, ks)
% Open each section in ks for hand-drawing; saves to <bird dir>\midlineCorrections.mat
% (replacing earlier entries for the same section). Broken outlines: draw a new
% outline first (click points, double-click the polygon to accept). Midline: drag
% the line ends (dorsal end first) and double-click the line to accept. A/B images
% start with a line on the medial side (A: image left, B: image right). Images
% with joined pieces are shown joined, so corrections are in joined pixels.
corr_file = fullfile(bird.dir, 'midlineCorrections.mat');
corrections = struct('k', {}, 'midLineCorr', {}, 'edge', {}, 'mirror', {});
if isfile(corr_file), corrections = load(corr_file).corrections; end
ds = bird.edge_scale;
for k = ks
    [img, full_size] = readSectionImage(bird.image_files{k});
    j = bird.join(k);
    if ~isempty(j.steps), img = joinPieces(img, j.lab, j.steps, exp(j.peaks(1))); end
    img = double(img);
    fig = figure('Name', sprintf('%s section %d', bird.name, k), 'Position', [100 100 1000 800]);
    imagesc([1 full_size(2)], [1 full_size(1)], img, prctile(img(:), [1 99.5])); colormap(gray); axis image; hold on
    [~, fname] = fileparts(bird.image_files{k});
    new_edge = [];
    e = bird.meta(k).edge;
    if sectionProblem(bird, k) == "outline missing or broken"
        if ~isempty(e), plot(ds*e(:, 2), ds*e(:, 1), 'r-'); end
        title(sprintf('%s (k=%d): click around the brain outline; double-click the polygon to accept', fname, k), 'Interpreter', 'none')
        roi = drawpolygon(gca, 'Color', 'c'); wait(roi);
        new_edge = roi.Position(:, [2 1]) / ds;          % [row col] at 1/ds, as metaData.edge
        E = roi.Position;
    else
        E = ds * e(:, [2 1]);
    end
    plot(E([1:end 1], 1), E([1:end 1], 2), 'c-')
    m = bird.meta(k).midLineCorr;
    if isequal(size(m), [2 2])
        plot(m(:, 1), m(:, 2), 'y--', 'LineWidth', 1); start = m;
    elseif bird.single_hemi(k)
        x_med = min(E(:, 1)) * (bird.ab(k) == "A") + max(E(:, 1)) * (bird.ab(k) == "B");
        start = [x_med min(E(:, 2)); x_med max(E(:, 2))];
    else
        start = [mean(E(:, 1)) min(E(:, 2)); mean(E(:, 1)) max(E(:, 2))];
    end
    title(sprintf('%s (k=%d): drag the midline (dorsal end first); double-click it to accept', fname, k), 'Interpreter', 'none')
    roi = drawline(gca, 'Position', start, 'Color', 'y'); wait(roi);
    old = corrections([corrections.k] == k);
    mirror = []; if ~isempty(old), mirror = old.mirror; end   % keep a hand-set mirror
    corrections([corrections.k] == k) = [];
    corrections(end+1) = struct('k', k, 'midLineCorr', roi.Position, 'edge', new_edge, 'mirror', mirror); %#ok<AGROW>
    save(corr_file, 'corrections');
    fprintf('  section %d saved to %s\n', k, corr_file);
end
end


function top = refSliceTop(ref, ap, rot)
% Dorsal surface (max DV per ML) of the reference, rotated rot deg nose down, in
% the coronal plane AP = ap of the rotated frame. Nearest-voxel sampling.
[Yq, Zq] = meshgrid(ref.yy, ref.zz);                       % zz x yy
[Xn, Zn] = rotateNoseDown(ap*ones(size(Zq)), Zq, -rot);          % back to the reference frame
dI = ref.xx(2) - ref.xx(1);
ix = round((Xn - ref.xx(1))/dI) + 1; iy = round((Yq - ref.yy(1))/(ref.yy(2) - ref.yy(1))) + 1; iz = round((Zn - ref.zz(1))/dI) + 1;
ok = ix >= 1 & ix <= numel(ref.xx) & iz >= 1 & iz <= numel(ref.zz);
S = false(size(Zq));
S(ok) = ref.Ibrain(sub2ind(size(ref.Ibrain), iy(ok), ix(ok), iz(ok)));
top = nan(1, numel(ref.yy));
has = any(S, 1);
[~, i_top] = max(flipud(S), [], 1);
top(has) = ref.zz(end + 1 - i_top(has));
end


function d = overlapWithRef(T, ref)
% Dice overlap of the template volume with the reference brain rotated into the
% template's frame (nearest voxel on T's grid), over the AP planes T covers
planes = find(squeeze(any(any(T.Ibrain, 1), 3)));
[Y, Z] = ndgrid(T.yy, T.zz);
rot = T.P.default_angle - T.P.ref_angle;
dI = [ref.xx(2) - ref.xx(1), ref.yy(2) - ref.yy(1), ref.zz(2) - ref.zz(1)];
n_and = 0; n_sum = 0;
for p = planes
    [Xn, Zn] = rotateNoseDown(T.xx(p)*ones(size(Z)), Z, -rot);
    ix = round((Xn - ref.xx(1))/dI(1)) + 1; iy = round((Y - ref.yy(1))/dI(2)) + 1; iz = round((Zn - ref.zz(1))/dI(3)) + 1;
    ok = ix >= 1 & ix <= numel(ref.xx) & iy >= 1 & iy <= numel(ref.yy) & iz >= 1 & iz <= numel(ref.zz);
    B = false(size(Y));
    B(ok) = ref.Ibrain(sub2ind(size(ref.Ibrain), iy(ok), ix(ok), iz(ok)));
    A = squeeze(T.Ibrain(:, p, :));
    n_and = n_and + nnz(A & B); n_sum = n_sum + nnz(A) + nnz(B);
end
d = 2 * n_and / n_sum;
end


function mask = fillPolygons(o, yy, zz)
% Mask on the grid (yy x zz) inside the polygons o = [ML DV] (several separated
% by NaN rows)
mask = false(numel(yy), numel(zz));
dI = yy(2) - yy(1);
brk = [0; find(isnan(o(:, 1))); size(o, 1) + 1];
for i = 1:numel(brk) - 1
    p = o(brk(i)+1:brk(i+1)-1, :);
    if size(p, 1) < 3, continue; end
    mask = mask | poly2mask((p(:, 2) - zz(1))/dI + 1, (p(:, 1) - yy(1))/dI + 1, numel(yy), numel(zz));
end
end


function o = maskOutline(mask, yy, zz)
% Outer boundaries of a mask on the grid (yy x zz) as closed polygons [ML DV],
% separated by NaN rows
B = bwboundaries(mask, 'noholes');
o = cell2mat(cellfun(@(b) [yy(b(:, 1))', zz(b(:, 2))'; NaN NaN], B(:), 'UniformOutput', false));
if ~isempty(o), o(end, :) = []; end
end


function [centers, top] = upperEnvelope(x, h, bin)
% Highest h in each x bin
edges = floor(min(x)/bin)*bin : bin : max(x) + bin;
[~, ~, b] = histcounts(x, edges);
ok = b > 0;
top = accumarray(b(ok), h(ok), [numel(edges) - 1, 1], @max, NaN);
centers = edges(1:end-1)' + bin/2;
end


function P = densify(P, step)
% Resample closed polygons (several separated by NaN rows) so consecutive points
% are at most step apart
if any(isnan(P(:, 1)))
    brk = [0; find(isnan(P(:, 1))); size(P, 1) + 1];
    parts = arrayfun(@(i) [densify(P(brk(i)+1:brk(i+1)-1, :), step); NaN NaN], 1:numel(brk) - 1, 'UniformOutput', false);
    P = vertcat(parts{:}); P(end, :) = [];
    return
end
Q = P([1:end 1], :);
seg = vecnorm(diff(Q), 2, 2);
out = cell(numel(seg), 1);
for i = 1:numel(seg)
    t = (0:max(1, ceil(seg(i)/step)) - 1)' / max(1, ceil(seg(i)/step));
    out{i} = Q(i, :) + t .* (Q(i + 1, :) - Q(i, :));
end
P = vertcat(out{:});
end


function tf = fitDryadTransform(T, bird, dry)
% tf = fitDryadTransform(T, bird, dry) fits the rotation and DV offset that map a
% bird's Dryad DV onto its DV in template T (built by buildBirdTemplate from the
% same bird's sections). The result is applied with dryadToTemplate, also to birds
% that have no section outlines.
%
% bird  from loadAnatomyBird; uses bird.labelled (per wavelength: .cells with
%       .pos/.label per section, .wavelength) and bird.section_ap
% dry   Dryad rows for this bird: AP, ML, DV, wavelength
%
% Per section, cells with label > 0 are mapped from pixels with
% pixelsToTemplate and compared with the Dryad cells at the same AP and
% wavelength by their mean DV (row order differs, so cells are not paired).
% Model: frame DV = DV rotated theta deg nose down (after AP is scaled and shifted
% like the sections, see dryadToTemplate) + c.
%
% tf: theta, c, ap_shift, ap_scale, ap_mid, rms, plus per-section table (AP, n, mean DV in frame,
% mean Dryad DV, residual, mean ML from pixels and from Dryad).
%
% Hannah Payne, 2026

rows = {};
for w = 1:numel(bird.labelled)
    L = bird.labelled(w).cells;
    dw = dry(string(dry.wavelength) == bird.labelled(w).wavelength, :);
    for k = 1:numel(L.pos)
        lab = L.label{k} > 0;
        if ~any(lab) || isnan(T.sections.shift(k)), continue; end
        q = abs(dw.AP - bird.section_ap(k)) < 1e-6;
        if ~any(q), continue; end
        [ml, dv] = pixelsToTemplate(L.pos{k}(lab, :), T.sections(k, :));
        rows(end+1, :) = {k, bird.labelled(w).wavelength, bird.section_ap(k), nnz(lab), ...
            mean(dv), mean(dw.DV(q)), mean(dw.AP(q)), mean(ml), mean(dw.ML(q))}; %#ok<AGROW>
    end
end
S = cell2table(rows, 'VariableNames', {'k', 'wavelength', 'AP', 'n', 'dv_frame', 'dv_dryad', 'ap_dryad', 'ml_pixels', 'ml_dryad'});

tf = struct('theta', 0, 'c', 0, 'ap_shift', T.P.ap_shift, ...
    'ap_scale', T.P.section_step / 0.1, 'ap_mid', T.P.ap_mid);
gap = @(th) S.dv_frame - modelDV(S.ap_dryad, S.dv_dryad, tf, th);
cost = @(th) sum(S.n .* (gap(th) - wmean(gap(th), S.n)).^2);
tf.theta = fminbnd(cost, -30, 30);
tf.c = wmean(gap(tf.theta), S.n);
S.residual = gap(tf.theta) - tf.c;
tf.rms = sqrt(wmean(S.residual.^2, S.n));
tf.per_section = S;
end


function dv = modelDV(AP, DV, tf, theta)
% Frame DV of Dryad (AP, DV) under tf with rotation theta and no DV offset
tf.theta = theta; tf.c = 0;
[~, ~, dv] = dryadToTemplate(AP, zeros(size(AP)), DV, tf);
end


function m = wmean(x, w)
m = sum(w .* x) / sum(w);
end
