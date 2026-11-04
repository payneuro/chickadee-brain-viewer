function brainViewer
%% Interactive viewer: chickadee brain templates and retrogradely labeled cells
% Tab 1: coronal, sagittal and horizontal slices plus a 3D view (top left) of
% a template brain and the selected cells. Sliders set the slice positions; the
% table picks injection sites and labeled regions. "Probe track": click a point in
% the brain on the coronal slice to plan a track to it at the set roll angle; the
% table below gives the target, its depth below the dorsal surface, the entrance
% point and the depth along the track, and the track is drawn in all four views.
% Tab 2: MC022 coronal sections, one at a time, as placed in the MC022 frame
% (detached pieces joined, A/B halves combined; renderSection), in one of four
% channels, with the outline and MC022's labelled cells. Slider, spinner or the
% scroll wheel step through sections; the section follows tab 1's AP. The
% region/site table and display switches are shared by both tabs.
% Measuring box: square button on the slice and section axes toolbars.
%
% Frames. Each template is stored at its own skull-groove angle (MC022.mat P):
%   MC029_fixed brain + hippocampus + DL/CDL, at 65 deg (P.ref_angle, nominal and
%               trusted), ML stretched by P.section_step/0.1
%   MC022       brain outline only (built by makeAnatomyTemplates), sectioned at 85 deg
%               (P.default_angle): MC029's angle + the tilt fitted to
%               MC022's sections (fitReferencePlacement: 19.9, rounded to 20). MC022 was injected at
%               75 deg; a skull groove of 85 shows its sections flat.
% Cells: Dryad AP/ML/DV mapped into the MC022 frame by dryadToTemplate,
% with the transform fitted on MC022 (stored in MC022.mat); MC022's own cells
% come from its section images instead (labelled regions only). Everything is
% rotated about lambda to the "Skull groove" angle for display (default
% display_angle, 75 deg, the injection angle of most anatomy birds including
% MC022; 65 deg = 30 deg bite bar).
%
% Data (in dataDir: .\data, the shared copy on Z:, or downloaded): MC022.mat,
% MC029_fixed.mat, MC029_orig.mat, dryad_applegate2023\anatomy.mat (the Dryad
% cells, read by loadDryadCells). Volumes, 3D surfaces and
% cell coordinates are cached in tempdir (loadViewerCache). Tab 2 reads
% MC022.mat's sections table and the 1/8 images in MC022_ds8\.
% Needs loadBrain3d, plotBrainRegion3d, cellDensity3d, loadDryadCells
% dryadToTemplate, rotateNoseDown, renderSection
% (pixelsToTemplateImage, readSectionImage, joinPieces).
%
% Hannah Payne, 2026

addpath(fullfile(fileparts(mfilename('fullpath')), 'src'));
data_dir = dataDir();

%% Data (precomputed once, see loadViewerCache)
C = loadViewerCache(data_dir);
tmpl = C.tmpl;                          % templates, each in its native frame
cells = C.cells;
ap_frame = C.ap; ml_frame = C.ml; dv_frame = C.dv;   % cells in the MC022 frame
frame_angle = tmpl(1).native_angle;     % skull-groove angle of the cell frame (MC022's)
display_angle = 75;                     % initial "Skull groove" (deg): the injection angle
tf = C.tf;
for it0 = 1:numel(tmpl)                 % one interpolant per volume, for tilted slices
    for nm = ["Ibrain" "Ihippo" "Idl"]
        if isfield(tmpl(it0).D, nm)
            tmpl(it0).F.(nm) = griddedInterpolant({tmpl(it0).D.yy, tmpl(it0).D.xx, tmpl(it0).D.zz}, ...
                tmpl(it0).D.(nm), 'linear', 'none');
        end
    end
end
it_cur = 1;                             % current template

% Display grid: the template grid, padded so rotated brains stay inside
G = tmpl(1).D; dI = G.xx(2) - G.xx(1);
Dd.xx = (G.xx(1) - 1) : dI : (G.xx(end) + 1);
Dd.yy = G.yy;
Dd.zz = (G.zz(1) - 1) : dI : (G.zz(end) + 2);

draw_key = rand(RandStream('twister', 'Seed', 0), numel(ap_frame), 1);   % fixed random subset for "Cells drawn" (own stream, global one untouched)
slider_step = 0.1;                      % sliders snap to this; cell AP sits on 100 um sections
% AP snaps to the anatomy birds' coronal sections (spacing and position as in
% dryadToTemplate), ML and DV to slider_step
step_dim = slider_step * [tf.ap_scale, 1, 1];   % AP, ML, DV
origin_dim = [tf.ap_mid*(1 - tf.ap_scale) + tf.ap_shift, 0, 0];   % frame AP of Dryad AP 0

% Injection sites, Table 1 of Applegate et al. 2023. For DL, Long is the
% topographically matched Hp location; AP is the actual injection coordinate.
% AP and ML (mm) from lambda, depth (mm) below the brain surface, at head_angle.
sites = table(["Hp";"Hp";"Hp";"Hp";"DL";"DL";"DL"], [-1;1;3;5;1;3;5], ...
    (1:7)', [-1;1;3;5;2;3.5;5.5], [1.9;0.4;0.4;0.4;2.3;1.7;0.7], [0.4;0.4;0.4;1.0;0.4;0.4;0.3], [90;75;75;75;75;75;75], ...
    'VariableNames', {'Transverse','Long','Loc','AP','ML','depth','head_angle'});
site_labels = compose('%s AP %.1f', upper(sites.Transverse), sites.AP);   % HP, DL (data say Hp)
[site_ap, site_dv] = sitePositions(sites, tf, tmpl(1).F.Ibrain, frame_angle, Dd.zz);   % MC022 frame
site_ap_d = site_ap; site_dv_d = site_dv;   % display frame (updateCellDisplay)
site_window = 0.5;   % mm: a slice or section shows the sites within this distance of its plane
[~, cell_site] = ismember(string(cells.Transverse) + "_" + cells.Long, ...
    sites.Transverse + "_" + sites.Long);

% Labeled regions, pallial first
pallial = ["DL/CDL","HD","NFL","NIL_M","NIL_L"];
roi_all = string(categories(categorical(cells.ROI)));
regions = [pallial(ismember(pallial, roi_all)), setdiff(roi_all, pallial, 'stable')'];
[~, cell_region] = ismember(string(cells.ROI), regions);

example_birds = ["MC017","MC025"];
[birds, ~, cell_bird] = unique(string(cells.bird));

% Colors
color_brain = .5*[1 1 1];
color_hippo = [190 30 45]/255;
color_dl    = [247 148 29]/255;
% Regions: pallial ones set by hand (DL/CDL red, the others cool); the rest
% spread along turbo without its dark ends (they vanish on black section images)
pallial_colors = [0.90 0.15 0.15; 0.55 0.25 0.80; 0.00 0.65 0.60; 0.10 0.20 0.65; 0.45 0.70 1.00];   % DL/CDL, HD, NFL, NIL_M, NIL_L
is_pal = ismember(regions, pallial);
region_colors = zeros(numel(regions), 3);
[~, i_pal] = ismember(regions(is_pal), pallial);
region_colors(is_pal, :) = pallial_colors(i_pal, :);
region_colors(~is_pal, :) = interp1(linspace(0, 1, 256), turbo(256), linspace(0.15, 0.92, nnz(~is_pal)));
% Injection sites: Hp warm, DL cool, each location its own hue, posterior to anterior
site_colors = [0.85 0.10 0.75; 0.90 0.20 0.15; 1.00 0.60 0.00; 0.85 0.75 0.00; ...   % Hp -1, 1, 3, 5: magenta, red, orange, gold
               0.15 0.20 0.75; 0.00 0.70 0.80; 0.20 0.70 0.25];                       % DL 2, 3.5, 5.5: blue, teal, green
blob_level = 0.1;   % isosurface level for density blobs

%% Layout
scr = get(groot, 'ScreenSize');                    % main screen; opens maximized, restores to at most 85% of it
fig_size = round(min([1300 860], 0.85*scr(3:4)));
fig = uifigure('Name', 'Chickadee brain viewer', 'Position', [scr(1:2) + round((scr(3:4) - fig_size)/2), fig_size]);
fig.WindowState = 'maximized';
set(fig, 'defaultAxesToolbarExpanded', 'off');      % axes toolbars collapsed until hover
tabs = uitabgroup(uigridlayout(fig, [1 1], 'Padding', 0));
tab1 = uitab(tabs, 'Title', 'Template and cells');
tab2 = uitab(tabs, 'Title', 'MC022 sections');

%% Tab 2: sections in the MC022 frame
S2 = load(fullfile(data_dir, 'MC022.mat'), 'sections', 'injections');
secs = S2.sections;
% Lambda shift, as the cells (loadViewerCache): section APs and outlines here,
% the image's DV where it is made (h_sec_img)
secs.AP = secs.AP - C.lambda(1);
has_out = ~cellfun(@isempty, secs.outline_smooth);
secs.outline_smooth(has_out) = cellfun(@(p) p - [0 C.lambda(2)], secs.outline_smooth(has_out), 'UniformOutput', false);
sec_orders = unique(secs.trueOrder)';
[sec_ML, sec_DV] = meshgrid(-7:0.02:7, -10:0.02:2);
% Image channels (file suffix) and their buttons: tracer wavelength with its
% injection site (MC022's injections table), DAPI. Each channel has a color
% near its emission; MC022's sites take their tracer's color in tab 2.
sec_channels = ["640" "555" "475" "395"];
sec_channel_names = ["647" "555" "488" "DAPI"];
sec_channel_colors = [0.95 0.15 0.15; 1 0.6 0; 0.15 0.8 0.2; 0.3 0.45 1];   % red, orange, green, blue
sec_site_colors = site_colors;
for ic = 1:numel(sec_channel_names)
    [~, i_inj] = ismember(sec_channel_names(ic), S2.injections.wavelength);
    if i_inj == 0, continue; end
    is_site = sites.Transverse == S2.injections.Transverse(i_inj) & sites.Long == S2.injections.Long(i_inj);
    if ~any(is_site), continue; end
    sec_channel_names(ic) = sec_channel_names(ic) + ", " + site_labels{is_site};
    sec_site_colors(is_site, :) = sec_channel_colors(ic, :);
end
sec_cache = cell(numel(sec_orders), numel(sec_channels));   % rendered images, kept once shown
sec_ap = arrayfun(@(o) secs.AP(find(secs.trueOrder == o, 1)), sec_orders);   % frame AP of each section
color_outline = [0 0.8 1];

main2 = uigridlayout(tab2, [1 2]);
main2.ColumnWidth = {390, '1x'};
% Rows: section slider, the section's AP, image and trueOrder, channel, outline
% checkbox, the controls shared with tab 1 (row shared_row(2), moved here by
% tabChanged), histogram, display range.
ctrl2 = uigridlayout(main2, [7 2]);
ctrl2.RowHeight = {42, 44, 52, 26, '1x', 160, 42};
ctrl2.ColumnWidth = {120, '1x'};
[sld_sec, spn_sec] = sliderRow(ctrl2, sprintf('Section (of %d):', numel(sec_orders)), [1 numel(sec_orders)], 1, 1, '%.0f');
sld_sec.MajorTicks = 1:10:numel(sec_orders);
sld_sec.MajorTickLabels = string(sec_orders(sld_sec.MajorTicks));   % trueOrder
sld_sec.MinorTicks = [];
lbl_sec = uilabel(ctrl2, 'Text', '', 'WordWrap', 'on', 'VerticalAlignment', 'top');
lbl_sec.Layout.Row = 2; lbl_sec.Layout.Column = [1 2];
lbl_channel = uilabel(ctrl2, 'Text', 'Channel:', 'FontWeight', 'bold');
lbl_channel.Layout.Row = 3; lbl_channel.Layout.Column = 1;
bg_channel = uibuttongroup(ctrl2, 'BorderType', 'none');
bg_channel.Layout.Row = 3; bg_channel.Layout.Column = 2;
rb_channel = gobjects(1, numel(sec_channels));
for ic = 1:numel(sec_channels)
    rb_channel(ic) = uiradiobutton(bg_channel, 'Text', sec_channel_names(ic), ...
        'FontColor', 0.75*sec_channel_colors(ic, :), ...   % darker, readable on the panel
        'Position', [mod(ic-1, 2)*120, 28 - 26*floor((ic-1)/2), 118, 22]);   % 2 x 2
end
% Checkbox text in a darker shade of the line color, readable on the white panel
chk_sec_outline = uicheckbox(ctrl2, 'Text', 'Outline', 'Value', true, 'FontColor', 0.7*color_outline);
chk_sec_outline.Layout.Row = 4; chk_sec_outline.Layout.Column = [1 2];
% Display range: histogram of the section's image with a min/max slider below,
% both on one fixed log scale for all sections: log10(1 + intensity), from 0
% to clim_top. Image units: 0 = background, 1 = a default maximum (renderSection).
% The range set stays from section to section (drawSection).
clim_top = 100;          % top of the histogram and slider (intensity)
toLog = @(I) log10(1 + max(I, 0));
clim_ticks = [0 0.3 1 3 10 30 100];
ax_hist = uiaxes(ctrl2);
ax_hist.Layout.Row = 6; ax_hist.Layout.Column = [1 2];
set(ax_hist, 'FontSize', 12, 'TickDir', 'out', 'YScale', 'log', 'Box', 'off', 'Color', 'none', ...
    'XLim', toLog([0 clim_top]), 'XTick', toLog(clim_ticks), 'XTickLabel', string(clim_ticks));
ax_hist.YAxis.Visible = 'off';   % counts: only the shape matters
hold(ax_hist, 'on')
h_hist = stairs(ax_hist, NaN, NaN, 'k-');
h_clim = xline(ax_hist, toLog([0 1]), 'r-');
xlabel(ax_hist, 'Intensity (0 = background, 1 = default max)')
sld_clim = uislider(ctrl2, 'range', 'Limits', toLog([0 clim_top]), 'Value', toLog([0 1]), ...
    'MajorTicks', toLog(clim_ticks), 'MajorTickLabels', string(clim_ticks), 'MinorTicks', []);
sld_clim.Layout.Row = 7; sld_clim.Layout.Column = [1 2];

ax_sec = uiaxes(main2);
set(ax_sec, 'FontSize', 12, 'TickDir', 'out');
hold(ax_sec, 'on')
h_sec_img = imagesc(ax_sec, sec_ML(1, :), sec_DV(:, 1) - C.lambda(2), zeros(size(sec_ML)), [0 1]);
colormap(ax_sec, gray); ax_sec.YDir = 'normal';
xline(ax_sec, 0, ':', 'Color', [1 1 0.6]);
h_sec_out = plot(ax_sec, NaN, NaN, '-', 'Color', color_outline, 'LineWidth', 1);
h_sec_sites = siteMarkers(ax_sec, height(sites));   % injection sites near the section (drawSection)
h_sec_groups = gobjects(0);    % MC022's cells or density contours per group (drawSection)
h_sec_keys = gobjects(0);      % their legend stand-ins
groups2 = struct('idx', {}, 'show', {}, 'color', {}, 'name', {}, 'label', {}, 'I', {});   % MC022-only groups (updateGroups)
axis(ax_sec, 'image'); xlim(ax_sec, sec_ML(1, [1 end])); ylim(ax_sec, sec_DV([1 end], 1));
set(ax_sec, 'XTick', -7:7, 'YTick', -10:2);   % 1 mm
xlabel(ax_sec, 'ML (mm)'); ylabel(ax_sec, 'DV (mm)')

main = uigridlayout(tab1, [1 2]);
main.ColumnWidth = {390, '1x'};

% Rows: 3 sliders, template, skull groove, 2 checkbox rows, probe track (2 rows),
% the controls shared with tab 2 (row shared_row(1)), birds, bird list
ctrl = uigridlayout(main, [12 2]);
ctrl.RowHeight = [repmat({42}, 1, 3), repmat({26}, 1, 5), {92, '1x', 26, 130}];   % sliders need room for tick labels
ctrl.RowSpacing = 6;
ctrl.ColumnWidth = {120, '1x'};

[sld_ap, spn_ap] = sliderRow(ctrl, 'AP (mm):', Dd.xx([1 end]), 4, step_dim(1), '%.2f');   % spinner steps one section
[sld_ml, spn_ml] = sliderRow(ctrl, 'ML (mm):', Dd.yy([1 end]), 2, step_dim(2), '%.2f');
[sld_dv, spn_dv] = sliderRow(ctrl, 'DV (mm):', Dd.zz([1 end]), -2, step_dim(3), '%.2f');
for sld = [sld_ap sld_ml sld_dv]
    sld.MajorTicks = ceil(sld.Limits(1)/2)*2 : 2 : sld.Limits(2);   % every 2 mm
    sld.MinorTicks = [];
end
uilabel(ctrl, 'Text', 'Template:', 'FontWeight', 'bold');
dd_tmpl = uidropdown(ctrl, 'Items', {tmpl.label}, 'ItemsData', 1:numel(tmpl), 'Value', it_cur);
uilabel(ctrl, 'Text', 'Skull groove (deg):', 'FontWeight', 'bold', 'Tooltip', 'Skull groove 65 deg above horizontal = beak bar 30 deg below horizontal');
spn_angle = uispinner(ctrl, 'Value', display_angle, 'Step', 5, 'Limits', [30 100], ...
    'Tooltip', 'Skull groove 65 deg above horizontal = beak bar 30 deg below horizontal');
% Checkboxes in two rows: outlines on the first, display options on the second
chk_grid_lay = uigridlayout(ctrl, [2 3], 'Padding', 0, 'RowSpacing', 6, ...
    'RowHeight', {26, 26}, 'ColumnWidth', {'fit', 'fit', 'fit'}, 'ColumnSpacing', 16);
chk_grid_lay.Layout.Row = [6 7]; chk_grid_lay.Layout.Column = [1 2];
chk_brain = uicheckbox(chk_grid_lay, 'Text', 'Brain', 'Value', true);
chk_hippo = uicheckbox(chk_grid_lay, 'Text', 'Hippocampus', 'Value', true, 'FontColor', color_hippo);
chk_dl    = uicheckbox(chk_grid_lay, 'Text', 'DL/CDL', 'Value', false, 'FontColor', color_dl);
chk_planes = uicheckbox(chk_grid_lay, 'Text', 'Slices in 3D', 'Value', true);
chk_grid  = uicheckbox(chk_grid_lay, 'Text', '1 mm grid', 'Value', true);
chk_minor = uicheckbox(chk_grid_lay, 'Text', '0.2 mm minor grid', 'Value', false);
% Probe track: checkbox and roll angle, then a one-row table of the planned
% coordinates (display frame, mm from lambda; drawProbe)
color_probe = [0.85 0 0.6];
probe_lay = uigridlayout(ctrl, [2 3], 'Padding', 0, 'RowSpacing', 6, 'ColumnSpacing', 6, ...
    'RowHeight', {26, 92}, 'ColumnWidth', {'1x', 'fit', 70});
probe_lay.Layout.Row = [8 9]; probe_lay.Layout.Column = [1 2];
chk_probe = uicheckbox(probe_lay, 'Text', 'Probe track (click coronal)', 'Value', false, ...
    'FontColor', 0.8*color_probe, 'Tooltip', 'Click a point in the brain on the coronal slice to plan a track to it');
uilabel(probe_lay, 'Text', 'Roll (deg):', 'FontWeight', 'bold');
spn_roll = uispinner(probe_lay, 'Value', 0, 'Step', 5, 'Limits', [-60 60], ...
    'Tooltip', 'Angle from vertical in the coronal plane; + aims toward the midline');
tbl_probe = uitable(probe_lay, 'Data', zeros(0, 7), 'RowName', {}, ...
    'ColumnName', {'AP', 'ML', 'DV', 'Roll', sprintf('Entr.\nAP'), sprintf('Entr.\nML'), 'Depth'}, ...
    'ColumnFormat', repmat({'bank'}, 1, 7), 'ColumnWidth', repmat({'1x'}, 1, 7), 'ColumnEditable', true(1, 7), ...
    'Tooltip', sprintf(['Target AP, ML (mm from lambda, at the skull groove when clicked); DV below the dorsal\n' ...
    'surface straight above it; roll (+ toward midline); entrance on the brain surface;\n' ...
    'depth along the track from the entrance. Edit any value: AP, ML or DV moves the target; roll, entrance or depth moves the target along the edited track.']));
tbl_probe.Layout.Row = 2; tbl_probe.Layout.Column = [1 3];
probe_ml = []; probe_ap = []; probe_dv = [];   % target in the MC022 frame, like the cells (probeClick)
probe_angle = [];   % skull groove the track was planned at; its table is in that frame
% Controls used by both tabs, in one block that tabChanged moves to the tab showing
shared_row = [10 5];   % its row in ctrl (tab 1) and ctrl2 (tab 2)
% Columns: label, switch's left option (right-aligned, so the switches line up),
% switch, right option
shared = uigridlayout(ctrl, [6 4], 'Padding', 0, 'RowSpacing', 6, 'ColumnSpacing', 6, ...
    'ColumnWidth', {'fit', 'fit', 45, '1x'}, 'RowHeight', {26, 26, 26, 26, 26, '1x'});
shared.Layout.Row = shared_row(1); shared.Layout.Column = [1 2];
sw_display  = switchRow(shared, 1, 'Display:', {'Cells', 'Density'}, 'Cells');
sw_fraction = switchRow(shared, 2, 'Cells drawn:', {'100%', '10%'}, '100%');
sw_color    = switchRow(shared, 3, 'Color by:', {'Region', 'Injection site'}, 'Region');
sw_sites    = switchRow(shared, 4, 'Show injection sites:', {'Off', 'On'}, 'Off');
sw_sites.Tooltip = sprintf('Stars at the sites selected in the table; slices show sites within %.1f mm', site_window);
% Regions: buttons set which rows of the table are checked, in the site columns
% in use (all sites if none are)
lbl_regions = uilabel(shared, 'Text', 'Regions:', 'FontWeight', 'bold');
lbl_regions.Layout.Row = 5; lbl_regions.Layout.Column = 1;
region_btns = uigridlayout(shared, [1 3], 'Padding', 0, 'ColumnSpacing', 4);
region_btns.Layout.Row = 5; region_btns.Layout.Column = [2 4];
btn_reg_all  = uibutton(region_btns, 'Text', 'All');
btn_reg_pal  = uibutton(region_btns, 'Text', 'Pallial', 'Tooltip', strjoin(pallial, ', '));
btn_reg_none = uibutton(region_btns, 'Text', 'None');
% Selection table: one row per labeled region, one column per injection site
% (headers: HP or DL / AP). Clicking a region name checks or clears
% its whole row; names of fully checked rows are shaded.
sel0 = false(numel(regions), height(sites));
sel0(regions == "HD", 2:4) = true;                     % start with HD, Hp locations 2-4
site_headers = cellstr(compose("%s\n%.1f", upper(sites.Transverse), sites.AP))';
tbl = uitable(shared, 'Data', [table(regions') array2table(sel0)], ...
    'ColumnName', [{sprintf('Region\nAP (mm)')}, site_headers], 'RowName', {}, ...
    'ColumnEditable', [false true(1, height(sites))], ...
    'ColumnWidth', [{80}, repmat({'1x'}, 1, height(sites))]);   % site columns share the rest of the width
tbl.Layout.Row = 6; tbl.Layout.Column = [1 4];
style_row_on = uistyle('BackgroundColor', [0.80 0.88 1.00], 'FontWeight', 'bold');
% Birds: multi-select list; counts (set in updateGroups) are cells in the table selection
lbl_birds = uilabel(ctrl, 'Text', 'Birds:', 'FontWeight', 'bold');
lbl_birds.Layout.Row = shared_row(1) + 1; lbl_birds.Layout.Column = 1;
bird_btns = uigridlayout(ctrl, [1 3], 'Padding', 0, 'ColumnSpacing', 4);
btn_all  = uibutton(bird_btns, 'Text', 'All');
btn_ex   = uibutton(bird_btns, 'Text', 'Examples', 'Tooltip', strjoin(example_birds, ', '));
btn_none = uibutton(bird_btns, 'Text', 'None');
lst_birds = uilistbox(ctrl, 'Multiselect', 'on', 'Items', cellstr(birds), ...
    'ItemsData', 1:numel(birds), 'Value', find(birds == "MC022"));   % other birds' Dryad placement is still off
lst_birds.Layout.Column = [1 2];

axes_grid = uigridlayout(main, [2 2]);
ax_cor = uiaxes(axes_grid);
ax_sag = uiaxes(axes_grid);
ax_hor = uiaxes(axes_grid);
ax_3d  = uiaxes(axes_grid);
% 3D top left, sagittal top right, coronal bottom left, horizontal bottom right
ax_3d.Layout.Row = 1;  ax_3d.Layout.Column = 1;
ax_sag.Layout.Row = 1; ax_sag.Layout.Column = 2;
ax_cor.Layout.Row = 2; ax_cor.Layout.Column = 1;
ax_hor.Layout.Row = 2; ax_hor.Layout.Column = 2;
set([ax_cor ax_sag ax_hor ax_3d], 'FontSize', 12, 'TickDir', 'out');   % 12 matches the uifigure controls
slice_ax = [ax_cor ax_sag ax_hor];
slice_lims = [];   % data limits the slice views must show (restoreAxes3d, scaleSlices)

%% 3D: template surfaces inside an hgtransform (rotated to the display angle)
ax_3d.YDir = 'reverse';   % +ML is the animal's right; with +AP anterior and +DV dorsal that is left-handed
h_tf = hgtransform(ax_3d);
h_brain = gobjects(1); h_hippo = gobjects(1); h_dl = gobjects(1);
h_planes = gobjects(3, 1);       % current coronal, sagittal, horizontal slices; drawSlices sets the data
for ip = 1:3                     % not kk: nested functions would share it
    h_planes(ip) = patch(ax_3d, NaN(4,1), NaN(4,1), NaN(4,1), .3*[1 1 1], 'FaceAlpha', 0.1, ...
        'EdgeColor', .3*[1 1 1], 'FaceLighting', 'none', 'HitTest', 'off');
end
h_groups3d = gobjects(0);
h_key3d = gobjects(0);       % legend stand-ins: cell markers drawn larger than in the plot
h_sites3d = siteMarkers(ax_3d, height(sites));   % injection sites (drawSites3d)
groups = struct('idx', {}, 'show', {}, 'color', {}, 'name', {}, 'label', {}, 'I', {});
cell_ap = []; cell_ml = []; cell_dv = [];   % cells in the display frame (updateCellDisplay)
label_margin = 1.3;       % 3D zoom: view this much wider than the box, for its tick labels (restoreAxes3d)
view_3d_default = [30 10];   % az, el: from the right side and a bit in front, slightly above, anterior right
view(ax_3d, view_3d_default);
drawTemplate3d();         % limits around the template brains (restoreAxes3d)

%% Slice views: graphics objects made once, then only their data is updated
outline_names = ["Ibrain" "Ihippo" "Idl"];
outline_colors = [color_brain; color_hippo; color_dl];
title_fmt = {'Coronal, AP %.2f mm', 'Sagittal, ML %.2f mm', 'Horizontal, DV %.2f mm'};
h_outline = gobjects(3, 3);     % slice x outline
h_cross = gobjects(3, 2);       % slice x (vertical, horizontal)
h_groups2d = {gobjects(0), gobjects(0), gobjects(0)};
h_sites2d = gobjects(3, height(sites));   % slice x site (setupSlices)
setupSlices();
% Probe track, entrance to target (circle), in each slice and in 3D (drawProbe)
h_probe = gobjects(4, 1);
probe_ax = [slice_ax ax_3d];
for ip = 1:4
    h_probe(ip) = line(probe_ax(ip),'XData', NaN, 'YData', NaN, 'ZData', NaN, ...
        'Color', color_probe, 'LineWidth', 2, 'Marker', 'o', 'MarkerIndices', 1, ...
        'MarkerFaceColor', color_probe, 'HitTest', 'off');
end

%% Measuring box: a button on each slice's and the section image's axes toolbar;
% drag out a box, labelled with its width and height in mm while it moves. One
% box at a time.
h_box = gobjects(0);
box_icon = ones(16, 16, 3); box_icon([3 14], 3:14, :) = 0; box_icon(3:14, [3 14], :) = 0;   % square outline
for ax_m = [ax_cor ax_sag ax_hor ax_sec]
    tb = axtoolbar(ax_m, {'export', 'datacursor', 'pan', 'zoomin', 'zoomout', 'restoreview'});
    axtoolbarbtn(tb, 'push', 'Icon', box_icon, 'Tooltip', 'Measure box (mm)', ...
        'ButtonPushedFcn', @(~, ~) measureBox(ax_m));
end

%% Callbacks
sld_ap.ValueChangingFcn = @(~, ev) sliderMoved(sld_ap, ev.Value, 1);   % AP moves the coronal slice
sld_ml.ValueChangingFcn = @(~, ev) sliderMoved(sld_ml, ev.Value, 2);
sld_dv.ValueChangingFcn = @(~, ev) sliderMoved(sld_dv, ev.Value, 3);
sld_ap.ValueChangedFcn = @(~, ev) sliderMoved(sld_ap, ev.Value, 1);
sld_ml.ValueChangedFcn = @(~, ev) sliderMoved(sld_ml, ev.Value, 2);
sld_dv.ValueChangedFcn = @(~, ev) sliderMoved(sld_dv, ev.Value, 3);
spn_ap.ValueChangingFcn = @(~, ev) sliderMoved(sld_ap, ev.Value, 1);   % spinners move their slider
spn_ml.ValueChangingFcn = @(~, ev) sliderMoved(sld_ml, ev.Value, 2);
spn_dv.ValueChangingFcn = @(~, ev) sliderMoved(sld_dv, ev.Value, 3);
spn_ap.ValueChangedFcn = @(~, ev) sliderMoved(sld_ap, ev.Value, 1);
spn_ml.ValueChangedFcn = @(~, ev) sliderMoved(sld_ml, ev.Value, 2);
spn_dv.ValueChangedFcn = @(~, ev) sliderMoved(sld_dv, ev.Value, 3);
dd_tmpl.ValueChangedFcn   = @(~, ~) templateChanged();
spn_angle.ValueChangedFcn = @(~, ~) angleChanged();
chk_brain.ValueChangedFcn = @(~, ~) toggleOutlines();
chk_hippo.ValueChangedFcn = @(~, ~) toggleOutlines();
chk_dl.ValueChangedFcn    = @(~, ~) toggleOutlines();
chk_grid.ValueChangedFcn  = @(~, ~) arrayfun(@applyGrid, [ax_cor ax_sag ax_hor ax_3d]);
chk_minor.ValueChangedFcn = @(~, ~) arrayfun(@applyGrid, [ax_cor ax_sag ax_hor ax_3d]);
chk_planes.ValueChangedFcn = @(~, ~) set(h_planes, 'Visible', chk_planes.Value);
sw_sites.ValueChangedFcn   = @(~, ~) updateGroups();
chk_probe.ValueChangedFcn  = @(~, ~) drawProbe();
spn_roll.ValueChangedFcn   = @(~, ~) drawProbe();
tbl_probe.CellEditCallback = @(~, ev) probeEdited(ev);
ax_cor.ButtonDownFcn       = @(~, ~) probeClick();   % on the axes: zoom/rotate modes clear fig.WindowButtonDownFcn
sw_display.ValueChangedFcn = @(~, ~) updateGroups();
sw_fraction.ValueChangedFcn = @(~, ~) updateGroups();
lst_birds.ValueChangedFcn  = @(~, ~) updateGroups();
btn_all.ButtonPushedFcn    = @(~, ~) selectBirds(1:numel(birds));
btn_ex.ButtonPushedFcn     = @(~, ~) selectBirds(find(ismember(birds, example_birds))');
btn_none.ButtonPushedFcn   = @(~, ~) selectBirds([]);
btn_reg_all.ButtonPushedFcn  = @(~, ~) selectRegions(true(numel(regions), 1));
btn_reg_pal.ButtonPushedFcn  = @(~, ~) selectRegions(ismember(regions, pallial)');
btn_reg_none.ButtonPushedFcn = @(~, ~) selectRegions(false(numel(regions), 1));
sw_color.ValueChangedFcn   = @(~, ~) updateGroups();
tbl.CellEditCallback       = @(~, ~) updateGroups();
tbl.CellSelectionCallback  = @(~, ev) tableClicked(ev);
sld_sec.ValueChangingFcn   = @(~, ev) sectionMoved(ev.Value);
sld_sec.ValueChangedFcn    = @(~, ev) sectionMoved(ev.Value);
spn_sec.ValueChangingFcn   = @(~, ev) sectionMoved(ev.Value);
spn_sec.ValueChangedFcn    = @(~, ev) sectionMoved(ev.Value);
chk_sec_outline.ValueChangedFcn = @(~, ~) drawSection();
bg_channel.SelectionChangedFcn = @(~, ~) drawSection();
sld_clim.ValueChangingFcn  = @(~, ev) setDisplayRange(ev.Value);
sld_clim.ValueChangedFcn   = @(~, ev) setDisplayRange(ev.Value);
tabs.SelectionChangedFcn   = @(~, ~) tabChanged();

% White background: the figure, its layout containers and the axes around the
% plots; the control columns a light grey
fig.Color = 'w';
is_container = {'Type', 'uigridlayout', '-or', 'Type', 'uitab', '-or', 'Type', 'uibuttongroup', '-or', 'Type', 'axes'};
set(findall(fig, is_container{:}), 'BackgroundColor', 'w');
set([findall(ctrl, is_container{:}); findall(ctrl2, is_container{:})], 'BackgroundColor', 0.97*[1 1 1]);

% Settings remembered between sessions: saved when the window goes, restored
% here (anything that no longer fits, e.g. after the data change, is skipped)
settings_file = fullfile(prefdir, 'brainViewerSettings.mat');
remembered = struct('ap', sld_ap, 'ml', sld_ml, 'dv', sld_dv, 'template', dd_tmpl, ...
    'angle', spn_angle, 'brain', chk_brain, 'hippo', chk_hippo, 'dl', chk_dl, ...
    'planes', chk_planes, 'grid', chk_grid, 'minor', chk_minor, 'probe', chk_probe, ...
    'roll', spn_roll, 'display', sw_display, 'fraction', sw_fraction, 'color', sw_color, ...
    'sites', sw_sites, 'birds', lst_birds, 'sec_outline', chk_sec_outline, 'clim', sld_clim);
saved = loadSettings();
fig.DeleteFcn = @(~, ~) saveSettings();   % on closing, and on delete()

sld_ap.Value = snapValue(sld_ap, sld_ap.Value, 1);   % AP on a section, shared with tab 2
sld_sec.Value = nearestSection(sld_ap.Value);
syncSpinners();
templateChanged();     % sets outline checkboxes for the template, then updateGroups
if isfield(saved, 'hippo')   % templateChanged turns hippocampus back on
    chk_hippo.Value = saved.hippo && chk_hippo.Enable;
    toggleOutlines();
end
fig.AutoResizeChildren = 'off';   % (the grid layouts resize anyway) so SizeChangedFcn runs
fig.SizeChangedFcn = @(~, ~) scaleSlices();
% 3D drag-rotate is inert after the setup above until re-enabled. (No click
% callbacks here or on the figure: they turn it off again.)
enableDefaultInteractivity(ax_3d);
% Right-click view presets for the 3D view, oriented like the slice views
views_3d = struct('name', {'Top (horizontal)', 'Side (sagittal)', 'Coronal', 'Default'}, ...
    'azel', {[0 90], [0 0], [-90 0], view_3d_default});
ax_3d.ContextMenu = uicontextmenu(fig);
for iv = 1:numel(views_3d)
    uimenu(ax_3d.ContextMenu, 'Text', views_3d(iv).name, 'MenuSelectedFcn', @(~, ~) view(ax_3d, views_3d(iv).azel));
end


    function rot = templateRotation()
        % deg nose down from the current template's native frame to the display frame
        rot = spn_angle.Value - tmpl(it_cur).native_angle;
    end

    function s = loadSettings()
        % Restore the controls, table selection, channel, probe target and 3D
        % view saved by saveSettings; returns what was saved (empty struct if none)
        s = struct();
        if ~isfile(settings_file), return; end
        try
            s = load(settings_file).s;
        catch
            return
        end
        for name = string(fieldnames(remembered))'
            if isfield(s, name)
                try remembered.(name).Value = s.(name); catch, end   % out of range now: keep the default
            end
        end
        if isfield(s, 'sel_regions') && size(s.sel, 2) == height(sites)   % matched by region name
            [has, row] = ismember(regions, s.sel_regions);
            sel = false(numel(regions), height(sites));
            sel(has, :) = s.sel(row(has), :);
            tbl.Data{:, 2:end} = sel;
        end
        if isfield(s, 'channel') && isscalar(s.channel) && s.channel <= numel(rb_channel)
            bg_channel.SelectedObject = rb_channel(s.channel);
        end
        if isfield(s, 'probe_target') && numel(s.probe_target) == 4
            probe_ml = s.probe_target(1); probe_ap = s.probe_target(2); probe_dv = s.probe_target(3);
            probe_angle = s.probe_target(4);
        end
        if isfield(s, 'view'), view(ax_3d, s.view); end
    end

    function saveSettings()
        % Save the settings for next time (when the window is deleted, however)
        s = struct();
        for name = string(fieldnames(remembered))'
            s.(name) = remembered.(name).Value;
        end
        s.sel = tbl.Data{:, 2:end}; s.sel_regions = regions;
        s.channel = find(rb_channel == bg_channel.SelectedObject);
        s.probe_target = [probe_ml probe_ap probe_dv probe_angle];
        [az, el] = view(ax_3d); s.view = [az el];
        try save(settings_file, 's'); catch, end   % never block closing
    end

    function updateCellDisplay()
        % Cells from the MC022 frame to the display frame
        [cell_ap, cell_dv] = rotateNoseDown(ap_frame, dv_frame, spn_angle.Value - frame_angle);
        cell_ml = ml_frame;
        [site_ap_d, site_dv_d] = rotateNoseDown(site_ap, site_dv, spn_angle.Value - frame_angle);
    end

    function show = sitesShown()
        % Injection sites to mark: the table's site columns with any region checked
        show = strcmp(sw_sites.Value, 'On') & any(tbl.Data{:, 2:end}, 1)';
    end

    function drawSites3d()
        uistack([h_sites3d; h_probe(4)], 'top');   % over the cells
        show = sitesShown();
        for kk = 1:height(sites)
            if show(kk)
                set(h_sites3d(kk), 'XData', site_ap_d(kk), 'YData', sites.ML(kk), 'ZData', site_dv_d(kk));
            else
                set(h_sites3d(kk), 'XData', NaN, 'YData', NaN, 'ZData', NaN);
            end
        end
    end

    function colorSites()
        % Injection site stars black (white edged), or in their site colors
        % (black edged) when cells are colored by injection site (tab 2: MC022's
        % channel colors)
        by_site = strcmp(sw_color.Value, 'Injection site');
        edge = [1 1 1]; if by_site, edge = [0 0 0]; end
        for kk = 1:height(sites)
            c = [0 0 0]; c2 = [0 0 0];
            if by_site, c = site_colors(kk, :); c2 = sec_site_colors(kk, :); end
            set([h_sites2d(:, kk); h_sites3d(kk)], 'MarkerFaceColor', c, 'MarkerEdgeColor', edge);
            set(h_sec_sites(kk), 'MarkerFaceColor', c2, 'MarkerEdgeColor', edge);
        end
    end

    function templateChanged()
        it_cur = dd_tmpl.Value;
        has = isfield(tmpl(it_cur).D, ["Ihippo" "Idl"]);
        % Outlines the template lacks are unchecked and greyed out; hippocampus
        % comes back on when a template has it, DL/CDL stays as it was
        chk_hippo.Value = has(1); chk_dl.Value = chk_dl.Value && has(2);
        chk_hippo.Enable = has(1); chk_dl.Enable = has(2);
        drawTemplate3d();
        updateGroups();
        drawProbe();             % the surface and entrance depend on the template
    end

    function angleChanged()
        setTemplateTransform();
        updateGroups();          % cells and densities move with the angle
        drawProbe();
    end

    function in = insideBrain(ml, ap, dv, angle)
        % Points inside the current template's brain outline, in the display
        % frame at skull groove angle (default: the current one)
        if nargin < 4, angle = spn_angle.Value; end
        [Xn, Zn] = rotateNoseDown(ap, dv, tmpl(it_cur).native_angle - angle);
        v = tmpl(it_cur).F.Ibrain(ml, Xn, Zn);
        in = ~isnan(v) & v > 0.5;
    end

    function surface = dorsalSurface(ap, ml, angle)
        % DV of the brain's top at ap, ml (frame at skull groove angle; empty if
        % none), 0.01 mm steps
        z = (Dd.zz(end) : -0.01 : Dd.zz(1))';
        surface = z(find(insideBrain(ml*ones(size(z)), ap*ones(size(z)), z, angle), 1));
    end

    function probeClick()
        % Probe track on: a click inside the brain on the coronal slice plans a
        % new track at the current skull groove, target rounded to 0.1 mm
        if ~chk_probe.Value || tabs.SelectedTab ~= tab1, return; end
        p = round([sld_ap.Value, ax_cor.CurrentPoint(1, 1:2)] * 10) / 10;   % AP, ML, DV
        if ~insideBrain(p(2), p(1), p(3)), return; end
        probe_angle = spn_angle.Value;
        setProbeTarget(p(1), p(2), p(3));
    end

    function probeEdited(ev)
        % Edited probe table (coordinates at the planning skull groove). The
        % target (AP, ML, DV below the dorsal surface) moves there; an edit to
        % the track (roll, entrance AP or ML, depth) keeps the other three and
        % moves the target to the end of the track. Then redraw; anything that
        % puts the target outside the brain reverts.
        v = tbl_probe.Data;
        if isnan(ev.NewData), drawProbe(); return; end
        if ev.Indices(2) >= 4
            roll = min(max(v(4), spn_roll.Limits(1)), spn_roll.Limits(2));
            ap = v(5);
            e_dv = dorsalSurface(ap, v(6), probe_angle);   % entrance on the surface
            if isempty(e_dv), drawProbe(); return; end
            side = sign(v(6)) + (v(6) == 0);               % toward the midline = toward ML 0
            ml = v(6) - side*sind(roll)*abs(v(7));
            dv = e_dv - cosd(roll)*abs(v(7));
            if ~insideBrain(ml, ap, dv, probe_angle), drawProbe(); return; end
            spn_roll.Value = roll;
            setProbeTarget(ap, ml, dv);
            return
        end
        ap = v(1); ml = v(2);
        surface = dorsalSurface(ap, ml, probe_angle);
        dv = surface - abs(v(3));       % DV below the surface, either sign typed
        if isempty(surface) || ~insideBrain(ml, ap, dv, probe_angle), drawProbe(); return; end   % outside the brain: revert
        setProbeTarget(ap, ml, dv);
    end

    function [e_ml, e_dv, depth] = probeEntrance(ap, ml, dv, roll)
        % Where a track to target ap, ml, dv (frame at probe_angle) at roll (deg,
        % + toward the midline) enters the brain: the last point in the brain
        % going up the track; NaN if the target is outside it
        t = (0 : 0.01 : 15)';                              % along the track, upward
        side = sign(ml) + (ml == 0);                       % toward the midline = toward ML 0
        t_ml = ml + side*sind(roll)*t; t_dv = dv + cosd(roll)*t;
        k = find(insideBrain(t_ml, ap*ones(size(t)), t_dv, probe_angle), 1, 'last');
        if isempty(k), e_ml = NaN; e_dv = NaN; depth = NaN; return; end
        e_ml = t_ml(k); e_dv = t_dv(k); depth = t(k);
    end

    function setProbeTarget(ap, ml, dv)
        % Target at ap, ml, dv in the frame at the planning skull groove
        % (probe_angle); slices move to where it is shown now (snapped)
        probe_ml = ml;
        [probe_ap, probe_dv] = rotateNoseDown(ap, dv, frame_angle - probe_angle);
        [ap_d, dv_d] = rotateNoseDown(ap, dv, spn_angle.Value - probe_angle);
        sliderMoved(sld_ap, ap_d, 1);
        sliderMoved(sld_ml, ml, 2);
        sliderMoved(sld_dv, dv_d, 3);
        drawProbe();
    end

    function drawProbe()
        % Track from the brain surface to the target at the roll angle (in the
        % coronal plane, + toward the midline), planned at skull groove
        % probe_angle: the table is in that frame. The track is fixed in the
        % brain, so at another skull groove it is drawn tilted, in all four views.
        set(h_probe, 'MarkerIndices', 1, 'XData', NaN, 'YData', NaN, 'ZData', NaN);   % marker: the target, once drawn
        tbl_probe.Data = zeros(0, 7);
        if ~chk_probe.Value || isempty(probe_ml), return; end
        [ap, dv] = rotateNoseDown(probe_ap, probe_dv, probe_angle - frame_angle);
        ml = probe_ml;
        roll = spn_roll.Value;
        surface = dorsalSurface(ap, ml, probe_angle);      % straight above the target
        [e_ml, e_dv, depth] = probeEntrance(ap, ml, dv, roll);
        if isempty(surface) || isnan(depth), return; end   % target outside this template's brain
        tbl_probe.Data = [ap, ml, dv - surface, roll, ap, e_ml, depth];
        P = [ap e_ml e_dv; ap ml dv];                      % entrance, target: AP, ML, DV
        [P(:, 1), P(:, 3)] = rotateNoseDown(P(:, 1), P(:, 3), spn_angle.Value - probe_angle);   % to the display frame
        axes_of = {[2 3], [1 3], [1 2]};                   % slice x, y (as in drawSlices)
        for dim = 1:3
            set(h_probe(dim), 'XData', P(:, axes_of{dim}(1)), 'YData', P(:, axes_of{dim}(2)), 'ZData', []);
        end
        set(h_probe(4), 'XData', P(:, 1), 'YData', P(:, 2), 'ZData', P(:, 3));
        set(h_probe, 'MarkerIndices', 2);
        arrayfun(@(h) uistack(h, 'top'), h_probe);
    end

    function setTemplateTransform()
        rot = templateRotation();
        c = cosd(rot); s = sind(rot);
        h_tf.Matrix = [c 0 s 0; 0 1 0 0; -s 0 c 0; 0 0 0 1];   % same rotation as rotateNoseDown
    end

    function lims = brainLimits3d()
        % 3D axes limits [x; y; z] around all templates' brains (MC029's is the
        % largest) at the display angle, out to the next whole mm
        lo = inf(1, 3); hi = -inf(1, 3);
        for kk = 1:numel(tmpl)
            v = tmpl(kk).surf.Ibrain.vertices;
            [v(:, 1), v(:, 3)] = rotateNoseDown(v(:, 1), v(:, 3), spn_angle.Value - tmpl(kk).native_angle);
            lo = min(lo, min(v)); hi = max(hi, max(v));
        end
        lims = [floor(lo); ceil(hi)]';
    end

    function drawTemplate3d()
        % Template surfaces (native frame) under h_tf; plotBrainRegion3d resets the view
        [az, el] = view(ax_3d);
        delete(h_tf.Children);
        T = tmpl(it_cur);
        h_brain = plotBrainRegion3d(T.D.xx, T.D.yy, T.D.zz, T.surf.Ibrain, 0.5, .7*[1 1 1], 0.15, ax_3d);
        h_brain.SpecularStrength = 0.1; h_brain.Parent = h_tf;
        h_hippo = gobjects(1); h_dl = gobjects(1);
        if isfield(T.surf, 'Ihippo')
            h_hippo = plotBrainRegion3d(T.D.xx, T.D.yy, T.D.zz, T.surf.Ihippo, 0.5, color_hippo, 0.4, ax_3d);
            h_hippo.Parent = h_tf;
        end
        if isfield(T.surf, 'Idl')
            h_dl = plotBrainRegion3d(T.D.xx, T.D.yy, T.D.zz, T.surf.Idl, 0.5, color_dl, 0.4, ax_3d);
            h_dl.Parent = h_tf;
        end
        setTemplateTransform();
        toggleOutlines3d();
        restoreAxes3d(az, el);
    end

    function restoreAxes3d(az, el)
        % Undo what plotBrainRegion3d changes: limits around the template brains
        % (fixed, so planes and cells don't rescale them; equal aspect), the view,
        % grid, and one headlight
        lims = brainLimits3d();
        daspect(ax_3d, [1 1 1]);              % equal aspect; axis equal would override the limits
        set(ax_3d, 'XLim', lims(1, :), 'YLim', lims(2, :), 'ZLim', lims(3, :));
        view(ax_3d, az, el);
        axis(ax_3d, 'vis3d');                 % hold the aspect while rotating
        % Zoom: seen from straight above, the box's AP x ML face fills the panel
        % like the slice views, with room for the tick labels (label_margin);
        % kept at every rotation (corners may be cut off)
        ax_3d.CameraTargetMode = 'auto'; ax_3d.CameraPositionMode = 'auto';
        span = max(diff(lims(1:2, :), 1, 2));
        d = norm(ax_3d.CameraPosition - ax_3d.CameraTarget);
        ax_3d.CameraViewAngle = 2*atand(span/2 * label_margin / d);
        grid(ax_3d, 'off')
        applyGrid(ax_3d);                     % ticks for the new limits
        slice_lims = lims;
        scaleSlices();
        delete(findall(ax_3d, 'Type', 'light')); camlight(ax_3d, 'headlight')
    end

    function scaleSlices()
        % Slice views: at least slice_lims (AP, ML, DV rows, from restoreAxes3d),
        % widened so all three show the same mm per pixel (corresponding axes
        % equally long on screen); again when the window is resized
        if isempty(slice_lims), return; end
        dims_of = {[2 3], [1 3], [1 2]};      % slice x, y (as in sliceSpec)
        mid = mean(slice_lims, 2); range = diff(slice_lims, 1, 2);
        for pass = 1:2                        % new limits change the tick labels, so the plot area: settle once
            drawnow
            px = zeros(3, 2);
            for dim = 1:3
                p = slice_ax(dim).InnerPosition; px(dim, :) = p(3:4);   % pixels for the plot
            end
            mm_per_px = max(range([dims_of{:}])' ./ reshape(px', 1, []));   % the tightest fit sets the scale
            for dim = 1:3
                d = dims_of{dim};
                half = px(dim, :) * mm_per_px / 2;
                set(slice_ax(dim), 'XLim', mid(d(1)) + [-1 1]*half(1), 'YLim', mid(d(2)) + [-1 1]*half(2));
                applyGrid(slice_ax(dim));
            end
        end
    end

    function tableClicked(ev)
        % Region name toggles its row: all on, or all off if it was full
        if isempty(ev.Indices) || ev.Indices(1, 2) ~= 1, return; end
        row = ev.Indices(1, 1);
        tbl.Data{row, 2:end} = ~all(tbl.Data{row, 2:end});
        tbl.Selection = [];   % so clicking the same name again fires again
        updateGroups();
    end

    function styleTable()
        % Region names in their plot colors (color by region) or site columns
        % tinted with theirs (color by injection site); names of fully checked
        % rows shaded
        removeStyle(tbl);
        full_rows = find(all(tbl.Data{:, 2:end}, 2));
        if ~isempty(full_rows)
            addStyle(tbl, style_row_on, 'cell', [full_rows ones(size(full_rows))]);
        end
        if strcmp(sw_color.Value, 'Region')
            for kk = 1:numel(regions)
                addStyle(tbl, uistyle('FontColor', region_colors(kk, :), 'FontWeight', 'bold'), 'cell', [kk 1]);
            end
        else
            for kk = 1:height(sites)
                addStyle(tbl, uistyle('BackgroundColor', 1 - 0.4*(1 - site_colors(kk, :))), 'column', kk + 1);
            end
        end
    end

    function selectRegions(rows)
        % Check the region rows (logical, one per region) in the site columns in
        % use (all sites if none are); clear everything else
        sel = tbl.Data{:, 2:end};
        cols = any(sel, 1);
        if ~any(cols), cols(:) = true; end
        tbl.Data{:, 2:end} = rows & cols;
        updateGroups();
    end

    function selectBirds(bird_idx)
        lst_birds.Value = bird_idx;
        updateGroups();
    end

    function updateBirdCounts(n)
        % List each bird's cell count for the table selection; grey out birds with none
        lst_birds.Items = cellstr(compose('%s (n=%d)', birds, n));
        removeStyle(lst_birds);
        if any(n == 0)
            addStyle(lst_birds, uistyle('FontColor', .6*[1 1 1]), 'item', find(n == 0));
        end
    end

    function updateGroups()
        % Collect the selected cells into color groups, then redraw everything
        fig_ptr = fig.Pointer; fig.Pointer = 'watch'; drawnow
        updateCellDisplay();
        colorSites();
        sel = tbl.Data{:, 2:end}';   % sites x regions
        styleTable();
        use_blobs = strcmp(sw_display.Value, 'Density');
        by_region = strcmp(sw_color.Value, 'Region');
        in_sel = false(size(cell_ap));
        [si, ri] = find(sel);
        for kk = 1:numel(si)
            in_sel = in_sel | (cell_site == si(kk) & cell_region == ri(kk));
        end
        in_sel = in_sel & ~isnan(cell_dv);
        updateBirdCounts(accumarray(cell_bird(in_sel), 1, [numel(birds) 1]));
        if by_region
            keys = unique(ri)'; key_of_cell = cell_region;
        else
            keys = unique(si)'; key_of_cell = cell_site;
        end

        if tabs.SelectedTab == tab2
            % Tab 2: MC022's own cells only, in the MC022 frame
            in_022 = in_sel & cell_bird == find(birds == "MC022");
            groups2 = makeGroups(in_022, keys, key_of_cell, by_region, use_blobs, ap_frame, ml_frame, dv_frame, sec_site_colors);
            drawSection();
        else
            in_sel = in_sel & ismember(cell_bird, lst_birds.Value);
            groups = makeGroups(in_sel, keys, key_of_cell, by_region, use_blobs, cell_ap, cell_ml, cell_dv, site_colors);
            draw3d(use_blobs);
            drawSites3d();
            makeSliceGroups(use_blobs);
            drawSlices();
        end
        fig.Pointer = fig_ptr;
    end

    function grp = makeGroups(in_sel, keys, key_of_cell, by_region, use_blobs, ap, ml, dv, site_cols)
        % One color group per region or injection site (colors site_cols) among
        % the cells in_sel; density volumes (use_blobs) from coordinates ap, ml,
        % dv on the display grid
        if by_region
            names = regions; colors = region_colors;
        else
            names = string(site_labels); colors = site_cols;
        end
        grp = struct('idx', {}, 'show', {}, 'color', {}, 'name', {}, 'label', {}, 'I', {});
        for key = keys
            idx = find(in_sel & key_of_cell == key);
            if isempty(idx), continue; end
            g.idx = idx;
            g.show = idx(draw_key(idx) < str2double(erase(sw_fraction.Value, '%'))/100);   % points drawn; blobs use all cells
            g.color = colors(key, :);
            g.name = names(key);
            g.label = sprintf('%s (n=%d)', g.name, numel(idx));
            g.I = [];
            if use_blobs
                g.I = cellDensity3d(ap(idx), ml(idx), dv(idx), Dd.xx, Dd.yy, Dd.zz);
            end
            grp(end+1) = g; %#ok<AGROW>
        end
    end

    function h_keys = groupLegend(ax, grp, h, labels, use_blobs)
        % Legend of the groups; for cells, large stand-in markers (returned) so
        % the colors are readable
        h_keys = gobjects(0);
        if isempty(grp), legend(ax, 'off'); return; end
        if ~use_blobs
            h_keys = gobjects(numel(grp), 1);
            for kk = 1:numel(grp)
                h_keys(kk) = line(ax, 'XData', NaN, 'YData', NaN, 'ZData', NaN, 'LineStyle', 'none', 'Marker', 'o', ...
                    'Color', grp(kk).color, 'MarkerFaceColor', grp(kk).color, 'MarkerSize', 10);
            end
            h = h_keys;
        end
        legend(ax, h, labels, 'Interpreter', 'none', 'Location', 'northeast', ...
            'AutoUpdate', 'off', 'FontSize', 12);   % objects drawn later (template, lights) stay out of it
    end

    function makeSliceGroups(use_blobs)
        % One empty line (cells) or contour (blobs) per group in each slice view
        for dim = 1:3
            delete(h_groups2d{dim});
            s = sliceSpec(dim);
            h = gobjects(numel(groups), 1);
            for kk = 1:numel(groups)
                if use_blobs
                    h(kk) = newContour(slice_ax(dim), s.hx, s.hy, groups(kk).color, blob_level, 2);
                else
                    h(kk) = plot(slice_ax(dim), NaN, NaN, '.', 'Color', groups(kk).color, 'MarkerSize', 6, 'HitTest', 'off');
                end
            end
            h_groups2d{dim} = h;
            uistack([h_sites2d(dim, :) h_probe(dim)], 'top');   % injection sites and probe over the cells
        end
    end

    function draw3d(use_blobs)
        [az, el] = view(ax_3d);             % plotBrainRegion3d resets the view
        delete(h_groups3d); delete(h_key3d);
        h_groups3d = gobjects(numel(groups), 1);
        for kk = 1:numel(groups)
            g = groups(kk);
            if use_blobs
                h_groups3d(kk) = plotBrainRegion3d(Dd.xx, Dd.yy, Dd.zz, g.I, blob_level, g.color, 1, ax_3d);
            else
                % line, not plot3: a group can have no cells drawn at 10% (plot3 then returns nothing)
                h_groups3d(kk) = line(ax_3d, 'XData', cell_ap(g.show), 'YData', cell_ml(g.show), 'ZData', cell_dv(g.show), ...
                    'LineStyle', 'none', 'Marker', '.', 'Color', g.color, 'MarkerSize', 4);
            end
        end
        restoreAxes3d(az, el);
        applyGrid(ax_3d);
        h_key3d = groupLegend(ax_3d, groups, h_groups3d, {groups.label}, use_blobs);
    end

    function toggleOutlines3d()
        h = [h_brain h_hippo h_dl];
        on = [chk_brain.Value, chk_hippo.Value && chk_hippo.Enable, chk_dl.Value && chk_dl.Enable];
        for kk = 1:3
            if isgraphics(h(kk)), h(kk).Visible = on(kk); end
        end
    end

    function toggleOutlines()
        toggleOutlines3d();
        drawSlices();
    end

    function value = snapValue(sld, value, dim)
        % Nearest slider position on dim's section grid, within the slider limits
        o = origin_dim(dim); step = step_dim(dim);
        value = min(max(o + round((value - o)/step)*step, sld.Limits(1)), sld.Limits(2));
    end

    function sliderMoved(sld, value, dim)
        sld.Value = snapValue(sld, value, dim);
        if dim == 1, sld_sec.Value = nearestSection(sld.Value); end   % tab 2 follows AP
        syncSpinners();
        drawSlices(dim);
    end

    function syncSpinners()
        % Spinners show their slider's (snapped) value
        spinners = [spn_ap spn_ml spn_dv spn_sec];
        sliders = [sld_ap sld_ml sld_dv sld_sec];
        for kk = 1:numel(spinners)
            if spinners(kk).Value ~= sliders(kk).Value, spinners(kk).Value = sliders(kk).Value; end
        end
    end

    function measureBox(ax)
        delete(h_box);
        % Leave pan/zoom/data tips first: while one is on, it takes the clicks
        pan(ax, 'off'); zoom(ax, 'off'); datacursormode(fig, 'off');
        h_box = drawrectangle(ax, 'Color', [0.85 0 0.6], 'FaceAlpha', 0, 'LabelVisible', 'on', ...
            'LabelTextColor', 'w');   % white on the magenta label: readable on black images and white slices
        labelBox(h_box);
        addlistener(h_box, 'MovingROI', @(roi, ~) labelBox(roi));
    end

    function labelBox(roi)
        if isempty(roi.Position), return; end
        roi.Label = sprintf('\x0394X %.2f, \x0394Y %.2f mm', roi.Position(3), roi.Position(4));
    end

    function si = nearestSection(ap)
        [~, si] = min(abs(sec_ap - ap));
    end

    function tabChanged()
        % Move the shared controls to the tab now showing and redraw it
        it = 1 + (tabs.SelectedTab == tab2);
        parents = [ctrl ctrl2];
        shared.Parent = parents(it);
        shared.Layout.Row = shared_row(it); shared.Layout.Column = [1 2];
        % Scroll wheel steps sections in tab 2 only: a WindowScrollWheelFcn turns
        % off the axes' built-in drag/zoom (tab 1's 3D rotation)
        fig.WindowScrollWheelFcn = '';
        if it == 2, fig.WindowScrollWheelFcn = @(~, ev) sectionScroll(ev); end
        updateGroups();
    end

    function s = sliceSpec(dim)
        % dim 1 coronal (ML x DV at AP), 2 sagittal (AP x DV at ML), 3 horizontal (AP x ML at DV),
        % all in the display frame. s.plane: display coordinates of the slice grid.
        switch dim
            case 1
                s.pos = sld_ap.Value; [~, k] = min(abs(Dd.xx - s.pos));
                s.hx = Dd.yy; s.hy = Dd.zz; s.index = @(V) squeeze(V(:, k, :))';
                [Y, Z] = meshgrid(Dd.yy, Dd.zz); s.plane = {s.pos*ones(size(Y)), Y, Z};
                s.cx = cell_ml; s.cy = cell_dv; s.c_depth = cell_ap;
                s.cross = [sld_ml.Value sld_dv.Value];
                s.xl = 'Lateral (mm)'; s.yl = 'Dorsal (mm)';
            case 2
                s.pos = sld_ml.Value; [~, k] = min(abs(Dd.yy - s.pos));
                s.hx = Dd.xx; s.hy = Dd.zz; s.index = @(V) squeeze(V(k, :, :))';
                [X, Z] = meshgrid(Dd.xx, Dd.zz); s.plane = {X, s.pos*ones(size(X)), Z};
                s.cx = cell_ap; s.cy = cell_dv; s.c_depth = cell_ml;
                s.cross = [sld_ap.Value sld_dv.Value];
                s.xl = 'Anterior (mm)'; s.yl = 'Dorsal (mm)';
            case 3
                s.pos = sld_dv.Value; [~, k] = min(abs(Dd.zz - s.pos));
                s.hx = Dd.xx; s.hy = Dd.yy; s.index = @(V) V(:, :, k);
                [X, Y] = meshgrid(Dd.xx, Dd.yy); s.plane = {X, Y, s.pos*ones(size(X))};
                s.cx = cell_ap; s.cy = cell_ml; s.c_depth = cell_dv;
                s.cross = [sld_ap.Value sld_ml.Value];
                s.xl = 'Anterior (mm)'; s.yl = 'Lateral (mm)';
        end
    end

    function S = sampleTemplate(name, s)
        % Template volume on the slice's display-frame plane: rotate the plane back
        % to the template's native frame and interpolate
        [Xn, Zn] = rotateNoseDown(s.plane{1}, s.plane{3}, -templateRotation());
        S = tmpl(it_cur).F.(name)(s.plane{2}, Xn, Zn);
        S(isnan(S)) = 0;
    end

    function setupSlices()
        % Axes settings, outline contours and crosshairs, made once
        for dim = 1:3
            ax = slice_ax(dim);
            s = sliceSpec(dim);
            hold(ax, 'on')
            for kk = 1:3
                h_outline(dim, kk) = newContour(ax, s.hx, s.hy, outline_colors(kk, :), 0.5, 1);
            end
            h_cross(dim, 1) = xline(ax, s.cross(1), ':', 'Color', .3*[1 1 1], 'HitTest', 'off');
            h_cross(dim, 2) = yline(ax, s.cross(2), ':', 'Color', .3*[1 1 1], 'HitTest', 'off');
            h_sites2d(dim, :) = siteMarkers(ax, height(sites));
            daspect(ax, [1 1 1]);
            if dim == 3, ax.YDir = 'reverse'; end     % horizontal: seen from above, anterior right, so right side down
            xlabel(ax, s.xl); ylabel(ax, s.yl);
            applyGrid(ax);
        end
    end

    function drawSlices(dims)
        % Update the slice data for dims (default all); crosshairs always
        if nargin < 1, dims = 1:3; end
        x = Dd.xx([1 end end 1]); y = Dd.yy([1 1 end end]); z = Dd.zz([1 1 end end]);
        set(h_planes(1), 'XData', sld_ap.Value*[1 1 1 1], 'YData', Dd.yy([1 end end 1]), 'ZData', z);   % y must cycle opposite to z
        set(h_planes(2), 'XData', x, 'YData', sld_ml.Value*[1 1 1 1], 'ZData', z);
        set(h_planes(3), 'XData', x, 'YData', y, 'ZData', sld_dv.Value*[1 1 1 1]);
        show_outline = [chk_brain.Value, chk_hippo.Value && chk_hippo.Enable, chk_dl.Value && chk_dl.Enable];
        for dim = 1:3
            s = sliceSpec(dim);
            % Depth window, one step (step_dim) wide. Coronal: centered on pos,
            % which snaps to a section. Sagittal: centered on |pos| (window
            % mirrored when ML < 0). Horizontal: symmetric about pos.
            depth_sign = 1;
            step = step_dim(dim);
            if dim == 1
                lo = s.pos - step/2;
            elseif dim == 2
                depth_sign = sign(s.pos) + (s.pos == 0);
                lo = abs(s.pos) - step/2;
            else
                lo = s.pos - step/2 - 1e-9;
            end
            hi = lo + step + 2e-9*(dim == 3);
            h_cross(dim, 1).Value = s.cross(1);
            h_cross(dim, 2).Value = s.cross(2);
            if ~ismember(dim, dims), continue; end
            for kk = 1:3
                if show_outline(kk)
                    setContour(h_outline(dim, kk), sampleTemplate(outline_names(kk), s));
                else
                    h_outline(dim, kk).Visible = 'off';
                end
            end
            for kk = 1:numel(groups)
                g = groups(kk);
                if ~isempty(g.I)
                    setContour(h_groups2d{dim}(kk), s.index(g.I));
                else
                    c = depth_sign*s.c_depth(g.show);
                    in = g.show(c >= lo & c < hi);
                    set(h_groups2d{dim}(kk), 'XData', s.cx(in), 'YData', s.cy(in));
                end
            end
            % Injection sites within site_window of the plane
            P = [site_ap_d, sites.ML, site_dv_d];
            axes_of = {[2 3 1], [1 3 2], [1 2 3]};   % slice x, slice y, depth
            P = P(:, axes_of{dim});
            near = sitesShown() & abs(depth_sign*P(:, 3) - (lo + hi)/2) <= site_window;
            for kk = 1:height(sites)
                if near(kk), set(h_sites2d(dim, kk), 'XData', P(kk, 1), 'YData', P(kk, 2));
                else, set(h_sites2d(dim, kk), 'XData', NaN, 'YData', NaN); end
            end
            slice_ax(dim).Title.String = sprintf(title_fmt{dim}, s.pos);
        end
        drawnow limitrate
    end

    function sectionMoved(value)
        % Tab 2 slider/spinner: snap to a whole section and redraw; tab 1's AP follows
        sld_sec.Value = min(max(round(value), 1), numel(sec_orders));
        sld_ap.Value = snapValue(sld_ap, sec_ap(sld_sec.Value), 1);
        syncSpinners();
        drawSection();
    end

    function setDisplayRange(lims)
        % Tab 2 min/max slider (log scale, toLog): image color limits, marked
        % on the histogram
        if diff(lims) <= 0, return; end
        ax_sec.CLim = 10.^lims - 1;
        h_clim(1).Value = lims(1); h_clim(2).Value = lims(2);
    end

    function sectionScroll(ev)
        % Scroll wheel over tab 2's image steps sections (elsewhere, e.g. the
        % table, it scrolls as usual)
        p = getpixelposition(ax_sec, true);
        xy = fig.CurrentPoint;
        if all(xy >= p(1:2) & xy <= p(1:2) + p(3:4))
            sectionMoved(sld_sec.Value + sign(ev.VerticalScrollCount));
        end
    end

    function drawSection()
        % Tab 2: the section's combined image in the selected channel (rendered
        % once, then kept), its outlines and cells
        si = round(sld_sec.Value);
        so = sec_orders(si);
        ich = find(rb_channel == bg_channel.SelectedObject);
        if isempty(sec_cache{si, ich})
            sec_cache{si, ich} = renderSection(secs, so, sec_ML, sec_DV, sec_channels(ich));
        end
        img = sec_cache{si, ich};
        h_sec_img.CData = img;
        % Histogram of this section's intensities on the fixed log scale (values
        % at or below 0 in the first bin, above clim_top in the last); the
        % display range stays as set
        vals = double(img(~isnan(img)));   % NaN: off the images
        vals = vals(vals > -0.5);          % far below the background: padding inside stitched images
        h_edges = linspace(0, toLog(clim_top), 141);
        n_hist = histcounts(min(toLog(vals), h_edges(end)), h_edges);
        set(h_hist, 'XData', h_edges, 'YData', max([n_hist n_hist(end)], 0.5));   % 0.5: empty bins on the log axis
        setDisplayRange(sld_clim.Value);
        ka = find(secs.trueOrder == so)';
        ks = ka(~isnan(secs.shift(ka)));
        if isempty(ks)
            set(h_sec_out, 'XData', NaN, 'YData', NaN);
        else
            setPolys(h_sec_out, secs.outline_smooth(ks(1)), chk_sec_outline.Value);   % A/B halves share it
        end
        % MC022's cells of each group on this section, or its density cut at the
        % section's AP (as the tab 1 slices)
        use_blobs = strcmp(sw_display.Value, 'Density');
        delete(h_sec_groups);
        h_sec_groups = gobjects(numel(groups2), 1);
        labels = cell(1, numel(groups2));
        [~, ix] = min(abs(Dd.xx - sec_ap(si)));
        for kk = 1:numel(groups2)
            g = groups2(kk);
            on_sec = g.idx(abs(ap_frame(g.idx) - sec_ap(si)) < step_dim(1)/2);
            labels{kk} = sprintf('%s (n=%d)', g.name, numel(on_sec));   % cells on this section
            if use_blobs
                h_sec_groups(kk) = newContour(ax_sec, Dd.yy, Dd.zz, g.color, blob_level, 2);
                setContour(h_sec_groups(kk), squeeze(g.I(:, ix, :))');
            else
                in = intersect(on_sec, g.show);
                h_sec_groups(kk) = line(ax_sec, 'XData', ml_frame(in), 'YData', dv_frame(in), ...
                    'LineStyle', 'none', 'Marker', '.', 'Color', g.color, 'MarkerSize', 6);
            end
        end
        delete(h_sec_keys);
        h_sec_keys = groupLegend(ax_sec, groups2, h_sec_groups, labels, use_blobs);
        near = sitesShown() & abs(site_ap - sec_ap(si)) <= site_window;   % section image: MC022 frame
        for kk = 1:height(sites)
            if near(kk), set(h_sec_sites(kk), 'XData', sites.ML(kk), 'YData', site_dv(kk));
            else, set(h_sec_sites(kk), 'XData', NaN, 'YData', NaN); end
        end
        uistack(h_sec_sites, 'top');   % over the cells
        [~, names] = fileparts(secs.image_file(ka));
        % secs.AP: placed AP (Dryad AP with the vibratome scale and ap_shift), same as tab 1's AP
        lbl_sec.Text = sprintf('AP %.2f mm, %s (trueOrder %d)', secs.AP(ka(1)), strjoin(names, ' + '), so);
        if isempty(ks), lbl_sec.Text = [lbl_sec.Text ' (not aligned)']; end
        drawnow limitrate
    end

    function applyGrid(ax)
        % Ticks every 1 mm, minor ticks every 0.2 mm; grid lines per the checkboxes
        rulers = {ax.XAxis, ax.YAxis, ax.ZAxis};
        for kk = 1:numel(rulers)
            lims = rulers{kk}.Limits;
            rulers{kk}.TickValues = ceil(lims(1)) : floor(lims(2));
            rulers{kk}.MinorTickValues = ceil(lims(1)/0.2)*0.2 : 0.2 : lims(2);
        end
        set(ax, 'XGrid', chk_grid.Value, 'YGrid', chk_grid.Value, 'ZGrid', chk_grid.Value, ...
            'XMinorGrid', chk_minor.Value, 'YMinorGrid', chk_minor.Value, 'ZMinorGrid', chk_minor.Value);
    end
end


function h = siteMarkers(ax, n)
% n hidden stars, one per injection site (colorSites sets their colors); data
% set when shown
h = gobjects(n, 1);
for k = 1:n
    h(k) = line(ax, 'XData', NaN, 'YData', NaN, 'ZData', NaN, 'LineStyle', 'none', 'Marker', 'p', ...
        'MarkerSize', 12, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'w', 'HitTest', 'off');
end
end


function [ap, dv] = sitePositions(sites, tf, F_brain, frame_angle, zz)
% Injection sites (Table 1: AP, ML from lambda; depth below the brain surface at
% head_angle) in the MC022 frame: AP mapped like the cells (dryadToTemplate),
% the surface found along the vertical of the head_angle frame, at +ML (the
% injected side). F_brain: MC022's brain interpolant {yy, xx, zz}.
ap = nan(height(sites), 1); dv = ap;
ap_frame = tf.ap_mid + (sites.AP - tf.ap_mid) * tf.ap_scale + tf.ap_shift;
for k = 1:height(sites)
    rot = sites.head_angle(k) - frame_angle;               % frame -> head_angle frame, deg nose down
    ap_h = ap_frame(k);                                    % Table 1 AP is in the head_angle frame
    z = zz(:);
    [x, zf] = rotateNoseDown(ap_h * ones(size(z)), z, -rot);   % vertical line back in the frame
    inside = F_brain(sites.ML(k) * ones(size(z)), x, zf) > 0.5;
    if ~any(inside), continue; end
    [ap(k), dv(k)] = rotateNoseDown(ap_h, max(z(inside)) - sites.depth(k), -rot);
end
end


function h = newContour(ax, hx, hy, color, level, width)
% Hidden placeholder contour at one level; setContour fills in the slice
Z = zeros(numel(hy), numel(hx)); Z(1) = 1;   % non-constant, so contour doesn't warn
[~, h] = contour(ax, hx, hy, Z, [level level], 'LineColor', color, 'LineWidth', width, 'Visible', 'off', 'HitTest', 'off');
end


function setContour(h, S)
% Swap in a new 2D slice; hidden if the slice never reaches the level
reaches = max(S(:)) >= h.LevelList(1);
if reaches, h.ZData = S; end
h.Visible = reaches;
end


function C = loadViewerCache(data_dir)
% Everything brainViewer derives from the data files: smoothed template volumes
% (loadBrain3d defaults), their 3D surfaces (isosurface at 0.5 with normals, native
% frame) and cell coordinates in the MC022 frame (Dryad through dryadToTemplate;
% MC022 from its section images). Saved in tempdir and rebuilt when any source
% file's date changes or cache_version does (bump it when the code that makes the
% cache changes: here, loadDryadCells, dryadToTemplate, rotateNoseDown,
% loadBrain3d). Lambda: the origin moves P.lambda_ap anterior along MC029's AP
% axis (MC029 is the ground truth; the same point, rotated, in MC022's frame),
% so the templates and cells move back by it; C.tf.ap_shift and C.lambda (AP,
% DV in MC022's frame) carry it to the injection sites, AP snapping and tab 2.
cache_version = 13;   % 13: ML mirror undone (the 3D view's ML axis is reversed instead)
cache_file = fullfile(tempdir, 'brainViewerCache.mat');
src = fullfile(data_dir, {'MC022.mat', 'MC029_fixed.mat', fullfile('dryad_applegate2023', 'anatomy.mat'), 'MC029_orig.mat'});
src_time = max(cellfun(@(f) dir(f).datenum, src));
if isfile(cache_file)
    stamp = load(cache_file, 'version', 'src_time');
    if isfield(stamp, 'version') && stamp.version == cache_version && stamp.src_time == src_time
        C = load(cache_file);
        return
    end
end

fprintf('Building brainViewer cache (once per data change)...\n');
P = load(src{1}, 'P').P;    % MC029's fitted angle and the vibratome scale (its ML: sagittal sections)
spec = struct('file', src(1:2), 'name', {'MC022', 'MC029_fixed'}, ...
    'label', {sprintf('MC022 (sectioned %g deg)', P.default_angle), sprintf('MC029 (%g deg)', P.ref_angle)}, ...
    'native_angle', {P.default_angle, P.ref_angle}, 'ml_scale', {1, P.section_step/0.1});
for it = 1:numel(spec)
    D = loadBrain3d(spec(it).file);
    D.yy = D.yy * spec(it).ml_scale;
    [dx, dz] = rotateNoseDown(P.lambda_ap, 0, spec(it).native_angle - P.ref_angle);   % lambda in this frame
    D.xx = D.xx - dx; D.zz = D.zz - dz;
    [X, Y, Z] = meshgrid(D.xx, D.yy, D.zz);
    t = rmfield(spec(it), 'file');
    t.D = struct('xx', D.xx, 'yy', D.yy, 'zz', D.zz);
    t.surf = struct();
    for name = intersect(["Ibrain" "Ihippo" "Idl"], string(fieldnames(D))', 'stable')
        s = isosurface(X, Y, Z, D.(name), 0.5);
        s.normals = isonormals(X, Y, Z, D.(name), s.vertices);
        t.surf.(name) = s;
        t.D.(name) = single(D.(name));   % for slice contours; single halves the file
    end
    C.tmpl(it) = t;
end
cells = loadDryadCells(data_dir, P);
C.tf = load(src{1}, 'tf').tf;
[C.ap, C.ml, C.dv] = dryadToTemplate(cells.AP, cells.ML, cells.DV, C.tf);
C.cells = cells(:, {'bird', 'ROI', 'Transverse', 'Long'});
% MC022: the cells mapped from its own section images (exact in the MC022
% frame) replace its Dryad cells. Only the regions labelled in its labelledCells
% files are available.
S = load(src{1}, 'sections', 'injections');
img = vertcat(S.sections.cells{:});
img_ap = repelem(S.sections.AP, cellfun(@height, S.sections.cells));
[~, i_inj] = ismember(img.wavelength, S.injections.wavelength);
has_inj = i_inj > 0;   % cells of a channel with no injection listed are left out
img = img(has_inj, :); img_ap = img_ap(has_inj); i_inj = i_inj(has_inj);
img_cells = table(repmat({'MC022'}, height(img), 1), categorical(img.roi), ...
    categorical(S.injections.Transverse(i_inj)), S.injections.Long(i_inj), ...
    'VariableNames', C.cells.Properties.VariableNames);
dry022 = strcmp(C.cells.bird, 'MC022');
C.cells = [C.cells(~dry022, :); img_cells];
C.ap = [C.ap(~dry022); img_ap]; C.ml = [C.ml(~dry022); img.ml]; C.dv = [C.dv(~dry022); img.dv];
[C.lambda(1), C.lambda(2)] = rotateNoseDown(P.lambda_ap, 0, P.default_angle - P.ref_angle);   % lambda in MC022's frame
C.ap = C.ap - C.lambda(1); C.dv = C.dv - C.lambda(2);
C.tf.ap_shift = C.tf.ap_shift - C.lambda(1);    % sites and AP snapping, via tf
C.cells.bird = categorical(C.cells.bird);               % text columns as categorical:
C.cells.Transverse = categorical(C.cells.Transverse);   % cellstr makes the file huge and slow
C.version = cache_version;
C.src_time = src_time;
save(cache_file, '-struct', 'C', '-v7');
end


function setPolys(h, polys, show)
% Line h through the polygons in the cell polys ([x y], NaN-separated pieces allowed)
P = cellfun(@(p) [p; NaN NaN], polys(:), 'UniformOutput', false);
P = vertcat(P{:}, [NaN NaN]);
set(h, 'XData', P(:, 1), 'YData', P(:, 2), 'Visible', show);
end


function sw = switchRow(parent, row, label, items, value)
% Row of parent (4 columns): label, items{1}, a two-way switch, items{2}. The
% switch shows no text of its own; its Value is the selected item.
lbl = uilabel(parent, 'Text', label, 'FontWeight', 'bold');
lbl.Layout.Row = row; lbl.Layout.Column = 1;
lbl = uilabel(parent, 'Text', items{1}, 'HorizontalAlignment', 'right');
lbl.Layout.Row = row; lbl.Layout.Column = 2;
sw = uiswitch(parent, 'Items', {'', ''}, 'ItemsData', items, 'Value', value);
sw.Layout.Row = row; sw.Layout.Column = 3;
lbl = uilabel(parent, 'Text', items{2});
lbl.Layout.Row = row; lbl.Layout.Column = 4;
end


function [sld, spn] = sliderRow(parent, label, limits, value, step, fmt)
% Label in one column; in the next, a slider and a spinner showing its value
% (the spinner's arrows step by step and keep stepping while held)
uilabel(parent, 'Text', label, 'FontWeight', 'bold');
g = uigridlayout(parent, [2 2], 'Padding', 0, 'RowSpacing', 0, 'ColumnSpacing', 6, ...
    'RowHeight', {24, '1x'}, 'ColumnWidth', {'1x', 80});
sld = uislider(g, 'Limits', limits, 'Value', value);
sld.Layout.Row = [1 2]; sld.Layout.Column = 1;
spn = uispinner(g, 'Limits', limits, 'Value', value, 'Step', step, 'ValueDisplayFormat', fmt);
spn.Layout.Row = 1; spn.Layout.Column = 2;
end
