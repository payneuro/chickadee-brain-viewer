%% Tissue outlines for every section image, with detached pieces joined
% For each section image (1/8 scale, readSectionImage): detectSectionOutline
% finds the tissue; if it is in several detached pieces, findJoin finds how to
% move them together. Not joined: A/B images (one hemisphere each) and images with two
% drawn midlines (two hemispheres mounted apart), which are joined along their
% midlines instead, and images without a usable midline yet.
%
% Output: data\<bird>_outlines.mat, struct array S (one per image, metaData
% order): image_file, n_pieces, joined, outline (cell of [row col] boundaries at
% 1/8 scale, after joining), lab (piece labels before joining, uint8; joined
% images only), steps (findJoin), peaks and threshold (detectSectionOutline),
% error ("" or why detection failed; the run continues).
% Joined images are not saved: joinPieces(readSectionImage(f), lab, steps, bg)
% rebuilds one in a fraction of a second.
% Takes ~3 min plus ~15-30 s per joined image; run with matlab -batch.
%
% Hannah Payne, 2026

bird = 'MC022';
no_join = 1;   % images never joined: k = 1 has two hemispheres mounted apart but one (wrong) midline, to be drawn
bird_dir = fullfile('Z:\Hannah\imaging', [bird ' sortedSlices']);
code_dir = fileparts(fileparts(mfilename('fullpath')));   % repository root
data_dir = fullfile(code_dir, 'data');

meta = load(fullfile(bird_dir, 'metaDataUpdated3.mat')).metaData;
corr_file = fullfile(bird_dir, 'midlineCorrections.mat');
if isfile(corr_file)                                     % hand-drawn midlines replace metaData's
    for c = load(corr_file).corrections
        if ~isempty(c.midLineCorr), meta(c.k).midLineCorr = c.midLineCorr; end
    end
end
files = dir(fullfile(bird_dir, 'images', 'slice_*_640.tif'));
S = struct('image_file', {}, 'n_pieces', {}, 'joined', {}, 'outline', {}, 'lab', {}, ...
    'steps', {}, 'peaks', {}, 'threshold', {}, 'error', {});
for k = 1:numel(files)
    tic
    f = fullfile(files(k).folder, files(k).name);
    img = readSectionImage(f);
    try
        [bw, info] = detectSectionOutline(img);
    catch err
        S(k) = struct('image_file', string(f), 'n_pieces', 0, 'joined', false, 'outline', {{}}, ...
            'lab', [], 'steps', [], 'peaks', [NaN NaN], 'threshold', NaN, 'error', string(err.message));
        fprintf('%d/%d %s: FAILED, %s\n', k, numel(files), files(k).name, err.message);
        continue
    end
    lab = bwlabel(bw);
    s = struct('image_file', string(f), 'n_pieces', max(lab(:)), 'joined', false, 'outline', [], ...
        'lab', [], 'steps', [], 'peaks', info.peaks, 'threshold', info.threshold, 'error', "");
    joinable = isempty(regexp(files(k).name, '\d[AB]_', 'once')) && size(meta(k).midLineCorr, 1) == 2 && ~ismember(k, no_join);
    if s.n_pieces > 1 && joinable
        [s.steps, ~, bw] = findJoin(bw);
        s.joined = true;
        s.lab = uint8(lab);
    end
    s.outline = bwboundaries(bw, 'noholes');
    S(k) = s;
    fprintf('%d/%d %s: %d pieces%s (%.0f s)\n', k, numel(files), files(k).name, s.n_pieces, ...
        repmat(', joined', 1, s.joined), toc);
end
save(fullfile(data_dir, [bird '_outlines.mat']), 'S', '-v7');
fprintf('Joined %d images: k = %s\n', nnz([S.joined]), mat2str(find([S.joined])));
fprintf('Failed %d images: k = %s\n', nnz([S.error] ~= ""), mat2str(find([S.error] ~= "")));
