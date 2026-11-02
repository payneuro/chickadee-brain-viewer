function [bw, info] = detectSectionOutline(img)
% [bw, info] = detectSectionOutline(img) finds the tissue in a section image at
% 1/8 scale (readSectionImage; stitched images are 0-padded).
%   1. Valley threshold on the smoothed log image: the lowest point of the
%      histogram between its two highest peaks (background and tissue), so small
%      bright modes such as saturated strips don't move it. Blobs touching the
%      image or padding edge are dropped (unless that removes most of the mask);
%      blobs > 20% of the largest are kept, holes filled.
%   2. Edge cut: the mask is cut along Canny edges of the log image; pieces that
%      touch a component's core (deeper than 200 um, or than half that
%      component's depth) and not the image or padding edge are kept (drops glow
%      and partner hemispheres); the cuts are closed, holes filled.
% bw: tissue mask (may hold several detached pieces, see findJoin).
% info: peaks (log intensity of background and tissue, for display contrast),
% threshold (log intensity).
%
% Hannah Payne, 2026

valid = imerode(img > 0, strel('disk', 3));
L = imgaussfilt(log(max(double(img), 1)), 4);
[info.threshold, info.peaks] = valleyThreshold(L(valid));
bw = L > info.threshold & valid;
bw = clearEdgePieces(bw, valid);
bw = imfill(keepLarge(bw, 0.2), 'holes');

% Edge cut
E = edge(L, 'canny') & valid;
pieces = bw & ~imdilate(E, strel('disk', 1));
D = bwdist(~bw);
comp = bwlabel(bw);
depth = accumarray(comp(comp > 0), D(comp > 0), [], @max);
core = comp > 0;
core(core) = D(core) > min(25, depth(comp(core)) / 2);    % 25 px = 200 um
lab = bwlabel(pieces);
bw = ismember(lab, unique(lab(core & lab > 0)));
bw = clearEdgePieces(bw, valid);
bw = imfill(imclose(keepLarge(bw, 0.02), strel('disk', 3)), 'holes');   % pieces already touch a core: only crumbs go
end


function [t, peaks] = valleyThreshold(v)
% Lowest point of the smoothed histogram between its two highest peaks. Peaks
% are found on log counts (prominence >= 1, a factor of e), so a tissue mode a
% few % the height of the background one (mostly empty slides) still counts.
% Empty bins pad both ends, so a peak at the edge of the range keeps its height.
edges = linspace(prctile(v, 0.1), prctile(v, 99.9), 257);
centers = (edges(1:end-1) + edges(2:end)) / 2;
n = smoothdata(histcounts(v, edges), 'gaussian', 9);
pk = find(islocalmax(log([0 n 0] + 1), 'MinProminence', 1)) - 1;
assert(numel(pk) >= 2, 'detectSectionOutline: intensity histogram has one peak');
[~, order] = sort(n(pk), 'descend'); pk = sort(pk(order(1:2)));
[~, i] = min(n(pk(1):pk(2)));
t = centers(pk(1) + i - 1);
peaks = centers(pk);
end


function bw = clearEdgePieces(bw, valid)
% Drop blobs touching the image or padding edge, unless that leaves < 20%
cleared = imclearborder(bw | ~valid) & valid;
if nnz(cleared) > 0.2 * nnz(bw), bw = cleared; end
end


function bw = keepLarge(bw, frac)
% Blobs larger than frac of the largest
a = regionprops(bw, 'Area');
if ~isempty(a), bw = bwareaopen(bw, round(frac * max([a.Area]))); end
end
