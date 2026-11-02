function [img_j, bw_j, masks] = joinPieces(img, lab, steps, bg)
% [img_j, bw_j, masks] = joinPieces(img, lab, steps, bg) applies the moves found
% by findJoin: each piece of lab (piece labels) is rotated and shifted by every
% step that moved it, in order, and its tissue pasted at its new place on a
% background of bg. img may be [] to move only the masks.
% img_j: joined image (class of img); bw_j: joined mask; masks: cell, each
% piece's moved mask.
%
% Hannah Payne, 2026

n = max(lab(:));
masks = cell(1, n);
img_j = [];
if ~isempty(img), img_j = bg * ones(size(img)); end
for p = 1:n
    m = double(lab == p);
    I = double(img);
    for s = steps
        if ~ismember(p, s.moved), continue; end
        m = transformImage(m, s, 'nearest');
        if ~isempty(img), I = transformImage(I, s, 'linear'); end
    end
    masks{p} = m > 0.5;
    if ~isempty(img), img_j(masks{p}) = I(masks{p}); end
end
bw_j = any(cat(3, masks{:}), 3);
if ~isempty(img), img_j = cast(img_j, 'like', img); end
end


function out = transformImage(in, s, method)
% Rotate by s.theta deg about s.ctr = [row col], then shift by s.shift; 0 outside
[cc, rr] = meshgrid(1:size(in, 2), 1:size(in, 1));
R = [cosd(s.theta) -sind(s.theta); sind(s.theta) cosd(s.theta)];
src = ([rr(:) cc(:)] - s.shift - s.ctr) * R + s.ctr;   % inverse: un-shift, rotate back
out = reshape(interp2(in, src(:, 2), src(:, 1), method, 0), size(in));
end
