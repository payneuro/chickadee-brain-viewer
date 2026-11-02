function [steps, lab, bw_j] = findJoin(bw, opts)
% [steps, lab, bw_j] = findJoin(bw, opts) finds how to move the detached
% tissue pieces of one section image (bw, from detectSectionOutline) together:
% the gaps between them are mounting artifacts. Closest pair first: the smaller
% group (by area) is fitted against the larger one, which stays put; joined
% pieces then move together. A fit tries each rotation in opts.angles (deg,
% about the moving group's centroid), moves the group straight toward the other
% until first contact, then slides it (within +-opts.reach px) to minimize the
% overall gap: the mean distance from its boundary to the other group, capped
% at opts.cap px so its far side doesn't count, never closer than 2 px.
% steps: struct per join with moved and onto (piece numbers), theta (deg),
% shift ([drow dcol] px) and ctr ([row col], rotation center); lab: piece labels
% of bw; bw_j: the joined mask. joinPieces moves images and masks by steps;
% joinPoints (in makeAnatomyTemplates) moves positions.
%
% Hannah Payne, 2026

if nargin < 2, opts = struct(); end
defaults = struct('angles', -10:10, 'reach', 125, 'cap', 60);   % reach 125 px = 1 mm at 1/8
for f = fieldnames(defaults)'
    if ~isfield(opts, f{1}), opts.(f{1}) = defaults.(f{1}); end
end
lab = bwlabel(bw);
n = max(lab(:));
group = 1:n;
steps = struct('moved', {}, 'onto', {}, 'theta', {}, 'shift', {}, 'ctr', {});
cur_mask = arrayfun(@(p) lab == p, 1:n, 'UniformOutput', false);
while numel(unique(group)) > 1
    gids = unique(group);
    G = arrayfun(@(g) any(cat(3, cur_mask{group == g}), 3), gids, 'UniformOutput', false);
    best = inf;
    for i = 1:numel(gids)
        D = bwdist(G{i});
        for j = i+1:numel(gids)
            d = min(D(G{j}));
            if d < best, best = d; pair = [i j]; end
        end
    end
    [~, small] = min([nnz(G{pair(1)}) nnz(G{pair(2)})]);
    mv = pair(small); fx = pair(3 - small);
    [rr, cc] = find(G{mv}); ctr = mean([rr cc], 1);
    move = fitGroup(G{mv}, ctr, G{fx}, opts);
    steps(end+1) = struct('moved', find(group == gids(mv)), 'onto', find(group == gids(fx)), ...
        'theta', move(1), 'shift', move(2:3), 'ctr', ctr); %#ok<AGROW>
    group(group == gids(mv)) = gids(fx);
    [~, bw_j, cur_mask] = joinPieces([], lab, steps);
end
if n <= 1, bw_j = bw; end
end


function move = fitGroup(piece, ctr, fixed, opts)
% Best [theta drow dcol] for piece against fixed (see findJoin)
[D, idx] = bwdist(fixed);
[nr, nc] = size(D);
[r, c] = find(bwperim(piece));
gapAt = @(P) sampleD(D, P, nr, nc);
best = struct('cost', inf, 'move', [0 0 0]);
for theta = opts.angles
    R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
    P = ([r c] - ctr) * R' + ctr;                          % rotated boundary, [row col]
    d = gapAt(P);
    if any(isnan(d)), continue; end
    [~, j] = min(d);
    [tr, tc] = ind2sub([nr nc], idx(sub2ind([nr nc], round(P(j, 1)), round(P(j, 2)))));
    v = [tr tc] - P(j, :);
    lo = 0; hi = 1;                                        % largest fraction of v without contact
    for it = 1:12
        t = (lo + hi) / 2;
        dt = gapAt(P + t * v);
        if any(isnan(dt) | dt < 2), hi = t; else, lo = t; end
    end
    s0 = lo * v;
    fit = struct('cost', mean(min(gapAt(P + s0), opts.cap)), 'ds', [0 0]);
    for step = [4 1]                                       % coarse grid, then fine around the best
        span = -opts.reach:step:opts.reach; if step == 1, span = -4:4; end
        center = fit.ds;
        for dr = center(1) + span
            for dc = center(2) + span
                dd = gapAt(P + s0 + [dr dc]);
                if any(isnan(dd) | dd < 2), continue; end
                cost = mean(min(dd, opts.cap));
                if cost < fit.cost, fit = struct('cost', cost, 'ds', [dr dc]); end
            end
        end
    end
    if fit.cost < best.cost, best = struct('cost', fit.cost, 'move', [theta, s0 + fit.ds]); end
end
move = best.move;
end


function d = sampleD(D, P, nr, nc)
% Distance map at the nearest pixel of each point; NaN off the image
ri = round(P(:, 1)); ci = round(P(:, 2));
d = nan(size(ri));
in = ri >= 1 & ri <= nr & ci >= 1 & ci <= nc;
d(in) = D(sub2ind([nr nc], ri(in), ci(in)));
end

