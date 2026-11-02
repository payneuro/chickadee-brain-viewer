function [canvas, ks] = renderSection(sections, order, ML, DV, channel)
% [canvas, ks] = renderSection(sections, order, ML, DV, [channel]) combines the
% images of one section (trueOrder = order: a whole section, or an A/B pair) in
% the template frame at points ML, DV (mm): each 1/8 image (readSectionImage),
% with its detached pieces joined (joinPieces, sections.join_lab/join_steps),
% scaled so its background is 0 and 1 is a default display maximum (not
% clipped, so the display range can be set later), placed by
% pixelsToTemplateImage, and the images merged by maximum. canvas is single, NaN
% off the images (so the background keeps its values below 0). sections: the
% template's sections table. ks: its rows for this section that are aligned
% (shift not NaN); canvas is all NaN if none.
% channel (optional): image channel, e.g. "395" (DAPI), "475", "555"; default the
% channel of sections.image_file (640, where outlines and joins were found).
% Scaling: for the outline channel, from sections.peaks (log intensity of the
% background and tissue histogram peaks): 0 = background peak, 1 = tissue peak
% plus 2x the gap between the peaks. For other channels, 0 = median of the
% image (mostly background), 1 = its 99.9th percentile.
%
% Hannah Payne, 2026

ds = 8;
contrast_spread = 2;   % display max = tissue peak + this x (tissue - background peak), log scale
ks = find(sections.trueOrder == order & ~isnan(sections.shift))';
canvas = NaN(size(ML), 'single');
for k = ks
    file = sections.image_file(k);
    own_channel = nargin < 5 || endsWith(file, "_" + channel + ".tif");
    if ~own_channel
        file = regexprep(file, '_\d+\.tif$', "_" + channel + ".tif");
    end
    img = readSectionImage(file);
    if own_channel
        pk = sections.peaks(k, :);
        lim = exp([pk(1), pk(2) + contrast_spread * diff(pk)]);
    else
        v = double(img(img > 0));   % 0: padding in stitched images
        lim = [median(v), prctile(v, 99.9)];
    end
    if ~isempty(sections.join_steps{k})
        img = joinPieces(img, sections.join_lab{k}, sections.join_steps{k}, lim(1));
    end
    img = (double(img) - lim(1)) / diff(lim);
    canvas = max(canvas, single(pixelsToTemplateImage(sections(k, :), img, ds, ML, DV)));   % max ignores NaN
end
end
