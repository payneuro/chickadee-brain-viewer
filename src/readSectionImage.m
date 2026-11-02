function [img, full_size] = readSectionImage(image_file)
% [img, full_size] = readSectionImage(image_file) returns a section image at 1/8
% of full resolution (8x8 block average, uint16) and the full-res [height width].
% Images are cached as 16-bit PNGs in <dataDir>\<bird>_ds8\ (bird from the
% "<bird> sortedSlices" folder), since the full-res TIFFs (~1 GB each) take
% 5-8 s to read from Z:. Once cached, the TIFFs are not needed.
%
% Hannah Payne, 2026

ds = 8;
[images_dir, name] = fileparts(char(image_file));
[~, bird_folder] = fileparts(fileparts(images_dir));
bird = strtrim(erase(bird_folder, 'sortedSlices'));
cache_dir = fullfile(dataDir(), [bird '_ds8']);
cache_file = fullfile(cache_dir, [name '.png']);
if isfile(cache_file)
    img = imread(cache_file);
    full_size = sscanf(imfinfo(cache_file).Description, 'full_size %d %d')';
    return
end
full = imread(image_file);
full_size = size(full, [1 2]);
img = imresize(full, 1/ds, 'box');   % 8x8 block average (every 8th pixel alone is speckled)
if ~isfolder(cache_dir), mkdir(cache_dir); end
imwrite(img, cache_file, 'Description', sprintf('full_size %d %d', full_size));
end
