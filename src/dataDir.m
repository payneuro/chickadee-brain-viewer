function data_dir = dataDir()
% data_dir = dataDir() returns the folder holding the data files (MC022.mat,
% MC029_fixed.mat, MC029_orig.mat, dryad_applegate2023\anatomy.mat, MC022_ds8\
% images), first of:
%   1. data\ in the repository (a local copy; the build scripts write here)
%   2. the shared copy on the Payne lab server (shared_dir below)
%   3. a fresh download of the shared Dropbox folder (dropbox_url) into data\,
%      after asking (~650 MB)
%
% Hannah Payne, 2026

shared_dir = 'Z:\Hannah\imaging\plotChickadeeBrain\data';
dropbox_url = 'https://www.dropbox.com/scl/fo/1kjz4tkh3fw4lckc4yrw2/AIGJRO2nt89LfSsVwOjWvlE?rlkey=ec833tj9nj63ctl9vs4zgunwh&dl=1';
required = {'MC022.mat', 'MC029_fixed.mat', 'MC029_orig.mat', fullfile('dryad_applegate2023', 'anatomy.mat')};

persistent found
if ~isempty(found), data_dir = found; return; end

local_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'data');
has_data = @(d) all(cellfun(@(f) isfile(fullfile(d, f)), required)) && isfolder(fullfile(d, 'MC022_ds8'));
if has_data(local_dir)
    data_dir = local_dir;
elseif has_data(shared_dir)
    data_dir = shared_dir;
else
    answer = input(sprintf(['No data in %s or %s.\n' ...
        'Download it from Dropbox into %s (~650 MB)? [y/n]'], local_dir, shared_dir, local_dir), 's');
    if ~strcmpi(strtrim(answer), 'y'), error('dataDir:noData', 'No data folder found.'); end
    zip_file = [tempname '.zip'];
    unzip_dir = tempname;
    cleanup = onCleanup(@() removeTemp(zip_file, unzip_dir));   % also after a failed download
    fprintf('Downloading...\n');
    websave(zip_file, dropbox_url, weboptions('Timeout', 600));
    unzip(zip_file, unzip_dir);
    hit = dir(fullfile(unzip_dir, '**', required{1}));   % the zip may or may not hold a top folder
    if ~isempty(hit)
        if ~isfolder(local_dir), mkdir(local_dir); end
        copyfile(fullfile(hit(1).folder, '*'), local_dir, 'f');   % merges into a partial data\
    end
    if ~has_data(local_dir), error('dataDir:noData', 'The download did not contain the data files.'); end
    data_dir = local_dir;
end
found = data_dir;
end


function removeTemp(zip_file, unzip_dir)
if isfile(zip_file), delete(zip_file); end
if isfolder(unzip_dir), rmdir(unzip_dir, 's'); end
end
