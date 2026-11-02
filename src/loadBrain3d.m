function D = loadBrain3d(filename, downsamp, smoothparams, brain_smooth_factor)
%% Load (and optionally downsample and smooth) pre-analyzed data from imageTracePlot.m
% D = loadBrain3d(filename, [downsamp], [smoothparams], [brain_smooth_factor])
%
% e.g. D = loadBrain3d('HZ09.mat', 2) to load HZ09 brain and downsample by 2
% run D:\hannah\Dropbox\alab\Code\imaging\imageTracePlot.m first to save
% processed brain
% run plotBrain3d(D) after to plot the results
%
% Hannah Payne, 2020


% Control smoothing here
if exist('smoothparams','var') && ~isempty(smoothparams)
    sigmasmooth = smoothparams(1);
    nsmooth = smoothparams(2);
else
    sigmasmooth = 3; 		% sigma of gaussian smoothing
    nsmooth = 7;  			% size of window
end

if ~exist('brain_smooth_factor','var') || isempty(brain_smooth_factor)
    brain_smooth_factor = 1;
end

D = load(filename, '-regexp', '^(I.*|xx|yy|zz|P)$');   % volumes and grid only (template files also hold section tables)

if ~exist('downsamp','var') || isempty(downsamp)
    downsamp = 2;
end

% Get field names to process
field_name = fieldnames(D);
Iind = find(strncmp(field_name,'I',1));

if downsamp
    D.P.dI_downsample = D.P.dI*downsamp;
    
    for ii = 1:length(Iind)
        D.(field_name{Iind(ii)}) = D.(field_name{Iind(ii)})([1 2:downsamp:end-1 end], [1 2:downsamp:end-1 end], [1 2:downsamp:end-1 end]);
    end
    
    D.xx = D.xx(2:downsamp:end-1);
    D.xx = [D.xx(1)-D.P.dI_downsample D.xx D.xx(end)+D.P.dI_downsample];
    D.yy = D.yy(2:downsamp:end-1);
    D.yy = [D.yy(1)-D.P.dI_downsample D.yy D.yy(end)+D.P.dI_downsample];
    D.zz = D.zz(2:downsamp:end-1);
    D.zz = [D.zz(1)-D.P.dI_downsample D.zz D.zz(end)+D.P.dI_downsample];
end

for ii = 1:length(Iind)
    if strcmp(field_name{Iind(ii)},'Ibrain')
        nsmooth_curr =  min(15,floor(nsmooth*brain_smooth_factor/2)*2+1);
        sigma_curr = sigmasmooth*brain_smooth_factor;
    else
        nsmooth_curr = nsmooth;
        sigma_curr = sigmasmooth;
    end
    
    Icurr = D.(field_name{Iind(ii)}) ;
    Icurr = smooth3(double(Icurr), 'gaussian',nsmooth_curr, sigma_curr);
    Icurr(:,:, [1 end]) = 0;
    Icurr(:,[1 end], :) = 0;
    Icurr([1 end], :, :) = 0;
    D.(field_name{Iind(ii)}) = Icurr;
    
end







