function I = cellDensity3d(ap, ml, dv, xx, yy, zz, smoothparams)
% I = cellDensity3d(ap, ml, dv, xx, yy, zz, [smoothparams])
% Smoothed 3D cell density for isosurface plotting: bin cells at (ap,ml,dv) onto
% grid xx/yy/zz, normalize by the median nonzero count, gaussian smooth. Output
% is [numel(yy) numel(xx) numel(zz)], matching loadBrain3d volumes.
% smoothparams = [sigma window], default [3 7].
%
% Hannah Payne, 2026

if nargin < 7 || isempty(smoothparams)
    smoothparams = [3 7];
end
sigma = smoothparams(1);
nsmooth = smoothparams(2);

xx = xx(:); yy = yy(:); zz = zz(:);
ex = [-inf; (xx(1:end-1)+xx(2:end))/2; +inf];
ey = [-inf; (yy(1:end-1)+yy(2:end))/2; +inf];
ez = [-inf; (zz(1:end-1)+zz(2:end))/2; +inf];

xBin = discretize(ap, ex);
yBin = discretize(ml, ey);
zBin = discretize(dv, ez);
valid = ~isnan(xBin) & ~isnan(yBin) & ~isnan(zBin);

counts = accumarray([yBin(valid), xBin(valid), zBin(valid)], 1, [numel(yy) numel(xx) numel(zz)]);
nonzero_counts = counts(counts>0);
I = counts ./ median(nonzero_counts);

% Smooth and zero the borders (so the isosurface closes), matching loadBrain3d
I = smooth3(I, 'gaussian', nsmooth, sigma);
I(:,:,[1 end]) = 0;
I(:,[1 end],:) = 0;
I([1 end],:,:) = 0;
