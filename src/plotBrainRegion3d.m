function h_patch = plotBrainRegion3d(xx, yy, zz, Iregion, isolevel, color, facealpha, ax)
% plotBrainRegion3d(xx,yy,zz, Iregion, isolevel, color, facealpha, [ax]) plots
% the isosurface of volume Iregion (or a precomputed surface struct with faces,
% vertices, normals) in ax (default: current axis)
%
% Example:
% figure;
% h_brain = plotBrainRegion3d(D.xx, D.yy, D.zz,D.Ibrain, 0.5, .7*[1 1 1], 0.3);
% h_hippo = plotBrainRegion3d(D.xx, D.yy, D.zz,D.Ihippo, 0.6, [0.5 0.5 1], 0.6);
%
% Hannah Payne, 2020

if isempty(isolevel)
    isolevel = 0.5;
end
if nargin < 8 || isempty(ax)
    ax = gca;
end

% Plot the region; Iregion is a volume, or a precomputed surface struct
% (faces, vertices, normals) whose isolevel is already applied
if isstruct(Iregion)
    h_patch = patch(ax, 'Faces', Iregion.faces, 'Vertices', Iregion.vertices, ...
        'VertexNormals', Iregion.normals); hold(ax, 'on')
else
    [X,Y,Z] = meshgrid(xx, yy, zz);
    iso_region = isosurface(X,Y,Z,Iregion,isolevel);
    h_patch = patch(ax, iso_region); hold(ax, 'on')
    isonormals(X,Y,Z,Iregion,h_patch)
end
h_patch.EdgeColor = 'none';
h_patch.FaceVertexCData = repmat(color, size(h_patch.Vertices,1),1);
h_patch.FaceColor = 'interp';
h_patch.FaceAlpha = facealpha;
h_patch.DiffuseStrength = 1;
h_patch.SpecularStrength = 1;
h_patch.SpecularExponent = 10;

xlabel(ax, 'Anterior (mm)','Horiz','center');
zlabel(ax, 'Dorsal (mm)','Horiz','center')
dt = 4;
set(ax,'XTick',floor(min(get(ax,'XLim'))/dt)*dt:dt:ceil(max(get(ax,'XLim'))/dt)*dt,...
    'YTick', floor(min(get(ax,'YLim'))/dt)*dt:dt:ceil(max(get(ax,'YLim'))/dt)*dt,...
    'ZTick', floor(min(get(ax,'ZLim'))/dt)*dt:dt:ceil(max(get(ax,'ZLim'))/dt)*dt)
grid(ax, 'on')

view(ax, [20 20]);

delete(findall(ax,'Type','light'))
camlight(ax, 'headlight')
h_patch.FaceLighting = 'gouraud';   % lighting(ax,...) opens a stray figure when ax is in a uifigure

axis(ax, 'equal');
axis(ax, 'vis3d');   % after equal, so hand-rotating keeps the aspect ratio
if ~isa(ax, 'matlab.ui.control.UIAxes')   % uiaxes rotate already; rotate3d costs ~1.5 s per call there
    rotate3d(ax, 'on')
end



