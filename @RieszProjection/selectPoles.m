function selectPoles(RP, radius)
%SELECTPOLES Select poles interactively with the mouse
%   SELECTPOLES(RP, radius) All Poles within a circle of the given radius 
%   and with its center at the point you clicked on will be added to a 
%   contour. To start a new contour press enter. When you are finished 
%   press enter twice. The default radius is the property of the same name.

% This file is part of the software RPExpand
% Copyright: 2021 Zuse Institute Berlin
% Authors: Fridtjof Betz, Felix Binkowski
% Updated: May-2021

disp(['Please select poles with right-clicks and confirm with enter. ' ...
    'This allows to sum up contributions of several modes. When you ' ...
    'are finished press enter twice.'])
if nargin==1
    radius = RP.radius;
end
figure(1)
clf(1)
hold on
plot(RP.expansionPoints,zeros(size(RP.expansionPoints)),'k.')
plot(RP.poles, 'bo','MarkerFaceColor','b');
axis equal
sPoles = cell(1,10);
index = 1;
while ishandle(1)
    try
        [x, y] = ginput;
    catch
        break
    end
    if isempty(x), break; end
    points = x + 1i*y;
    [~, indices] = find(abs(points-RP.poles)<radius);
    if isempty(indices), continue; end
    sPoles{index} = indices;
    plot(RP.poles(indices), 'ro','MarkerFaceColor','r')
    index = index+1;
end
sPoles = sPoles(~cellfun('isempty',sPoles));
RP.selectedPoles = sPoles;
end
