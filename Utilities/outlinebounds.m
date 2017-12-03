function hnew = outlinebounds(hl, hp)
%OUTLINEBOUNDS Outline the patch of a boundedline
%
% hnew = outlinebounds(hl, hp)
%
% This function adds an outline to the patch objects created by
% boundedline, matching the color of the central line associated with each
% patch.
%
% Input variables:
%
%   hl:     handles to line objects from boundedline
%
%   hp:     handles to patch objects from boundedline
%
% Output variables:
%
%   hnew:   handle to new line objects

% Copyright 2012 Kelly Kearney


hnew = zeros(size(hl));
for il = 1:length(hp)
    col = get(hl(il), 'color');
    xy = get(hp(il), {'xdata','ydata'});
    ax = ancestor(hl(il), 'axes');
    
    hnew(il) = line(xy{1}, xy{2}, 'parent', ax, 'linestyle', '-', 'color', col);
end
    
