function cMap = green_magenta_Cmap(nColors)
% create green to magenta colormap, with black in the middle.
% INPUT:
%   nColors - a fi

arguments
    nColors (1,1) double {mustBeFinite,mustBeInteger,mustBeGreaterThan(nColors,1)} = 256
end

switch nColors
    case 9
        cMap = [1,     0,     1;
                0.5,   0,     0.5;
                0.25,  0,     0.25;
                0.125, 0,     0.125;
                0,     0,     0;
                0,     0.125, 0;
                0,     0.25,  0;
                0,     0.5,   0;
                0,     1,     0];
%     case 8
%          cMap = [1,     0,     1;
%                 0.5,   0,     0.5;
%                 0.25,  0,     0.25;
%                 0.125, 0,     0.125;
%                 0,     0.125, 0;
%                 0,     0.25,  0;
%                 0,     0.5,   0;
%                 0,     1,     0];
    case 7
        cMap = [1,    0,    1;
                0.5,  0,    0.5;
                0.25, 0,    0.25;
                0,    0,    0;
                0,    0.25, 0;
                0,    0.5,  0;
                0,    1,    0];
    case 6
        cMap = [1,    0,    1;
                0.5,  0,    0.5;
                0.25, 0,    0.25;
                0,    0.25, 0;
                0,    0.5,  0;
                0,    1,    0];
    case 5
        cMap = [1,   0,   1;
                0.5, 0,   0.5;
                0,   0,   0;
                0,   0.5, 0;
                0,   1,   0];
    case 4
        cMap = [1,   0,   1;
                0.5, 0,   0.5;
                0,   0.5, 0;
                0,   1,   0];
    case 3
        cMap = [1,   0,   1;
                0  , 0,   0;
                0,   1,   0];
    case 2
        cMap = [1,   0,   1;
                0,   1,   0];
    otherwise
        half_map = nColors/2;
        if floor(half_map) == half_map
            cMap = [ [linspace(1,0,half_map)' zeros(half_map,1) linspace(1,0,half_map)'];...
                [zeros(half_map,1) linspace(0,1,half_map)' zeros(half_map,1)] ];
        else
            half_map = floor(half_map);
            cMap = [ [linspace(1,0,half_map)' zeros(half_map,1) linspace(1,0,half_map)'];...
                zeros(1,3);...
                [zeros(half_map,1) linspace(0,1,half_map)' zeros(half_map,1)] ];
        end
end

% EOF