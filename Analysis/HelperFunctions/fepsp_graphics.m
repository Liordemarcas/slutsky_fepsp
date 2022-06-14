function hndl = fepsp_graphics(hndl)

% gets a figure / axes handle and applyes default graphics without changing
% the root 

hndlType = get(hndl, 'type');

if strcmp(hndlType, 'axes')
    set(hndl, 'FontSize', 16)
    set(hndl, 'FontUnits', 'points')
    set(hndl, 'box', 'off')
    set(hndl, 'TickLength', [0 0])
    set(hndl, 'TitleFontSizeMultiplier', 1.4)
    set(hndl, 'LabelFontSizeMultiplier', 1.2)
    set(hndl, 'FontName', 'Arial')
    
elseif strcmp(hndlType, 'figure')
    set(hndl, 'units', 'normalized')
    set(hndl, 'position', [0.15 0.15 0.7 0.7])
    set(hndl, 'Color', 'w')

end
    
end

% EOF
