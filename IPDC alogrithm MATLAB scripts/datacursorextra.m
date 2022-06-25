%% <https://stackoverflow.com/questions/25477401/matlab-permanently-set-text-update-function-for-data-cursor-in-figures/39330765#39330765>
function datacursorextra(fig)
    % Use current figure as default
    if nargin<1
        fig = gcf;
    end

    % Get the figure's datacursormode, and set the update function
    h = datacursormode(fig);
    set(h,'UpdateFcn',@myupdatefcn)

    % The actual update function
    function txt = myupdatefcn(~,event)
        % Short-hand to write X, Y and if available Z, with 10 digit precision:
        lbl = 'XYZ';
        txt = arrayfun(@(s,g)sprintf('%s: %.10g',s,g), lbl(1:length(event.Position)), event.Position,'uniformoutput',false);

        % If a DataIndex is available, show that also:
        info = getCursorInfo(h);
        if isfield(info,'DataIndex')
            txt{end+1} = sprintf('Index: %d', info.DataIndex);
        end
    end
end