function stop = FOkinPlotObjectiveModel(results, state)

if (size(results.VariableDescriptions, 2) == 1)
    % 1D problem, execute the original function
    stop = plotObjectiveModel(results, state);
else
    % 2D problem, draw a vertical line across the position of the minimum.
    h_model_fig = findobj('Tag', 'bayesopt.ObjectiveModel');
    if ~isempty(h_model_fig)
        hlin = findobj(h_model_fig, 'DisplayName', 'Model minimum feasible');
        if ~isempty(hlin)
            delete(hlin)
        end
    end
    
    stop = plotObjectiveModel(results, state);
    
    if ~strcmp(state, 'initial')
        h_model_fig = findobj('Tag', 'bayesopt.ObjectiveModel');
        hlin = findobj(h_model_fig, 'DisplayName', 'Model minimum feasible');
        hax = hlin.Parent;
        x = hlin.XData;
        y = hlin.YData;
        delete(hlin)
        line(hax, [x, x], [y, y], hax.ZLim, 'Color','c',...
            'LineStyle', ':',  'LineWidth', 3,...
            'DisplayName', 'Model minimum feasible' );
        drawnow
    end
end
end