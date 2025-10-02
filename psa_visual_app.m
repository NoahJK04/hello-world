function app = psa_visual_app()
%PSA_VISUAL_APP Interactive PSA nitrogen purification visualiser.
%   Launches a UI that allows real-time tuning of key process parameters and
%   visualises the resulting nitrogen concentration, pressure trajectory, and
%   cycle metrics.
%
%   See also PSA_SIMULATION, PSA_DEFAULT_PARAMETERS, PSA_DEFAULT_SCHEDULE.

app = struct();
app.state = struct('autoUpdate', true, ...
                   'isRunning', false, ...
                   'needsRerun', false, ...
                   'lastResults', []);

app.fig = uifigure('Name', 'PSA Nitrogen Visual Simulator', ...
                   'Color', 'w', ...
                   'Position', [100, 100, 1180, 680]);
mainGrid = uigridlayout(app.fig, [1, 2]);
mainGrid.ColumnWidth = {320, '1x'};
mainGrid.RowHeight = {'1x'};

% Control panel ------------------------------------------------------------
controlPanel = uipanel(mainGrid, 'Title', 'Process controls', ...
    'FontWeight', 'bold');
controlPanel.Scrollable = 'on';
controlGrid = uigridlayout(controlPanel, [14, 2]);
controlGrid.RowHeight = repmat({'fit'}, 1, 14);
controlGrid.ColumnWidth = {'1x', 90};

row = 1;
addSpinner('Feed pressure (bar)', [3, 12], 6.5, 0.1, 'feedPressure');
addSpinner('Low pressure (bar)', [0.8, 3], 1.2, 0.05, 'lowPressure');
addSpinner('Pressurisation time (s)', [20, 120], 40, 2, 'pressurisationTime');
addSpinner('Adsorption time (s)', [60, 320], 150, 5, 'adsorptionTime');
addSpinner('Blowdown time (s)', [30, 160], 70, 5, 'blowdownTime');
addSpinner('Purge time (s)', [60, 220], 110, 5, 'purgeTime');
addSpinner('Feed flow (mol/s)', [0.4, 1.8], 1.1, 0.05, 'feedFlow');
addSpinner('Purge flow (mol/s)', [0.1, 0.8], 0.35, 0.02, 'purgeFlow');
addSpinner('Mass transfer coeff (1/s)', [0.05, 0.3], 0.12, 0.005, 'massTransfer');
addSpinner('Cycles to simulate', [2, 12], 8, 1, 'cycles');

% Auto update toggle and buttons
uilabel(controlGrid, 'Text', 'Auto update', 'FontWeight', 'bold', ...
    'Layout', struct('Row', row, 'Column', 1));
app.controls.autoSwitch = uiswitch(controlGrid, 'slider');
app.controls.autoSwitch.Items = {'Manual', 'Auto'};
app.controls.autoSwitch.Value = 'Auto';
app.controls.autoSwitch.Tooltip = 'Auto-runs the simulation when parameters change.';
app.controls.autoSwitch.Layout.Row = row;
app.controls.autoSwitch.Layout.Column = 2;
app.controls.autoSwitch.ValueChangedFcn = @toggleAutoUpdate;
row = row + 1;

app.controls.runButton = uibutton(controlGrid, 'Text', 'Run simulation', ...
    'ButtonPushedFcn', @(src, evt) triggerSimulation());
app.controls.runButton.Layout.Row = row;
app.controls.runButton.Layout.Column = [1 2];
row = row + 1;

app.controls.resetButton = uibutton(controlGrid, 'Text', 'Reset defaults', ...
    'ButtonPushedFcn', @(src, evt) resetDefaults());
app.controls.resetButton.Layout.Row = row;
app.controls.resetButton.Layout.Column = [1 2];
row = row + 1;

app.outputs.statusLabel = uilabel(controlGrid, 'Text', 'Ready', ...
    'FontAngle', 'italic', 'FontColor', [0.25 0.25 0.25]);
app.outputs.statusLabel.Layout.Row = row;
app.outputs.statusLabel.Layout.Column = [1 2];
row = row + 1;

app.legendPanel = uipanel(controlGrid, 'Title', 'Cycle steps');
app.legendPanel.Layout.Row = row;
app.legendPanel.Layout.Column = [1 2];
app.legendGrid = uigridlayout(app.legendPanel, [1, 1]);
app.legendGrid.RowHeight = {'fit'};
app.legendGrid.ColumnWidth = {'1x'};
row = row + 1;

% Visualisation side -------------------------------------------------------
visualPanel = uipanel(mainGrid, 'Title', 'Simulation outputs', 'FontWeight', 'bold');
visualGrid = uigridlayout(visualPanel, [3, 1]);
visualGrid.RowHeight = {'fit', '1x', '1x'};
visualGrid.Padding = [10 10 10 10];

summaryPanel = uipanel(visualGrid, 'Title', 'Cycle metrics');
summaryGrid = uigridlayout(summaryPanel, [2, 4]);
summaryGrid.RowHeight = {'fit', 'fit'};
summaryGrid.ColumnWidth = {'1x', '1x', '1x', '1x'};

app.outputs.purityLabel = addSummaryLabel('Product purity', '–', 1);
app.outputs.recoveryLabel = addSummaryLabel('N_2 recovery', '–', 2);
app.outputs.cycleTimeLabel = addSummaryLabel('Cycle time (min)', '–', 3);
app.outputs.productFlowLabel = addSummaryLabel('Average product flow (mol/s)', '–', 4);

app.axes.pressure = uiaxes(visualGrid);
app.axes.pressure.Layout.Row = 2;
app.axes.pressure.XLabel.String = 'Time (min)';
app.axes.pressure.YLabel.String = 'Pressure (bar)';
app.axes.pressure.Title.String = 'Column pressure profile';
app.axes.pressure.Toolbar.Visible = 'off';
app.axes.pressure.Interactions = [];

axtitlecolor = [0.1 0.1 0.1];
app.axes.pressure.Title.Color = axtitlecolor;
app.axes.pressure.XLabel.FontWeight = 'bold';
app.axes.pressure.YLabel.FontWeight = 'bold';
app.axes.pressure.FontSize = 11;

axt = uiaxes(visualGrid);
axt.Layout.Row = 3;
axt.XLabel.String = 'Time (min)';
axt.YLabel.String = 'Gas-phase mole fraction';
axt.Title.String = 'Gas composition and product window';
axt.Toolbar.Visible = 'off';
axt.Interactions = [];
axt.Title.Color = axtitlecolor;
axt.XLabel.FontWeight = 'bold';
axt.YLabel.FontWeight = 'bold';
axt.FontSize = 11;
app.axes.composition = axt;

updateStatus('Initialising simulation...');
triggerSimulation();

if nargout > 0
    app.fig.UserData = app;
end

    function lbl = addSummaryLabel(titleText, valueText, position)
        lblTitle = uilabel(summaryGrid, 'Text', titleText, ...
            'FontWeight', 'bold', 'HorizontalAlignment', 'center');
        lblTitle.Layout.Row = 1;
        lblTitle.Layout.Column = position;
        lbl = uilabel(summaryGrid, 'Text', valueText, ...
            'FontSize', 18, 'FontWeight', 'bold', ...
            'HorizontalAlignment', 'center');
        lbl.Layout.Row = 2;
        lbl.Layout.Column = position;
    end

    function addSpinner(labelText, limits, defaultValue, step, field)
        lbl = uilabel(controlGrid, 'Text', labelText, 'FontWeight', 'bold');
        lbl.Layout.Row = row;
        lbl.Layout.Column = 1;
        spinner = uispinner(controlGrid, 'Limits', limits, ...
            'Value', defaultValue, 'Step', step);
        spinner.Layout.Row = row;
        spinner.Layout.Column = 2;
        spinner.ValueChangedFcn = @(src, evt) parameterChanged();
        spinner.ValueDisplayFormat = '%.3g';
        spinner.Tooltip = labelText;
        app.controls.(field) = spinner;
        row = row + 1;
    end

    function toggleAutoUpdate(~, ~)
        app.state.autoUpdate = strcmp(app.controls.autoSwitch.Value, 'Auto');
        if app.state.autoUpdate
            updateStatus('Auto update enabled.');
            triggerSimulation();
        else
            updateStatus('Auto update disabled – click Run simulation to update.');
        end
    end

    function resetDefaults()
        app.controls.feedPressure.Value = 6.5;
        app.controls.lowPressure.Value = 1.2;
        app.controls.pressurisationTime.Value = 40;
        app.controls.adsorptionTime.Value = 150;
        app.controls.blowdownTime.Value = 70;
        app.controls.purgeTime.Value = 110;
        app.controls.feedFlow.Value = 1.1;
        app.controls.purgeFlow.Value = 0.35;
        app.controls.massTransfer.Value = 0.12;
        app.controls.cycles.Value = 8;
        updateStatus('Defaults restored.');
        triggerSimulation();
    end

    function parameterChanged()
        if app.state.autoUpdate
            triggerSimulation();
        else
            updateStatus('Parameters changed – click Run simulation.');
        end
    end

    function triggerSimulation()
        if app.state.isRunning
            app.state.needsRerun = true;
            return;
        end
        app.state.isRunning = true;
        updateStatus('Running simulation...');
        drawnow limitrate;
        try
            results = runSimulation();
            app.state.lastResults = results;
            updateVisuals(results);
            updateStatus(sprintf('Updated: purity %.1f%%, recovery %.1f%%.', ...
                100*results.metrics.product_purity, 100*results.metrics.n2_recovery));
        catch simulationError
            warning(simulationError.identifier, '%s', simulationError.message);
            uialert(app.fig, simulationError.message, 'Simulation error');
            updateStatus('Simulation failed – see command window for details.');
        end
        app.state.isRunning = false;
        if app.state.needsRerun
            app.state.needsRerun = false;
            triggerSimulation();
        end
    end

    function results = runSimulation()
        settings = readSettings();
        overrides = struct('feed_pressure', settings.feedPressurePa, ...
                           'low_pressure', settings.lowPressurePa, ...
                           'mass_transfer_coeff', settings.massTransfer);
        params = psa_default_parameters(overrides);
        schedule = psa_default_schedule(params);
        schedule = tuneSchedule(schedule, params, settings);
        results = psa_simulation('plot', false, ...
                                 'cycles', settings.cycles, ...
                                 'params', overrides, ...
                                 'schedule', schedule);
        results.schedule = schedule; %#ok<STRNU>
        results.settings = settings;
    end

    function settings = readSettings()
        settings = struct();
        settings.feedPressureBar = app.controls.feedPressure.Value;
        settings.lowPressureBar = app.controls.lowPressure.Value;
        settings.feedPressurePa = settings.feedPressureBar * 1e5;
        settings.lowPressurePa = settings.lowPressureBar * 1e5;
        settings.pressurisationTime = app.controls.pressurisationTime.Value;
        settings.adsorptionTime = app.controls.adsorptionTime.Value;
        settings.blowdownTime = app.controls.blowdownTime.Value;
        settings.purgeTime = app.controls.purgeTime.Value;
        settings.feedFlow = app.controls.feedFlow.Value;
        settings.purgeFlow = app.controls.purgeFlow.Value;
        settings.massTransfer = app.controls.massTransfer.Value;
        settings.cycles = round(app.controls.cycles.Value);
        app.controls.cycles.Value = settings.cycles;
    end

    function schedule = tuneSchedule(schedule, params, settings)
        schedule(1).duration = settings.pressurisationTime;
        schedule(1).P_target = params.feed_pressure;
        schedule(1).F_in = max(settings.feedFlow * 0.65, 0.05);
        schedule(1).F_out = 0.05;

        schedule(2).duration = settings.adsorptionTime;
        schedule(2).P_target = params.feed_pressure;
        schedule(2).F_in = settings.feedFlow;
        schedule(2).F_out = settings.feedFlow;

        schedule(3).duration = settings.blowdownTime;
        schedule(3).P_target = params.low_pressure;
        schedule(3).F_out = max(settings.feedFlow * 1.2, 0.3);

        schedule(4).duration = settings.purgeTime;
        schedule(4).P_target = params.low_pressure * 1.2;
        schedule(4).F_in = settings.purgeFlow;
        schedule(4).F_out = max(settings.purgeFlow * 1.1, settings.purgeFlow);
    end

    function updateVisuals(results)
        schedule = results.schedule;
        stepColours = lines(numel(schedule));
        updateLegend(schedule, stepColours);

        timeMin = results.time / 60;
        pressureBar = results.states(:, 5) / 1e5;
        nN2 = results.states(:, 1);
        nO2 = results.states(:, 2);
        yN2 = nN2 ./ max(nN2 + nO2, eps);
        yO2 = nO2 ./ max(nN2 + nO2, eps);
        stepIdx = results.step_indices;

        plotWithStepBackground(app.axes.pressure, timeMin, pressureBar, stepIdx, stepColours);
        app.axes.pressure.YLimMode = 'auto';

        plotComposition(app.axes.composition, timeMin, yN2, yO2, stepIdx, schedule, stepColours);

        avgProductFlow = results.metrics.product_total / max(results.metrics.cycle_time * results.settings.cycles, eps);
        app.outputs.purityLabel.Text = sprintf('%.1f%%', 100 * results.metrics.product_purity);
        app.outputs.recoveryLabel.Text = sprintf('%.1f%%', 100 * results.metrics.n2_recovery);
        app.outputs.cycleTimeLabel.Text = sprintf('%.2f', results.metrics.cycle_time / 60);
        app.outputs.productFlowLabel.Text = sprintf('%.3f', avgProductFlow);
    end

    function plotWithStepBackground(ax, x, y, stepIdx, colours)
        cla(ax);
        hold(ax, 'on');
        ypad = max(range(y), eps) * 0.08;
        ylo = min(y) - ypad;
        yhi = max(y) + ypad;
        segments = segmentIndices(stepIdx);
        for s = 1:size(segments, 1)
            idxRange = segments(s, :);
            colour = colours(stepIdx(idxRange(1)), :);
            patch(ax, [x(idxRange(1)), x(idxRange(2)), x(idxRange(2)), x(idxRange(1))], ...
                [ylo, ylo, yhi, yhi], colour, 'FaceAlpha', 0.08, 'EdgeColor', 'none');
        end
        plot(ax, x, y, 'LineWidth', 1.8, 'Color', [0.05 0.3 0.7]);
        ax.Box = 'on';
        ax.XGrid = 'on';
        ax.YGrid = 'on';
        hold(ax, 'off');
    end

    function plotComposition(ax, x, yN2, yO2, stepIdx, schedule, colours)
        cla(ax);
        hold(ax, 'on');
        ypad = 0.05;
        ylo = max(min([yN2; yO2]) - ypad, 0);
        yhi = min(max([yN2; yO2]) + ypad, 1);
        segments = segmentIndices(stepIdx);
        for s = 1:size(segments, 1)
            idxRange = segments(s, :);
            colour = colours(stepIdx(idxRange(1)), :);
            patch(ax, [x(idxRange(1)), x(idxRange(2)), x(idxRange(2)), x(idxRange(1))], ...
                [ylo, ylo, yhi, yhi], colour, 'FaceAlpha', 0.05, 'EdgeColor', 'none');
        end
        productMask = arrayfun(@(idx) schedule(idx).product_stream, stepIdx);
        productSegments = contiguousSegments(productMask);
        for p = 1:size(productSegments, 1)
            idxRange = productSegments(p, :);
            patch(ax, [x(idxRange(1)), x(idxRange(2)), x(idxRange(2)), x(idxRange(1))], ...
                [ylo, ylo, yhi, yhi], [0.05 0.5 0.2], 'FaceAlpha', 0.12, 'EdgeColor', 'none');
        end
        plot(ax, x, yN2, 'LineWidth', 1.8, 'Color', [0 0.45 0.74]);
        plot(ax, x, yO2, '--', 'LineWidth', 1.2, 'Color', [0.85 0.33 0.1]);
        legend(ax, {'y_{N_2}', 'y_{O_2}'}, 'Location', 'best');
        ax.Box = 'on';
        ax.XGrid = 'on';
        ax.YGrid = 'on';
        ax.YLim = [ylo, yhi];
        hold(ax, 'off');
    end

    function seg = segmentIndices(stepIdx)
        changeIdx = [1; find(diff(stepIdx)) + 1; numel(stepIdx) + 1];
        seg = [changeIdx(1:end-1), changeIdx(2:end) - 1];
    end

    function seg = contiguousSegments(mask)
        mask = mask(:) > 0;
        changeIdx = [1; find(diff(mask)) + 1; numel(mask) + 1];
        seg = zeros(0, 2);
        for k = 1:numel(changeIdx) - 1
            if mask(changeIdx(k))
                seg(end+1, :) = [changeIdx(k), changeIdx(k+1) - 1]; %#ok<AGROW>
            end
        end
    end

    function updateLegend(schedule, colours)
        delete(app.legendGrid.Children);
        nSteps = numel(schedule);
        app.legendGrid.RowHeight = repmat({'fit'}, 1, nSteps);
        for k = 1:nSteps
            lbl = uilabel(app.legendGrid, 'Text', schedule(k).name, ...
                'FontWeight', 'bold', 'HorizontalAlignment', 'left');
            lbl.Layout.Row = k;
            lbl.Layout.Column = 1;
            lbl.BackgroundColor = colours(k, :);
            lbl.FontColor = contrastColour(colours(k, :));
            lbl.Padding = [6 3 6 3];
        end
    end

    function c = contrastColour(baseColour)
        if numel(baseColour) ~= 3
            c = [1 1 1];
            return;
        end
        brightness = 0.299 * baseColour(1) + 0.587 * baseColour(2) + 0.114 * baseColour(3);
        if brightness < 0.5
            c = [1 1 1];
        else
            c = [0 0 0];
        end
    end

    function updateStatus(msg)
        app.outputs.statusLabel.Text = msg;
    end

end
