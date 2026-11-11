function [progress, fig] = parforWorkerProgress(nIterations)
%PARFORWORKERPROGRESS Display per-worker progress for a parfor loop.
%
% Usage:
%   N = 200;
%   [progress, fig] = parforWorkerProgress(N);
%
%   parfor i = 1:N
%       % Your work here
%       pause(0.05);
%       progress.step(1);
%   end
%
%   close(fig)
%
% Each worker bar represents that worker's completed iterations as a
% fraction of the total number of loop iterations.

    arguments
        nIterations (1,1) {mustBePositive, mustBeInteger}
    end

    pool = gcp();
    nWorkers = pool.NumWorkers;

    % Queue for progress messages sent from workers to the MATLAB client.
    queue = parallel.pool.DataQueue;

    % Creates one stable identifier for each worker.
    workerID = parallel.pool.Constant(@localWorkerID);

    % ---- Create one display containing stacked bars --------------------
    fig = figure( ...
        'Name', 'parfor worker progress', ...
        'NumberTitle', 'off', ...
        'Color', 'w', ...
        'Position', [300 300 720 max(180, 100 + 55*nWorkers)]);

    ax = axes(fig, ...
        'Position', [0.20 0.13 0.73 0.76], ...
        'Box', 'on', ...
        'XLim', [0 1], ...
        'YLim', [0.5 nWorkers + 0.5], ...
        'YDir', 'reverse', ...
        'XTick', 0:0.1:1, ...
        'XTickLabel', compose('%d%%', 0:10:100), ...
        'YTick', 1:nWorkers, ...
        'YTickLabel', compose('Worker %d', 1:nWorkers));

    xlabel(ax, 'Share of total iterations completed');
    title(ax, sprintf('parfor progress: 0 / %d iterations', nIterations));

    % Background and foreground rectangles for each worker.
    backgroundBars = gobjects(nWorkers, 1);
    progressBars = gobjects(nWorkers, 1);
    labels = gobjects(nWorkers, 1);

    for k = 1:nWorkers
        backgroundBars(k) = rectangle(ax, ...
            'Position', [0, k - 0.32, 1, 0.64], ...
            'FaceColor', [0.92 0.92 0.92], ...
            'EdgeColor', [0.75 0.75 0.75]);

        progressBars(k) = rectangle(ax, ...
            'Position', [0, k - 0.32, 0, 0.64], ...
            'FaceColor', [0.15 0.55 0.85], ...
            'EdgeColor', 'none');

        labels(k) = text(ax, 0.01, k, sprintf('0 / %d', nIterations), ...
            'VerticalAlignment', 'middle', ...
            'HorizontalAlignment', 'left', ...
            'Color', 'k', ...
            'FontWeight', 'bold');
    end

    % State is maintained on the MATLAB client.
    completedByWorker = zeros(nWorkers, 1);
    knownWorkerIDs = nan(nWorkers, 1);
    nextFreeBar = 1;

    % Run updateProgress on the client for every queue message.
    afterEach(queue, @updateProgress);

    % This structure is safe to use inside parfor.
    progress.queue = queue;
    progress.workerID = workerID;
    progress.step = @(amount) sendProgress(queue, workerID, amount);

    function updateProgress(message)
        % message = [workerID, numberOfIterationsCompleted]
        if ~isgraphics(fig)
            return
        end

        id = message(1);
        amount = message(2);

        % Associate the worker ID with one of the displayed rows.
        barIndex = find(knownWorkerIDs == id, 1);

        if isempty(barIndex)
            barIndex = nextFreeBar;
            knownWorkerIDs(barIndex) = id;
            nextFreeBar = nextFreeBar + 1;
        end

        completedByWorker(barIndex) = completedByWorker(barIndex) + amount;

        fraction = min(completedByWorker(barIndex) / nIterations, 1);

        % Resize this worker's foreground bar.
        pos = progressBars(barIndex).Position;
        pos(3) = fraction;
        progressBars(barIndex).Position = pos;

        % Move the label just after the bar where possible.
        labelX = min(fraction + 0.015, 0.88);
        labels(barIndex).Position(1) = labelX;
        labels(barIndex).String = sprintf('%d / %d', ...
            completedByWorker(barIndex), nIterations);

        % Update overall progress in the title.
        totalCompleted = sum(completedByWorker);
        title(ax, sprintf('parfor progress: %d / %d iterations', ...
            totalCompleted, nIterations));

        drawnow limitrate
    end
end

function sendProgress(queue, workerIDConstant, amount)
%SENDPROGRESS Send one progress update from a parfor worker.
    send(queue, [workerIDConstant.Value, amount]);
end

function id = localWorkerID()
%LOCALWORKERID Return the identifier of the current parallel worker.
    task = getCurrentTask();

    if isempty(task)
        id = 1;
    else
        id = task.ID;
    end
end