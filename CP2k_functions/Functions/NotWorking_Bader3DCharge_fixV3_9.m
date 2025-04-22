function Bader3DCharge_fixV3_9(TrajXYZ, TrajABC, TrajQmc, varargin)
% Bader3DCharge(TrajXYZ, TrajABC, TrajQmc, SelectedIndices, ...)
% Plots 3D charge distribution with complete simulation box and view control
%
% Parameters:
%   'PlayVideo',      false      - Auto-play animation
%   'SaveVideo',      false      - Export as MP4 video
%   'FrameRate',      24         - Frames per second
%   'OutputFileName', 'output.mp4' - Video filename
%   'Interactive',    true       - Show UI controls
%   'ViewAxis',       'X'        - Main axis facing viewer ('X' or 'Y')
%   'VideoView',      []         - Custom view angles for video [azimuth,elevation]

% Input parsing
p = inputParser;
addRequired(p, 'TrajXYZ', @(x) isnumeric(x) || iscell(x));
addRequired(p, 'TrajABC', @(x) isnumeric(x) || iscell(x));
addRequired(p, 'TrajQmc', @(x) isnumeric(x) || iscell(x));
addOptional(p, 'SelectedIndices', [], @(x) isempty(x) || isvector(x));
addParameter(p, 'PlayVideo', false);
addParameter(p, 'SaveVideo', false);
addParameter(p, 'FrameRate', 24);
addParameter(p, 'OutputFileName', 'charge_video.mp4');
addParameter(p, 'Interactive', true);
addParameter(p, 'ViewAxis', 'X', @(x) ismember(x,{'X','Y'}));
addParameter(p, 'VideoView', [], @isnumeric);
addParameter(p, 'LoopVideo', false, @islogical); % Add loop parameter
parse(p, TrajXYZ, TrajABC, TrajQmc, varargin{:});


params = p.Results;
if params.PlayVideo
    params.Interactive = false;  % Disable UI controls for auto-play
end
SelectedIndices = params.SelectedIndices;

% Convert inputs to cell arrays
% Convert inputs to cell arrays FIRST
[TrajXYZ, TrajABC, TrajQmc] = convert_to_cells(TrajXYZ, TrajABC, TrajQmc); % Fixed line
% [TrajXYZ, TrajABC, TrajQmc] = convert_to_cells(TrajXYZ, TrajABC, TrajQmc);
numFrames = numel(TrajXYZ);
Natoms = size(TrajXYZ{1}, 1);

% Set default selected indices
if isempty(SelectedIndices)
    SelectedIndices = 1:Natoms;
end

% Initialize figure and axes
fig = figure('Color','w', 'Name','3D Charge Viewer',...
    'Position',[100 100 1200 800], 'NumberTitle','off');
ax = axes('Parent', fig, 'Position',[0.1 0.1 0.7 0.8]);
hold(ax, 'on');
axis(ax, 'equal');
grid(ax, 'on');

% Set axis labels
xlabel(ax, 'X-axis');
ylabel(ax, 'Y-axis');
zlabel(ax, 'Z-axis');

% Set initial view
if strcmpi(params.ViewAxis, 'X')
    view(ax, [0 0]); % X facing viewer
else
    view(ax, [90 0]); % Y facing viewer
end
set(ax, 'CameraUpVector', [0 0 1]); % Keep Z up

colorbar(ax); % Add colorbar

% Sphere template
[xs,ys,zs] = sphere(8);
Radius = 1;

% Video writer initialization
if params.SaveVideo
    [outputDir, ~, ~] = fileparts(params.OutputFileName);
    if ~isempty(outputDir) && ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end
    
    writerObj = VideoWriter(params.OutputFileName, 'MPEG-4');
    writerObj.FrameRate = params.FrameRate;
    open(writerObj);
    
    % Store original view if using custom video view
    originalView = ax.View;
    if ~isempty(params.VideoView)
        view(ax, params.VideoView);
    end
end

% UI controls (positioned at right side)
if params.Interactive
    uiX = 0; % Left-aligned controls

    uicontrol('Style','togglebutton', 'String','Play',...
        'Position',[uiX 700 80 30], 'Tag','playBtn',...
        'Callback',@(src,~) playCallback(src));

    uicontrol('Style','slider',...
        'Position',[uiX 750 250 30],...
        'Min',1, 'Max',numFrames, 'Value',1,...
        'SliderStep',[1/numFrames 10/numFrames],...
        'Tag','frameSlider',...
        'Callback',@(src,~) sliderCallback(src));

    % uicontrol('Style','text', 'String','Frame:',...
    %     'Position',[uiX 670 100 20]);
end

% Initialize animation variables
currentFrame = 1;
hAtoms = [];
hVectors = [];

% Initial plot
[hAtoms, hVectors] = updatePlot(ax, TrajXYZ{1}, TrajABC{1}, TrajQmc{1},...
    xs, ys, zs, Radius, SelectedIndices, currentFrame, numFrames, hAtoms, hVectors);

% Video export routine
if params.SaveVideo
    for frame = 1:numFrames
        [hAtoms, hVectors] = updatePlot(ax, TrajXYZ{frame}, TrajABC{frame},...
            TrajQmc{frame}, xs, ys, zs, Radius, SelectedIndices, frame, numFrames, hAtoms, hVectors);
        writeVideo(writerObj, getframe(fig));
        drawnow('limitrate');
        if params.Interactive
            set(findobj(fig,'Tag','frameSlider'), 'Value', frame);
        end
    end
    close(writerObj);
    
    % Restore original view if changed
    if ~isempty(params.VideoView)
        view(ax, originalView);
    end
end

% Animation controls
% In the video playback section, modify the updatePlot call to include ALL parameters:
if params.PlayVideo
    writerObj = VideoWriter(params.OutputFileName, 'MPEG-4');
    writerObj.FrameRate = params.FrameRate;
    open(writerObj);
    
    
    % Initialize handles
    hAtoms = [];
    hVectors = [];
    
    % Get all required parameters
    xs = xs; % From sphere template
    ys = ys;
    zs = zs;
    Radius = 1; % Your sphere radius
    
    try
        frameCount = 0;
        while true % Infinite loop for continuous playback
            for frame = 1:numFrames
                % PROPERLY PASS ALL ARGUMENTS TO updatePlot
                [hAtoms, hVectors] = updatePlot(...
                    ax,...
                    TrajXYZ{frame},...
                    TrajABC{frame},...
                    TrajQmc{frame},...
                    xs, ys, zs,...
                    Radius,...
                    SelectedIndices,...
                    frame,...
                    numFrames,...
                    hAtoms,...
                    hVectors...
                );
                
                writeVideo(writerObj, getframe(fig));
                drawnow('limitrate');
                
                % Break loop if not in continuous mode
                if ~params.LoopVideo && frame == numFrames
                    break;
                end
            end
            
            % Exit loop if not in continuous mode
            if ~params.LoopVideo
                break;
            end
        end
    catch ME
        close(writerObj);
        rethrow(ME);
    end
    close(writerObj);
end

% Callback functions
% Modified playCallback function
    function playCallback(src)
        persistent isPlaying
        if isempty(isPlaying), isPlaying = false; end
        
        % Toggle state
        isPlaying = ~isPlaying;
        set(src, 'Value', double(isPlaying));
        
        while loopActive && ishandle(src)
            % Update frame number with wrap-around
            if currentFrame < numFrames
                currentFrame = currentFrame + 1;
            else
                if params.LoopVideo
                    currentFrame = 1;
                else
                    loopActive = false;
                end
            end
            
            % Update visualization
            [hAtoms, hVectors] = updatePlot(ax, TrajXYZ{currentFrame},...
                TrajABC{currentFrame}, TrajQmc{currentFrame},...
                xs, ys, zs, Radius, SelectedIndices, currentFrame, numFrames, hAtoms, hVectors);
            
            % Update slider position
            set(findobj(fig,'Tag','frameSlider'), 'Value', currentFrame);
            
            % Process ALL pending events
            drawnow;
            
            % Check actual button state
            currentButtonState = get(src,'Value');
            if ~currentButtonState
                loopActive = false;
                break;
            end
            
            % Timing control with event processing
            tStart = tic;
            while toc(tStart) < 1/params.FrameRate && loopActive
                drawnow;
                % Check for stop condition
            if ~get(src, 'Value')
                isPlaying = false;
                break;
            end
            
            % Add small pause for UI responsiveness
            pause(0.01);
        end
        set(src, 'Value', 0);
    end
end

    function sliderCallback(src)
        frame = round(get(src,'Value'));
        if frame ~= currentFrame
            currentFrame = frame;
            [hAtoms, hVectors] = updatePlot(ax, TrajXYZ{currentFrame},...
                TrajABC{currentFrame}, TrajQmc{currentFrame},...
                xs, ys, zs, Radius, SelectedIndices, currentFrame, numFrames, hAtoms, hVectors);
        end
    end

% Modified auto-play for continuous support
    function autoPlayAnimation()
        playCount = 0;
        while params.PlayVideo && (params.LoopVideo || playCount == 0)
            for frame = 1:numFrames
                % Break if interrupted
                if ~params.PlayVideo && ~params.LoopVideo
                    break;
                end
                
                % Update visualization
                [hAtoms, hVectors] = updatePlot(ax, TrajXYZ{frame},...
                    TrajABC{frame}, TrajQmc{frame},...
                    xs, ys, zs, Radius, SelectedIndices, frame, numFrames, hAtoms, hVectors);
                
                % Update UI elements
                if params.Interactive
                    set(findobj(fig,'Tag','frameSlider'), 'Value', frame);
                end
                
                % Maintain frame rate
                pause(1/params.FrameRate);
            end
            playCount = playCount + 1;
        end
    end

hold(ax, 'off');
end

%% Helper functions
function [TrajXYZ, TrajABC, TrajQmc] = convert_to_cells(TrajXYZ, TrajABC, TrajQmc)
% Input conversion and validation (same as previous versions)
% ... (identical content from previous implementation) ...
% numFrames = max([size(TrajXYZ, 1), 1]);
% Proper cell conversion for all trajectory inputs

% Convert XYZ trajectories
if ~iscell(TrajXYZ)
    if ndims(TrajXYZ) == 3
        % Correct 3D array handling: [frames x atoms x coordinates]
        TrajXYZ = arrayfun(@(i) squeeze(TrajXYZ(i,:,:)), 1:size(TrajXYZ,1), 'UniformOutput', false)';
    else
        % Wrap 2D array in cell
        TrajXYZ = {TrajXYZ};
    end
end

% Convert ABC trajectories (same logic)
if ~iscell(TrajABC)
    if ndims(TrajABC) == 3
        TrajABC = arrayfun(@(i) squeeze(TrajABC(i,:,:)), 1:size(TrajABC,1), 'UniformOutput', false)';
    else
        TrajABC = {TrajABC};
    end
end

% Convert Qmc trajectories 
if ~iscell(TrajQmc)
    if ndims(TrajQmc) == 3
        TrajQmc = arrayfun(@(i) squeeze(TrajQmc(i,:)), 1:size(TrajQmc,1), 'UniformOutput', false)';
    else
        TrajQmc = {TrajQmc};
    end
end

% Validate frame counts
numFrames = numel(TrajXYZ);
if numel(TrajABC) == 1 && numFrames > 1
    TrajABC = repmat(TrajABC, numFrames, 1);
end
if numel(TrajQmc) == 1 && numFrames > 1
    TrajQmc = repmat(TrajQmc, numFrames, 1);
end

assert(numel(TrajABC) == numFrames, 'ABC frame count mismatch');
assert(numel(TrajQmc) == numFrames, 'Qmc frame count mismatch');
end
% end

function [hAtoms, hVectors] = updatePlot(ax, XYZ, ABC, Qmc, xs, ys, zs, Radius, SelectedIndices, currentFrame, numFrames, hAtoms, hVectors)
% Update plot elements with complete simulation box

% Delete old frame text
delete(findobj(ax,'Type','text','Tag','FrameText'));

% Create color map
[cmap, minQ, maxQ] = createColormap(Qmc(SelectedIndices));
colormap(ax, cmap);
caxis(ax, [minQ, maxQ]);

% Update atoms
if isempty(hAtoms) || ~isvalid(hAtoms(1))
    hAtoms = createAtoms(ax, XYZ, xs, ys, zs, Radius, Qmc, SelectedIndices, cmap, minQ, maxQ);
else
    updateAtoms(hAtoms, XYZ, xs, ys, zs, Radius, Qmc, SelectedIndices, cmap, minQ, maxQ);
end

% Update unit cell vectors
Vec = diag(ABC);
if isempty(hVectors) || ~isvalid(hVectors(1))
    hVectors = createVectors(ax, Vec);
else
    updateVectors(hVectors, Vec);
end

% Add frame counter
text(ax, 0.05, 0.95, 0.9, sprintf('Frame: %d/%d', currentFrame, numFrames),...
    'Units','normalized', 'Tag','FrameText');
end

function hAtoms = createAtoms(ax, XYZ, xs, ys, zs, Radius, Qmc, SelectedIndices, cmap, minQ, maxQ)
% Create new atom surfaces
hAtoms = gobjects(size(XYZ,1), 1);
for i = 1:size(XYZ,1)
    color = [0.7 0.7 0.7];
    alpha = 0.3;
    
    if ismember(i, SelectedIndices)
        if maxQ ~= minQ
            cidx = round((Qmc(i) - minQ)/(maxQ - minQ) * (size(cmap,1)-1)) + 1;
        else
            cidx = 1;
        end
        color = cmap(cidx,:);
        alpha = 1;
    end
    
    xAt = xs*Radius + XYZ(i,1);
    yAt = ys*Radius + XYZ(i,2);
    zAt = zs*Radius + XYZ(i,3);
    
    hAtoms(i) = surf(ax, xAt, yAt, zAt,...
        'FaceColor', color,...
        'EdgeColor', 'none',...
        'FaceAlpha', alpha);
end
end

function updateAtoms(hAtoms, XYZ, xs, ys, zs, Radius, Qmc, SelectedIndices, cmap, minQ, maxQ)
% Update existing atoms
for i = 1:size(XYZ,1)
    color = [0.7 0.7 0.7];
    alpha = 0.3;
    if ismember(i, SelectedIndices)
        if maxQ ~= minQ
            cidx = round((Qmc(i) - minQ)/(maxQ - minQ) * (size(cmap,1)-1)) + 1;
        else
            cidx = 1;
        end
        color = cmap(cidx,:);
        alpha = 1;
    end
    xAt = xs*Radius + XYZ(i,1);
    yAt = ys*Radius + XYZ(i,2);
    zAt = zs*Radius + XYZ(i,3);
    set(hAtoms(i), 'XData', xAt, 'YData', yAt, 'ZData', zAt,...
        'FaceColor', color, 'FaceAlpha', alpha);
end
end

function hVectors = createVectors(ax, Vec)
% Create all 12 simulation box edges
hVectors = gobjects(12,1);

% Base vectors
hVectors(1) = plot3(ax, [0 Vec(1,1)],[0 Vec(1,2)],[0 Vec(1,3)],'-k','LineWidth',2);
hVectors(2) = plot3(ax, [0 Vec(2,1)],[0 Vec(2,2)],[0 Vec(2,3)],'-k','LineWidth',2);
hVectors(3) = plot3(ax, [0 Vec(3,1)],[0 Vec(3,2)],[0 Vec(3,3)],'-k','LineWidth',2);

% Bottom connections
hVectors(4) = plot3(ax, [Vec(1,1) Vec(1,1)+Vec(2,1)],...
                     [Vec(1,2) Vec(1,2)+Vec(2,2)],...
                     [Vec(1,3) Vec(1,3)+Vec(2,3)],'-k','LineWidth',2);
hVectors(5) = plot3(ax, [Vec(2,1) Vec(1,1)+Vec(2,1)],...
                     [Vec(2,2) Vec(1,2)+Vec(2,2)],...
                     [Vec(2,3) Vec(1,3)+Vec(2,3)],'-k','LineWidth',2);

% Vertical connections
hVectors(6) = plot3(ax, [Vec(1,1) Vec(1,1)+Vec(3,1)],...
                     [Vec(1,2) Vec(1,2)+Vec(3,2)],...
                     [Vec(1,3) Vec(1,3)+Vec(3,3)],'-k','LineWidth',2);
hVectors(7) = plot3(ax, [Vec(2,1) Vec(2,1)+Vec(3,1)],...
                     [Vec(2,2) Vec(2,2)+Vec(3,2)],...
                     [Vec(2,3) Vec(2,3)+Vec(3,3)],'-k','LineWidth',2);
hVectors(8) = plot3(ax, [Vec(1,1)+Vec(2,1) Vec(1,1)+Vec(2,1)+Vec(3,1)],...
                     [Vec(1,2)+Vec(2,2) Vec(1,2)+Vec(2,2)+Vec(3,2)],...
                     [Vec(1,3)+Vec(2,3) Vec(1,3)+Vec(2,3)+Vec(3,3)],'-k','LineWidth',2);

% Top connections
hVectors(9) = plot3(ax, [Vec(3,1) Vec(3,1)+Vec(1,1)],...
                     [Vec(3,2) Vec(3,2)+Vec(1,2)],...
                     [Vec(3,3) Vec(3,3)+Vec(1,3)],'-k','LineWidth',2);
hVectors(10) = plot3(ax, [Vec(3,1) Vec(3,1)+Vec(2,1)],...
                      [Vec(3,2) Vec(3,2)+Vec(2,2)],...
                      [Vec(3,3) Vec(3,3)+Vec(2,3)],'-k','LineWidth',2);
hVectors(11) = plot3(ax, [Vec(3,1)+Vec(1,1) Vec(3,1)+Vec(1,1)+Vec(2,1)],...
                      [Vec(3,2)+Vec(1,2) Vec(3,2)+Vec(1,2)+Vec(2,2)],...
                      [Vec(3,3)+Vec(1,3) Vec(3,3)+Vec(1,3)+Vec(2,3)],'-k','LineWidth',2);
hVectors(12) = plot3(ax, [Vec(3,1)+Vec(2,1) Vec(3,1)+Vec(2,1)+Vec(1,1)],...
                      [Vec(3,2)+Vec(2,2) Vec(3,2)+Vec(2,2)+Vec(1,2)],...
                      [Vec(3,3)+Vec(2,3) Vec(3,3)+Vec(2,3)+Vec(1,3)],'-k','LineWidth',2);
end

function updateVectors(hVectors, Vec)
% Update all 12 simulation box edges
set(hVectors(1), 'XData', [0 Vec(1,1)], 'YData', [0 Vec(1,2)], 'ZData', [0 Vec(1,3)]);
set(hVectors(2), 'XData', [0 Vec(2,1)], 'YData', [0 Vec(2,2)], 'ZData', [0 Vec(2,3)]);
set(hVectors(3), 'XData', [0 Vec(3,1)], 'YData', [0 Vec(3,2)], 'ZData', [0 Vec(3,3)]);

set(hVectors(4), 'XData', [Vec(1,1) Vec(1,1)+Vec(2,1)],...
                'YData', [Vec(1,2) Vec(1,2)+Vec(2,2)],...
                'ZData', [Vec(1,3) Vec(1,3)+Vec(2,3)]);
set(hVectors(5), 'XData', [Vec(2,1) Vec(1,1)+Vec(2,1)],...
                'YData', [Vec(2,2) Vec(1,2)+Vec(2,2)],...
                'ZData', [Vec(2,3) Vec(1,3)+Vec(2,3)]);

set(hVectors(6), 'XData', [Vec(1,1) Vec(1,1)+Vec(3,1)],...
                'YData', [Vec(1,2) Vec(1,2)+Vec(3,2)],...
                'ZData', [Vec(1,3) Vec(1,3)+Vec(3,3)]);
set(hVectors(7), 'XData', [Vec(2,1) Vec(2,1)+Vec(3,1)],...
                'YData', [Vec(2,2) Vec(2,2)+Vec(3,2)],...
                'ZData', [Vec(2,3) Vec(2,3)+Vec(3,3)]);
set(hVectors(8), 'XData', [Vec(1,1)+Vec(2,1) Vec(1,1)+Vec(2,1)+Vec(3,1)],...
                'YData', [Vec(1,2)+Vec(2,2) Vec(1,2)+Vec(2,2)+Vec(3,2)],...
                'ZData', [Vec(1,3)+Vec(2,3) Vec(1,3)+Vec(2,3)+Vec(3,3)]);

set(hVectors(9), 'XData', [Vec(3,1) Vec(3,1)+Vec(1,1)],...
                'YData', [Vec(3,2) Vec(3,2)+Vec(1,2)],...
                'ZData', [Vec(3,3) Vec(3,3)+Vec(1,3)]);
set(hVectors(10), 'XData', [Vec(3,1) Vec(3,1)+Vec(2,1)],...
                 'YData', [Vec(3,2) Vec(3,2)+Vec(2,2)],...
                 'ZData', [Vec(3,3) Vec(3,3)+Vec(2,3)]);
set(hVectors(11), 'XData', [Vec(3,1)+Vec(1,1) Vec(3,1)+Vec(1,1)+Vec(2,1)],...
                 'YData', [Vec(3,2)+Vec(1,2) Vec(3,2)+Vec(1,2)+Vec(2,2)],...
                 'ZData', [Vec(3,3)+Vec(1,3) Vec(3,3)+Vec(1,3)+Vec(2,3)]);
set(hVectors(12), 'XData', [Vec(3,1)+Vec(2,1) Vec(3,1)+Vec(2,1)+Vec(1,1)],...
                 'YData', [Vec(3,2)+Vec(2,2) Vec(3,2)+Vec(2,2)+Vec(1,2)],...
                 'ZData', [Vec(3,3)+Vec(2,3) Vec(3,3)+Vec(2,3)+Vec(1,3)]);
end

function [cmap, minQ, maxQ] = createColormap(Q)
% Create signed colormap
minQ = min(Q);
maxQ = max(Q);

if minQ == maxQ
    cmap = [0 0 1];
    return
end

xi = linspace(minQ, maxQ, 256);
if any(Q < 0) && any(Q > 0)
    cmap = interp1([minQ, 0, maxQ], [1 0 0; 0 0.8 0; 0 0 1], xi);
elseif all(Q >= 0)
    cmap = interp1([0, maxQ], [0 0.8 0; 0 0 1], xi);
else
    cmap = interp1([minQ, 0], [1 0 0; 0 0.8 0], xi);
end
end