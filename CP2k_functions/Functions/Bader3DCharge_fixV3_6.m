function Bader3DCharge_fixV3_6(TrajXYZ, TrajABC, TrajQmc, varargin)
% Bader3DCharge(TrajXYZ, TrajABC, TrajQmc, SelectedIndices, ...)
% Plots charge distribution for trajectories with interactive controls
% 
% Optional parameters (name-value pairs):
%   'PlayVideo',      false      - Auto-play animation
%   'SaveVideo',      false      - Export as MP4 video
%   'FrameRate',      24         - Frames per second
%   'OutputFileName', 'output.mp4' - Video filename
%   'Interactive',    true       - Show UI controls

% Bader3DCharge(TrajXYZ, TrajABC, TrajQmc, [SelectedIndices], ...)
% Features:
% - 3D charge visualization with sign-aware colormap
% - Trajectory animation support
% - Interactive controls (play/pause, frame slider)
% - Video export capability
% - Complete unit cell vector plotting
% - Context atoms with faint visualization

%% For trajectory visualization with video export
% Bader3DCharge(trajXYZ, trajABC, trajQmc, [1,3,5],...
%     'SaveVideo', true,...
%     'OutputFileName', 'results/charge_animation.mp4',...
%     'FrameRate', 30);
%%

%% V3.6 usage example
%Bader3DCharge_fixV3_6(..., 'SaveVideo', true,...
    % 'VideoView', [45 30]);  % Azimuth=45°, Elevation=30°
%% 


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
% [...] (keep previous input parsing section, add these parameters)
% Modified input parser section:
addParameter(p, 'ViewAxis', 'X', @(x) ismember(x,{'X','Y'}));  % New parameter
addParameter(p, 'VideoView', [], @isnumeric);  % Custom view for video
% [...] (previous code until axes initialization)
parse(p, TrajXYZ, TrajABC, TrajQmc, varargin{:});

params = p.Results;
SelectedIndices = params.SelectedIndices;

% Convert inputs to cell arrays
[TrajXYZ, TrajABC, TrajQmc] = convert_to_cells(TrajXYZ, TrajABC, TrajQmc);
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
% Set axis labels BEFORE setting view
xlabel(ax, 'X-axis');
ylabel(ax, 'Y-axis');
zlabel(ax, 'Z-axis');

% Configure view based on parameters
if strcmpi(params.ViewAxis, 'X')
    view(ax, [0 0]);  % X towards viewer
else
    view(ax, [90 0]); % Y towards viewer
end
set(ax, 'CameraUpVector', [0 0 1]); % Keep Z-axis up

% [...] (rest of previous code until video writer section)

view(ax, [0 0]); % Set view with X towards viewer, Y right, Z up
set(ax, 'CameraUpVector', [0 0 1]); % Ensure Z-axis is up
colorbar(ax); % Add colorbar

% Sphere template
[xs,ys,zs] = sphere(8);
Radius = 1;

% Video writer initialization
% Modified video recording section:
if params.SaveVideo
    % Store original view if using custom video view
    if ~isempty(params.VideoView)
        originalView = ax.View;
        view(ax, params.VideoView);
    end
    
    % [...] (previous video setup code)
    [outputDir, ~, ~] = fileparts(params.OutputFileName);
    if ~isempty(outputDir) && ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end
    
    writerObj = VideoWriter(params.OutputFileName, 'MPEG-4');
    writerObj.FrameRate = params.FrameRate;
    open(writerObj);
end

% UI controls
if params.Interactive
    uicontrol('Style','togglebutton', 'String','Play',...
        'Position',[820 750 100 30], 'Tag','playBtn',...
        'Callback',@(src,~) playCallback(src));
    
    uicontrol('Style','slider',...
        'Position',[820 700 350 30],...
        'Min',1, 'Max',numFrames, 'Value',1,...
        'SliderStep',[1/numFrames 10/numFrames],...
        'Tag','frameSlider',...
        'Callback',@(src,~) sliderCallback(src));
    
    uicontrol('Style','text', 'String','Frame:',...
        'Position',[820 670 100 20]);
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
    % Store original view if using custom video view
    if ~isempty(params.VideoView)
        originalView = ax.View;
        view(ax, params.VideoView);
    end

    % [...] (previous video setup code)
    
    for frame = 1:numFrames
        [hAtoms, hVectors] = updatePlot(ax, TrajXYZ{frame}, TrajABC{frame},...
            TrajQmc{frame}, xs, ys, zs, Radius, SelectedIndices, frame, numFrames, hAtoms, hVectors);
        % writeVideo(writerObj, getframe(fig));
        writeVideo(writerObj, getframe(fig));
        drawnow('limitrate');  % Improved rendering speed
        if params.Interactive
            set(findobj(fig,'Tag','frameSlider'), 'Value', frame);
        end
    end
    % Restore original view if changed
    if ~isempty(params.VideoView)
        view(ax, originalView);
    end
    % [...] (rest of video code)
    close(writerObj);
end

% Animation controls
if params.PlayVideo && params.Interactive
    autoPlayAnimation();
end

% Callback functions
    function playCallback(src)
        isPlaying = get(src,'Value');
        while isPlaying && currentFrame < numFrames
            currentFrame = currentFrame + 1;
            [hAtoms, hVectors] = updatePlot(ax, TrajXYZ{currentFrame},...
                TrajABC{currentFrame}, TrajQmc{currentFrame},...
                xs, ys, zs, Radius, SelectedIndices, currentFrame, numFrames, hAtoms, hVectors);
            set(findobj(fig,'Tag','frameSlider'), 'Value', currentFrame);
            pause(1/params.FrameRate);
            isPlaying = get(src,'Value');
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

    function autoPlayAnimation()
        for frame = 1:numFrames
            [hAtoms, hVectors] = updatePlot(ax, TrajXYZ{frame},...
                TrajABC{frame}, TrajQmc{frame},...
                xs, ys, zs, Radius, SelectedIndices, frame, numFrames, hAtoms, hVectors);
            if params.Interactive
                set(findobj(fig,'Tag','frameSlider'), 'Value', frame);
            end
            pause(1/params.FrameRate);
        end
    end

hold(ax, 'off');
end

%% Helper functions
function [TrajXYZ, TrajABC, TrajQmc] = convert_to_cells(TrajXYZ, TrajABC, TrajQmc)
% Handle input conversions and ABC replication
numFrames = max([size(TrajXYZ, 1), 1]);

% Convert XYZ
if ~iscell(TrajXYZ)
    if ndims(TrajXYZ) == 3
        TrajXYZ = arrayfun(@(i) squeeze(TrajXYZ(i,:,:)), 1:size(TrajXYZ,1), 'UniformOutput', false)';
    else
        TrajXYZ = {TrajXYZ};
    end
end
numFrames = numel(TrajXYZ);

% Convert ABC with replication
if ~iscell(TrajABC)
    if ndims(TrajABC) == 3
        TrajABC = arrayfun(@(i) squeeze(TrajABC(i,:,:)), 1:size(TrajABC,1), 'UniformOutput', false)';
    else
        TrajABC = {TrajABC};
    end
end
if numel(TrajABC) == 1 && numFrames > 1
    TrajABC = repmat(TrajABC, numFrames, 1);
end

% Convert Qmc
if ~iscell(TrajQmc)
    if ndims(TrajQmc) == 3
        TrajQmc = arrayfun(@(i) squeeze(TrajQmc(i,:)), 1:size(TrajQmc,1), 'UniformOutput', false)';
    else
        TrajQmc = {TrajQmc};
    end
end
if numel(TrajQmc) == 1 && numFrames > 1
    TrajQmc = repmat(TrajQmc, numFrames, 1);
end

% Validate frame counts
assert(numel(TrajABC) == numFrames, 'ABC frame count mismatch');
assert(numel(TrajQmc) == numFrames, 'Qmc frame count mismatch');
end

function [hAtoms, hVectors] = updatePlot(ax, XYZ, ABC, Qmc, xs, ys, zs, Radius, SelectedIndices, currentFrame, numFrames, hAtoms, hVectors)
% Update plot elements by reusing existing graphics objects

% Create color map
[cmap, minQ, maxQ] = createColormap(Qmc(SelectedIndices));
colormap(ax, cmap);
caxis(ax, [minQ, maxQ]);

% Check if handles are provided and valid
if isempty(hAtoms) || ~isvalid(hAtoms(1))
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
        hAtoms(i) = surf(ax, xAt, yAt, zAt, 'FaceColor', color, 'EdgeColor', 'none', 'FaceAlpha', alpha);
    end
else
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
        set(hAtoms(i), 'XData', xAt, 'YData', yAt, 'ZData', zAt, 'FaceColor', color, 'FaceAlpha', alpha);
    end
end

% Handle vectors
Vec = diag(ABC);
if isempty(hVectors) || ~isvalid(hVectors(1))
    % Create new vectors
    hVectors = gobjects(0);
    % Base vectors
    hVectors(end+1) = plot3(ax, [0 Vec(1,1)], [0 Vec(1,2)], [0 Vec(1,3)], '-k', 'LineWidth', 2);
    hVectors(end+1) = plot3(ax, [0 Vec(2,1)], [0 Vec(2,2)], [0 Vec(2,3)], '-k', 'LineWidth', 2);
    hVectors(end+1) = plot3(ax, [0 Vec(3,1)], [0 Vec(3,2)], [0 Vec(3,3)], '-k', 'LineWidth', 2);
    % Additional edges
    hVectors(end+1) = plot3(ax, Vec(1,1)+[0 Vec(3,1)], Vec(1,2)+[0 Vec(3,2)], Vec(1,3)+[0 Vec(3,3)], '-k', 'LineWidth', 2);
    hVectors(end+1) = plot3(ax, Vec(2,1)+[0 Vec(3,1)], Vec(2,2)+[0 Vec(3,2)], Vec(2,3)+[0 Vec(3,3)], '-k', 'LineWidth', 2);
    hVectors(end+1) = plot3(ax, Vec(1,1)+Vec(2,1)+[0 Vec(3,1)], Vec(1,2)+Vec(2,2)+[0 Vec(3,2)], Vec(1,3)+Vec(2,3)+[0 Vec(3,3)], '-k', 'LineWidth', 2);
    hVectors(end+1) = plot3(ax, Vec(3,1)+[0 Vec(1,1)], Vec(3,2)+[0 Vec(1,2)], Vec(3,3)+[0 Vec(1,3)], '-k', 'LineWidth', 2);
    hVectors(end+1) = plot3(ax, Vec(3,1)+[0 Vec(2,1)], Vec(3,2)+[0 Vec(2,2)], Vec(3,3)+[0 Vec(2,3)], '-k', 'LineWidth', 2);
    hVectors(end+1) = plot3(ax, Vec(3,1)+Vec(1,1)+[0 Vec(2,1)], Vec(3,2)+Vec(1,2)+[0 Vec(2,2)], Vec(3,3)+Vec(1,3)+[0 Vec(2,3)], '-k', 'LineWidth', 2);
    hVectors(end+1) = plot3(ax, Vec(3,1)+Vec(2,1)+[0 Vec(1,1)], Vec(3,2)+Vec(2,2)+[0 Vec(1,2)], Vec(3,3)+Vec(2,3)+[0 Vec(1,3)], '-k', 'LineWidth', 2);
else
    % Update existing vectors
    set(hVectors(1), 'XData', [0 Vec(1,1)], 'YData', [0 Vec(1,2)], 'ZData', [0 Vec(1,3)]);
    set(hVectors(2), 'XData', [0 Vec(2,1)], 'YData', [0 Vec(2,2)], 'ZData', [0 Vec(2,3)]);
    set(hVectors(3), 'XData', [0 Vec(3,1)], 'YData', [0 Vec(3,2)], 'ZData', [0 Vec(3,3)]);
    set(hVectors(4), 'XData', Vec(1,1)+[0 Vec(3,1)], 'YData', Vec(1,2)+[0 Vec(3,2)], 'ZData', Vec(1,3)+[0 Vec(3,3)]);
    set(hVectors(5), 'XData', Vec(2,1)+[0 Vec(3,1)], 'YData', Vec(2,2)+[0 Vec(3,2)], 'ZData', Vec(2,3)+[0 Vec(3,3)]);
    set(hVectors(6), 'XData', Vec(1,1)+Vec(2,1)+[0 Vec(3,1)], 'YData', Vec(1,2)+Vec(2,2)+[0 Vec(3,2)], 'ZData', Vec(1,3)+Vec(2,3)+[0 Vec(3,3)]);
    set(hVectors(7), 'XData', Vec(3,1)+[0 Vec(1,1)], 'YData', Vec(3,2)+[0 Vec(1,2)], 'ZData', Vec(3,3)+[0 Vec(1,3)]);
    set(hVectors(8), 'XData', Vec(3,1)+[0 Vec(2,1)], 'YData', Vec(3,2)+[0 Vec(2,2)], 'ZData', Vec(3,3)+[0 Vec(2,3)]);
    set(hVectors(9), 'XData', Vec(3,1)+Vec(1,1)+[0 Vec(2,1)], 'YData', Vec(3,2)+Vec(1,2)+[0 Vec(2,2)], 'ZData', Vec(3,3)+Vec(1,3)+[0 Vec(2,3)]);
    set(hVectors(10), 'XData', Vec(3,1)+Vec(2,1)+[0 Vec(1,1)], 'YData', Vec(3,2)+Vec(2,2)+[0 Vec(1,2)], 'ZData', Vec(3,3)+Vec(2,3)+[0 Vec(1,3)]);
end

title(ax, sprintf('Frame: %d/%d', currentFrame, numFrames));
end

function [cmap, minQ, maxQ] = createColormap(Q)
% Create red-green-blue colormap based on charge values
minQ = min(Q);
maxQ = max(Q);

if minQ == maxQ
    cmap = [0 0 1]; % Single color
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