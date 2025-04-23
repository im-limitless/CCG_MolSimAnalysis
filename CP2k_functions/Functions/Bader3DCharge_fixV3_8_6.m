function Bader3DCharge_fixV3_8_6(TrajXYZ, TrajABC, TrajQmc, varargin)
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

% Bader3DCharge(TrajXYZ, TrajABC, TrajQmc, SelectedIndices, ...)
% Plots 3D charge distribution with multiple selection sets
%
% New Parameters:
%   'AdditionalSelectionSets', {} - Cell array of additional index sets
%                                   Each set can be:
%                                   - Numeric vector (static for all frames)
%                                   - Cell array of vectors (per-frame indices)

% Input parsing
p = inputParser;
addRequired(p, 'TrajXYZ', @(x) isnumeric(x) || iscell(x));
addRequired(p, 'TrajABC', @(x) isnumeric(x) || iscell(x));
addRequired(p, 'TrajQmc', @(x) isnumeric(x) || iscell(x));
addOptional(p, 'SelectedIndices', [], @(x) isempty(x) || isvector(x) || (iscell(x) && all(cellfun(@isvector, x))));
addParameter(p, 'AdditionalSelectionSets', {}, @iscell);
addParameter(p, 'PlayVideo', false);
addParameter(p, 'SaveVideo', false);
addParameter(p, 'FrameRate', 24);
addParameter(p, 'OutputFileName', 'charge_video.mp4');
addParameter(p, 'Interactive', true);
addParameter(p, 'ViewAxis', 'X', @(x) ismember(x,{'X','Y'}));
addParameter(p, 'VideoView', [], @isnumeric);
addParameter(p, 'LoopVideo', false, @islogical);
parse(p, TrajXYZ, TrajABC, TrajQmc, varargin{:});

params = p.Results;
SelectedIndices = params.SelectedIndices;
AdditionalSets = params.AdditionalSelectionSets;

% Convert trajectory inputs to cells
[TrajXYZ, TrajABC, TrajQmc] = convert_to_cells(TrajXYZ, TrajABC, TrajQmc);
numFrames = numel(TrajXYZ);
Natoms = size(TrajXYZ{1}, 1);

% Process all selection sets =============================================
allSets = [{SelectedIndices}, AdditionalSets];
validSets = allSets(~cellfun(@isempty, allSets));

% Convert each set to per-frame cell array
processedSets = cell(numel(validSets), 1);
for sIdx = 1:numel(validSets)
    currSet = validSets{sIdx};
    
    if iscell(currSet)
        % Validate cell array format
        if numel(currSet) == 1
            % Single cell - repeat for all frames
            processedSets{sIdx} = repmat(currSet(1), numFrames, 1);
        elseif numel(currSet) == numFrames
            % Already per-frame format
            processedSets{sIdx} = currSet(:);
        else
            error('Selection set %d has %d elements (expected 1 or %d)',...
                  sIdx, numel(currSet), numFrames);
        end
    else
        % Numeric vector - expand to all frames
        processedSets{sIdx} = repmat({currSet(:)}, numFrames, 1);
    end
end

% Combine indices from all sets per frame
CombinedIndices = cell(numFrames, 1);
for f = 1:numFrames
    frameIndices = [];
    for sIdx = 1:numel(processedSets)
        if ~isempty(processedSets{sIdx}{f})
            frameIndices = union(frameIndices, processedSets{sIdx}{f}(:));
        end
    end
    
    % Default to all atoms if no selections
    if isempty(frameIndices) && ~isempty(validSets)
        frameIndices = (1:Natoms)';
    end
    CombinedIndices{f} = frameIndices;
end

% Initialize figure and axes (unchanged)
% ... [Rest of the original code] ...

% Initialize figure and axes
fig = figure('Color','w', 'Name','3D Charge Viewer',...
    'Position',[100 100 1200 800], 'NumberTitle','off');
ax = axes('Parent', fig, 'Position',[0.1 0.1 0.7 0.8]);
hold(ax, 'on');
axis(ax, 'equal');
grid(ax, 'on');

xlabel(ax, 'X-axis');
ylabel(ax, 'Y-axis');
zlabel(ax, 'Z-axis');

if strcmpi(params.ViewAxis, 'X')
    view(ax, [0 0]);
else
    view(ax, [90 0]);
end
set(ax, 'CameraUpVector', [0 0 1]);

colorbar(ax);

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
    
    originalView = ax.View;
    if ~isempty(params.VideoView)
        view(ax, params.VideoView);
    end
end

% UI controls
if params.Interactive
    uiX = 0;
    
    uicontrol('Style','togglebutton', 'String','Play',...
        'Position',[uiX 700 80 30], 'Tag','playBtn',...
        'Callback',@(src,~) playCallback(src));
    
    uicontrol('Style','slider',...
        'Position',[uiX 750 250 30],...
        'Min',1, 'Max',numFrames, 'Value',1,...
        'SliderStep',[1/numFrames 10/numFrames],...
        'Tag','frameSlider',...
        'Callback',@(src,~) sliderCallback(src));
end

% Initialize animation variables
currentFrame = 1;
hAtoms = [];
hVectors = [];

% Initial plot
[hAtoms, hVectors] = updatePlot(ax, TrajXYZ{1}, TrajABC{1}, TrajQmc{1},...
    xs, ys, zs, Radius, CombinedIndices{1}, currentFrame, numFrames, hAtoms, hVectors);

% Update all calls to updatePlot to use CombinedIndices
% Example in video export:

% Video export
if params.SaveVideo
    for frame = 1:numFrames
        [hAtoms, hVectors] = updatePlot(ax, TrajXYZ{frame}, TrajABC{frame},...
            TrajQmc{frame}, xs, ys, zs, Radius, CombinedIndices{frame}, frame, numFrames, hAtoms, hVectors);
        % ... [Remaining code] ...
        writeVideo(writerObj, getframe(fig));
        drawnow('limitrate');
        if params.Interactive
            set(findobj(fig,'Tag','frameSlider'), 'Value', frame);
        end
    end
    close(writerObj);
    
    if ~isempty(params.VideoView)
        view(ax, originalView);
    end
end

% Animation controls
if params.PlayVideo && params.Interactive
    autoPlayAnimation();
end

% Callback functions
    function playCallback(src)
        persistent loopActive
        if isempty(loopActive)
            loopActive = false;
        end
        
        if loopActive
            loopActive = false;
            return
        end
        
        loopActive = true;
        set(src, 'Value', 1);
        initialFrame = currentFrame;
        
        while loopActive && ishandle(src)
            if currentFrame < numFrames
                currentFrame = currentFrame + 1;
            else
                if params.LoopVideo
                    currentFrame = 1;
                else
                    loopActive = false;
                end
            end
            
            [hAtoms, hVectors] = updatePlot(ax, TrajXYZ{currentFrame},...
                TrajABC{currentFrame}, TrajQmc{currentFrame},...
                xs, ys, zs, Radius, CombinedIndices{currentFrame}, currentFrame, numFrames, hAtoms, hVectors);
            
            set(findobj(fig,'Tag','frameSlider'), 'Value', currentFrame);
            
            drawnow;
            
            currentButtonState = get(src,'Value');
            if ~currentButtonState
                loopActive = false;
                break;
            end
            
            tStart = tic;
            while toc(tStart) < 1/params.FrameRate && loopActive
                drawnow;
                if get(src,'Value') == 0
                    loopActive = false;
                    break;
                end
            end
        end
        
        set(src, 'Value', 0);
        loopActive = false;
        
        if currentFrame > numFrames
            currentFrame = numFrames;
            updatePlot(ax, TrajXYZ{currentFrame},...
                TrajABC{currentFrame}, TrajQmc{currentFrame},...
                xs, ys, zs, Radius, CombinedIndices{currentFrame}, currentFrame, numFrames, hAtoms, hVectors);
        end
    end

    function sliderCallback(src)
        frame = round(get(src,'Value'));
        if frame ~= currentFrame
            currentFrame = frame;
            [hAtoms, hVectors] = updatePlot(ax, TrajXYZ{currentFrame},...
                TrajABC{currentFrame}, TrajQmc{currentFrame},...
                xs, ys, zs, Radius, CombinedIndices{currentFrame}, currentFrame, numFrames, hAtoms, hVectors);
        end
    end

    function autoPlayAnimation()
        playCount = 0;
        while params.PlayVideo && (params.LoopVideo || playCount == 0)
            for frame = 1:numFrames
                if ~params.PlayVideo && ~params.LoopVideo
                    break;
                end
                
                [hAtoms, hVectors] = updatePlot(ax, TrajXYZ{frame},...
                    TrajABC{frame}, TrajQmc{frame},...
                    xs, ys, zs, Radius, CombinedIndices{frame}, frame, numFrames, hAtoms, hVectors);
                
                if params.Interactive
                    set(findobj(fig,'Tag','frameSlider'), 'Value', frame);
                end
                
                pause(1/params.FrameRate);
            end
            playCount = playCount + 1;
        end
    end

hold(ax, 'off');
end

%% Helper functions
function [TrajXYZ, TrajABC, TrajQmc] = convert_to_cells(TrajXYZ, TrajABC, TrajQmc)
% Convert inputs to cell arrays

% XYZ conversion
if ~iscell(TrajXYZ)
    if ndims(TrajXYZ) == 3
        TrajXYZ = arrayfun(@(i) squeeze(TrajXYZ(i,:,:)), 1:size(TrajXYZ,1), 'UniformOutput', false)';
    else
        TrajXYZ = {TrajXYZ};
    end
end

% ABC conversion
if ~iscell(TrajABC)
    if ndims(TrajABC) == 3
        TrajABC = arrayfun(@(i) squeeze(TrajABC(i,:,:)), 1:size(TrajABC,1), 'UniformOutput', false)';
    else
        TrajABC = {TrajABC};
    end
end

% Qmc conversion
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

function [hAtoms, hVectors] = updatePlot(ax, XYZ, ABC, Qmc, xs, ys, zs, Radius, CombinedIndices, currentFrame, numFrames, hAtoms, hVectors)
% Update plot elements

delete(findobj(ax,'Type','text','Tag','FrameText'));

[cmap, minQ, maxQ] = createColormap(Qmc(CombinedIndices));
colormap(ax, cmap);
caxis(ax, [minQ, maxQ]);

if isempty(hAtoms) || ~isvalid(hAtoms(1))
    hAtoms = createAtoms(ax, XYZ, xs, ys, zs, Radius, Qmc, CombinedIndices, cmap, minQ, maxQ);
else
    updateAtoms(hAtoms, XYZ, xs, ys, zs, Radius, Qmc, CombinedIndices, cmap, minQ, maxQ);
end

Vec = diag(ABC);
if isempty(hVectors) || ~isvalid(hVectors(1))
    hVectors = createVectors(ax, Vec);
else
    updateVectors(hVectors, Vec);
end

text(ax, 0.05, 0.95, 0.9, sprintf('Frame: %d/%d', currentFrame, numFrames),...
    'Units','normalized', 'Tag','FrameText');
end

function hAtoms = createAtoms(ax, XYZ, xs, ys, zs, Radius, Qmc, CombinedIndices, cmap, minQ, maxQ)
hAtoms = gobjects(size(XYZ,1), 1);
for i = 1:size(XYZ,1)
    color = [0.7 0.7 0.7];
    alpha = 0.3;
    
    if ismember(i, CombinedIndices)
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

function updateAtoms(hAtoms, XYZ, xs, ys, zs, Radius, Qmc, CombinedIndices, cmap, minQ, maxQ)
for i = 1:size(XYZ,1)
    color = [0.7 0.7 0.7];
    alpha = 0.3;
    if ismember(i, CombinedIndices)
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
hVectors = gobjects(12,1);

hVectors(1) = plot3(ax, [0 Vec(1,1)],[0 Vec(1,2)],[0 Vec(1,3)],'-k','LineWidth',2);
hVectors(2) = plot3(ax, [0 Vec(2,1)],[0 Vec(2,2)],[0 Vec(2,3)],'-k','LineWidth',2);
hVectors(3) = plot3(ax, [0 Vec(3,1)],[0 Vec(3,2)],[0 Vec(3,3)],'-k','LineWidth',2);

hVectors(4) = plot3(ax, [Vec(1,1) Vec(1,1)+Vec(2,1)],...
                     [Vec(1,2) Vec(1,2)+Vec(2,2)],...
                     [Vec(1,3) Vec(1,3)+Vec(2,3)],'-k','LineWidth',2);
hVectors(5) = plot3(ax, [Vec(2,1) Vec(1,1)+Vec(2,1)],...
                     [Vec(2,2) Vec(1,2)+Vec(2,2)],...
                     [Vec(2,3) Vec(1,3)+Vec(2,3)],'-k','LineWidth',2);

hVectors(6) = plot3(ax, [Vec(1,1) Vec(1,1)+Vec(3,1)],...
                     [Vec(1,2) Vec(1,2)+Vec(3,2)],...
                     [Vec(1,3) Vec(1,3)+Vec(3,3)],'-k','LineWidth',2);
hVectors(7) = plot3(ax, [Vec(2,1) Vec(2,1)+Vec(3,1)],...
                     [Vec(2,2) Vec(2,2)+Vec(3,2)],...
                     [Vec(2,3) Vec(2,3)+Vec(3,3)],'-k','LineWidth',2);
hVectors(8) = plot3(ax, [Vec(1,1)+Vec(2,1) Vec(1,1)+Vec(2,1)+Vec(3,1)],...
                     [Vec(1,2)+Vec(2,2) Vec(1,2)+Vec(2,2)+Vec(3,2)],...
                     [Vec(1,3)+Vec(2,3) Vec(1,3)+Vec(2,3)+Vec(3,3)],'-k','LineWidth',2);

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
if isempty(Q)
    cmap = [0 0 1];
    minQ = 0;
    maxQ = 1;
    return;
end

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