%% Updated Main Function
function density_analyzer_rashid(xyz_file, restart_file)
    % Global constants
    global AVOGADRO ANGSTROM3_TO_CM3 ATOMIC_MASSES GROUPS
    
    AVOGADRO = 6.02214076e23;
    ANGSTROM3_TO_CM3 = 1e-24;
    ATOMIC_MASSES = containers.Map(...
        {'H', 'O', 'C', 'N', 'S', 'Cl', 'Na', 'K', 'Al', 'Pt', 'Ag', 'Au', 'Fe', 'Mg', 'Si', 'P'}, ...
        [1.008, 16.00, 12.01, 14.01, 32.07, 35.45, 22.99, 39.10, 26.98, 195.08, 107.87, 196.97, 55.85, 24.31, 28.09, 30.97]);
    
    GROUPS = struct(...
        'Water', {{'O', 'H'}}, ...
        'Oxygen', {{'O'}}, ...
        'Hydrogen', {{'H'}});

    % Check inputs
    if nargin ~= 2
        error('Usage: density_analyzer <trajectory.xyz> <restart_file.restart>');
    end
    
    try
        [Lx, Ly, Lz] = parse_restart_box(restart_file);
        box_size = [Lx, Ly, Lz];
        n_z = 100;
        densities = calculate_densities(xyz_file, box_size, n_z);
        
        % Define smoothing regions
        smoothing_regions = {...
            struct('z_range', [0, 10], 'window', 3, 'polyorder', 2), ...
            struct('z_range', [10, 30], 'window', 15, 'polyorder', 3), ...
            struct('z_range', [30, 50], 'window', 3, 'polyorder', 2)};
        
        plot_results(densities, box_size, 5, 0.1, smoothing_regions);
    catch ME
        fprintf('Error: %s\n', ME.message);
        exit(1);
    end
end

%% Parse Restart File
function [Lx, Ly, Lz] = parse_restart_box(filename)
    content = fileread(filename);
    cell_section = regexp(content, '&CELL\s*([^&]*)&END CELL', 'tokens', 'dotexceptnewline');
    if isempty(cell_section)
        error('CELL section not found');
    end
    
    cell_lines = strsplit(strtrim(cell_section{1}{1}), '\n');
    vectors = struct('A', [], 'B', [], 'C', []);
    
    for i = 1:length(cell_lines)
        line = strtrim(cell_lines{i});
        if startsWith(line, {'A', 'B', 'C'})
            parts = strsplit(line);
            vec_type = parts{1};
            coords = str2double(parts(2:4));
            vectors.(vec_type) = coords;
        end
    end
    
    if isempty(vectors.A) || isempty(vectors.B) || isempty(vectors.C)
        error('Missing cell vectors');
    end
    
    Lx = norm(vectors.A);
    Ly = norm(vectors.B);
    Lz = norm(vectors.C);
end

%% Calculate Densities
function avg_densities = calculate_densities(xyz_file, box_size, n_z)
    global AVOGADRO ANGSTROM3_TO_CM3 ATOMIC_MASSES GROUPS
    
    Lx = box_size(1);
    Ly = box_size(2);
    Lz = box_size(3);
    dz = Lz / n_z;
    volume_per_bin = Lx * Ly * dz * ANGSTROM3_TO_CM3;
    
    groups = fieldnames(GROUPS);
    sum_densities = struct();
    for i = 1:length(groups)
        sum_densities.(groups{i}) = zeros(n_z, 1);
    end
    frame_count = 0;
    
    fid = fopen(xyz_file, 'r');
    while true
        header = fgetl(fid);
        if header == -1
            break;
        end
        n_atoms = str2double(header);
        fgetl(fid); % Skip comment line
        
        current_mass = struct();
        for i = 1:length(groups)
            current_mass.(groups{i}) = zeros(n_z, 1);
        end
        
        for j = 1:n_atoms
            line = fgetl(fid);
            parts = strsplit(line);
            symbol = parts{1};
            x = str2double(parts{2});
            y = str2double(parts{3});
            z = str2double(parts{4});
            
            % Periodic boundary condition
            z_wrapped = mod(z, Lz);
            bin_idx = floor(z_wrapped / dz) + 1;
            if bin_idx > n_z
                bin_idx = n_z;
            end
            
            % Assign mass to groups
            for k = 1:length(groups)
                group = groups{k};
                elements = GROUPS.(group);
                if ismember(symbol, elements)
                    mass = ATOMIC_MASSES(symbol);
                    current_mass.(group)(bin_idx) = current_mass.(group)(bin_idx) + mass;
                end
            end
        end
        
        % Accumulate densities
        frame_count = frame_count + 1;
        for k = 1:length(groups)
            group = groups{k};
            mass_g = current_mass.(group) / AVOGADRO;
            density = mass_g / volume_per_bin;
            sum_densities.(group) = sum_densities.(group) + density;
        end
    end
    fclose(fid);
    
    % Average densities
    avg_densities = struct();
    for k = 1:length(groups)
        group = groups{k};
        avg_densities.(group) = sum_densities.(group) / frame_count;
    end
end

%% Plot Results
function plot_results(avg_densities, box_size, smooth_window, min_prominence, smoothing_regions)
    global GROUPS
    
    Lz = box_size(3);
    groups = fieldnames(GROUPS);
    n_z = length(avg_densities.(groups{1}));
    dz = Lz / n_z;
    z_centers = (0.5:1:n_z) * dz;
    
    figure('Position', [100, 100, 1200, 600]);
    
    for k = 1:length(groups)
        group = groups{k};
        densities = avg_densities.(group);
        smoothed = zeros(size(densities));
        processed = false(n_z, 1);
        
        % Apply region-specific smoothing
        if ~isempty(smoothing_regions)
            for r = 1:length(smoothing_regions)
                region = smoothing_regions{r};
                z_min = region.z_range(1);
                z_max = region.z_range(2);
                window = region.window;
                polyorder = region.polyorder;
                
                mask = (z_centers >= z_min) & (z_centers <= z_max);
                indices = find(mask);
                
                if isempty(indices)
                    continue;
                end
                
                subset = densities(indices);
                window_adj = min(window, length(subset));
                window_adj = max(window_adj, 3);
                if mod(window_adj, 2) == 0
                    window_adj = window_adj - 1;
                end
                
                if window_adj >= 3
                    subset_smoothed = sgolayfilt(subset, polyorder, window_adj);
                    smoothed(indices) = subset_smoothed;
                    processed(indices) = true;
                end
            end
        end
        
        % Apply default smoothing to unprocessed regions
        unprocessed = find(~processed);
        if ~isempty(unprocessed)
            subset = densities(unprocessed);
            window_adj = min(smooth_window, length(subset));
            window_adj = max(window_adj, 3);
            if mod(window_adj, 2) == 0
                window_adj = window_adj - 1;
            end
            if window_adj >= 3
                subset_smoothed = sgolayfilt(subset, 3, window_adj);
                smoothed(unprocessed) = subset_smoothed;
            end
        end
        
        % Clip negative values
        smoothed = max(smoothed, 0);
        
        % Find minima
        [~, locs] = findpeaks(-smoothed, 'MinPeakProminence', min_prominence);
        if ~isempty(locs)
            min_z = z_centers(locs);
            fprintf('\n%s minima at:\n', group);
            fprintf('  - %.2f Å\n', min_z);
        else
            fprintf('\n%s: No clear minima found\n', group);
        end
        
        % Plot
        plot(z_centers, smoothed, 'DisplayName', group);
        hold on;
    end
    
    xlabel('Z-coordinate (Å)');
    ylabel('Density (g/cm³)');
    title('Density Profile with Region-Specific Smoothing');
    legend('show');
    grid on;
    hold off;
end