classdef GeoRasterGrid < matlab.mixin.Copyable
%GEORASTERGRID Manage gridded geospatial raster data.
%
%   obj = GEORASTERGRID(___) manages a tiled dataset of geospatial information
%   from multiple georaster files.  It provides optimized access to raster values and
%   can scale up to datasets too large to fit in memory.  It is primarily intended to
%   be used with sets of GeoTIFF, .dt0/1/2, or other georaster data files.
%
%   Methods:
%
%       Constructor:
%       <a href="matlab:help GeoRasterGrid.GeoRasterGrid">GeoRasterGrid</a>
%
%       Public methods:
%       <a href="matlab:help GeoRasterGrid.get">get</a> (values at specific lat/lon)
%       <a href="matlab:help GeoRasterGrid.roi">roi</a> (values in a region of interest)
%       <a href="matlab:help GeoRasterGrid.clear">clear</a> (drop all tiles in memory)
%       <a href="matlab:help GeoRasterGrid.plot">plot</a>

%   Author:     Austin Fite
%   Contact:    akfite@gmail.com
%   Date:       10-2022

    properties
        capacity(1,1) uint16 = 36 % max number of tiles to store
    end

    properties (SetAccess = private)
        % current tiles loaded in memory
        tiles(:,1) GeoRasterTile
        index(:,1) double % index to each tile that is currently-loaded
    end

    properties (SetAccess = private)
        % metadata about ALL tiles
        raster_files(:,1) string % filepath to each raster
        lat_extents(:,2) double % bounds of each raster; degrees
        lon_extents(:,2) double % bounds of each raster; degrees
    end

    properties (Access = private)
        % currently-displayed axis
        ax

        % metadata to enable performance optimization
        lat_intervals(:,2) double
        lon_intervals(:,2) double
        grid_idx_map(:,:) double
    end

    properties (Constant, Hidden)
        % list of supported file extensions
        supported_filetypes = [
            ".tif",".tiff",".adf",".asc",".grd",".flt",".dt0",".dt1",".dt2",...
            ".ddf",".dem",".ers",".dat",".img",".grc",".hgt"
            ]';
    end

    %% Constructor
    methods
        function this = GeoRasterGrid(files, read_limits)
            %GEORASTERGRID Constructor.
            %
            %   Usage:
            %
            %       obj = GEORASTERGRID(files)
            %       obj = GEORASTERGRID(files, read_limits)
            %
            %   Inputs:
            %
            %       files <1xN char or Nx1 string/cellstr>
            %           - the location of the georaster files
            %           - folder path inputs will be recursively searched for
            %             files with the following extensions:
            %
            %                   .tif, .tiff
            %                   .adf
            %                   .asc, .grd
            %                   .flt
            %                   .dt0, .dt1, .dt2
            %                   .ddf
            %                   .dem
            %                   .ers
            %                   .dat
            %                   .img
            %                   .grc
            %                   .hgt
            %
            %           - file path inputs will be interpreted as files that can
            %             be read directly by readgeoraster()
            %           - while multiple file types are supported, they are NOT
            %             simultaneously supported.  if more than one file type is
            %             found, the most frequently-occuring filetype will be kept
            %             and the others discarded
            %
            %       read_limits (=@GeoRasterTile.read_limits) <1x1 function_handle>
            %           - function handle to read the latitude & longitude limits of
            %             the raster given its filepath
            %           - by default it will support files readable by readgeoraster,
            %             but if the data is in some custom format (e.g. JPEG2000 with
            %             geospatial info encoded in the filepath) the function can be
            %             replaced with a user-defined method
            %           - user-defined functions must accept 1 input (the filepath to
            %             a raster file) and return two outputs: [min_lat max_lat], and
            %             [min_lon max_lon] for that file
            %
            %   For more methods, see <a href="matlab:help GeoRasterGrid">GeoRasterGrid</a>
            
            if nargin < 2
                read_limits = @GeoRasterTile.read_limits;
            end

            if nargin < 1
                files = "D:\data\blue-marble\dem";
            end

            validateattributes(files, {'string','char','cell'}, {});
            validateattributes(read_limits, {'function_handle'}, {'scalar'});

            % enforce common format
            if ~isstring(files)
                files = string(files);
            end

            raster_files = files(isfile(files)); % take filepaths as-is
            files(isfile(files)) = []; % do not search under filepaths

            % accumulate files for backup search so that we don't crawl through
            % the filesystem twice
            all_files = string.empty;

            % search remaining folders for raster files
            for i = 1:numel(files)
                found_files = dir(fullfile(files{i}, '**\*.*'));
                found_files = string(fullfile({found_files.folder}, {found_files.name})');

                [~,~,ext] = fileparts(found_files);
                isgeoraster = contains(ext, this.supported_filetypes, 'ignorecase', true);

                if any(isgeoraster)
                    % found expected file type; add to list
                    raster_files = vertcat(raster_files, found_files(isgeoraster)); %#ok<*AGROW>
                else
                    % accumulate for backup search of non-standard image types
                    all_files = vertcat(all_files, found_files);
                end
            end

            % last-ditch effort: check for other image file types if user provided a
            % custom read function and some directory to search
            if isempty(raster_files) && ~isequal(read_limits, @GeoRasterTile.read_limits)
                [~,~,ext] = fileparts(all_files);

                isgeoraster = contains(ext, ...
                    [".jp2",".jpeg2000",".png",".jpg",".jpeg",".bmp"], ...
                    'ignorecase', true);

                if any(isgeoraster)
                    raster_files = all_files(isgeoraster);
                end
            end

            assert(~isempty(raster_files),'GeoRasterGrid:none_found', ...
                'Failed to find any raster files.');

            % check that all files have the same extension
            [~,~,ext] = fileparts(raster_files);
            assert(numel(unique(lower(ext))) == 1, 'GeoRasterGrid:mixed_types', ...
                'All files must have the same extension (no mixing of filetypes allowed)');

            % call user-configurable function to parse tile boundaries
            [lat_lim, lon_lim] = cellfun(read_limits, raster_files, 'uniform', false);

            % assign to object
            this.raster_files = raster_files;
            this.lat_extents = vertcat(lat_lim{:});
            this.lon_extents = vertcat(lon_lim{:});

            % performance optimization for uniformly-gridded data will improve the time it
            % takes to lookup the correct tile for each lat/lon from O(n_tiles) to O(1).
            % when tested on a dataset with 10k+ tiles, lookup speed was >700x faster
            ok = compute_dense_grid_mapping(this);

            if ~ok
                warning('GeoRasterGrid:optimization_failed',...
                    'Raster grid edges are not aligned; performance optimizations disabled.');
            end
        end
    end

    %% Public interface
    methods
        function value = get(this, lat, lon, idx)
            %GEORASTERGRID/GET Interpolate the value at specific lat/lon point(s).
            %
            %   Usage (assuming map is a 1x1 GeoRasterGrid):
            %
            %       value = map.get(lat, lon)
            %       value = map.get(lat, lon, idx)
            %
            %   Inputs:
            %
            %       lat <double matrix>
            %           - geodetic latitude (degrees)
            %           - positive North from equator
            %
            %       lon <double matrix>
            %           - geodetic longitude (degrees)
            %           - positive East from Greenwich meridian
            %
            %       idx <double matrix>
            %           - optional; purely a performance optimzation
            %           - this is the index into the tile that each lat/lon pair
            %             will be sourced from (in case the caller can more
            %             efficiently calculate it)
            %
            %   Outputs:
            %
            %       value <Nx1xC double>
            %           - the raster value at each lat/lon point
            %           - 3rd dimension reserved for multi-channel raster values
            %             (i.e. an RGB raster will return Nx1x3)
            %
            %   For more methods, see <a href="matlab:help GeoRasterGrid">GeoRasterGrid</a>

            assert(all(size(lat) == size(lon)), ...
                'GeoRasterGrid:size_mismatch', ...
                'Latitude & longitude must have the same size.');
            assert(ismatrix(lat) && ismatrix(lon), ...
                'GeoRasterGrid:not_matrix', ...
                'Both latitude & longitude must be matrices (2d arrays)');

            orig_sz = size(lat);

            lat = lat(:);
            lon = lon(:);

            if nargin < 4
                idx = this.latlon2tileindex(lat, lon);
            end

            idx = idx(:);

            % get index to all tiles we'll need to access
            unique_tiles = unique(idx);
            unique_tiles = unique_tiles(~isnan(unique_tiles));

            % process cached (already-loaded) tiles first by
            % putting the cached indices at the front of the array
            icached = ismember(unique_tiles, this.index);
            unique_tiles = [unique_tiles(icached); unique_tiles(~icached)];
            n_cached = sum(icached);

            % pre-allocate output
            value = nan(size(idx));

            % process in batches, tile-by-tile (starting with those already in memory)
            for i = 1:numel(unique_tiles)
                tile_idx = unique_tiles(i); % index of tile to load
                iout = find(idx == tile_idx); % index into output array

                if i <= n_cached
                    % tile is already in memory
                    tile = this.tiles(this.index == tile_idx);
                else
                    % not in memory; load from disk
                    tile = this.load_tile(tile_idx);
                end

                % GeoRasterTile object manages value lookup
                raster_values = tile.get(lat(iout), lon(iout));
                
                % if the raster is multi-channel, we need to adjust our pre-allocated array
                if size(raster_values,3) ~= size(value,3)
                    value = repmat(value, [1 1 size(raster_values,3)]);
                end

                value(iout,1,:) = raster_values;
            end

            value = reshape(value, [orig_sz size(value,3)]);
        end
        
        function ax = show(this, ax)
            %GEORASTERGRID/SHOW Display the current state of the map.
            %
            %   Usage:
            %
            %       ax = obj.show()
            %       ax = obj.show(ax)
            %
            %   Description:
            %
            %       Displays the state of the map into an axis.  If no axis is provided,
            %       a new figure will be created to host the axis.
            %
            %   For more methods, see <a href="matlab:help GeoRasterGrid">GeoRasterGrid</a>

            N = numel(this.raster_files);

            if nargin < 2
                ax = this.ax;
            end

            if isempty(ax) || ~isvalid(ax)
                hfig = figure('Name',class(this),'NumberTitle','off');
                colordef(hfig,'black'); %#ok<*COLORDEF>
                ax = axes('nextplot','add','parent',hfig);
                colormap(ax,'gray');
            end

            % first-time setup
            if isempty(ax.Children)
                hold(ax,'on');

                % overlay on map of the Earth
                earth_texture = imread('world_map.jpg');
                image(...
                    linspace(-180,180, size(earth_texture,2)), ...
                    linspace(90,-90, size(earth_texture,1)), ...
                    earth_texture, ...
                    'parent', ax)
    
                % show the outline of all the possible tiles in gray
                bx = nan(6,N);
                by = nan(6,N);
                for i = 1:N
                    bx(1:5,i) = this.lon_extents(i,[1 1 2 2 1]);
                    by(1:5,i) = this.lat_extents(i,[1 2 2 1 1]);
                end
                plot(ax,bx(:),by(:),':',...
                    'Color',[0.75 0.75 0.75],... gray
                    'HitTest','off');

                box(ax,'on');
                xlim(ax,[-180 180]);
                ylim(ax,[-90 90]);
                xlabel(ax,'longitude (degrees)');
                ylabel(ax,'latitude (degrees)');
                axis(ax,'equal');
                axis(ax,'tight');
                zoom(ax,'on');
            end

            delete(findobj(ax,'tag','active-tiles'));

            % refresh overlay of tiles currently in memory
            if ~isempty(this.tiles)
                shapes = polyshape(this.tiles);
                plot(shapes, ...
                    'FaceColor','g',...
                    'FaceAlpha',0.125,...
                    'parent',ax,...
                    'EdgeColor','g',...
                    'HitTest','off',...
                    'Tag','active-tiles');
            end

            title(ax, sprintf('%s: %d of %d grid tiles in memory (max capacity = %d)', ...
                class(this), ...
                numel(this.tiles), ...
                numel(this.raster_files), ...
                this.capacity));

            this.ax = ax;
            drawnow;
        end
    end

    %% Private helper methods
    methods (Access = private)
        function tile = load_tile(this, index)
            %GEORASTERGRID/LOAD_TILE Load the tile at a given index.
            %
            %   Usage:
            %
            %       tile = obj.LOAD_TILE(index)
            %
            %   Inputs:
            %
            %       index <1x1 double>
            %           - the index to the tile to load
            %           - must be valid on 1:numel(obj.raster_files)
            %
            %   Outputs:
            %
            %       tile <1x1 GeoRasterTile>
            %           - the map tile that was loaded
            %           - this tile will also be inserted into obj.tiles
            %           - if obj.capacity has been reached, this method will remove
            %             obj.tiles(1) to make room for the new tile

            tile = GeoRasterTile(...
                this.raster_files{index}, ...
                this.lat_extents(index,:), ...
                this.lon_extents(index,:));

            if numel(this.tiles) == this.capacity
                % delete the first to make room (first in, first out)
                delete(this.tiles(1));
                this.tiles(1) = [];
                this.index(1) = [];
            end

            this.tiles(end+1) = tile;
            this.index(end+1) = index;

            if ~isempty(this.ax) && isvalid(this.ax)
                show(this, this.ax); % refresh plot if already active
            end
        end

        function idx = latlon2tileindex(this, lat, lon)
            %GEORASTERGRID/LATLON2TILEINDEX Fast tile lookup for many lat-lon points.
            %
            %   Usage:
            %
            %       idx = obj.LATLON2TILEINDEX(lat, lon)
            %
            %   Inputs:
            %
            %       lat, lon <double matrix>
            %           - geodetic latitude & longitude
            %
            %   Outputs:
            %
            %       idx <double matrix>
            %           - the index to the tile that contains each lat-lon pair
            %           - NaN if no tile contains a point
            %           - matches the size of the input arrays
            %
            %   Notes:
            %
            %       This method attempts to use an O(1) solution, but it requires that
            %       the tiles are aligned to a uniform grid.  If the tiles are not aligned,
            %       a brute-force O(n_tiles) method will be employed instead.

            if ~isempty(this.grid_idx_map)
                idx = local_optimized_lookup(lat, lon);
            else
                idx = local_bruteforce_lookup(lat, lon);
            end

            function index = local_optimized_lookup(lat, lon)
                % O(1) solution
                ix = findinterval(this.lon_intervals(:,1), this.lon_intervals(:,2), lon);
                iy = findinterval(this.lat_intervals(:,1), this.lat_intervals(:,2), lat);
                imap = sub2ind(size(this.grid_idx_map), iy, ix);
                index = this.grid_idx_map(imap);
            end

            function index = local_bruteforce_lookup(lat, lon)
                % O(n_tiles) solution
                [i,j] = find(... "i" is the tile index, and "j" is the query index
                    lat(:)' >= this.lat_extents(:,1) & lat(:)' <= this.lat_extents(:,2) ...
                    & ...
                    lon(:)' >= this.lon_extents(:,1) & lon(:)' < this.lon_extents(:,2));
    
                % allocate output & only fill out the entries that were found
                index = nan(size(lat));
                index(j) = i;
            end
        end

        function is_optimized = compute_dense_grid_mapping(this)
            %GEORASTERGRID/COMPUTE_DENSE_GRID_MAPPING Enable performance optimization.
            %
            %   Usage:
            %
            %       is_optimized = obj.COMPUTE_DENSE_GRID_MAPPING()
            %
            %   Inputs:
            %
            %       none
            %
            %   Outputs:
            %
            %       is_optimized <1x1 logical>
            %           - true if the optimization was successful
            %
            %   Description:
            %
            %       This method converts a sparsely-connected grid of tiles into
            %       a dense grid so that tile lookups can occur in O(1) time (as
            %       opposed to having to check the bounds of every tile--O(n_tiles)).
            %       However, it can only be done if the gridded tiles are all aligned
            %       to lines of longitude & latitude.  In other words, if you were to
            %       extend each of the 4 lines that make up the sides of any tile to
            %       +Inf and -Inf, those lines should only be coincident with the EDGES 
            %       of other tiles.  The sides of any tile (when extended to infinity)
            %       can NEVER intersect the body of any other tile, else this optimization
            %       will be automatically disabled and the method will return false.
            %
            %       One way of thinking about this method is that it is "filling in the gaps"
            %       for a set of tiles.  For example, the tiles may form some irregular
            %       shape as they trace out land mass and have gaps over the oceans.  But to
            %       calculate an index directly, the shape must be a dense matrix.  This method
            %       creates "ghost" tiles to fill in the rectangle that encloses the sparse
            %       set of real tiles.  Tile index lookups happen against this dense grid in
            %       constant time, and the values of the dense grid map back to indices in
            %       the sparse grid to be used directly with obj.load_tile().

            is_optimized = false;

            % check that all tiles are aligned to the same grid lines
            checked = false(size(this.raster_files));
            cy = mean(this.lat_extents,2);
            y_intervals = [];

            for i = 1:numel(this.raster_files)
                if checked(i), continue; end

                start_lat = this.lat_extents(i,1);
                end_lat = this.lat_extents(i,2);
                
                % find all tiles that have a center point inside this interval
                icheck = cy >= start_lat & cy <= end_lat;
                grid_start = uniquetol(this.lat_extents(icheck,1), 1e-12);
                grid_end = uniquetol(this.lat_extents(icheck,2), 1e-12);

                if numel(grid_start) > 1 || numel(grid_end) > 1
                    % grids are not aligned; abort
                    return
                else
                    % don't check these tiles again
                    checked(icheck) = true;
                    y_intervals = vertcat(y_intervals, [grid_start grid_end]);
                end
            end

            checked = false(size(this.raster_files));
            cx = mean(this.lon_extents,2);
            x_intervals = [];

            for i = 1:numel(this.raster_files)
                if checked(i), continue; end
                
                start_lon = this.lon_extents(i,1);
                end_lon = this.lon_extents(i,2);
                
                % find all tiles that have a center point inside this interval
                icheck = cx >= start_lon & cx <= end_lon;
                grid_start = uniquetol(this.lon_extents(icheck,1), 1e-12);
                grid_end = uniquetol(this.lon_extents(icheck,2), 1e-12);

                if numel(grid_start) > 1 || numel(grid_end) > 1
                    % grids are not aligned; abort
                    return
                else
                    % don't check these tiles again
                    checked(icheck) = true;
                    x_intervals = vertcat(x_intervals, [grid_start grid_end]);
                end
            end

            % OK, tiles are aligned to the same grid--we can optimize this! start
            % by sorting the intervals
            x_intervals = sortrows(x_intervals);
            y_intervals = sortrows(y_intervals);

            % check that intervals never overlap
            if any(x_intervals(2:end,1) - x_intervals(1:end-1,2) < 0)
                return
            end
            if any(y_intervals(2:end,1) - y_intervals(1:end-1,2) < 0)
                return
            end

            % define a dense grid and map indices from the dense grid to the sparse grid
            center_lon = mean(x_intervals,2);
            center_lat = mean(y_intervals,2);
            [center_lat, center_lon] = ndgrid(center_lat, center_lon);

            this.grid_idx_map = latlon2tileindex(this, center_lat, center_lon);
            this.lon_intervals = x_intervals;
            this.lat_intervals = y_intervals;

            is_optimized = true;
        end
    end

end
