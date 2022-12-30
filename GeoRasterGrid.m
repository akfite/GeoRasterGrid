classdef GeoRasterGrid < matlab.mixin.Copyable
%GEORASTERGRID Fast access to tiled geospatial data
%
%   obj = GEORASTERGRID(___) manages a tiled dataset of geospatial information
%   such as a collection of GeoTIFF images.  It provides O(1) tile lookup, a
%   data cache to reuse recently-accessed data, and the ability to operate without
%   any toolbox licenses.
%
%   Methods:
%
%       Constructor:
%       <a href="matlab:help GeoRasterGrid.GeoRasterGrid">GeoRasterGrid</a>
%
%       Public methods:
%       <a href="matlab:help GeoRasterGrid.get">get</a> (values at specific lat/lon)
%       <a href="matlab:help GeoRasterGrid.roi">roi</a> (get rectangular region of interest)
%       <a href="matlab:help GeoRasterGrid.latlon2tileidx">latlon2tileidx</a> (lookup tile index)
%       <a href="matlab:help GeoRasterGrid.get_tile">get_tile</a> (cached tile accessor)
%       <a href="matlab:help GeoRasterGrid.store_tile">store_tile</a> (add tile to obj state)
%       <a href="matlab:help GeoRasterGrid.clear">clear</a> (drop tiles from obj state)
%       <a href="matlab:help GeoRasterGrid.show">show</a> (display the current state)

%   Author:     Austin Fite
%   Contact:    akfite@gmail.com
%   Date:       10-2022

    properties (Dependent)
        grid_optimized(1,1) logical
    end

    properties
        capacity(1,1) uint16 = 16 % max number of tiles to store
    end

    properties (SetAccess = private, Transient)
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
        lat_intervals(:,2) double % lower & upper bound of each unique interval
        lon_intervals(:,2) double % lower & upper bound of each unique interval
        grid_idx_map(:,:) double % each element is an index into obj.raster_files
    end

    properties (Constant, Hidden)
        georaster_filetypes = [
            ".tif",".tiff",".adf",".asc",".grd",".flt",".dt0",".dt1",".dt2",...
            ".ddf",".dem",".ers",".dat",".img",".grc",".hgt"
            ]';
    end

    %% Events
    events (NotifyAccess = protected)
        TileAdded
        TileRemoved
    end

    %% Constructor
    methods
        function this = GeoRasterGrid(files, varargin)
            %GEORASTERGRID Constructor.
            %
            %   Usage:
            %
            %       obj = GEORASTERGRID(files)
            %       obj = GEORASTERGRID(files, read_limits)
            %       obj = GEORASTERGRID(files, read_limits, capacity)
            %       obj = GEORASTERGRID(files, capacity, read_limits) % order invariant
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
            %             but if the data is in some custom format (e.g. JPEG with
            %             geospatial info encoded in the filepath) the function can be
            %             replaced with a user-defined method
            %           - user-defined functions must accept 1 input (the filepath to
            %             a raster file) and return two outputs: [min_lat max_lat], and
            %             [min_lon max_lon] for that file
            %
            %       capacity (=16) <1x1 integer>
            %           - the number of tiles that can be stored in memory at any time
            %
            %   For more methods, see <a href="matlab:help GeoRasterGrid">GeoRasterGrid</a>
            
            read_limits = @GeoRasterTile.read_limits;
            capacity = 16;

            for i = 1:numel(varargin)
                if ischar(varargin{i})
                    read_limits = str2func(varargin{i});
                elseif isa(varargin{i}, 'function_handle')
                    read_limits = varargin{i};
                elseif isnumeric(varargin{i}) && isscalar(varargin{i})
                    capacity = varargin{i};
                else
                    error('Unexpected input format.');
                end
            end

            validateattributes(files, {'string','char','cell'}, {});
            validateattributes(read_limits, {'function_handle'}, {'scalar'});
            validateattributes(capacity,{'numeric'},{'integer','positive','scalar'});

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
                found_files = dir(fullfile(files{i}, ['**' filesep '*.*']));
                found_files = string(fullfile({found_files.folder}, {found_files.name})');

                [~,~,ext] = fileparts(found_files);
                isgeoraster = contains(ext, this.georaster_filetypes, 'ignorecase', true);

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
            this.capacity = min(capacity, numel(raster_files));

            % performance optimization for uniformly-gridded data will improve the time it
            % takes to lookup the correct tile for each lat/lon from O(n_tiles) to something
            % between O(1) and O(sqrt(n_tiles)).  when tested on a dataset with 15k tiles, 
            % lookup speed was >700x faster
            if ~compute_dense_grid_mapping(this)
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
            %   Usage:
            %
            %       value = obj.get(lat, lon)
            %       value = obj.get(lat, lon, idx)
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
            %           - if you get this wrong, it will completely invalidate
            %             the result!  only use this if you are 100% sure 
            %
            %   Outputs:
            %
            %       value <Nx1xC single or double>
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
            
            % prefer to work with reduced-precision integers for speed
            if numel(this.raster_files) < intmax('uint8')
                idx = uint8(idx);
            elseif numel(this.raster_files) < intmax('uint16')
                idx = uint16(idx);
            else
                % leave as double
            end

            OPTIMIZE_FOR_SINGLE_TILE_ACCESS = isscalar(idx) && idx > 0;

            if OPTIMIZE_FOR_SINGLE_TILE_ACCESS
                % fastest code path (used particularly by method roi())
                unique_tiles = idx;
                n_cached = sum(this.index == idx);
                idx = repmat(idx, size(lat));
            else
                % get index to all tiles we'll need to access
                unique_tiles = unique(idx);
                if isinteger(idx)
                    unique_tiles = unique_tiles(unique_tiles > 0);
                else
                    unique_tiles = unique_tiles(~isnan(unique_tiles));
                end
    
                % process cached (already-loaded) tiles first by
                % putting the cached indices at the front of the array
                icached = ismember(unique_tiles, this.index);
                unique_tiles = [unique_tiles(icached); unique_tiles(~icached)];
                n_cached = sum(icached);
            end

            % process in batches, tile-by-tile (starting with those already in memory)
            for i = 1:numel(unique_tiles)
                tile_idx = unique_tiles(i); % index of tile to load

                if i <= n_cached
                    % tile is already in memory
                    tile = this.tiles(this.index == tile_idx);
                else
                    % not in memory; load from disk
                    tile = this.load_tile(tile_idx);
                    this.store_tile(tile);
                end

                if OPTIMIZE_FOR_SINGLE_TILE_ACCESS
                    value = tile.get(lat, lon);
                else
                    iout = find(idx == tile_idx); % index into output array
                    raster_values = tile.get(lat(iout), lon(iout));

                    if i == 1
                        % pre-allocate output now that type is known
                        value = nan(size(idx),'like',raster_values);
                    end
                
                    % if the raster is multi-channel, we need to adjust our pre-allocated array
                    if size(raster_values,3) ~= size(value,3)
                        value = repmat(value, [1 1 size(raster_values,3)]);
                    end

                    value(iout,1,:) = raster_values;
                end
            end

            value = reshape(value, [orig_sz size(value,3)]);
        end

        function [value, lat, lon] = roi(this, lat_lim, lon_lim, res)
            %GEORASTERGRID/ROI Get the values in a rectangular Region Of Interest (ROI).
            %
            %   Usage:
            %
            %       [value, lat, lon] = obj.roi(lat_lim, lon_lim)
            %       [value, lat, lon] = obj.roi(lat_lim, lon_lim, res)
            %
            %   Inputs:
            %
            %       lat_lim, lon_lim <1x2 double>
            %           - the min & max values of the rectangular region to sample
            %           - degrees
            %
            %       res (=1024) <1x1 double>
            %           - the number of sample points to get along the longer side of 
            %             the ROI rectangle.  the shorter side will be adjusted to
            %             keep the sampling aspect ratio square
            %           - if you wish to sample at the native resolution of the dataset,
            %             you may set RES to -1.  HOWEVER, it is very easy to run out of
            %             memory if you request a large area or your rasters are at very
            %             high resolution, so please be cautious and know what you're doing!
            %
            %   Outputs:
            %
            %       value <NxMxC single or double>
            %           - the raster value at each lat/lon point
            %           - 3rd dimension reserved for multi-channel raster values
            %             (i.e. an RGB raster will return NxMx3)
            %
            %       lat, lon <1xN and 1xM double>
            %           - the latitude & longitude vectors on which the ROI was sampled
            %
            %   For more methods, see <a href="matlab:help GeoRasterGrid">GeoRasterGrid</a>

            if nargin < 4
                res = 1024;
            end

            % it's easy to mess up the input order--let's make it hard to mess up
            lat_lim = sort(lat_lim);
            lon_lim = sort(lon_lim);

            if res == -1 % sample at native resolution
                % first step is to figure out which tiles are involved
                dx = min(abs(diff(lon_lim)), min(diff(this.lon_extents,[],2)));
                dy = min(abs(diff(lat_lim)), min(diff(this.lat_extents,[],2)));

                lat_grid = lat_lim(1):(dy/3):lat_lim(2);
                lon_grid = lon_lim(1):(dx/3):lon_lim(2);
                [lat_grid,lon_grid] = ndgrid(lat_grid, lon_grid);
                idx = this.latlon2tileindex(lat_grid, lon_grid);
                unique_tiles = unique(idx(:));
                
                assert(any(~isnan(unique_tiles)), 'GeoRasterGrid:invalid_roi', ...
                    'Failed to find any valid data in the region of interest.');

                % load first tile to kick things off (we'll assume other tiles have same
                % post spacings)
                tile = this.get_tile(unique_tiles(1));

                % initialize lat/lon to the tile grid points
                lat = tile.lat;
                lon = tile.lon;

                % expand to meet lat/lon limits (while staying on native grid points!)
                dy = abs(diff(tile.lat(1:2)));
                dx = abs(diff(tile.lon(1:2)));

                if lat_lim(1) < lat(1)
                    lat = [fliplr((lat(1)-dy):-dy:lat_lim(1)), lat];
                end
                if lat_lim(2) > lat(end)
                    lat = [lat, lat(end)+dy:dy:lat_lim(2)];
                end
                if lon_lim(1) < lon(1)
                    lon = [fliplr((lon(1)-dx):-dx:lon_lim(1)), lon];
                end
                if lon_lim(2) > lon(end)
                    lon = [lon, (lon(end)+dx):dx:lon_lim(2)];
                end

                % trim
                lat = lat(lat >= lat_lim(1) & lat <= lat_lim(2));
                lon = lon(lon >= lon_lim(1) & lon <= lon_lim(2));

                % sanity check memory usage...
                mem_usage = numel(lat)*numel(lon)*size(tile.raster,3)*4/1e9;

                if mem_usage > 4
                    warning('GeoRasterGrid:high_memory_usage',...
                        'Creating a %d x %d (%.1f GB) matrix; keep an eye on your RAM...', ...
                        numel(lat), numel(lon), mem_usage);
                end

                [lat_grid, lon_grid] = ndgrid(lat, lon);

                if all(unique_tiles == unique_tiles(1))
                    value = this.get(lat_grid, lon_grid, unique_tiles(1));
                else
                    value = this.get(lat_grid, lon_grid);
                end
            else
                % sample at custom resolution (interpolate) 
                validateattributes(res, {'numeric'},{'scalar','positive','real'});

                height = abs(diff(lat_lim));
                width = abs(diff(lon_lim));
    
                if height/width <= 1
                    % shorter in latitude
                    res = round([res*(height/width) res]);
                else
                    % shorter in longitude
                    res = round([res res*(width/height)]);
                end
    
                lat = linspace(lat_lim(1), lat_lim(2), res(1));
                lon = linspace(lon_lim(1), lon_lim(2), res(2));
    
                % first try a shortcut: if diagonal corners are in the same tile, it
                % follows that every other point is in the same tile
                corner_tile = this.latlon2tileindex(lat_lim(:), lon_lim(:));
    
                [lat_grid, lon_grid] = ndgrid(lat, lon);
                
                if all(corner_tile == corner_tile(1))
                    % all data is in the same tile; use slightly faster access path
                    value = this.get(lat_grid, lon_grid, corner_tile(1));
                else
                    % data spans multiple tiles
                    value = this.get(lat_grid, lon_grid);
                end
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
            %       This method attempts to use an optimized solution, but it requires that
            %       the tiles are aligned to a uniform grid.  If the tiles are not aligned,
            %       a brute-force O(n_tiles) method will be employed instead.
            %
            %   For more methods, see <a href="matlab:help GeoRasterGrid">GeoRasterGrid</a>

            if this.grid_optimized
                idx = local_optimized_lookup(lat, lon); return
            else
                idx = local_bruteforce_lookup(lat, lon); return
            end

            % LOCAL (NESTED) FUNCTION (% something like 300x-1000x faster than brute force)
            function index = local_optimized_lookup(lat, lon)
                ix = findinterval(this.lon_intervals(:,1), this.lon_intervals(:,2), lon);
                iy = findinterval(this.lat_intervals(:,1), this.lat_intervals(:,2), lat);

                % we'll just assume all values are in range and if we error we'll go
                % back and correct it.  this is way faster for like 99% of use-cases
                try
                    imap = sub2ind(size(this.grid_idx_map), iy, ix);
                    index = this.grid_idx_map(imap);
                catch mexcept
                    if contains(mexcept.identifier, 'IndexOutOfRange')
                        % some values were outside the valid map area...
                        index = nan(size(lat));
                        valid = iy > 0 & ix > 0;
                        imap = sub2ind(size(this.grid_idx_map), iy(valid), ix(valid));
                        index(valid) = this.grid_idx_map(imap);
                    else
                        rethrow(mexcept);
                    end
                end
            end

            % LOCAL (NESTED) FUNCTION  (brute force O(n_tiles) solution)
            function index = local_bruteforce_lookup(lat, lon)
                [i,j] = find(... "i" is the tile index, and "j" is the query index
                    lat(:)' >= this.lat_extents(:,1) & lat(:)' <= this.lat_extents(:,2) ...
                    & ...
                    lon(:)' >= this.lon_extents(:,1) & lon(:)' < this.lon_extents(:,2));
    
                % allocate output & only fill out the entries that were found
                index = nan(size(lat));
                index(j) = i;
            end
        end

        function tile = get_tile(this, idx)
            %GEORASTERGRID/GET_TILE Cached tile accessor.
            %
            %   Usage:
            %
            %       tile = obj.get_tile(index)
            %
            %   Description:
            %
            %       Gets a tile by index either from the cache (if it has already
            %       been stored) or loads it from disk, if required.
            %
            %   Inputs:
            %
            %       index <1x1 numeric>
            %           - the index of the tile to retrieve
            %           - must be a member of 1:numel(obj.raster_files)
            %
            %   Outputs:
            %
            %       tile <1x1 GeoRasterTile>
            %           - the tile either pulled from object state or loaded from disk
            %
            %   For more methods, see <a href="matlab:help GeoRasterGrid">GeoRasterGrid</a>

            itile = this.index == idx;

            if any(itile)
                tile = this.tiles(itile);
            else
                tile = this.load_tile(idx);
            end
        end

        function store_tile(this, tile)
            %GEORASTERGRID/STORE_TILE Save a tile to the object state.
            %
            %   Usage:
            %
            %       obj.store_tile(tile)
            %
            %   Description:
            %
            %       Caches a GeoRasterTile into the buffer of a GeoRasterGrid so that
            %       future lookups do not require loading from disk.  Subject to the
            %       memory constraint imposed by obj.capacity.  The buffer operates
            %       in a first-in, first-out (FIFO) manner.
            %
            %   Inputs:
            %
            %       tile <1x1 GeoRasterTile>
            %           - the tile to store
            %
            %   Events:
            %
            %       Triggers event 'TileAdded' after the tile has been added to the
            %       cache.  The tile added is always the last in the list: obj.tiles(end).
            %
            %   For more methods, see <a href="matlab:help GeoRasterGrid">GeoRasterGrid</a>

            validateattributes(tile, {'GeoRasterTile'},{'scalar'});

            idx = find(strcmp(this.raster_files, tile.file),1,'first');
            assert(~isempty(idx), 'Map has no record of raster file "%s"', tile.file);

            if any(idx == this.index)
                return % tile is already cached!
            end

            if numel(this.tiles) == this.capacity
                % delete the first to make room (first in, first out)
                this.clear(1);
            end

            this.tiles(end+1) = tile;
            this.index(end+1) = idx;

            notify(this, 'TileAdded');

            if ~isempty(this.ax) && isvalid(this.ax)
                show(this, this.ax); % refresh plot if already active
            end
        end

        function this = clear(this, idx)
            %GEORASTERGRID/CLEAR Delete saved tiles from memory.
            %
            %   Usage:
            %
            %       obj.clear()
            %       obj.clear(idx)
            %
            %   Description:
            %
            %       Deletes one or more tiles from object state.
            %
            %   Inputs:
            %
            %       idx (=1:numel(obj.tiles)) <numeric integer>
            %           - the indices (in obj.tiles) to remove
            %           - defaults to all tiles
            %
            %   Events:
            %
            %       Triggers event 'TileRemoved' after deletion.
            %
            %   For more methods, see <a href="matlab:help GeoRasterGrid">GeoRasterGrid</a>

            if nargin < 2
                idx = 1:numel(this.tiles);
            end

            % keep values in range
            idx = idx(ismember(idx, 1:numel(this.tiles)));

            if isempty(idx)
                return
            end

            % remove entries
            delete(this.tiles(idx));
            this.tiles(idx) = [];
            this.index(idx) = [];

            notify(this, 'TileRemoved');

            if ~isempty(this.ax) && isvalid(this.ax)
                show(this, this.ax); % refresh plot if already active
            end
        end
        
        function ax = show(this, ax)
            %GEORASTERGRID/SHOW Display the current state of the obj.
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
                xlabel(ax,'longitude (degrees)');
                ylabel(ax,'latitude (degrees)');
                axis(ax,'equal');
                axis(ax,'tight');
                xlim(ax,[-180 180]);
                ylim(ax,[-90 90]);
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
            %       tile = obj.load_tile(index)
            %
            %   Description:
            %
            %       Loads a tile from disk.
            %
            %   Inputs:
            %
            %       index <1x1 double>
            %           - the index of the tile to load
            %           - must be a member of 1:numel(obj.raster_files)
            %
            %   Outputs:
            %
            %       tile <1x1 GeoRasterTile>
            %           - the tile that was loaded from disk
            %
            %   For more methods, see <a href="matlab:help GeoRasterGrid">GeoRasterGrid</a>

            validateattributes(index, ...
                {'numeric'},...
                {'scalar','integer','>=',1,'<=',numel(this.raster_files)});

            tile = GeoRasterTile(this.raster_files{index}, ...
                this.lat_extents(index,:), ...
                this.lon_extents(index,:));
        end

        function grid_optimized = compute_dense_grid_mapping(this)
            %GEORASTERGRID/COMPUTE_DENSE_GRID_MAPPING Enable performance optimization.
            %
            %   Usage:
            %
            %       grid_optimized = obj.COMPUTE_DENSE_GRID_MAPPING()
            %
            %   Inputs:
            %
            %       none
            %
            %   Outputs:
            %
            %       grid_optimized <1x1 logical>
            %           - true if the optimization was successful
            %
            %   Description:
            %
            %       This method represents a sparsely-connected grid of tiles as
            %       a dense grid so that tile lookups can occur in near-constant time
            %       (as opposed to checking the bounds of every tile... very slow).
            %       However, this can only be done if the gridded tiles are all aligned
            %       to lines of longitude & latitude.  In other words, if you were to
            %       extend each of the 4 lines that make up the sides of any tile to
            %       +Inf and -Inf, those lines should only be coincident with the EDGES 
            %       of other tiles.  The edges of any tile (when extended to infinity)
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
            %
            %   For more methods, see <a href="matlab:help GeoRasterGrid">GeoRasterGrid</a>

            grid_optimized = false;

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
                grid_start = this.lat_extents(icheck,1);
                grid_end = this.lat_extents(icheck,2);

                grid_start_u = uniquetol(grid_start, eps(single(pi)));
                grid_end_u = uniquetol(grid_end, eps(single(pi)));

                if numel(grid_start_u) > 1 || numel(grid_end_u) > 1
                    % grids are not aligned; abort
                    return
                else
                    % they really need to be aligned to machine precision... if
                    % that's not the case we'll force it to be true.
                    tol = eps(max(this.lat_extents(:)));
                    if any(abs(grid_start - grid_start(1)) > tol)
                        grid_start(:) = mean(grid_start);
                    end
                    if any(abs(grid_end - grid_end(1)) > tol)
                        grid_end(:) = mean(grid_end);
                    end

                    % don't check these tiles again
                    checked(icheck) = true;
                    y_intervals = vertcat(y_intervals, [grid_start(1) grid_end(1)]);
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
                grid_start = this.lon_extents(icheck,1);
                grid_end = this.lon_extents(icheck,2);

                grid_start_u = uniquetol(grid_start, eps(single(pi)));
                grid_end_u = uniquetol(grid_end, eps(single(pi)));

                if numel(grid_start_u) > 1 || numel(grid_end_u) > 1
                    % grids are not aligned; abort
                    return
                else
                    % they really need to be aligned to machine precision... if
                    % that's not the case we'll force it to be true.
                    tol = eps(max(this.lon_extents(:)));
                    if any(abs(grid_start - grid_start(1)) > tol)
                        grid_start(:) = mean(grid_start);
                    end
                    if any(abs(grid_end - grid_end(1)) > tol)
                        grid_end(:) = mean(grid_end);
                    end

                    % don't check these tiles again
                    checked(icheck) = true;
                    x_intervals = vertcat(x_intervals, [grid_start(1) grid_end(1)]);
                end
            end

            % OK, tiles are aligned to the same grid--we can optimize this! start
            % by sorting the intervals
            x_intervals = sortrows(x_intervals);
            y_intervals = sortrows(y_intervals);

            % define a dense grid and map indices from the dense grid to the sparse grid
            center_lon = mean(x_intervals,2);
            center_lat = mean(y_intervals,2);
            [center_lat, center_lon] = ndgrid(center_lat, center_lon);

            this.grid_idx_map = latlon2tileindex(this, center_lat, center_lon);
            this.lon_intervals = x_intervals;
            this.lat_intervals = y_intervals;

            grid_optimized = true;
        end
    end

    %% Accessors
    methods
        function ok = get.grid_optimized(this)
            ok = ~isempty(this.grid_idx_map);
        end
    end

end
