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
        hfig % currently-displayed figure
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
        end
    end

    %% Public interface
    methods
        function value = get(this, lat, lon, idx)
            %GEORASTERGRID/get Interpolate the value at specific lat/lon point(s).
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
        
        function ax = show(this)
            %GEORASTERGRID/SHOW Display the current state of the map.
            %
            %   Usage:
            %
            %       ax = obj.show()
            %
            %   Description:
            %
            %       Displays the state of the map on a 
            %
            %   For more methods, see <a href="matlab:help GeoRasterGrid">GeoRasterGrid</a>

            N = numel(this.raster_files);

            if isempty(this.hfig) || ~isvalid(this.hfig)
                this.hfig = figure('Name',class(this),'NumberTitle','off');
                colordef(this.hfig,'black'); %#ok<*COLORDEF>
                ax = axes('nextplot','add','parent',this.hfig);
                colormap(ax,'gray');

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
                    'Color',[0.5 0.5 0.5],... gray
                    'HitTest','off');

                box(ax,'on');
                xlim(ax,[-180 180]);
                ylim(ax,[-90 90]);
                xlabel(ax,'longitude (degrees)');
                ylabel(ax,'latitude (degrees)');
                axis(ax,'equal');
                axis(ax,'tight');
            else
                ax = findobj(this.hfig,'type','axes');
            end

            delete(findobj(ax,'tag','active-tiles'));

            % refresh overlay of tiles currently in memory
            if ~isempty(this.tiles)
                shapes = polyshape(this.tiles);
                plot(shapes, ...
                    'FaceColor','g',...
                    'FaceAlpha',0.25,...
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

            drawnow;
        end
    end

    %% Private helper methods
    methods (Access = private)
        function tile = load_tile(this, index)
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

            if ~isempty(this.hfig) && isvalid(this.hfig)
                show(this); % refresh plot if already active
            end
        end

        function idx = latlon2tileindex(this, lat, lon)
            % "i" is the tile index, and "j" is the query index
            [i,j] = find(...
                lat(:)' >= this.lat_extents(:,1) & lat(:)' <= this.lat_extents(:,2) ...
                & ...
                lon(:)' >= this.lon_extents(:,1) & lon(:)' <= this.lon_extents(:,2));

            % pre-allocate output to the correct size in case some entries
            % were not found
            idx = nan(size(lat));
            idx(j) = i;
        end
    end

end
