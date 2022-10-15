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
            %           - function handle to create a GeoRasterTile given the
            %             filepath to a raster
            %           - by default it will support files readable by readgeoraster,
            %             but if the data is in some custom format (e.g. JPEG2000 with
            %             geospatial info encoded in the filepath) the function can be
            %             replaced with a wrapper that calls GeoRasterTile with 2+ inputs
            %             after parsing info from the filepath (or wherever it is stored)
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

            % search remaining folders for raster files
            for i = 1:numel(files)
                found_files = dir(fullfile(files{i}, '**\*.*'));
                found_files = string(fullfile({found_files.folder}, {found_files.name})');

                [~,~,ext] = fileparts(found_files);
                isgeoraster = contains(ext, this.supported_filetypes, 'ignorecase', true);

                if any(isgeoraster)
                    raster_files = vertcat(raster_files, found_files(isgeoraster)); %#ok<*AGROW> 
                end
            end

            assert(~isempty(raster_files),'GeoRasterGrid:none_found',...
                'Failed to find any georaster files.');

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
            %       Displays the tiles that the object currently has cached on
            %       a world map.
            %
            %   For more methods, see <a href="matlab:help GeoRasterGrid">GeoRasterGrid</a>

            hfig = figure;
            colordef(hfig,'black'); %#ok<*COLORDEF>
            ax = axes('nextplot','add');
            colormap(ax,'gray');

            dx = GeoRasterGrid.dx; %#ok<*PROP>
            dy = GeoRasterGrid.dy;

            % overlay on map of the Earth
            earth_texture = imread('landOcean.jpg');
            image(...
                linspace(-180,180, size(earth_texture,2)), ...
                fliplr(linspace(-90,90, size(earth_texture,1))), ...
                earth_texture, ...
                'parent', ax)

            % overlay the currently-loaded tiles
            for i = 1:numel(this.tiles)
                lat = fliplr(linspace(this.lat(i,1)+dy/2, this.lat(i,2)-dy/2, this.height));
                lon = linspace(this.lon(i,1)+dx/2, this.lon(i,2)-dx/2, this.width);
                imagesc(lon,lat,this.tiles{i},'parent',ax);
            end

            % show the outline of all the possible tiles in gray
            for i = 1:numel(this.paths)
                [bx,by] = griddata2d(this.lon_extents(i,:), this.lat_extents(i,:));
                plot(ax,bx,by,':','color',[0.5 0.5 0.5]);
            end

            % show the outline of all the cached tiles in red
            for i = 1:numel(this.tiles)
                [bx,by] = griddata2d(this.lon(i,:), this.lat(i,:));
                plot(ax,bx,by,'r-');
            end

            box(ax,'on');

            xlim(ax,[-180 180]);
            ylim(ax,[-90 90]);

            title(ax, sprintf('%s: %d/%d capacity', ...
                class(this), numel(this.tiles), this.capacity));
            xlabel(ax,'longitude (degrees)');
            ylabel(ax,'latitude (degrees)');
        end
    end

    %% Private helper methods
    methods (Access = private)
        function tile = load_tile(this, index)
            tile = GeoRasterTile(this.raster_files{index});

            if numel(this.tiles) == this.capacity
                % delete the first to make room (first in, first out)
                delete(this.tiles(1));
                this.tiles(1) = [];
                this.index(1) = [];
            end

            this.tiles(end+1) = tile;
            this.index(end+1) = index;
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

        function [row,col] = latlon2rowcol(this, lat, lon, index)
            dx = GeoRasterGrid.dx;
            dy = GeoRasterGrid.dy;

            % x0 and y0 are the pixel edge, so dx/dy will step us across whole
            % pixels (as opposed to starting at pixel center, which would step
            % across two half-pixels)
            x0 = this.lon(index,1);
            y0 = this.lat(index,1);

            x = lon;
            y = lat;

            % y-axis is inverted, x-axis is normal
            row = this.height - ceil((y - y0)./dy) + 1;
            col = ceil((x - x0)./dx);
        end
    end

    %% Static methods
    methods (Static)
        function write_json_metadata(files, read_metadata)
            if nargin < 2
                read_metadata = @GeoRasterGrid.read_limits;
            end


        end
    end

end
