classdef GeoRasterMap < matlab.mixin.Copyable
%GEORASTERMAP Manage tiled geospatial information.
%
%   obj = GEORASTERMAP(___) manages a tiled dataset of geospatial information
%   from one or more georaster files.
%
%   Methods:
%
%       Constructor:
%       <a href="matlab:help GeoRasterMap.GeoRasterMap">GeoRasterMap</a>
%
%       Public methods:
%       <a href="matlab:help GeoRasterMap.get">get</a> (values at specific lat/lon)
%       <a href="matlab:help GeoRasterMap.roi">roi</a> (values in a region of interest)
%       <a href="matlab:help GeoRasterMap.clear">clear</a> (drop all tiles in memory)
%       <a href="matlab:help GeoRasterMap.plot">plot</a>

%   Author:     Austin Fite
%   Contact:    akfite@gmail.com
%   Date:       10-2022

    properties
        capacity(1,1) uint16 = 16 % max number of tiles to store
    end

    properties (SetAccess = private)
        % current tiles loaded in memory
        tiles(:,1) GeoRasterTile

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
        function this = GeoRasterMap(files, capacity)
            %GEORASTERMAP Constructor.
            %
            %   Usage:
            %
            %       obj = GEORASTERMAP(files)
            %       obj = GEORASTERMAP(files, capacity)
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
            %                   .dt0, .dt1, dt2
            %                   .DDF
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
            %       capacity (=16) <1x1 unsigned integer>
            %           - the maximum number of raster tiles that may be in memory
            %             at any one time
            %           - this value can also be set via the public class property
            %             at any time after construction
            %
            %   For more methods, see <a href="matlab:help GeoRasterMap">GeoRasterMap</a>

            if nargin < 2
                capacity = 16;
            end
            if nargin < 1
                files = "D:\data\blue-marble\dem";
            end

            validateattributes(files, {'string','char','cell'}, {});
            validateattributes(capacity, {'numeric'}, {'scalar','positive','integer'});

            % enforce common format
            if ~isstring(files)
                files = string(files);
            end

            raster_files = string.empty;

            % remove filepaths from directories provided
            raster_files = vertcat(raster_files, files(isfile(files)));

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

            assert(~isempty(raster_files),'GeoRasterMap:none_found',...
                'Failed to find any georaster files.');

            % import metadata (first pass: try to find a .json file that encodes metadata
            % without requiring the mapping toolbox).  note that we loop so that we can
            % break early without checking every file, which can be slow over slow networks
            json_files = raster_files + ".json";
            json_exists = true;

            for i = 1:numel(json_files)
                json_exists = json_exists && exist(json_files{i},'file');
                if ~json_exists, break; end
            end

            if json_exists % load metadata without mapping toolbox
                error('@TODO: not yet implemented')
            else % import metadata using mapping toolbox
                [ok, errmsg] = license('checkout','map_toolbox');
                assert(ok, 'GeoRasterMap:license_failure', errmsg);

                for i = numel(raster_files):-1:1
                    info = georasterinfo(raster_files{i});
                    R = info.RasterReference;

                    % assign to class props
                    this.raster_files(i) = raster_files(i);

                    switch R.AngleUnit
                        case 'degree'
                            this.lat_extents(i,:) = R.LatitudeLimits;
                            this.lon_extents(i,:) = R.LongitudeLimits;
                        case {'radian', 'radians'} % actually not sure what MATLAB uses...
                            this.lat_extents(i,:) = R.LatitudeLimits * pi/180;
                            this.lon_extents(i,:) = R.LongitudeLimits * pi/180;
                        otherwise
                            validatestring(R.AngleUnit,{'degree','radian','radians'});
                    end
                end
            end
        end
    end

    %% Public interface
    methods
        function value = get(this, varargin)
            %GEORASTERMAP/get Check the type of a point on the Earth.
            %
            %   Usage (assuming map is a 1x1 GeoRasterMap):
            %
            %       value = map.get(lla <Nx2 or Nx3>)
            %       value = map.get(lat, lon)
            %
            %   Inputs:
            %
            %       lat <vector of latitudes>
            %           - Geodetic latitude (degrees)
            %           - Positive North from equator
            %
            %       lon <vector of longitudes>
            %           - Geodetic longitude (degrees)
            %           - Positive East from Greenwich meridian
            %
            %   Outputs:
            %
            %       value <Nx1 BackgroundType>
            %           - the classification of each input point
            %
            %   For more methods, see <a href="matlab:help GeoRasterMap">GeoRasterMap</a>

            if nargin == 2
                lat = rad2deg(varargin{1}(:,1));
                lon = rad2deg(varargin{1}(:,2));
            else
                lat = rad2deg(varargin{1}(:));
                lon = rad2deg(varargin{2}(:));
            end

            idx = this.latlon2tileindex(lat, lon);

            % get unique, non-nan indices.  these are the tiles we
            % need to access
            unique_tiles = unique(idx);
            unique_tiles = unique_tiles(~isnan(unique_tiles));

            % process cached (already-loaded) tiles first (do this by
            % putting cached indices at the front of the array)
            icached = ismember(unique_tiles, this.index);
            unique_tiles = [unique_tiles(icached); unique_tiles(~icached)];
            n_cached = sum(icached);

            % pre-allocate output
            value = repmat(BackgroundType.INVALID, size(idx));

            % process in batches, tile-by-tile
            for i = 1:numel(unique_tiles)
                tile = unique_tiles(i);

                if i > n_cached
                    this.load_tile(tile);
                end

                % get index into output array
                iout = find(idx == tile);

                % get index into cached tiles
                icached = this.index == tile;

                % convert lat/lon to tile (image) coordinates
                [row,col] = this.latlon2rowcol(lat(iout), lon(iout), icached);

                % get image value at each pixel--images are valued 1 where land
                % mass exists, and 0 where ocean exists
                ind = sub2indq([this.height this.width], row, col);
                value(iout) = this.tiles{icached}(ind);
            end
        end
        
        function ax = show(this)
            %GEORASTERMAP/SHOW Display the current state of the map.
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
            %   For more methods, see <a href="matlab:help GeoRasterMap">GeoRasterMap</a>

            hfig = figure;
            colordef(hfig,'black'); %#ok<*COLORDEF>
            ax = axes('nextplot','add');
            colormap(ax,'gray');

            dx = GeoRasterMap.dx; %#ok<*PROP>
            dy = GeoRasterMap.dy;

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
        function load_tile(this, index)
            data = load(this.paths{index});

            if numel(this.tiles) == this.capacity
                % delete the first to make room (first in, first out)
                this.tiles(1) = [];
                this.lat(1,:) = [];
                this.lon(1,:) = [];
                this.index(1) = [];
            end

            this.tiles{end+1} = data.landmass;
            this.lat(end+1,:) = this.lat_extents(index,:);
            this.lon(end+1,:) = this.lon_extents(index,:);
            this.index(end+1) = index;
        end

        function idx = latlon2tileindex(this, lat, lon)
            [i,j] = find(...
                lat(:)' >= this.lat_extents(:,1) & lat(:)' < this.lat_extents(:,2) ...
                & ...
                lon(:)' >= this.lon_extents(:,1) & lon(:)' < this.lon_extents(:,2));

            % pre-allocate output to the correct size in case some entries
            % were not found
            idx = nan(size(lat));

            % we made a 96 x N matrix, where N is numel(lat) or numel(lon).  so
            % "i" is the tile index, and "j" is the query index
            idx(j) = i;
        end

        function [row,col] = latlon2rowcol(this, lat, lon, index)
            dx = GeoRasterMap.dx;
            dy = GeoRasterMap.dy;

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

end
