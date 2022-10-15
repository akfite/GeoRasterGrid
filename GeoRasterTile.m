classdef GeoRasterTile < matlab.mixin.Copyable
%GEORASTERTILE A rectangular section of geospatial data.
%
%

%   Author:     Austin Fite
%   Contact:    akfite@gmail.com
%   Date:       10-2022

    properties (SetAccess = immutable)
        file(1,:) char
    end
    
    properties (Dependent)
        raster double
        lat(1,:) double
        lon(1,:) double
    end

    properties
        interp(1,:) char = 'cubic' % setting this prop updates interpolant
    end

    properties (SetAccess = private)
    end

    properties (Access = private)
        interpolant % griddedInterpolant for each channel
    end
    
    %% Constructor
    methods
        function this = GeoRasterTile(raster_file, varargin)
            %GEORASTERTILE Construct an instance of this class.
            %
            %   Usage:
            %
            %       obj = GEORASTERTILE(raster)
            %       obj = GEORASTERTILE(raster, lat_extents, lon_extents)
            %
            %   Inputs:
            %
            %       raster <1xN char or MxNxC numeric>
            %           - ideally this is the filepath to a geo raster file
            %             (readable by readgeoraster), but it can also be an
            %             image or a path to an image if the latitude & longitude
            %             extents of that image are defined in the 2nd & 3rd args
            %           - note that rows & columns must be aligned to parallels of
            %             latitude & meridians of longitude
            %
            %       lat_extents, lon_extents <1x2 double>
            %           - only necessary if geo information is not encoded in the raster
            %           - values of latitude & longitude at the edges of the raster image
            %             (be mindful to use the true pixel EDGE--not the pixel center)
            %           - degrees
            %
            %   For more methods, see <a href="matlab:help GeoRasterTile">GeoRasterTile</a>

            narginchk(1,3);
            validateattributes(raster_file, {'char','string','numeric'},{});

            if ischar(raster_file) || isstring(raster_file)
                try
                    [img,R] = readgeoraster(raster_file);
                catch
                    assert(nargin > 1, ...
                        ['File must be readable by readgeoraster() when not ' ...
                        'providing spatial referencing info (nargin == 1).']);

                    try
                        img = imread(raster_file);
                    catch
                        img = readmatrix(raster_file);
                    end
                end
            else
                % provided an image directly
                validateattributes(raster_file,{'numeric'},{});
                img = raster_file;
                raster_file = '';
            end

            if nargin == 1 || nargin == 2
                % using mapping toolbox
                lat_extents = R.LatitudeLimits;
                lon_extents = R.LongitudeLimits;
            elseif nargin == 3
                % let user provide info explicitly (e.g. when map toolbox not avail)
                lat_extents = varargin{1};
                lon_extents = varargin{2};
            else
                error('Unexpected number of input arguments.');
            end

            validateattributes(lat_extents, {'numeric'},{'numel',2});
            validateattributes(lon_extents, {'numeric'},{'numel',2});

            % need to offset by 1/2 pixel at edges to account for the grid cell postings
            dx = diff(lon_extents) / size(img,2);
            dy = diff(lat_extents) / size(img,1);
            longrid = linspace(min(lon_extents)+dx/2, max(lon_extents)-dx/2, size(img,2));
            latgrid = linspace(min(lat_extents)+dy/2, max(lat_extents)-dy/2, size(img,1));

            % TODO: assume standard image coords, but use info from R if possible
%             img = flipud(img);

            [latgrid, longrid] = ndgrid(latgrid, longrid);

            this.interpolant = griddedInterpolant(latgrid, longrid, double(img), this.interp);
            this.file = raster_file;
        end
    end

    %% Public
    methods
        function values = get(this, lat, lon)
            %GEORASTERTILE/GET Access the value at lat/lon pairs.
            %
            %   Usage:
            %
            %       values = obj.get(lat, lon)
            %
            %   Inputs:
            %
            %       lat, lon <double matrix>
            %           - the query points
            %           - degrees on the interval [-90,90] and [-180, 180]
            %           - must be the same size
            %
            %   Outputs:
            %
            %       values <double array>
            %           - the raster value(s)
            %           - will always have type double (regardless of the native
            %             type of the raster)
            %           - first two dimensions will match the size of lat & lon,
            %             and the third dimension corresponds to the color channels
            %             of the raster (e.g. red, green, blue)
            %
            %   For more methods, see <a href="matlab:help GeoRasterTile">GeoRasterTile</a>

            assert(all(size(lat) == size(lon)), ...
                'GeoRasterTile:size_mismatch', ...
                'Latitude & longitude must have the same size.');
            assert(ismatrix(lat) && ismatrix(lon), ...
                'GeoRasterTile:not_matrix', ...
                'Both latitude & longitude must be matrices (2d arrays)');

            values = this.interpolant(lat, lon);
        end

        function [values, lat, lon] = roi(this, lat_bounds, lon_bounds, res)
            %GEORASTERTILE/ROI Get values in a region of interest (ROI).
            %
            %   Usage:
            %
            %       [values, lat, lon] = obj.roi(lat_bounds, lon_bounds)
            %       [values, lat, lon] = obj.roi(lat_bounds, lon_bounds, res)
            %
            %   Inputs:
            %
            %       lat_bounds, lon_bounds <1x2 double>
            %           - the min & max values of the region's rectangle
            %           - degrees on [-180, 180]
            %
            %       res <1x1 or 1x2 double>
            %           - sampling resolution in [lat lon] as fractional degrees
            %           - scalar inputs will be used in both dimensions
            %           - e.g. 1/120 to sample at 30 arcsec resolution
            %           - by default, samples at native resolution
            %
            %   Outputs:
            %
            %       values <NxMxC double>
            %           - the grid of raster values
            %           - sampled at the native resolution of the raster
            %             (at the exact post spacing the raster is defined on)
            %
            %       lat, lon <1xN and 1xM double>
            %           - the geodetic coordinates that the returned values
            %             are sampled on
            %
            %   For more methods, see <a href="matlab:help GeoRasterTile">GeoRasterTile</a>

            if nargin < 3
                % sample at native resolution
                ilat = find(this.lat >= lat_bounds(1) & this.lat <= lat_bounds(2));
                ilon = find(this.lon >= lon_bounds(1) & this.lon <= lon_bounds(2));
    
                % index directly (no interpolation)
                lat = this.lat(ilat);
                lon = this.lon(ilon);
                values = this.raster(ilat, ilon, :);
            else 
                % interpolate to custom resolution
                if isscalar(res)
                    res = [res res];
                end

                lat = lat_bounds(1):res(1):lat_bounds(2);
                lon = lon_bounds(1):res(2):lon_bounds(2);

                ilat = find(lat >= this.lat(1) & lat <= this.lat(end));
                ilon = find(lon >= lon_bounds(1) & lon <= lon_bounds(2));

                lat = lat(ilat);
                lon = lon(ilon);

                if isempty(ilat) || isempty(ilon)
                    values = []; return
                end

                [lat_grid, lon_grid] = ndgrid(lat, lon);
                values = this.get(lat_grid, lon_grid);
            end
        end

        function flag = contains(this, lat, lon)
            %GEORASTERTILE/CONTAINS Check if the tile contains some point(s).
            %
            %

            assert(isscalar(this),'GeoRasterTile:scalar_only',...
                'GeoRasterTile/contains cannot operate on more than one object at once.');
            
            lat_grid = this.lat;
            lon_grid = this.lon;
            dx = diff(lon_grid(1:2));
            dy = diff(lat_grid(1:2));

            flag = lat >= (lat_grid(1) - dy/2) ...
                & lat <= (lat_grid(end) + dy/2) ...
                & lon >= (lon_grid(1) - dx/2) ...
                & lon <= (lon_grid(end) + dx/2);
        end

        function [bx, by] = bounding_box(this)
            %GEORASTERTILE/BOUNDING_BOX Get bounding box coordinates in degrees.
            %
            %   Usage:
            %
            %       [bx, by] = obj.bounding_box();
            %
            %   Inputs:
            %
            %       obj <Nx1 GeoRasterTile>
            %           - one or more tiles
            %
            %   Outputs:
            %
            %       bx, by <Nx6 double>
            %           - x and y coordinates of bounding boxes (one row per tile)
            %           - use directly with plot() by flushing the result to column
            %             vectors: plot(bx(:), by(:))
            %           - latitude & longitude, in degrees
            %
            %   For more methods, see <a href="matlab:help GeoRasterTile">GeoRasterTile</a>

            if isscalar(this)
                bx = [this.lon([1 1 end end 1]) NaN];
                by = [this.lat([1 end end 1 1]) NaN];
                
                % need to adjust by 1/2 pixel at edges to account for grid cell postings
                % (coords are recorded at pixel center; we want the true edge of grid)
                dx = diff(this.lon(1:2));
                dy = diff(this.lat(1:2));
                bx = bx + ([-dx -dx dx dx -dx NaN] * 0.5);
                by = by + ([-dy dy dy -dy -dy NaN] * 0.5);
            else
                bx = nan(numel(this), 6);
                by = nan(numel(this), 6);

                for i = 1:numel(this)
                    [bx(i,:), by(i,:)] = this(i).bounding_box(); % call w/ scalar obj
                end
            end
        end

        function shape = polyshape(this)
            %GEORASTERTILE/POLYSHAPE Create a polyshape with raster boundaries in degrees.
            %
            %   Usage:
            %
            %       shape = obj.polyshape();
            %
            %   Inputs:
            %
            %       obj <Nx1 GeoRasterTile>
            %           - one or more tiles
            %
            %   Outputs:
            %
            %       shape <Nx1 polyshape>
            %           - polyshape objects for the raster boundaries
            %           - latitude & longitude, in degrees
            %
            %   For more methods, see <a href="matlab:help GeoRasterTile">GeoRasterTile</a>

            [bx,by] = this.bounding_box();
            for i = 1:size(bx,1)
                shape(i,1) = polyshape(bx(1:4), by(1:4), 'simplify', false); %#ok<AGROW> 
            end
        end
    end

    %% Property Setters/Getters
    methods
        function set.interp(this, value)
            %GEORASTERTILE/SET.INTERP Update the interpolation type.
            %
            %   This method re-creates the interpolant (rather than updating the
            %   existing interpolant's Method) because as of R2020b there is a
            %   bug with MATLAB's code that causes performance to massively degrade--
            %   on the order of 300x slowdown.

            validatestring(value, {'linear','nearest','spline','cubic','makima'});
            [latgrid, longrid] = ndgrid(this.lat, this.lon); %#ok<*MCSUP>
            this.interpolant = griddedInterpolant(latgrid, longrid, this.raster, value);
            this.interp = value;
        end

        function v = get.raster(this)
            v = this.interpolant.Values;
        end

        function lat = get.lat(this)
            lat = this.interpolant.GridVectors{2};
        end

        function lon = get.lon(this)
            lon = this.interpolant.GridVectors{1};
        end
    end

    %% Static methods
    methods (Static)
        function [lat_lim, lon_lim] = read_limits(file)
            [folder,name] = fileparts(file);

            if exist(fullfile(folder, name, '.json'),'file')
                json_txt = readascii(fullfile(folder, name, '.json'));
                R = jsondecode(json_txt);

                switch R.AngleUnit
                    case 'degree'
                        lat_lim = R.LatitudeLimits;
                        lon_lim = R.LongitudeLimits;
                    case {'radian', 'radians'} % actually not sure what MATLAB uses...
                        lat_lim = R.LatitudeLimits * pi/180;
                        lon_lim = R.LongitudeLimits * pi/180;
                    otherwise
                        validatestring(R.AngleUnit,{'degree','radian','radians'});
                end

            elseif license('checkout','map_toolbox')
                info = georasterinfo(file);
                R = info.RasterReference;

                % raster coordinates must be parallel to lat/lon lines
                validateattributes(R,...
                    {...
                        'map.rasterref.GeographicCellsReference', ...
                        'map.rasterref.GeographicPostingsReference' ...
                    }, ...
                    {'scalar'});

                switch R.AngleUnit
                    case 'degree'
                        lat_lim = R.LatitudeLimits;
                        lon_lim = R.LongitudeLimits;
                    case {'radian', 'radians'} % actually not sure what MATLAB uses...
                        lat_lim = R.LatitudeLimits * pi/180;
                        lon_lim = R.LongitudeLimits * pi/180;
                    otherwise
                        validatestring(R.AngleUnit,{'degree','radian','radians'});
                end
            else
                error('GeoRasterTile:read_limits:parse_fail',...
                    ['Failed to parse latitude & longitude limits from file.  You can define ' ...
                    'your own function to parse the limits if this default doesn''t work.']);
            end
        end
    end

end

