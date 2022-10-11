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
        interp(1,:) char = 'linear' % setting this prop updates interpolant
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

            % create grid vectors
            validateattributes(lat_extents, {'numeric'},{'numel',2});
            validateattributes(lon_extents, {'numeric'},{'numel',2});

            % need to offset by 1/2 pixel at edges to account for the grid cell postings
            dx = diff(lon_extents) / size(img,2);
            dy = diff(lat_extents) / size(img,1);

            longrid = linspace(min(lon_extents)+dx/2, max(lon_extents)-dx/2, size(img,2));
            latgrid = linspace(min(lat_extents)+dy/2, max(lat_extents)-dy/2, size(img,1));

            [latgrid, longrid] = ndgrid(latgrid, longrid);

            this.interpolant = griddedInterpolant(latgrid, longrid, double(img), this.interp);
            this.file = raster_file;
        end
    end

    %% Public
    methods
        function values = get(this, lat, lon, is_grid)
            %GEORASTERTILE/GET Access the value at lat/lon pairs.
            %
            %   Usage:
            %
            %       value = obj.get(lat, lon)
            %
            %   Inputs:
            %
            %       lat, lon <double arrays>
            %           - 
            assert(all(size(lat) == size(lon)), ...
                'GeoRasterTile:size_mismatch', ...
                'Latitude & longitude must have the same size.');

            values = this.interpolant(lat, lon);
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

    %% Property Accessors
    methods
        function v = get.raster(this)
            v = this.interpolant.Values;
        end

        function lat = get.lat(this)
            lat = this.interpolant.GridVectors{2};
        end

        function lon = get.lon(this)
            lon = this.interpolant.GridVectors{1};
        end

        function m = get.interp(this)
            m = this.interpolant.Method;
        end
    end

end

