classdef test_GeoRasterTile < matlab.unittest.TestCase

    properties (Constant)
        tol = single(0.0001); % meters vertical
    end

    properties (TestParameter)
        map_file(:,1) cell = get_blue_marble_paths()
    end

    %% Test cases
    methods (Test)
        function test_get(this, map_file)
            tile = GeoRasterTile(map_file);
            [img, R] = readgeoraster(map_file);

            % verify that accessing data at exact cell postings produces exact result,
            % so to do that we'll create a grid at integer pixel centers
            [xi, yi] = meshgrid(0:400:size(img,2), 0:400:size(img,1));
            xi(:,1) = 1;
            yi(1,:) = 1;
            [lat, lon] = R.intrinsicToGeographic(xi, yi);

            % verify that GeoRasterTile & direct indexing produce the same result
            ind = sub2ind(size(img), yi, xi);
            actual = tile.get(lat, lon);
            expected = single(img(ind));
            this.verifyEqual(actual, expected, 'AbsTol', this.tol);

            % create a sub-pixel grid in lat-lon coordinates and verify that interpolation
            % produces the same result
            lat = linspace(R.LatitudeLimits(1), R.LatitudeLimits(2), 29);
            lon = linspace(R.LongitudeLimits(1), R.LongitudeLimits(2), 29);
            [lat, lon] = ndgrid(lat, lon);
            [xi, yi] = R.geographicToIntrinsic(lat(:), lon(:));

            % we should get the same results using interp2
            actual = tile.get(lat(:), lon(:));
            expected = interp2(single(img), xi, yi, tile.interp, NaN);
            this.verifyEqual(actual, expected, 'AbsTol', this.tol);
        end

        function test_roi(this, map_file)
            tile = GeoRasterTile(map_file);
            [img, R] = readgeoraster(map_file);
            this.verifyTrue(true);

            % select a region at the center of the tile
            lat_bounds = R.LatitudeLimits + [30 -30];
            lon_bounds = R.LongitudeLimits + [30 -30];

            for version = 1:2
                switch version
                    case 1
                        % native resolution
                        args = {lat_bounds, lon_bounds};
                    case 2
                        % manually set the resolution (to an irrational #)
                        args = {lat_bounds, lon_bounds, pi*R.CellExtentInLatitude};
                end

                % use roi() to get a grid of points
                t = tic;
                [actual, lat, lon] = tile.roi(args{:});
                dt_tile = toc(t);
    
                % latitude should be decreasing, longitude increasing
                this.verifyTrue(issorted(lat), 'descend');
                this.verifyTrue(issorted(lon, 'ascend'));
    
                % back-solve for the same grid using native MATLAB/mapping toolbox
                [lat, lon] = ndgrid(lat, lon);
                [xi, yi] = R.geographicToIntrinsic(lat(:), lon(:));
    
                t = tic;
                expected = interp2(single(img), xi, yi, tile.interp, NaN);
                expected = reshape(expected, size(actual));
                dt_interp = toc(t);
                
                % we should get the same result less time than using interp2
                this.verifyEqual(actual, expected, 'AbsTol', this.tol);
                this.verifyLessThan(dt_tile, dt_interp);
            end
        end
    end

end

function paths = get_blue_marble_paths()
    root = fullfile(fileparts(mfilename('fullpath')),'..','data','blue-marble');
    files = dir(fullfile(root,'*.tif'));
    paths = fullfile(root, {files.name});
end
