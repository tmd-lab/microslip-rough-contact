% Function for Testing the Correctness of Asperity Ellipsoid Fitting
% Routine

clear;

addpath('../SURFACE/ASP_ID/')

%% Ellipse 

ax = 0.5;
ay = 0.3;
az = 1.2;

cx = 2;
cy = 3;
cz = 5;

th = linspace(pi/4, pi/2, 20);
phi = linspace(0, 2*pi, 20); 
[th, phi] = meshgrid(th, phi);

x = ax*cos(th(:)).*cos(phi(:)) + cx ;
y = ay*cos(th(:)).*sin(phi(:)) + cy ;

% z = cz + real(az*sqrt(1-((x/ax).^2+(y/ay).^2)));  % Ellipse formula 
z = az*sin(th(:)) + cz ;


%%%% Debugging Block to look at some example points
% tmp = load('../SURFACE/TMP_DebugPoints', 'XYZ');
% x = tmp.XYZ(:, 1);
% y = tmp.XYZ(:, 2);
% z = tmp.XYZ(:, 3);

figure(1)
clf()
plot3(x, y, z, 'o')
grid on 

% Routine Call 
xyz = [x y z];

[~, ~, ~, ~, center, axesq, csys, Serr, Bmat, Lvec, Csca, mxyz] = FIT_ELLIPSOID3D(xyz);

[xx, yy] = meshgrid(linspace(min(xyz(:,1)), max(xyz(:,1)), 100), ...
                linspace(min(xyz(:,2)), max(xyz(:,2)), 100));
            
a = Bmat(3,3);
b = Lvec(3);
c = diag([xx(:)-mxyz(1), yy(:)-mxyz(2)]*Bmat(1:2,1:2)*[xx(:)-mxyz(1) yy(:)-mxyz(2)]') + [xx(:)-mxyz(1), yy(:)-mxyz(2)]*Lvec(1:2) + Csca;

zz1 = (-b+sqrt(b^2-4*a*c))/(2*a)+mxyz(3);
zz2 = (-b-sqrt(b^2-4*a*c))/(2*a)+mxyz(3);


figure(1)
clf()
plot3(x, y, z, 'ko', 'MarkerFaceColor', 'k'); hold on 
surf(xx, yy, reshape(real(zz1), 100, 100), 'EdgeColor', 'none', 'FaceAlpha', 0.5);
grid on 
xlabel('X Coordinate')
ylabel('Y Coordinate')
zlabel('Z Coordinate')

%% Test Several Inputs and Angles

rng(42); % For repeatability

% Test Parameters
Ntest = 100;
radRange = [0.1, 10]; % min, max
locRange = [-10, 10]; % min, max

alpha_all = linspace(-2*pi, 2*pi, Ntest);

add_noise = [0, 0.001];

% Randomly Generated Test Parameters
axes_all = rand(Ntest, 3)*range(radRange) + radRange(1);
centers_all = rand(Ntest, 3)*range(locRange) + locRange(1);

% Surface parameterization
th = linspace(pi/4, pi/2, 20);
phi = linspace(0, 2*pi, 20); 
[th, phi] = meshgrid(th, phi);

% Tolerances:
localEllipsoidTol = 1e-14*radRange(2);
locationTol       = 1e-12*range(locRange);
radTol            = 1e-12*radRange(2)^2/radRange(1);
alphaTol          = 1e-12*pi/2;
noiseTolFactor    = 100;

fprintf('\n\n')
fprintf('Starting checking:\n');
fprintf('Alpha ranges from %f to %f \n', min(alpha_all), max(alpha_all));
fprintf('Axis lengths ranges from %f to %f \n', radRange(1), radRange(2));
fprintf('Locations ranges from %f to %f \n', locRange(1), locRange(2));
fprintf('Note: adding noise is expected to cause some failures due not appropriately modifying the radius tolerance\n');

for noise = add_noise
    fprintf('Testing with noise level %f, adding %f to all tolerances \n', noise, noise*noiseTolFactor);
for ii = 1:Ntest
    
    alpha = alpha_all(ii);
    axesxyz = axes_all(ii, :);
    
    axesxyz(1:2) = sort(axesxyz(1:2));
    
    xyzLocal = [axesxyz(1)*cos(th(:)).*cos(phi(:)), ...
                axesxyz(2)*cos(th(:)).*sin(phi(:)), ...
                axesxyz(3)*sin(th(:))];
    
    % Verify that the local ellipse is valid
    ellipsoidError = max(abs( xyzLocal(:, 1).^2/axesxyz(1)^2 ...
                                + xyzLocal(:, 2).^2/axesxyz(2)^2 ...
                                + xyzLocal(:, 3).^2/axesxyz(3)^2 ...
                                -1 ));

    if(ellipsoidError > localEllipsoidTol)
        warning('Local points differ from proper equation by %s', ellipsoidError);
    end
    
    % Transform
    Ploc_to_glob = [cos(alpha), -sin(alpha), 0;
                    sin(alpha),  cos(alpha), 0;
                             0,           0, 1];
    
    xyz = xyzLocal * Ploc_to_glob' + centers_all(ii, :);
    
    xyz = xyz + rand(size(xyz))*2*noise - noise;
    
    % Actually do the fitting
    [radii_out, xyzloc_out, alpha_out, valid, center, axesq, csys, ...
        Serr, Bmat, Lvec, Csca] = FIT_ELLIPSOID3D(xyz);

    % Error Checking on Results
    loc_error = max(abs( xyzloc_out - centers_all(ii, :) - [0, 0, axesxyz(3)] ));

    if(loc_error > locationTol + noise*noiseTolFactor)
        warning('Peak location differed from expected by %s', loc_error);
    end
    
    rad_error = max(abs( axesxyz(1:2).^2/axesxyz(3) - radii_out ));
    
    if(rad_error > radTol + noise*noiseTolFactor)
        
        rad_error_perc = abs( axesxyz(1:2).^2/axesxyz(3) - radii_out )./radii_out;
    
        warning('Principal radii differed from expected by %s or about %f percent', rad_error, max(rad_error_perc)*100);
    end
    
    alpha_error = alpha - alpha_out;
    alpha_error = alpha_error - pi*round(alpha_error / pi);
    
    if(alpha_error > alphaTol + noise*noiseTolFactor)
        warning('Rotation Angle differed from expected by %s', alpha_error);
    end
    
end
end

fprintf('Checking complete.\n');

%% Verify In Domain Validity Checks

ax = 0.5;
ay = 0.3;
az = 1.2;

cx = 2;
cy = 3;
cz = 5;

phiRange = [0, 2*pi;
            -0.9*pi/2, 0.9*pi/2;
            0.1*pi/2, 0.95*pi;
            1.1*pi/2, 2*pi-1.1*pi/2;
            1.05*pi, 2*pi-.1*pi/2];
            
expectedRes = [true; false; false; false; false];

fprintf('Checking domain validity:\n');

for ii = 1:size(phiRange, 1)

    th = linspace(pi/4, 0.9*pi/2, 20);
    phi = linspace(phiRange(ii, 1), phiRange(ii, 2), 20); 
    [th, phi] = meshgrid(th, phi);

    x = ax*cos(th(:)).*cos(phi(:)) + cx ;
    y = ay*cos(th(:)).*sin(phi(:)) + cy ;

    z = az*sin(th(:)) + cz ;

    % Routine Call 
    xyz = [x y z];
    [radii, xyzloc, alpha, valid, center, axesq, csys, ...
        Serr, Bmat, Lvec, Csca] = FIT_ELLIPSOID3D(xyz);
    
    if( valid ~= expectedRes(ii) )
        warning('Test %u unexpectedly returned valid = %u', ii, valid);
    end

end

fprintf('Finished domain checks.\n');

%% Verify Number of Points Validity Checks

ax = 0.5;
ay = 0.3;
az = 1.2;

cx = 2;
cy = 3;
cz = 5;

%%%% Test the following two point Ranges:
pointRangeTest = {1:9, 1:8, 1:7};

fprintf('Starting test: number of points used \n');

for ii = 1:length(pointRangeTest)

    pointRange = pointRangeTest{ii};

    th = linspace(pi/4, 0.9*pi/2, 3);
    phi = linspace(0, 2*pi, 4); 
    [th, phi] = meshgrid(th, phi);

    x = ax*cos(th(:)).*cos(phi(:)) + cx ;
    y = ay*cos(th(:)).*sin(phi(:)) + cy ;
    z = az*sin(th(:)) + cz ;

    % Routine Call 
    xyz = [x(pointRange) y(pointRange) z(pointRange)];
    [radii, xyzloc, alpha, valid, center, axesq, csys, Serr, Bmat, Lvec, Csca] = FIT_ELLIPSOID3D(xyz);

    if( valid ~= ( length(pointRange) >= 9) )
        warning('Unexpectedly returned valid = %u with %u points', valid, length(pointRange) );
    end
end

fprintf('Finished Test \n');


%% Verify Only Returns local maxima

ax = 0.5;
ay = 0.3;
az = 1.2;

cx = 2;
cy = 3;
cz = 5;

%%%% Test the following two point Ranges:
zmult = [-1, 1];


fprintf('Starting test: local maxima \n');

for ii = 1:length(zmult)


    th = linspace(pi/4, 0.9*pi/2, 3);
    phi = linspace(0, 2*pi, 4); 
    [th, phi] = meshgrid(th, phi);

    x = ax*cos(th(:)).*cos(phi(:)) + cx ;
    y = ay*cos(th(:)).*sin(phi(:)) + cy ;
    z = zmult(ii)*az*sin(th(:)) + cz ;

    % Routine Call 
    xyz = [x y z];
    [radii, xyzloc, alpha, valid, center, axesq, csys, Serr, Bmat, Lvec, Csca] = FIT_ELLIPSOID3D(xyz);

    if( valid ~= (zmult(ii) > 0) )
        warning('Unexpectedly returned valid = %u, when z multiplied by %u', valid, zmult(ii));
    end
end

fprintf('Finished Test \n');
