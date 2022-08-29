% Example of fitting an asperity for inclusion in paper / presentations. 

clear;

set(groot, 'defaultAxesTickLabelInterpreter','default');  %Tex
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');

%% Load Data


addpath('../SURFACE/ASP_ID/')

load('../SURFACE/OUT/watershed_toprofile_14sep21.mat', 'Xmat_tmp', 'Ymat_tmp', 'Zmat_tmp', 'peakMat_tmp');
input_valid = load('../SURFACE/OUT/asperities_14sep21_R1.mat');

%% Find a Valid Example

start_ii = 100000;
num_check = 500;


valid = false(1, start_ii+num_check);
radii = zeros(2, start_ii+num_check);
xyzloc = zeros(start_ii+num_check, 3);
alpha = zeros(1, start_ii+num_check);

%     warning('Fitting asperities in serial');
parfor ii = start_ii:(start_ii+num_check)

    [radii(:, ii), xyzloc(ii, :), alpha(ii), valid(ii), ...
        ~, ~, ~, ~, ~, ~, ~] = FIT_ELLIPSOID3D([Xmat_tmp(peakMat_tmp == ii), ...
                                                Ymat_tmp(peakMat_tmp == ii), ...
                                                Zmat_tmp(peakMat_tmp == ii)]);
end



%% Profile Asperity 
inds = find(valid);
index = 3;
ii = inds(end-index);

    
[radii_ii, xyzloc_ii, alpha_ii, valid_ii, center, axesq, ...
    csys, Serr, Bmat, Lvec, Csca, mxyz, Tm] ...
        = FIT_ELLIPSOID3D([Xmat_tmp(peakMat_tmp == ii), ...
                            Ymat_tmp(peakMat_tmp == ii), ...
                            Zmat_tmp(peakMat_tmp == ii)]);

assert(valid_ii, 'Do not plot a bad asperity for your example!');

%% 

% Plot Data
X = Xmat_tmp(peakMat_tmp == ii);
Y = Ymat_tmp(peakMat_tmp == ii);
Z = Zmat_tmp(peakMat_tmp == ii);

% Plot Surface Points
Nplot = 100;
offset_factor = 0.2;
Xrange = [min(X)-offset_factor*range(X), max(X)+offset_factor*range(X)];
Yrange = [min(Y)-offset_factor*range(Y), max(Y)+offset_factor*range(Y)];

Xsurf = linspace(Xrange(1), Xrange(2), Nplot);
Ysurf = linspace(Yrange(1), Yrange(2), Nplot);
[Xsurf, Ysurf] = meshgrid(Xsurf, Ysurf);

% XYsurf = Tm(1:2, 1:2)\[Xsurf(:) - mxyz(1), Ysurf(:) - mxyz(2)]';
XYsurf = [Xsurf(:) - mxyz(1), Ysurf(:) - mxyz(2)]';

% Quadratic coefficients: a z^2 + bz + c
a = Bmat(3,3);
b = Lvec(3);
c = sum(XYsurf.*(Bmat(1:2, 1:2)*XYsurf)) + Lvec(1:2)'*XYsurf + Csca;

Zsurf = (-b + sqrt(b.^2 - 4*a.*c) ) ./ (2 * a);

Zsurf(imag(Zsurf) ~=0) = NaN;

XYZsurf = [Xsurf(:), Ysurf(:), Zsurf(:)+mxyz(3)];

% Actual Plot

unit_convert = 10^6;

figure('Position', [250, 250, 700, 500]);
plot3((X- center(1))*unit_convert, ...
        (Y- center(2))*unit_convert, ...
        (Z- center(3))*unit_convert, ...
        'ko', 'MarkerFaceColor', 'k'); 
hold on 
surf(reshape(XYZsurf(:,1)- center(1), Nplot, Nplot)*unit_convert, ...
        reshape(XYZsurf(:,2)- center(2), Nplot, Nplot)*unit_convert, ...
        reshape(XYZsurf(:,3)- center(3), Nplot, Nplot)*unit_convert, ...
        'EdgeColor', 'none', 'FaceAlpha', 0.6);
grid on;
xlim((Xrange-center(1))*unit_convert)
ylim((Yrange-center(2))*unit_convert);
Zrange = [min(Z)-0.1*range(Z), max(max(Z), max(XYZsurf(:, 3)))];
zlim((Zrange-center(3))*unit_convert);

pbaspect([range(Xrange), range(Yrange), min(range(Xrange), 0.5*range(Yrange))]);


xlabel('$x$ [$\mu$m]')
ylabel('$y$ [$\mu$m]')
zlabel('$z$ [$\mu$m]')

view([-63, 21]);

font_size = 14;
set(gca, 'FontSize', font_size);

set(gcf, 'Renderer', 'painters');
drawnow;

% print('./OUTPUTS/paper_asperity_fit.eps', '-depsc', '-r400');
% print('./OUTPUTS/paper_asperity_fit.png', '-dpng', '-r600');
