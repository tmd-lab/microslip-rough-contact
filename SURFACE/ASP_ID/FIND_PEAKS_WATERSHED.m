function [peakMat] = FIND_PEAKS_WATERSHED(Zmat, Holes, settings)
% Function for identifying points that are on asperities
%
% Inputs:
%   Zmat - matrix of height data
%   Holes - true/false matrix with trues where the holes are. The Hole
%           areas are ignored for most parts and are excluded from being in
%           peakMat output
%   settings - structure with settings to be used
%       .cleanLength - number of points to include in imerode and imdilate
%                       for initial cleaning.
%
% Outputs:
%   peakMat - labeled segmented matrix for potential asperities.
% 
% Implementation of Algorithm from: 
%       Wen, Tan, Zhou, Li, 2020,  A Reconstruction and Contact Analysis 
%       Method of Three-Dimensional Rough Surface Based on Ellipsoidal 
%       Asperity
%
%       Implementation provided by Nidish Balaji from his PhD Dissertation.

    %% Fill NaNs from holes to prevent errors
    
    % Set holes below any scan data points
    Zmat(Holes) = min(Zmat(:));

    %% Morphological open for cleaning up the data
    se = strel('square', settings.cleanLength);
    zqo = imerode(Zmat, se);
    zqo = imdilate(zqo, se);

    [zqe,~,Gv,Gh] = edge(zqo);  % Edge command for Sobel Gradients

    Gradz = sqrt(Gv.^2+Gh.^2);

    %% Surface Parameters
    zqs = sort(zqo(:), 'descend');
    Sz = mean(zqs(1:5)-zqs(end:-1:end-4));  % Sz(ISO): 10-point height

    %% Apply watershed segmentation
    peakMat = watershed(Gradz);
    
    
    %% Remove Holes:
    
    peakMat(Holes) = 0;
    
    
end