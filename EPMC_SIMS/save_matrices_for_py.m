% Script for Resaving System Matrices for Python Simulations
% 
% Notes: 
%   1. Matrices are saved in a new .mat file same as previous save
%   2. Mesh information is saved just as matrices without a structure 
%       since structures can be harder to read when loading in python
%   3. Quadrature mapping matrices are saved so MATLAB scripts for those
%       do not need to be recreated.
%   4. This script does not resave all of the mapping matrices since some
%       are very large and are not needed in the planned analyses

clear;
%% Routines locations

addpath('../ROUTINES/')
addpath('../ROUTINES/FEM/')

%% File To Load

input_name = '../FJSIMS/ROMS/ROM_U_232ELS';
output_name = 'ROM_U_232ELS4py.mat';
output_dir = '../FJSIMS/ROMS/for_py';

%% Load Baseline Matrices

load('../FJSIMS/ROMS/ROM_U_232ELS');

%% Generate Quadrature Mapping

MESH.dpn = 3;

% %% ZTE Quadrature Matrices
No   = 1;                  % Number of GLQ integration points in each direction
[Qm,Tm] = ZTE_ND2QP(MESH,No);

%% Pull out Mesh info

% Each row represents an element. First column is element index. Next 3-4
% are node numbers
tri_elem_nodes = MESH.Tri; 
quad_elem_nodes = MESH.Quad;

node_coords = MESH.Nds;

%% Resave mat File

mkdir(output_dir)

output_loc = fullfile(output_dir, output_name);

save(output_loc, 'M', 'K', 'R', 'Fv', 'L', 'Qm', 'Tm', 'tri_elem_nodes', ...
                'quad_elem_nodes', 'node_coords', '-v7')






