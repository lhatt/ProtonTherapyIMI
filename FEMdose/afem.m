%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adaptive FEM in 2d for the elliptic problem
%  -Delta u = f   in Omega
%       u = gD  in Gamma_D
%   du/dn = gN  in Gamma_N
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clear all
%close all

%Change the parameters here for your prior sample
%Divide the domain into 3 chunks of tissue
% diffusion over region 1 (skin)
d1 = 1;
% diffusion over region 2 (bone)
d2 = 0.1;
% diffusion over region 3 (tumour)
d3 = 0.5;

% we first read all the initial data from the file init_data.m
init_data
% at this point we have filled structures
%       prob_data  and  adapt
% that are also global variables, and thus
% visible in all functions

% we now define three more global variables for the mesh to save memory
global  mesh  uh  fh

%% initialize mesh and finite element functions
% read the mesh from 'domain'
mesh.elem_vertices      = load([domain '/elem_vertices.txt']);
mesh.elem_neighbours    = load([domain '/elem_neighbours.txt']);
mesh.elem_boundaries    = load([domain '/elem_boundaries.txt']);
mesh.vertex_coordinates = load([domain '/vertex_coordinates.txt']);
mesh.n_elem     = size(mesh.elem_vertices, 1);
mesh.n_vertices = size(mesh.vertex_coordinates, 1);
uh = zeros(mesh.n_vertices, 1);
fh = zeros(mesh.n_vertices, 1);

% define the exact solution. If known it's stored through the boundary conditions
u_exact = prob_data.gD;

if (global_refinements)
    mesh.mark = global_refinements*ones(mesh.n_elem,1);
    refine_mesh;
end

% init error and estimator
err_H1_old = -1.0;
est_old = -1.0;
iter_counter = 1;

% Assemble and solve the system matrix
[t] = assemble_and_solve(d1,d2,d3);

% Calls the output script
output

XYZ=[mesh.vertex_coordinates(:,1), mesh.vertex_coordinates(:,2), ...
uh];
index=find(XYZ(:,2)==0.5);
XYZ_det=XYZ(index,:);

x=0:.01:1;
y=interp1(XYZ_det(:,1),XYZ_det(:,3),x);