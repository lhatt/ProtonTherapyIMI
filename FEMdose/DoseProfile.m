function [x,y]=DoseProfile(d1,d2,d3)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adaptive FEM in 2d for the elliptic problem
%  -Delta u = f   in Omega
%       u = gD  in Gamma_D
%   du/dn = gN  in Gamma_N
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clear all
%close all
global max_iterations
%Change the parameters here for your prior sample
%Divide the domain into 3 chunks of tissue
% diffusion over region 1 (skin)
%d1 = 1;
% diffusion over region 2 (bone)
%d2 = 0.1;
% diffusion over region 3 (tumour)
%d3 = 0.5;

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

while (1)
    
    % Assemble and solve the system matrix
    [t] = assemble_and_solve(d1,d2,d3);
        
    est(iter_counter) = estimate(prob_data, adapt);

    err_H1(iter_counter) = H1_err(mesh.elem_vertices, mesh.vertex_coordinates, uh, prob_data.grd_u_exact);

    ndof(iter_counter) = mesh.n_vertices;
    
    if err_H1_old < 0
        EOC_H1(iter_counter) = 0;
        EOC_est(iter_counter) = 0;
    else
        EOC_H1(iter_counter) = log(err_H1_old / err_H1(iter_counter) ) / log(2);
        EOC_est(iter_counter) = log(est_old / est(iter_counter) ) / log(2);
    end
    
    %fprintf(1,'n_dofs: %5d\n',mesh.n_vertices);
    %fprintf(1,'n_elements: %5d\n',mesh.n_elem);
    %fprintf(1,'Assemble/solve computational time[s]: %.2f \n',t);

    err_H1_old = err_H1(iter_counter);
    est_old = est(iter_counter);
    
    if ((iter_counter >= adapt.max_iterations) || (est(iter_counter) < adapt.tolerance))
        break;
    else
        mark_elements(adapt);
        refine_mesh;
    end
    
    iter_counter = iter_counter + 1;
    
end

% Calls the output script
%output

XYZ=[mesh.vertex_coordinates(:,1), mesh.vertex_coordinates(:,2), ...
uh];
index=find(XYZ(:,2)==0.5);
XYZ_det=XYZ(index,:);

x=0:.0001:1;
y=interp1(XYZ_det(:,1),XYZ_det(:,3),x);