function [t] = assemble_and_solve(d1,d2,d3)
% function assemble_and_solve
%   assemble the discrete system and solve it
%   all the data is in the global variables
%   mesh  prob_data
%   the right-hand side is stored in 
%   the global vector  fh
%   and the discrete solution is stored in
%   the global vector  uh

global mesh prob_data uh fh
% In order to simplify, we create the dirichlet and
% neumann variables as we did in the fixed mesh case.
% That is:
%   dirichlet is a vector containing the Dirichlet vertices
%   neumann   is a matrix containing the Neumann   segments

[dirichlet, neumann, interior] = get_dirichlet_neumann;

n_vertices = mesh.n_vertices;
n_elem = mesh.n_elem;

n_interior = n_vertices - length(dirichlet);

n_bound = length(dirichlet);

if length(interior) ~= n_interior
    fprintf(1,'Not assigning interior dofs correctly!');
    pause;
end

S = sparse(n_vertices, n_vertices);

fh = zeros(n_vertices, 1);

g = zeros(n_bound, 1);

% gradients of the basis functions in the reference element
grd_bas_fcts = [ -1 -1 ; 1 0 ; 0 1 ]' ;

% We loop through the elements of the mesh,
% and add the contributions of each element to the matrix A
% and the right-hand side fh

% At each element we use the cuadrature formula which uses 
% the function values at the midpoint of each side:
% \int_T  f  \approx  |T| ( f(m12) + f(m23) + f(m31) ) / 3.
% This formula is exact for quadratic polynomials
tic




for el = 1 : n_elem

    v_elem = mesh.elem_vertices( el, : );

    
    v1 = mesh.vertex_coordinates( v_elem(1), :)' ; % coords. of 1st vertex of elem
    v2 = mesh.vertex_coordinates( v_elem(2), :)' ; % coords. of 2nd vertex of elem
    v3 = mesh.vertex_coordinates( v_elem(3), :)' ; % coords. of 3rd vertex of elem
    
    m12 = (v1 + v2) / 2; % midpoint of side 1-2
    m23 = (v2 + v3) / 2; % midpoint of side 2-3
    m31 = (v3 + v1) / 2; % midpoint of side 3-1

    % derivative of the affine transformation from the reference
    % element onto the current element
    B = [ v2-v1  v3-v1 ];
    
    % element area
    el_area = abs(det(B)) * 0.5;
    
    % evaluation of f at the quadrature points

    f12 = feval(prob_data.f,m12);
    f23 = feval(prob_data.f,m23);
    f31 = feval(prob_data.f,m31);
    
    % computation of the element load vector
    f_el = [ (f12+f31)*0.5 ; (f12+f23)*0.5 ; (f23+f31)*0.5 ] * (el_area/3);
    
    % contributions added to the global load vector
    fh( v_elem ) = fh( v_elem ) + f_el;

    Binv = inv(B);
    bas_fcts = [0.5 0.5 0.5];
    % computation of the element matrix
    el_mat = el_area * ( grd_bas_fcts' * ( Binv*Binv') * grd_bas_fcts);
        
    %  barycentre of the element
    mp = 1/3*(m12+m23+m31);
    % assign diffusion as a piecewise contant
    if mp(1) < 1/3
        diff = d1*prob_data.epsilon;
    elseif mp(1) < 2/3
        diff = d2*prob_data.epsilon;
    else
        diff = d3*prob_data.epsilon;
    end
    
    ad_mat = el_area * (prob_data.advection*(grd_bas_fcts'*Binv)')'*bas_fcts;
    
    % contributions added to the global matrix
    S( v_elem, v_elem ) = S( v_elem, v_elem ) + diff*el_mat + ad_mat;
    
end

% ensure the Dirichlet boundary conditions are set
for i = 1:length(dirichlet)
    g(i) = feval(prob_data.gD, mesh.vertex_coordinates(dirichlet(i), :) );
    S(dirichlet(i),:) = 0;
    S(dirichlet(i),dirichlet(i)) = 1;
    fh(dirichlet(i)) = g(i);
    
end

% solve the system
uh = S\fh;

t = toc;

end