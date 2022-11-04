%%  Advanced Topics in Convex Optimization
%%  Nikolaos Paraskakis - A.M.: 2018030027



function proj = p1_proj_1(x,A,b)
    
    
    
    %% Description
    %
    %  This is a function that computes the proj_1 operator
    %  for FDPG that corresponds to the problem of exercise 12b.
    
    
    
    %% Input
    %
    %  x -> Vector of size nx1.
    %  A -> Matrix of size mxn.
    %  b -> Vector of size mx1.
    
    
    
    %% Output
    %
    %  proximal -> Proximal operator needed for FDPG algorithm.
    
    
    
    %% Function's body
    
    proj = projection_onto_Ax_Eq_b(x,A,b);
    
    
    
end