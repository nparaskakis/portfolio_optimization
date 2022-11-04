%%  Advanced Topics in Convex Optimization
%%  Nikolaos Paraskakis - A.M.: 2018030027



function proj = p1_proj_2(x,l,u)
    
    
    
    %% Description
    %
    %  This is a function that computes the proj_2 operator
    %  for FDPG that corresponds to the problem of exercise 12b.
    
    
    
    %% Input
    %
    %  x -> Vector of size nx1.
    %  l -> Lower bound of x.
    %  u -> Upper bound of x.
    
    
    
    %% Output
    %
    %  proximal -> Proximal operator needed for FDPG algorithm.
    
    
    
    %% Function's body
    
    proj = projection_onto_box(x,l,u);
    
    
    
end