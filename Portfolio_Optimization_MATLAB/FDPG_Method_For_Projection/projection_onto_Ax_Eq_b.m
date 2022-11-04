%% Advanced Topics in Convex Optimization
%% Nikolaos Paraskakis - A.M.: 2018030027



function projection = projection_onto_Ax_Eq_b(x,A,b)
    
    
    
    %% Description
    %
    %  This is a function that computes the projection
    %  of a point x onto set S = {x : Ax = b}.
    
    
    
    %% Input
    %
    %  x -> Vector of size nx1, that is going to be projected onto set S.
    %  A -> Matrix of size mxn.
    %  b -> Vector of size mx1.
    
    
    
    %% Output
    %
    %  projection -> The projection of input point x onto set S.
    
    
    
    %% Function's body
    
    projection = x - A'*((A*A')\(A*x-b));
    
    
    
end