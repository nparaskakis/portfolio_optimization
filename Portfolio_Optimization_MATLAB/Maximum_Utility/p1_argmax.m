%%  Advanced Topics in Convex Optimization
%%  Nikolaos Paraskakis - A.M.: 2018030027



function argmax = p1_argmax(w,d)
    
    
    
    %% Description
    %
    %  This is a function that computes the argmax operator
    %  for FDPG that corresponds to the problem of exercise 12b.
    
    
    
    %% Input
    %
    %  w -> Vector of size nxm.
    %  d -> Vector of size nx1.
    
    
    
    %% Output
    %
    %  argmax -> Argmax operator needed for FDPG algorithm.
    
    
    
    %% Function's body
    
    argmax = sum(w,2) + d;
    
    
    
end