%% Advanced Topics in Convex Optimization
%% Nikolaos Paraskakis - A.M.: 2018030027



function y = projection_onto_box(x,l,u)
    
    
    
    %% Description
    %
    %  This is a function that computes the
    %  projection of a point x onto Box[l,u].
    
    
    
    %% Input
    %
    %  x -> The vector whose projection we are looking for.
    %  l -> Lower bound of the box.
    %  u -> Upper bound of the box.
    
    
    
    %% Output
    %
    %  projection -> The projection of input vector x
    %                onto the Box[l,u].
    
    
    
    %% Function's body
    
    % Length of vector x.
    n = length(x);
    
    % Create vector y in which the projection of x will be saved.
    y = zeros(size(x));
    
    % Compute the projection of point x onto Box[l,u] component-wise.
    for i = (1:1:n)
        y(i) = min([max([x(i) l(i)]) u(i)]);
    end
    
    
    
end