%%  Advanced Topics in Convex Optimization
%%  Nikolaos Paraskakis - A.M.: 2018030027



function str = print_vector_inline(x)
    
    
    
    %% Description
    %
    %  This is a function that puts a vector x
    %  with its components in a string.
    
    
    
    %% Input
    %
    %  x -> The vector x that is going to be converted to string.
    
    
    
    %% Output
    %
    %  str -> The string with the vector x in a line.
    
    
    
    %% Function's body
    
    if (length(x) > 1)
        str = sprintf('[');
        str = strcat(str,sprintf('%.4f, ',x(1:end-1)));
        str = strcat(str,sprintf('%.4f]',x(end)));
    else
        str = sprintf('%.4f',x);
    end
    
    
    
end