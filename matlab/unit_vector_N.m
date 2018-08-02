function [ u_N ] = unit_vector_N( u_BC, u_AB )
%Calculates unit normal vector which is perpendicular to plane ABC

    cross_product = cross( u_BC , u_AB );
    u_N= cross_product / norm(cross_product);
    
end

