function [DD]= dof(node_number,ic,ndof)
DD=(node_number-1)*ndof+ic;
end