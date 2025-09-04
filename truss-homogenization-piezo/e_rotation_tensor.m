function [e_mat_rot]=e_rotation_tensor(e,R)
% % % % e = [6 X 3] piezoelectric matrix of the material in original basis
% % 
% % % % R = rotation matrix containing dircetion cosines for the rotation
% % 
% % 
% % % % e_mat_rot = piezoelectric matrix in the rotated basis
[e_base_tensor]=voigt_to_tensor_e(e);  % converting 6X3 e matrix to 3X3X3 tensor
e_tensor_rot=zeros(3,3,3);

% % Rotating the e tensor (e_i_j_k = R_i_p * R_j_q * R_k_r * e_p_q_r)
for i=1:3
    for j=1:3
        for k=1:3
            for p=1:3
                for q=1:3
                    for r=1:3
                        
                    e_tensor_rot(i,j,k)=e_tensor_rot(i,j,k)+R(i,p)*R(j,q)*R(k,r)*e_base_tensor(p,q,r);
                    
                    end
                end
            end
        end
    end
end

% %  converting 3X3X3 e tensor back to 6X3 e matrix
[e_mat_rot]=e_tensor_to_voigt(e_tensor_rot);

end

function [e_tensor]=voigt_to_tensor_e(e_voigt)
map=[1 6 5
    6 2 4
    5 4 3];
e_tensor=zeros(3,3,3);
for i=1:3
    for j=1:3
        for k=1:3
        c=map(i,j);
        e_tensor(i,j,k)=e_voigt(c,k);
        end
    end
end

end


function [e_mat]=e_tensor_to_voigt(e_tensor)

inv_map=[1,1;2,2;3,3;2,3;1,3;1,2];
e_mat=zeros(6,3);
for i=1:6
    for j=1:3
        [v]=inv_map(i,:);
        e_mat(i,j)=e_tensor(v(1),v(2),j);
    end
end

end


