function vec_map=fill_vec_map(DD,val,vec_map)
if vec_map(DD)==-1
    vec_map(DD)=val;
else
    disp('error : coupling already assigned to the node')
    return

end