function vec_map=fnalize_vec_map(vec_map)

for i=1:length(vec_map)

    if vec_map(i)==-1
        vec_map(i)=i;
    end
    
end

end