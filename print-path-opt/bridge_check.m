function [is_bridge]=bridge_check(edges,v)
current_degree__v=sum(edges==v,"all");

if current_degree__v==0
    is_bridge=1;
else
    is_bridge=0;

end
end