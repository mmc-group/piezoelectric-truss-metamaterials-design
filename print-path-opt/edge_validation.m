function [groups]=edge_validation(coords,edges_global, d_nozzle)
ind=0;
group_id=1;
edges=edges_global;
edge_ids=[1:size(edges,1)];
while ind==0
    current_group=[];
    for i=1:length(edges)
        edge_A=edges(i,:);
        A=[coords(edge_A(1),:);coords(edge_A(2),:);];

        x1=A(1,1);y1=A(1,2);z1=A(1,3);x2=A(2,1);y2=A(2,2);z2=A(2,3);

        for j=1:length(edges)
            if j==i
                valid=true;
            else
                edge_B=edges(j,:);
                B=[coords(edge_B(1),:);coords(edge_B(2),:);];

                x3=B(1,1);y3=B(1,2);z3=B(1,3);x4=B(2,1);y4=B(2,2);z4=B(2,3);
                angle_AB=Angle_bw_lines(A,B);
                Par_mat=[(x1-x2) (x3-x4)
                    (y1-y2) (y3-y4)];
                if abs(det(Par_mat))<1e-8
                    [num]=endpointsofAonB(A,B);
                    if (num)==0
                        valid=true;
                    elseif num==1
                        valid=true;
                    else
                        if (((z1+z2)/2)<((z3+z4)/2))
                            valid=true;
                        else
                            valid=false;
                        end
                    end

                else

                    [Px,Py] = linexline(A(:,1), A(:,2), B(:,1), B(:,2), 0);
                    %     P1=det([x1 y1;x2,y2]);P2=det([x1 1;x2 1]);P3=det([x3 y3;x4 y4]);
                    %     P4=det([x3 1;x4 1]);P5=det([y1 1;y2 1]);P6=det([y3 1;y4 1]);
                    %     P7=det([x3 1;x1 1]);
                    %     Px=(det([P1 P2;P3 P4]))/(det([P2 P5;P4 P6]));
                    %     Py=(det([P1 P5;P3 P6]))/(det([P2 P5;P4 P6]));
                    [checkPt_A, onEnd_A] = checkPointOnSegment([x1,y1;x2,y2], [Px,Py], 0);
                    [checkPt_B, onEnd_B] = checkPointOnSegment([x3,y3;x4,y4], [Px,Py], 0);
                    if checkPt_A==true && checkPt_B==true
                        if (x1-x3)==0 && (y1-y3)==0 || (x2-x3)==0 && (y2-y3)==0||...
                                (x1-x4)==0 && (y1-y4)==0 || (x2-x4)==0 && (y2-y4)==0
                            if angle_AB<30
                                if (((z1+z2)/2)<((z3+z4)/2))
                                    valid=true;
                                else
                                    valid=false;
                                end
                            else
                                valid=true;
                            end

                        else

                            [Az]=intersection_point(A,Px,Py);
                            [Bz]=intersection_point(B,Px,Py);

                            if Az<Bz
                                valid=true;
                            else
                                valid=false;
                            end

                        end
                    else
                        % if checkPt_A==true && checkPt_B==false
                        %     valid=true;
                        % elseif checkPt_A==false && checkPt_B==true
                        %     valid=true;
                        % else
                        [dist,~,~]=DistBetween2Segment(A(1,:),A(2,:),B(1,:),B(2,:));

                        if dist<d_nozzle
                            valid=false;
                        else
                            valid=true;
                        end
                        % end
                    end

                end
            end
            if valid==false
                break
            else
            end
            if j~=i

            end
        end
        % % if valid==false
        % plot3(A(:,1),A(:,2),A(:,3),'b')
        % hold on
        % plot3(B(:,1),B(:,2),B(:,3),'r')
        % % pause
        % %     end
        if valid==true
            current_group(end+1,1)=i;
        end
    end

    %%% split the current group into two, in case there is a disconnected 
    %%% floating sub path (not considered in the alogorithm by Weeks et al)
    if isempty(current_group)==false
        [current_group_split]=split_current_group(current_group, edges,coords);
        current_group=current_group_split{1,1};
        global_group=edge_ids(current_group);
        groups{group_id}=global_group;
        edges(current_group,:)=[];
        edge_ids(current_group)=[];
        group_id=group_id+1;
    else
        if i==length(edges)
            disp("error: no edges can be printed")
            %         return
        end
    end
    if size(edges,1)==0
        ind=1;
    end

end
end