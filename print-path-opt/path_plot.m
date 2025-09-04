function path_plot (Coords,global_path)
v=VideoWriter('path.avi');
v.FrameRate=4;
open(v)
figure('units','pixels','position',[0 0 1440 1080])
travel_motion=-1;
for gg=1:length(global_path)
    Connectivity=global_path{gg};
    travel_motion=travel_motion+1;
    for ii=1:size(Connectivity,1)
        plot3(Coords((Connectivity(ii,:)),1),Coords(Connectivity(ii,:),2),...
            Coords(Connectivity(ii,:),3),'Color',[0.2941,0.4902,0.6784],'LineWidth',3)
        if ii>1
            if Connectivity(ii-1,2)~=Connectivity(ii,1)
                travel_motion=travel_motion+1;
            else
            end
        end

        title("Path "+travel_motion,'FontSize',16)
        set(gcf,'color','w');
        %%% displaying the node numbers
        node1=Connectivity(ii,1);
        n1=num2str(node1);
        node2=Connectivity(ii,2);
        n2=num2str(node2);
        %     text(Coords(node1,1),Coords(node1,2),Coords(node1,3),n1,'FontSize',8)
        %     text(Coords(node2,1),Coords(node2,2),Coords(node2,3),n2,'FontSize',8)
        hold on
        xlabel('X')
        ylabel('Y')
        zlabel('Z')
        axis equal
        grid off
        offset = 0.1;
        xlim([min(Coords(:,1))-offset,max(Coords(:,1))+offset])
        ylim([min(Coords(:,2))-offset,max(Coords(:,2))+offset])
        zlim([min(Coords(:,3))-offset,max(Coords(:,3))+offset])
        set(gca,'View',[ -29.5000   28.4000])
        if travel_motion==1 || travel_motion==2
            if ii==1
                pause(2)
            else
                pause(0.1)
            end
        else
            pause (0.1)
        end
        frame=getframe(gcf);
        writeVideo(v,frame)
    end

end
close(v)
end