function truss_plotting (Coords,Connectivity)

for ii=1:size(Connectivity,1)
% 
    plot3(Coords((Connectivity(ii,:)),1),Coords(Connectivity(ii,:),2),...
        Coords(Connectivity(ii,:),3),'Color',[0.2941,0.4902,0.6784],'LineWidth',2)
    set(gcf,'color','w');
    node1=Connectivity(ii,1);
    n1=num2str(node1);
    node2=Connectivity(ii,2);
    n2=num2str(node2);
%     text(Coords(node1,1),Coords(node1,2),Coords(node1,3),n1,'FontSize',14)
%     text(Coords(node2,1),Coords(node2,2),Coords(node2,3),n2,'FontSize',14)
    hold on
    xlabel('X')
    ylabel('Y')
    zlabel('Z')

%         pause (0.2)
end
upper=[max(Coords(:,1)), max(Coords(:,2)), max(Coords(:,3))];
lower=[min(Coords(:,1)), min(Coords(:,2)), min(Coords(:,3))];
% plotcube([upper(1)-lower(1),upper(2)-lower(2),upper(3)-lower(3)],lower,0.1,[0.7,0.7,0.7])
axis equal
grid off
offset = 0.1;
xlim([min(Coords(:,1))-offset,max(Coords(:,1))+offset])
ylim([min(Coords(:,2))-offset,max(Coords(:,2))+offset])
zlim([min(Coords(:,3))-offset,max(Coords(:,3))+offset])
set(gca,'View',[ -29.5000   28.4000])

end