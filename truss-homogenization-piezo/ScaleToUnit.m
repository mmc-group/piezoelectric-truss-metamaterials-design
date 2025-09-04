function [Coords]=ScaleToUnit(Coords,connect)

xmin=min(Coords(:,1));
ymin=min(Coords(:,2));
zmin=min(Coords(:,3));

% % % translate to [0,0,0]

Coords(:,1)=Coords(:,1)-xmin;
Coords(:,2)=Coords(:,2)-ymin;
Coords(:,3)=Coords(:,3)-zmin;

% % % scale to Unit Cube


xmax=max(Coords(:,1));
ymax=max(Coords(:,2));
zmax=max(Coords(:,3));

Coords(:,1)=Coords(:,1)./xmax;
Coords(:,2)=Coords(:,2)./ymax;
Coords(:,3)=Coords(:,3)./zmax;

end