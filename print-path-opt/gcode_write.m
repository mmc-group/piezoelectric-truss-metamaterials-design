function gcode_write(path,lattice_coords)
%%% If writing from xlsx file
% lattice_coords=xlsread('outputs/e31e32e33_3X.xlsx','Coords');
% path=xlsread('outputs/e31e32e33_3X.xlsx','Path');
lattice_size=15; %mm
lattice_coords= lattice_coords*lattice_size; %mm
extmultiplier=2.75*0.0035;
nozdia=0.50;
ext_factor=1;     % for increasing extrusion at nodes
F_factor=0.1;         % for reducing speed at the nodes
travel_factor=0;
lattice_height=max(lattice_coords(:,3))-min(lattice_coords(:,3));
ind=1;
travel_move=0;
for i=1:length(path)
    if i>1
        x=lattice_coords(path(i),1)-lattice_coords(path(i-1),1);
        y=lattice_coords(path(i),2)-lattice_coords(path(i-1),2);
        z=lattice_coords(path(i),3)-lattice_coords(path(i-1),3);
        L=sqrt(x^2+y^2+z^2);
    else
        L=0;
    end
    E_extr=extmultiplier*L*pi/4*nozdia^2;
    E_travel=travel_factor*L;
    F_extr=50;
    F_travel=100;
    F_nodes=1;
    if i>2
        if path(i-1)==path(i)

            % %  it means the repeated node, i.e. the first node of the current strut is same as the
            % %  last node of the previous one.
            % %  Do Nothing
        else
            if path(i-1)~=path(i) && path(i-2)==path(i-1)
                % %      The travel path is being traversed without break. Add coordinates
                % %      to G-code
                path_coords(ind,1:3)=lattice_coords(path(i),:);
                path_coords(ind,4)=path_coords(ind-1,4)+E_extr;
                path_coords(ind,5)=F_extr;
                ind=ind+1;
            elseif path(i-1)~=path(i) && path(i-2)~=path(i-1)
                if travel_move==1
                    path_coords(ind,1:3)=lattice_coords(path(i),:);
                    path_coords(ind,4)=path_coords(ind-1,4)+E_extr;
                    path_coords(ind,5)=F_extr;
                    ind=ind+1;
                    travel_move=0;
                else
                    % %         travel motion required since the current strut has a different
                    % %         first node than the second node of previous strut
                    % %         moving upwards above the height of the lattice
                    path_coords(ind,1:3)=[lattice_coords(path(i-1),1:2),2*lattice_height]';
                    if abs(path_coords(ind-1,1)-max(lattice_coords(:,1)))<1e-5
                        path_coords(ind,1)=max(lattice_coords(:,1))+0.1*lattice_height;
                    elseif abs(path_coords(ind-1,1)-min(lattice_coords(:,:)))<1e-5
                        path_coords(ind,1)=min(lattice_coords(:,1))-0.1*lattice_height;
                    end
                    if abs(path_coords(ind-1,2)-max(lattice_coords(:,2)))<1e-5
                        path_coords(ind,2)=max(lattice_coords(:,2))+0.1*lattice_height;
                    elseif abs(path_coords(ind-1,2)-min(lattice_coords(:,2)))<1e-5
                        path_coords(ind,2)=min(lattice_coords(:,2))-0.1*lattice_height;
                    end
                    path_coords(ind,4)=path_coords(ind-1,4)+E_travel;
                    path_coords(ind,5)=F_travel;
                    ind=ind+1;
                    % %         moving to the beigining point of the next travel path, at an
                    % %         elevation 1.2*lattice height
                    path_coords(ind,1:3)=[lattice_coords(path(i),1:2),2*lattice_height]';
                    path_coords(ind,4)=path_coords(ind-1,4)+E_travel;
                    path_coords(ind,5)=F_travel;
                    ind=ind+1;
                    % %         Moving downward to the location of first node of the next
                    % %         travel path
                    path_coords(ind,1:3)=[lattice_coords(path(i),:)]';
                    path_coords(ind,4)=path_coords(ind-1,4)+E_travel;
                    path_coords(ind,5)=F_travel;
                    ind=ind+1;
                    travel_move=1;
                end
            end
        end
    elseif i==1
        path_coords(ind,1:3)=lattice_coords(path(i),:);
        path_coords(ind,4)=E_travel;
        path_coords(ind,5)=F_travel;
        ind=ind+1;
    elseif i==2
        path_coords(ind,1:3)=lattice_coords(path(i),:);
        path_coords(ind,4)=path_coords(ind-1,4)+E_extr;
        path_coords(ind,5)=F_extr;
        ind=ind+1;

    end

end

% % For adding a pause at vertices
path_coords1(1,:)=path_coords(1,:);
ind=1;
for ii=1:size(path_coords,1)-1
    ind=ind+1;

    path_coords1(ind,1)=(path_coords(ii,1)*0.95)+(path_coords(ii+1,1)*0.05);
    path_coords1(ind,2)=(path_coords(ii,2)*0.95)+(path_coords(ii+1,2)*0.05);
    path_coords1(ind,3)=(path_coords(ii,3)*0.95)+(path_coords(ii+1,3)*0.05);
    if abs(path_coords(ii+1,5)-F_travel)<1e-5
        x=path_coords1(ind,1)-path_coords1(ind-1,1);
        y=path_coords1(ind,2)-path_coords1(ind-1,2);
        z=path_coords1(ind,3)-path_coords1(ind-1,3);
        L=sqrt(x^2+y^2+z^2);
        E_extr=ext_factor*travel_factor*L*pi/4*nozdia^2;
        path_coords1(ind,4)=path_coords1(ind-1,4)+E_extr;
        path_coords1(ind,5)=F_travel;
    else
        x=path_coords1(ind,1)-path_coords1(ind-1,1);
        y=path_coords1(ind,2)-path_coords1(ind-1,2);
        z=path_coords1(ind,3)-path_coords1(ind-1,3);
        L=sqrt(x^2+y^2+z^2);
        E_extr=ext_factor*extmultiplier*L*pi/4*nozdia^2;
        path_coords1(ind,4)=path_coords1(ind-1,4)+E_extr;
        path_coords1(ind,5)=F_factor*F_extr;
    end


    ind=ind+1;
    path_coords1(ind,1)=(path_coords(ii,1)*0.05)+(path_coords(ii+1,1)*0.95);
    path_coords1(ind,2)=(path_coords(ii,2)*0.05)+(path_coords(ii+1,2)*0.95);
    path_coords1(ind,3)=(path_coords(ii,3)*0.05)+(path_coords(ii+1,3)*0.95);
    if abs(path_coords(ii+1,5)-F_travel)<1e-5
        x=path_coords1(ind,1)-path_coords1(ind-1,1);
        y=path_coords1(ind,2)-path_coords1(ind-1,2);
        z=path_coords1(ind,3)-path_coords1(ind-1,3);
        L=sqrt(x^2+y^2+z^2);
        E_extr=travel_factor*L*pi/4*nozdia^2;
        path_coords1(ind,4)=path_coords1(ind-1,4)+E_extr;
        path_coords1(ind,5)=F_travel;
    else
        x=path_coords1(ind,1)-path_coords1(ind-1,1);
        y=path_coords1(ind,2)-path_coords1(ind-1,2);
        z=path_coords1(ind,3)-path_coords1(ind-1,3);
        L=sqrt(x^2+y^2+z^2);
        E_extr=extmultiplier*L*pi/4*nozdia^2;
        path_coords1(ind,4)=path_coords1(ind-1,4)+E_extr;
        path_coords1(ind,5)=F_extr;
    end


    ind=ind+1;
    path_coords1(ind,1)=path_coords(ii+1,1);
    path_coords1(ind,2)=path_coords(ii+1,2);
    path_coords1(ind,3)=path_coords(ii+1,3);
    if abs(path_coords(ii+1,5)-F_travel)<1e-5
        x=path_coords1(ind,1)-path_coords1(ind-1,1);
        y=path_coords1(ind,2)-path_coords1(ind-1,2);
        z=path_coords1(ind,3)-path_coords1(ind-1,3);
        L=sqrt(x^2+y^2+z^2);
        E_extr=ext_factor*travel_factor*L*pi/4*nozdia^2;
        path_coords1(ind,4)=path_coords1(ind-1,4)+E_extr;
        path_coords1(ind,5)=F_travel;
    else
        x=path_coords1(ind,1)-path_coords1(ind-1,1);
        y=path_coords1(ind,2)-path_coords1(ind-1,2);
        z=path_coords1(ind,3)-path_coords1(ind-1,3);
        L=sqrt(x^2+y^2+z^2);
        E_extr=ext_factor*extmultiplier*L*pi/4*nozdia^2;
        path_coords1(ind,4)=path_coords1(ind-1,4)+E_extr;
        path_coords1(ind,5)=F_factor*F_extr;
    end

end

%%% plotting travel and printing paths
for ii=1:size(path_coords1,1)-1
    connect=[ii,ii+1];
    if abs(path_coords1(ii+1,5)-F_travel)<1e-5
        plot3(path_coords1(connect,1),path_coords1(connect,2),...
            path_coords1(connect,3),'*--','Color',[0,0,0],'LineWidth',2)
    else
        plot3(path_coords1(connect,1),path_coords1(connect,2),...
            path_coords1(connect,3),'*--','Color',[1,0,0],'LineWidth',2)
        % elseif ii
    end
%     title("Extrusion Rate"+path_coords1(connect(2),4))
% %     pause
hold on
% figure (2)
% plot(path_coords1(connect(2),4))
end
set(gca,'View',[ -29.5,28.4])

%%% Adjust for printer initial position and priming
path_coords1(:,1)=path_coords1(:,1)+73;
path_coords1(:,2)=path_coords1(:,2)+70;
path_coords1(:,4)=path_coords1(:,4)+0.002;   %for priming, priming value roughly about 0.002

%%Checking for negative cooridantes
xmin=min(path_coords1(:,1));
ymin=min(path_coords1(:,2));
zmin=min(path_coords1(:,3));

%%% moving travel motions at the boundaries to the outside of the lattice
%%% to avoid dripping of the ink and printing extra material
if xmin<0
    path_coords1(:,1)=path_coords1(:,1)-xmin;
end
if ymin<0
    path_coords1(:,2)=path_coords1(:,2)-ymin;
end
if zmin<0
    path_coords1(:,3)=path_coords1(:,3)-zmin;
end

%%% writing the g-code to a text file
formatSpec = 'G1 X%f Y%f Z%f E%f F%f\n';
filename='outputs/BCC_v05MAr_E2.75x.txt';
fileID = fopen(filename,'w');
fprintf(fileID,formatSpec,path_coords1');
fclose(fileID);
end