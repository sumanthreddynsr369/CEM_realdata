% Plot the voltage mesh and the conductivity mesh

load VoltageMesh_12

nt = size(Topol,1);
figure(1)
for j = 1:nt
    x = Vertices(1,[Topol(j,1),Topol(j,2),Topol(j,3),Topol(j,1)]);
    y = Vertices(2,[Topol(j,1),Topol(j,2),Topol(j,3),Topol(j,1)]);
    plot(x,y,'k-')
    hold on
end
for ell = 1:16
    x = 1.01*Vertices(1,ElectrodeNodes(ell,:));
    y = 1.01*Vertices(2,ElectrodeNodes(ell,:));
    plot(x,y,'r-','LineWidth',4)
end


load ConductivityMesh_120

nt = size(CondTopol,1);
for j = 1:nt
    x = 2.1+CondVertices(1,[CondTopol(j,1),CondTopol(j,2),CondTopol(j,3),CondTopol(j,1)]);
    y = CondVertices(2,[CondTopol(j,1),CondTopol(j,2),CondTopol(j,3),CondTopol(j,1)]);
    plot(x,y,'k-')
    hold on
end
hold off

axis('equal')
axis('off')



