%This function is used to create the graph Laplacian matrix on
  %prior discretization sture

%Input:
 %Vertices---a cell structure of prior discretization points in all the
             %rings
 %Vert_all---all the prior discretization points including the center            
 %n_th---a vector containing number of discretization points in each ring
 %n_neigh---number of neighbors in the graph Laplacian
 
 %Output:
 %L_2---graph Laplacian matrix.



function L_2=graph_laplace(Vertices,Vert_all,n_neigh,visualize)


L_2=zeros(size(Vertices,2));%allocate space for Laplacian matrix
   
%create a unit circle to visualize the nearest neghbors 
        r=1;%radius
        x=0;%center
        y=0;
        th = 0:pi/200:2*pi;%angular discretization on the boundary of unit disc
        xunit = r * cos(th) + x;
        yunit = r * sin(th) + y;



   for i=1:size(Vert_all,2)
       [~,I]=sort(sqrt((Vert_all(1,i)-Vert_all(1,:)).^2 + (Vert_all(2,i)-Vert_all(2,:)).^2));
       %find the distance of all the remaining points from current
       if sqrt((Vert_all(1,i))^2 + (Vert_all(2,i))^2) >1-1e-6
           ind=I(2:4);%take n_neigh nearest neighbors
       else
            ind=I(2:n_neigh+1);%take n_neigh nearest neighbors
       end
        %plot the neighbors   
           if strcmp(visualize,'yes')
           figure(5)
     
          plot(xunit, yunit,'k','LineWidth',2);
          hold on
           scatter(Vert_all(1,i),Vert_all(2,i),'*r')
           sz = 100;
            scatter(Vert_all(1,i),Vert_all(2,i),sz,'.b')
         xlim([-1 1])
            ylim([-1 1])
           axis equal
           hold on
           scatter(Vert_all(1,ind),Vert_all(2,ind),sz,'.r')
           xlim([-1 1])
            ylim([-1 1])
           axis equal
           hold off
            
           pause(0.005)
           end
           %graph Laplacian
           L_2(i,ind)= -1/length(ind);%averaging at the index of neighboring points
           L_2(i,i)=1;%main diagonal
   end

end
   