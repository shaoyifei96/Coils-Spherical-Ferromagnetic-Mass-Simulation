grid= 0:0.1:5;
[X,Y,Z] = meshgrid(grid);
i=1;j=1;k=1;
dot([X(i,j,k) Y(i,j,k) Z(i,j,k)],Location)