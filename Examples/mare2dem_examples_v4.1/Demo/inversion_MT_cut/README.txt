This is an example of adding a regularization penalty cut to the inversion. Load the .resistivity file into Mamba2D
to see the cuts, which are shown as white colored parameter boundaries.

When you want to make you own penalty cuts, it's helpful to add them at the start when building the model grid. 
This inversion grid was made by importing the file  topSalt.txt into Mamba, setting a penalty cut for 
its segments, and then generating the inversion parameter grid using the grid options (quads or triangles).
 