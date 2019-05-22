# cudaexpirement

Little something I'm working on that uses the gpu to do static-point sampling computations of a fluid simulation (in this case tectonics).

# (known) issues

The simulation is mostly intuitive and not physically accurate. As a consequence it has a bug that causes the values to stop changing (sort of like when a NN converges).

The client projection is awful. In addition I would like to add viewport controls.

# Using
(Probably as a note to myself in the future, this is too barebones of a project to be useful to other people.)

You can compile a new type of simulation by editing the `osmate()` device function in `sim.cu` and running `nvcc sim.cu`.

Afterwards you can run `python display.py <RES> <ALT>`. See console log for resolution and altitude. I highly recommend changing these values in `sim.cu`. (In the future I could invest some minimal effort to make these variables arguments to the binary call, but I'm way too lazy for that.)
