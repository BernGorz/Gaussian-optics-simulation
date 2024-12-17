Code by Bernard Gorzkowski\
Contact: bernard.gorzkowski@gmail.com

### Description ###

This code allows to simulate the behavior of a gaussian beam passing through a series of lenses, with their positions adjustable in real time. You can probe the beam for its parameters by clicking at any point on the graph.

### Parameters ###

The parameters of the simulation have to be modified directly in the code, in the "Parameters" section.

```
L = 3 * m             # Simulation range
w0 = 500 * um            # Initial beam waist
z0 = 0 * mm             # Initial beam position
wavelength = 810 * nm   # Wavelength

lens_f = [150*mm, 300*mm, 600*mm]         # List of focal lengths of your lenses
lens_z = [150*mm, 600*mm, 1500*mm]       # List of positions of your lenses

nb_points = 1024                # Number of points for graph
```

You can add or modify any number of lenses in the given list parameters.