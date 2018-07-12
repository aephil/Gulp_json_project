# Gulp_json_project

Phonon website
==============

Visualize phonon vibrational modes.

This project aims to provide a simple way to visualize the lattice vibrations of different materials. The temperature of a material is related to the agitation of its atoms. The atoms can move in any of the three cartesian directions. Combining the different possible ways the atoms can vibrate we obtain the eigenvectors. Each mode has associated a frequency of vibration that is related with the forces between the atoms.

How to use?
===========

In the phonon section you can click on any point in the dispersion relation and see an animation of how the atoms vibrate in that particular mode.
By default you can visualize the phonon dispersion of some materials we calculated.
If you want to see your own calculations done in Gulp, you can use this package.

Gulp
----------------
To read a Quantum Espresso calculation you need three files `<prefix>.gin` and `<prefix>.eig` and `<prefix>.disp`. The first one is the input file for gulp, the second and third one can be generated with `output phon <prefix>`, and `output eig <prefix>` in `<prefix>.gin` file. 
After installing the python scripts by:
	$ python setup.py develop


 you can obtain the `.json` files as:

    $ read_gulp.py prefix

You can then select the resulting `.json` file with the `Choose files` button.