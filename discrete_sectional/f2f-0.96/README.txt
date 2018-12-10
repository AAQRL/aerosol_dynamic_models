# f2f - Fortran 77 to Fortran 90/95 source code conversion

f2f is a [Perl](http://en.wikipedia.org/wiki/Perl) script which does much of the tedious work of converting Fortran 77 source code into [Fortran 90/95 form](http://en.wikipedia.org/wiki/Fortran). There seems to be a lot of Fortran hate in the world, and I think this comes from people who have been forced to use Fortran 77 at some time or another. Hopefully, this program will make you a less hateful person.

Download [f2f](http://bitbucket.org/lemonlab/f2f/downloads), release 0.92.

## USAGE:
> f2f [inputfile [outputfile]]

 e.g.:
> f2f legacycode.f legacycode.f90

I wrote this program for my own needs, and I have successfully used it many times. Since it might be helpful to others, I am making it available under the [GNU](http://www.gnu.org/) GPL. It has worked quite well for me with standard Fortran 77 source code, but can sometimes give problems on mixed 77/90 code. Give it a try, and see how it does on your code. I can, of course, make no guarantees that it will spit out code that suits your aesthetic tastes (it indents code according to my tastes), nor can I guarantee that it will generate code that compiles (you may need to make an edit or two first). In some cases it can generate really wacky code, especially if you feed it really wacky code.

## Bugs?
If you find a piece of standard Fortran 77 that my code munges into bad Fortran 90, I'd appreciate it if you could either send me your F77 code or a small section of the code (or at worst a detailed description of the few lines that cause the problems) so that I can try to fix f2f.

## Known Bug
f2f does not handle Windows line endings well (Carriage return/Line Feed). Perhaps I can fix this when I find time. In the meantime, you may want to convert Windows format text files to UNIX format before you use f2f. This can be done with various utilities or text editors.

### Alternative Options:

[f77tof90](http://www.soton.ac.uk/~fortran/tools/f77tof90/f77tof90.html), a Bourne shell script for dealing with non-standard F77 code. It converts Record and Structure F77 extensions into F90 Type statements. It also converts C-preprocessed #include statements into Fortran INCLUDE statements. Written by Jack Scheible.

[convert.f90](ftp://ftp.numerical.rl.ac.uk/pub/MandR/convert.f90), written by Mike Metcalf.

[to_f90.f90](http://users.bigpond.net.au/amiller/to_f90.f90), written by Alan Miller.

[ftof90.c](ftp://ftp.ifremer.fr/ifremer/ditigo/fortran90/ftof90.c.gz), for converting comments and continuation lines, written by Michael Olagnon.
