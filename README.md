Motiftoolbox
------------------------------------------

This is Motiftoolbox, a set of versatile tools to explore the dynamics of network motifs.




-------------------------------------------
INSTALL
-------------------------------------------

LINUX
-------------------------------------------

Motiftoolbox consists of a set of "./Tools".  To activate the tools, perform

$ cd Tools

$ make all

To use predefined programs, enter one of the program directories (e.g.
Leech_3-cell) and perform

$ ./run.py



MAC OSX
-------------------------------------------

 1. Download SourceTree application from http://www.sourcetreeapp.com/.
 2. Copy the link https://github.com/jusjusjus/Motiftoolbox.git.
 3. Open SourceTree, click on 'Clone Repository'.
 4. Paste the code into Source Path / URL and add .../Motitoolbox to the destination path
 5. Click clone. Let it finish and close SourceTree
 6. Open a terminal and type the following 
 

$ cd Motiftoolbox/Tools

$ make all

enter one of the exixting directories (e.g. cd Fitzhugh_3-cell )


$ ./run.py



Enjoy.



-------------------------------------------
Web Server
-------------------------------------------
One time installation :
sudo apt-get install python-pip
sudo pip install Flask
sudo pip install mpld3

The webserver looks for a comma separated file named '/usr/sbin/MotiftoolboxPwd'.

>username,password,
>EOF

Please, create one.


To start the server:
Enter Fitzhugh_3-cell directory
$ ./webrun.py

-------------------------------------------
Noise
-------------------------------------------
The program uses gsl libraries for ansi-c noise computation.  First you need to install the gsl package, and then you need to make sure that your gcc will link the libgsl, and libgslcblas.  Then, you can enable noise with the command.

make noise

Libraries using noise will be automatically included, and can be used.  Note, that Fitzhugh-n_cell uses the scipy noise generators.

-------------------------------------------
AUTO
-------------------------------------------
Tools/orbit.py uses the python module of AUTO-07p to compute periodic orbits exactly.  If not installed, the orbit is approximated through a long-enough forward integration.

Download AUTO-07p from http://sourceforge.net/projects/auto-07p/?source=typ_redirect and follow the directions.  Make sure to but in your .bashrc

> export PYTHONPATH=$PYTHONPATH:$HOME/auto/07p/python:
> source $HOME/auto/07p/cmds/auto.env.sh


-------------------------------------------
CONTRIBUTORS
-------------------------------------------

Justus Schwabedal

Drake Knapper

Krishna Pusuluri

Deniz Alacam

