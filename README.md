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

The webserver looks for a comma separated file named '~/web_password'.

>username,password,
>EOF

Please, create one.


To start the server:
Enter Fitzhugh_3-cell directory
$ ./webrun.py


-------------------------------------------
CONTRIBUTORS
-------------------------------------------

Justus Schwabedal

Drake Knapper

Krishna Pusuluri

Deniz Alacam

