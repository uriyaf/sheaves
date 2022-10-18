# sheaves

Welcome!

The Python code in this project was written by Uriya First and was used to test the modification process 
introduced in Section 12 this paper: https://arxiv.org/abs/2208.01778. It includes an implementation of 
sheaves in python and methods to compute their cohomology. The code also implements
some special simplicial complexes, e.g., triangulations of tori and some Ramanujan complexes. 

The main file is 'main.py'. It is meant to be 
run via the command line. It applies the modification process to a sheaf of your choice from the file 
'sheaf_examples.py'. Run 'main.py -h' for help and information about the process. 

The folder 'outputs' also includes a list of outputs, some of which took days to produce. The parameters used to
produce each output are listed on the first line of the file. Unfortunately, the stages in the output are
shifted by one. That is, the so-called "i-th stages" correspond to the (i+1)-th step in the modification
process considered in the paper, and dim(E'i) is dim(E'_{i+1}) in the paper.  

Quite a lot of running time is needed to produce sheaves and compute their cohomology. By default, the 
code will save the sheaves, and their cohomology bases, and it allows you to (manually) upload them in order to save time. The reason I did not include these extra files is because they are very big (over 
100MB in some cases). 

I did my best to test the code thoroughly. Most of the testing I did can be found in the file 
'tests.py'. To report bugs, please contact Uriya First. 

Disclaimer:
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND.

This code could have been better written. While I would value bug reports and extensions of the code, please do not 
write to me with suggestions about improving the code's design. It is possible that using SAGE would 
help in making the code more efficient, but I presently have no plans of exploring this.

What is still missing:
The are no functions for computing cup products. There used to be some in an older version, but I did 
not get the time to implement them in the present version. In the case of tori, the cup product can be 
computed by theoretical means.
