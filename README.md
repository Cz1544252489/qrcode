# code_in_my_article

This .txt file give all settings used in the article.
Settings includes 8 parameters:
    n_1, n_2, r, p, primarydata, intial, choice, way_fsk, (color); 
And they signify: 
    n_1: the row of primary matrix; 
    n_2: the column of primary matrix;
    r: the rank of primary matrix;
    p: the frequency of the index set Omega;
    primarydata: the parameter to choose which type of primary matrix;
    intial: the parameter to choose which type of intialization;
    choice: the parameter to choose which type of setting in the codes;
    way_fsk: the parameter to choose which type of line search;
    (color: the parameter to change the choice of color;)
And the choice of these parameters show on the function ‘maincode.m’;

Then how to use it? here is MATLAB code: 

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
maincode([1000 1000 36 0.2 1 1 1 1],'result.txt');
<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

The simplest way is just copy an array below to the input of the function.
if  we just code:

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
maincode();
<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
The default setting will be used.


part 1 Primary data generation

[1000 1000 36 0.2 1 1 1 1]
[1000 1000 36 0.2 0 1 1 1]
[1000 1000 36 0.2 1 1 2 0]
[1000 1000 36 0.2 0 1 2 0]

part 2 size and frequency of Omega

set=[1000 1000 0 0.1 1 1 1 1 0;
     1000 1000 0 0.2 1 1 1 1 1;
     1000 1000 0 0.1 1 1 3 0 1;
     1000 1000 0 0.2 1 1 3 0 2;
     2000 1000 0 0.1 1 1 1 1 0;
     2000 1000 0 0.2 1 1 1 1 1;
     2000 1000 0 0.1 1 1 3 0 1;
     2000 1000 0 0.2 1 1 3 0 2;
     2000 2000 0 0.1 1 1 1 1 0;
     2000 2000 0 0.2 1 1 1 1 1;
     2000 2000 0 0.1 1 1 3 0 1;
     2000 2000 0 0.2 1 1 3 0 2;
];









