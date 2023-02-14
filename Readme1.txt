There are three m-file for my experiment, "code4randomdata.m", "code4movielens.m" and "ratio4speedup.m".


Usage:
----------code4randomdata.m---------------------
Run the m-file code4randomdata.m

And you will get 8 xls-file to all the figures on random data part.

They name as 'p1beta6fig1.xls','p1beta6fig2.xls','p2beta6fig1.xls','p2beta6fig2.xls',
	    'p3beta6fig1.xls','p3beta6fig2.xls','p4beta6fig1.xls','p4beta6fig2.xls';

Each one of them contain five sheets 'set','a','b','c','d'. And sheet 'set' shows the setting of this experiment. The others contain a table with 250 rows, 8 columns, where row means the iteration times and every two columns mean the x-axis and y-axis of a line on the subfigure sorted by the sequence shown as legend.

It cost about 4 hours on my machine with intel i7-10700F CPU@2.9GHz  and 32GB memory for all.

Running them part by part will be a good idea.


-----------code4movielens.m------------------
run the m-file code4movielens and you will get a xls-file.

It names "MovieLens100k.xls", and contains eight sheets 'seta', 'a', 'setb', 'b', 'setc', 'c', 'setd', 'd'. 

Every two sheets represent a subfigure on experiment of real-world part and the sheet named with 'set' shows the setting or parameter of this experiment.

it contains the result of expriments on the dataset movielens100k.

------------ratio4speedup.m-------------------
run the m-file code4movielens and it will show the variable:

ran_SU_all: list contains ratio of speedup for each experiment on random data;
ran_ratio_of_good: mean of ratio of good speedup for each experiment on random data;
ran_ratio_of_vastly_good: proportion of vastly good speedup for experiments on random data.

real_SU_all: list contains ratio of speedup for each experiment on real data;
real_ratio_of_good: mean of ratio of good speedup for each experiment on real data;
real_ratio_of_vastly_good: proportion of vastly good speedup for experiments on real data.














