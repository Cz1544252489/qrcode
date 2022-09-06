clc,clear,close


%%  part 1
filename={'fig1a.xls','fig1b.xls','fig1c.xls','fig1d.xls'};
SET=[	1000 1000 36 0.2 1 1 1;	%respect to (a)
	1000 1000 36 0.2 0 1 1;	%respect to (b)
	1000 1000 36 0.2 1 1 0;	%respect to (c)
	1000 1000 36 0.2 0 1 0;];	%respect to (d)
for i=1:4
maincode1(SET(i,:),filename{i});
end

plot(fig1d(1,:),log(fig1d(2,:)),'r-',fig1d(3,:),log(fig1d(4,:)),'r--');
hold on;
plot(fig1d(5,:),log(fig1d(6,:)),'b-',fig1d(7,:),log(fig1d(8,:)),'b--');


%%   Part 2

filename1={'fig2a1.xls','fig2a2.xls','fig2b1.xls','fig2b2.xls','fig2c1.xls','fig2c2.xls','fig2d1.xls','fig2d2.xls'};
filename={'fig2a.xls','fig2b.xls','fig2c.xls','fig2d.xls'};
SET=[	2000 1000 24 0.1 1 1 1;2000 1000 48 0.2 1 1 1;	%respect to(a)
	2000 1000 24 0.1 0 1 1;2000 1000 48 0.2 0 1 1;	%respect to(b)
	2000 1000 24 0.1 1 1 1;2000 1000 48 0.2 1 1 1;	%respect to(c)
	2000 1000 24 0.1 0 0 0;2000 1000 48 0.2 0 0 0;];	%respect to(d)
for i=1:8
maincode2(SET(i,:),filename1{i});
end

for i=1:4
A=xlsread(filename{2*i-1});B=xlsread(filename{2*i});
writematrix([A;B],filename{i});
end
