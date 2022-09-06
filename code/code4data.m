clear,close

%%  part 1 

%%   数据生成
filename1={'fig1a.xls','fig1b.xls','fig1c.xls','fig1d.xls'};
SET1=[	1000 1000 36 0.2 1 1 1;	    %对应(a)
	    1000 1000 36 0.2 0 1 1;	    %对应(b)
	    1000 1000 36 0.2 1 1 0;	    %对应(c)
	    1000 1000 36 0.2 0 1 0;];	%对应(d)
% for i=1:4
%     maincode1(SET1(i,:),filename1{i});
% end

%%   图片绘制
figure
for i=1:4
    A=readmatrix(filename1{i});
    subplot(2,2,i);
    plot(A(1,:),log(A(2,:)),'r-',A(3,:),log(A(4,:)),'r--');
    hold on
    plot(A(5,:),log(A(6,:)),'b-',A(7,:),log(A(8,:)),'b--');
end


%% part 2
%%   数据生成
filename={'fig2a1.xls','fig2a2.xls','fig2b1.xls','fig2b2.xls','fig2c1.xls','fig2c2.xls','fig2d1.xls','fig2d2.xls'};
filename2={'fig2a.xls','fig2b.xls','fig2c.xls','fig2d.xls'};
SET2=[	2000 1000 24 0.1 1 1 1;2000 1000 48 0.2 1 1 1;	    %对应(a)
	    2000 1000 24 0.1 0 1 1;2000 1000 48 0.2 0 1 1;	    %对应(b)
	    2000 1000 24 0.1 1 1 1;2000 1000 48 0.2 1 1 1;	    %对应(c)
	    2000 1000 24 0.1 0 0 0;2000 1000 48 0.2 0 0 0;];	%对应(d)
% for i=1:8
% maincode2(SET2(i,:),filename{i});
% end
% 
% for i=1:4
% A=xlsread(filename{2*i-1});B=xlsread(filename{2*i});
% writematrix([A;B],filename2{i});
% end

%%   图片绘制
figure
for i=1:4
    A=readmatrix(filename2{i});
    subplot(2,2,i);
    plot(A(1,:),log(A(2,:)),'r-',A(3,:),log(A(4,:)),'r--');
    hold on
    plot(A(5,:),log(A(6,:)),'b-',A(7,:),log(A(8,:)),'b--');
end
