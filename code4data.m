function code4data(beta_type)
tic
%%  part 1
%beta_type = 1;

%%   数据生成
charset={'a','b','c','d'};
filename=cell(1,8);
file={strcat(['beta',num2str(beta_type),'fig1','.xls']),strcat(['beta',num2str(beta_type),'fig2','.xls'])};
for i=1:8
    filename{i}=strcat(['fig',num2str(floor((i-1)/4)+1),charset{mod((i-1),4)+1},'.xls']);
end

% r=floor(11/30*n_1*n_2*p/(n_1+n_2));

SET=[	1000 1000 36 0.2 1 1 1;	    %对应(1a)
        1000 1000 36 0.2 0 1 1;	    %对应(1b)
        1000 1000 36 0.2 1 1 0;	    %对应(1c)
        1000 1000 36 0.2 0 1 0; 	%对应(1d)
        1500 1000 44 0.2 1 1 1;     %对应(2a)
        1500 1000 44 0.2 0 1 1;     %对应(2b)
        1500 1000 44 0.2 1 1 0;     %对应(2c)
        1500 1000 44 0.2 0 1 0;     %对应(2d)
        ];
for i=1:8
    maincode(SET(i,:),file{floor((i-1)/4)+1},charset{mod((i-1),4)+1},beta_type);
end
% for i=1:4
%     A=readmatrix(filename1{i});
%     writematrix(A','fig1.xls','Sheet',charset{i},'Range','A1');
% end

%%   图片绘制
% figure
% for i=1:4
%     A=readmatrix(filename1{i});
%     subplot(2,2,i);
%     plot(A(1,:),log(A(2,:)),'r-',A(3,:),log(A(4,:)),'r--');
%     hold on
%     plot(A(5,:),log(A(6,:)),'b-',A(7,:),log(A(8,:)),'b--');
% end




%%   删除过渡文件
% for i=1:4
%     delete(filename1{i});
%     delete(filename2{i});
% end
% for i=1:8
%     delete(filename{i});
% end

toc
end