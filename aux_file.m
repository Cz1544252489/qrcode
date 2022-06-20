clear,clc,close


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


% r \in (mnp/3/(m+n),mnp/2.5/(m+n)) and r=floor(11*mnp/30/(m+n))

for i=1:size(set,1)
    set(i,3)=floor(11*set(i,1)*set(i,2)*set(i,4)/30/(set(i,1)+set(i,2)));
end

filename='table1.txt';
for j=1
    for i=2*j-1:2*j
        fid=fopen(filename,'a');
        OS=set(i,1)*set(i,2)*set(i,4)/(set(i,1)+set(i,2)-set(i,3))/set(i,3);
        if(mod(i,2)==1)
            fprintf(fid,'===========================\n');
        end
        fprintf(fid,'m=%d,n=%d,p=%f,r=%d,OS=%f\n',set(i,1),set(i,2),set(i,4),set(i,3),OS);
        fprintf(fid,'iteration\t RMSE\t time\n');
        fclose(fid);
        if(mod(i,2)==1)
            figure(floor(i/2)+1);
        end
        maincode(set(i,1:end),filename);
        hold on;
        if(mod(i,2)==0)
            legend('0.1 QR','0.1 not','0.2 QR','0.2 not');
        end
    end
end
