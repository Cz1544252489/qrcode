function [X]=maincode(SET,filename)
%%  This function can plot the figure and ouput the performance
%   An example shows below: 
%   maincode([1000 1000 36 0.2 1 1 1 1]);
%   The setting see file "ReadMeFirst.txt" for more.
%% parameter set needed to change
if nargin<2
    filename='table.txt';
    if nargin<1
        %     n_1  n_2  r  p   primarydata  initial choice way_fsk color
        SET=[1000 1000 36 0.2      1           1       1      0      1 ];
    end
end

%             max_iter  epsilon  delta  theta
SET_constant=[ 200      1e-15    1e-4     0.01];
if(length(SET)==8)
    SET=[SET 0];
end
%% definite global variable
global n_1;
global n_2;
global r;
global data;
global delta;
global Omega;
global theta;
global orth_type;
%% primary parameter
n_1=SET(1);n_2=SET(2);r=SET(3);
p=SET(4);%rates of Omega
M=rankmanQR(n_1,n_2,r);
max_iter=SET_constant(1);%max iteration times
epsilon=SET_constant(2);%stop parameter
delta=SET_constant(3);
theta=SET_constant(4);%parameter to justify weather it owns good orthogonality or not
% orth_type=1 %if we use QR method {1 use 0 don`t}
% way_fsk=0;%apply which line search {1  exact  0 inexact}
% method=0;%apply which descent direction  {1  negetive gradient  0  conjugate direction}
% if(~method)
%     beta_type=8;%which conjugate parameter will be choosed in conjugate direction
% end
%% produce primary data
ra=RandStream('mt19937ar','Seed',sum(100*clock));
RandStream.setGlobalStream(ra);
data_=M.gSetData(SET(5));
% generate Omega
if(1)   % sync. data or real data
    Omega=M.gOmega(p);
    data=data_.*Omega;%get new data  E(data)=0.01
else
    dataset_name='u1.base';
    [data,Omega]=real_network(dataset_name);
    [n_1,n_2]=size(data);r=20;
    M=rankmanQR(n_1,n_2,r);
    disp(SET(9));
end
fprintf('E(M)=%f\n',sum(sum(data))/nnz(data));
fprintf('OS=%f\n',nnz(Omega)/M.dim());
intial=SET(6);   %1 means spectral initialization  0 means random initial
if(intial)
    X0=M.gX0(data);
else
    X0=M.rands();
end
%%
para_set={(1),(1:8),(1),[1 0]};%four kinds of parameter: way_fsk , beta_type , orth_type
choice=SET(7);
for iter_var=para_set{choice}
    for orth_type=para_set{4}
        switch choice
            case 1      % for QRRGD and RGD
                way_fsk=1;
                method=1;
            case 2      % for beta in QRRCG and RCG
                method=0;
                way_fsk=SET(8);
                beta_type=iter_var;
                disp(beta_type);
            case 3      % for QRRCG and RCG 
                method=0;
                way_fsk=0;
                beta_type=8;
        end
        %% some important data record
        f0=zeros(1,max_iter);
        s=zeros(1,max_iter);
        RMSE=zeros(1,max_iter);
        
        X=X0;
        t=1;flag=1e5;
        time0=clock;
        if(~method)
            gl=M.grad(X);
            dl=zeros(n_1+n_2,r);
        end
        while(t<max_iter)
            %function value
            f0(t)=M.f(X);
            %compute graident
            g=M.grad(X);
            if(method)
                d=-g;
            else
                beta=M.fbeta(X,gl,g,dl,beta_type);
                d=-g+beta*dl;
            end
            %if stop or not
            if(flag<epsilon)
                break;
            end
            %find stepsize
            if(way_fsk)
                s(t)=M.fsk(X,-d);
            else
                s(t)=M.ufsk(X,g,d);
            end
            %update iterate
            X=M.next(X,s(t),d);
            %calculate RMSE
            H=Omega.*(X(1:n_1,:)*X(n_1+1:end,:)'-data);
            RMSE(t)=sqrt(sum(sum(H.^2))/nnz(Omega));
            %disp(RMSE(t));
            %stop parameter
            flag=M.metric(d,d,X);
            %update iterate parameter
            t=t+1;
            if(~method)
                dl=d;
                gl=g;
            end
        end
        time1=clock;
        time=(time1-time0)*[0;0;43200;3600;60;1];
        if(flag>epsilon)
            disp('over the max_iter!!');
        end
        if(choice==1)
            M.outplot(RMSE,time,[iter_var+1+SET(9),orth_type+1],filename);
        else
            M.outplot(RMSE,time,[iter_var+SET(9),orth_type+1],filename);
        end
        disp(flag);
    end
end
% if(choice==1)
%     legend('fsk,QR','fsk,not','ufsk,QR','ufsk,not');
% end
if(choice==2)
    legend('HS+','','FP','','PRP','','CD','','LS','','DY','','RMIL','','MRMIL','');
end