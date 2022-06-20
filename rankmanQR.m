function [ M ] = rankmanQR( n_1 , n_2 , r )
%   This funtion divides nine parts and the subfuntion in these:
%       1, common part: name, dim, size0
%       2, objective function: f, metric
%       3, metric and gradient: grad
%       4, QR decomposition: qr
%       5, randoms: rands, gX0, gOmega, gSetData
%       6, find stepsize: fsk, ufsk
%       7, update iteration point: next
%       8, conjugate method: fbeta
%       9, output the figure: outplot
    %%   common part
    M.name = @() sprintf('Manifold of %d x %d matrices of rank %d', n_1, n_2, r);
    
    M.dim = @() (n_1+n_2-r)*r;
    
    M.size0 = [n_1,n_2,r];
     
    %% objective function
    M.f=@f;
    function z=f(X)
        global data;
        global Omega;
        s=Omega.*(X(1:n_1,:)*X(n_1+1:n_1+n_2,:)'-data);
        z=0.5*trace(s'*s);
    end
    %% metric and gradient
    %metric
    M.metric = @metric;
    function z =metric( A , B , X )
        global delta;
        global orth_type;
        if nargin <3
            X=zeros(n_1+n_2,r);
        end
        Q=X(1:n_1,:);
        R=X(n_1+1:end,:)';
        if(orth_type)
            z=trace(A(1:n_1,:)'*B(1:n_1,:)*(R*R'+delta*eye(r)))+...
                trace(A(n_1+1:end,:)'*B(n_1+1:end,:)*(1+delta));
        else
            z=trace(A(1:n_1,:)'*B(1:n_1,:)*(R*R'+delta*eye(r)))+...
                trace(A(n_1+1:end,:)'*B(n_1+1:end,:)*(Q'*Q+delta*eye(r)));
        end
    end
    
    %gradient
    M.grad = @grad;
    function Z=grad( X )
        global Omega;
        global data;
        global delta;
        global orth_type;
        Q=X(1:n_1,:);R=X(n_1+1:end,:)';
        S=Omega.*(Q*R-data);
        Z=[S*R';S'*Q];
        if(orth_type)
            Z=[Z(1:n_1,:)/(R*R'+delta*eye(r));Z(n_1+1:end,:)/(1+delta)];
        else
            Z=[Z(1:n_1,:)/(R*R'+delta*eye(r));Z(n_1+1:end,:)/(Q'*Q+delta*eye(r))]; 
        end
    end
    
    %% QR decomposition
    % QR decomposition by MGS
    M.qr = @qr_MGS;
    function [Q,R]=qr_MGS(A)
        [s1,s2]=size(A);
        Q=zeros(s1,r);
        R=zeros(r,s2);
        for k=1:r
            R(k,k)=sqrt(A(:,k)'*A(:,k));
            Q(:,k)=A(:,k)/R(k,k);
            if (k<s2)
                for j=k+1:s2
                    R(k,j)=Q(:,k)'*A(:,j);
                    A(:,j)=A(:,j)-Q(:,k)*R(k,j);
                end
            end
        end
    end
    %% randoms
    % random initialization
    M.rands = @rands;
    function Z=rands(l)
        if nargin<1
            l=0;
        end
        temp=rand(n_1,r+l)*rand(r+l,n_2);
        [Q,R]=M.qr(temp);
        Z=[Q(:,1:r);R(1:r,:)'];
    end

    % generate X for figure process
    M.gXf = @gXf;
    function Z=gXf(data)
        [U,S,V]=svd(data);
        Z=[U(:,1:r)*sqrt(S(1:r,1:r));V(:,1:r)*(sqrt(S(1:r,1:r)))'];
    end

    M.randf = @randf;
    function Z=randf()
        temp=rand(n_1,r)*rand(r,n_2);
        temp=temp/max(temp(:))*255;
        [Q,R]=M.qr(temp);
        Z=[Q(:,1:r);R(1:r,:)'];
    end

    % spectral initialization
    M.gX0 = @gX0;
    function Z=gX0(data)
        [U,S,V]=svd(data);
        Z=[U(:,1:r)*sqrt(S(1:r,1:r));V(:,1:r)*(sqrt(S(1:r,1:r)))'];
    end
    
    % generate Omega
    M.gOmega = @gOmega;
    function Z=gOmega(p)
        Z=rand(n_1,n_2);
        Z(Z>1-p)=1;
        Z(Z<=1-p)=0;
    end
    
    % primary data generation
    M.gSetData = @gSetData;
    function Z=gSetData(input)
        l=0;
        Z=rand(n_1,r+l)*rand(r+l,n_2);  
        if(input)
            Z=Z/(max(Z(:))+max(-Z(:)));
        end
    end
    
    %% find stepsize
    % by exact line search
    M.fsk=@fsk;
    function [z,ss]=fsk(X,g)
       global Omega;
       global data;
       Q=X(1:n_1,:);R=X(n_1+1:n_1+n_2,:)';
       g_Q=g(1:n_1,:);g_R=g(n_1+1:n_1+n_2,:)';
       A=Omega.*(g_Q*g_R);B=Omega.*(Q*g_R+g_Q*R);O=Omega.*(Q*R-data);
       ss=[2*trace(A'*A) 3*trace(A'*B) 2*trace(A'*O)+trace(B'*B) trace(B'*O)];
       z=min(abs(roots(ss)));
    end
    
    % by inexact line search
    M.ufsk=@ufsk;
    function [z,im]=ufsk(X,g,d)
        z=1;
        %z=M.metric(g,d,X)/M.metric(d,d,X);
        im=0;im_max=20;
        while(1)
            z=z*0.5^im;
            if(M.f(X)-M.f(X+z*d)-0.001*z*M.metric(-g,d,X)>=0)
                break;
            else
                im=im+1;
            end
            if(im>=im_max)
                break;
            end
        end
    end
    
    %% update iteration point
    M.next=@next;
    function [Z]=next(X,s,d)
        global theta;
        global orth_type;
        if isempty(theta)
            disp('don`t set theta default 0.01 '); 
            theta=0.01;
        end
        Z=X+s*d;
        if(orth_type)
            h=abs(trace(Z(1:n_1,:)'*Z(1:n_1,:))-r)/r;
            if(h>theta)
                [Xq,Xr]=M.qr(Z(1:n_1,:));
                Z=[Xq;Z(n_1+1:end,:)*Xr'];
            end
        end
    end

    %% conjugate method
    % find conjugate parameter beta
    M.fbeta=@fbeta;
    function z=fbeta(X,gl,g,dl,type) % gl denotes gradient of last��dl denotes direction of last
        upper_border=1e3;
        switch type
            case 1  %HS+
                z=min(upper_border,max(0,M.metric(g-gl,g,X)/M.metric(g-gl,dl,X)));
            case 2  %FP
                z=min(upper_border,max(0,M.metric(g,g,X)/M.metric(gl,gl,X)));
            case 3  %PRP
                z=min(upper_border,max(0,M.metric(g-gl,g,X)/M.metric(gl,gl,X)));
            case 4  %CD
                z=min(upper_border,max(0,M.metric(g,g,X)/M.metric(-gl,dl,X)));
            case 5  %LS
                z=min(upper_border,max(0,M.metric(g-gl,g,X)/M.metric(-gl,dl,X)));
            case 6  %DY
                z=min(upper_border,max(0,M.metric(g,g,X)/M.metric(gl,g-gl,X)));
            case 7  %RMIL
                z=min(upper_border,max(0,M.metric(g,g-gl,X)/M.metric(dl,dl,X)));
            case 8  %MRMIL
                z=min(upper_border,max(0,M.metric(g,g-gl-dl,X)/M.metric(dl,dl,X)));
                
        end
    end

%% output the figure
    
    M.outplot=@outplot;
    function []=outplot(RMSE,time,type,filename)
        if nargin<3
            type=[1,1];
        end
        color={[0 0 0],[0 0 1],[1 1 0],[1 0 0],[1 0 1],[0 1 0],[0 1 1],[0.5,0.5,0.5]};
        linestyle={'-','--','-.','*'};
        % colors are black blue yellow red purple green cyan grey 
        RMSE(RMSE==0)=[];
        count=length(RMSE);
        %  plot RMSE
        x=zeros(1,count);
        for j=1:count
           x(j)=time/count*j;
        end
        plot(x,log(RMSE(1:count)),'Color',color{type(1)},...
            'LineStyle',linestyle{type(2)},'LineWidth',2);
        hold on;
%         axis([0,20,-25,2]);
        xlabel('time');ylabel('log(RMSE)');
        title('time to log(RMSE)');
        
        fid=fopen(filename,'a');
        fprintf(fid,'%3d\t%e\t%f\t',count,RMSE(end),time);
        if(type(2)-1)
            fprintf(fid,'QR method\n');
        else
            fprintf(fid,'not\n');
        end
        fclose(fid);
        
    end

    M.outplot1=@outplot1;
    function [v_RMSE]=outplot1(RMSE,time,nograd,type)
        if nargin<4
            type=[1,1];
        end
        tt={'b-','g-','r-','c-','m-','y-','k-';
            'b--','g--','r--','c--','m--','y--','k--';
            'bo','go','ro','co','mo','yo','ko';
            'bo-','go-','ro-','co-','mo-','yo-','ko-';};
        %b 蓝色  g 绿色  r 红色  c 青色  m 紫色  y 黄色  k 黑色
        RMSE(RMSE==0)=[];
        nograd(nograd==0)=[];
        count=length(RMSE);
        %%  画 RMSE
        x=zeros(1,count);
        for j=1:count
           x(j)=time/count*j;
        end
        subplot(2,3,4);
        plot(x,log(RMSE(1:count)),tt{type(1),type(2)},'LineWidth',2);
        axis([0,20,-25,2]);
        title('time to log(RMSE)');
        xlabel('time');ylabel('log(RMSE)');
        %%  画 v_RMSE
        v_RMSE=zeros(1,count-2);
        y=zeros(1,count-2);
        for i=1:count-2
            v_RMSE(i)=RMSE(i+1)/RMSE(i)-RMSE(i+2)/RMSE(i+1);
        end
        for j=1:count-2
           y(j)=time/(count-2)*j;
        end
        subplot(2,3,5);
        plot(y,v_RMSE(1:count-2),tt{type(1),type(2)},'LineWidth',2);
        %axis([0,20,-25,2]);
        title('time to v_{RMSE}');
        xlabel('time');ylabel('v_{RMSE}');
        %%   画梯度的范数
%         x=zeros(1,count);
%         for j=1:count
%            x(j)=time/count*j;
%         end
        subplot(2,3,3);
        plot(x,log(nograd(1:count)),tt{type(1),type(2)},'LineWidth',2);
        %axis([0,20,-25,2]);
        title('time to log(norm of grad)');
        xlabel('time');ylabel('log(norm of grad)');
    end

end

