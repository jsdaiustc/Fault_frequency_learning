function [res_x,f_sample] = fault_frequency_learning(y,f_sample,Fs)

N=length(f_sample);
[T,M]=size(y);
norm_y=norm(y,'fro')/sqrt(T*M);  y=y/norm_y;    % normalization


T_all=[0:1:T-1]';
A= exp(1i*2*pi*T_all/Fs*f_sample)/sqrt(N);   % dictionary matrix
B=(1i*2*pi*T_all/Fs).*A;                     % derivative of dictionary matrix
reslu=f_sample(2)-f_sample(1);

%% Initialization
maxiter=100;
beta0=1;
delta=ones(N,1);
a=1e-10;b=1e-10;
converged = false;
iter = 0;
tol=1e-20;
etc=20;

truncated_active=false;
AA=A'*A;


while ~converged && iter<maxiter
    %% update q(x)      Equation (23)
    Sigma=inv(beta0* AA  + diag(delta));
    mu = Sigma * (beta0 * (A' * y)   );
    
    
    %% update q(alpha)  Equation (24)
    resid=y-A*mu;
    term2=sum(diag( Sigma*AA));
    beta0=( T*M + a )/( b +  norm(resid(:), 'fro')^2+   M*real(term2) );
    
    
    %% update q(delta)  Equation (25)
    Exx = sum( mu.*conj(mu),2) + M*real(diag(Sigma));
    delta_last=delta;
    c_k=M+a;
    d_k=b+Exx;
    delta=c_k ./ d_k;
    
    
    %% off-grid  Equation (27)
    if ~truncated_active
        Pm=sum( mu.*conj(mu), 2);
        [~,sort_ind]=sort(Pm, 'descend');
        idx=sort_ind(1:etc);
        BHB = B(:,idx)' * B(:,idx);
        P2= M * Sigma(idx,idx);
        P = real( conj(BHB) .* ((mu(idx,:) * mu(idx,:)') +   P2   )  );
        v2= M * real(diag(B(:,idx)' * A * Sigma(:,idx)));
        v = sum( real(conj(mu(idx,:)) .* (B(:,idx)' * (y - A * mu))),2) -   v2;
        temp_grid=v./diag(P);
        temp_grid=temp_grid';
        theld=reslu/20*0.95^(iter);
        ind_small=find(abs(temp_grid)<theld);
        temp_grid(ind_small)=sign(temp_grid(ind_small))*theld;
        ind_unchang=find (abs(temp_grid)>reslu);
        temp_grid(ind_unchang)=sign(temp_grid(ind_unchang)) * reslu/20;
        f_sample(idx)=f_sample(idx) + temp_grid;
        F_active=exp(1i*2*pi*T_all/Fs*f_sample(idx))/sqrt(N);
        A(:,idx)=F_active;
        B(:,idx)=(1i*2*pi*T_all/Fs).*A(:,idx);
        AA(idx,:)=F_active'*A;
        AA(:,idx)=A'*F_active;
    else
        %%     Equation (29)
        BHB = B' * B;
        P2= M * Sigma;
        P = real( conj(BHB) .* ((mu * mu') +   P2   )  );
        v2= M * real(diag(B' * A * Sigma));
        v = sum( real(conj(mu) .* (B' * (y - A * mu))),2) -   v2;
        vect1=[1:N]';
        P=vect1'*P*vect1;
        temp_grid=(vect1'*v)/P;
        theld=reslu/20*0.95^(iter);
        ind_small=find(abs(temp_grid)<theld);
        temp_grid(ind_small)=sign(temp_grid(ind_small))*theld;
        ind_unchang=find (abs(temp_grid)>reslu);
        temp_grid(ind_unchang)=sign(temp_grid(ind_unchang)) * reslu/20;
        f_sample=f_sample + temp_grid*vect1';
        A=exp(1i*2*pi*T_all/Fs*f_sample)/sqrt(N);
        B=(1i*2*pi*T_all/Fs).*A;
        AA=A'*A;
    end
    
    %%  reset the new grid and Judge whether there is fault impulse signal or not
    if  iter==50 ||iter==80
        Pm=1./delta;
        [fn,unactive]=search_Pm(Pm,f_sample);    %% search the periodic peaks
        fn_all= fn*[0.5:0.5:20];                 %% reset the new grid
        if iter==50
            fn_all_temp= fn*[1:1:10];
            Aa= exp(1i*2*pi*T_all/Fs*fn_all_temp)/sqrt(N);
            ind_unactive= norm( y-Aa*(pinv(Aa)*y), 'fro' )^2/ norm(y,'fro')^2;
            threshold=length(fn_all_temp)/T;
            if  ind_unactive>(1-threshold*2) || unactive==1 
                delta=ones(N,1)*1e10;
                break;
            end          
        else
            fn_all= fn*[1:1:10];
            truncated_active=true;
        end
        f_sample=fn_all;
        A= exp(1i*2*pi*T_all/Fs*f_sample)/sqrt(N);
        B=(1i*2*pi*T_all/Fs).*A;
        AA=A'*A;
        N=length(f_sample);
        delta=ones(N,1)*1;
        etc=min(etc,N);
        delta_last=100;
    end
    
    
    %% prune out the irrelevant frequency component
    if iter<50
        theld2=20;
        delta_inv=1./abs(delta);
        ind_remove= find (delta_inv<(max(delta_inv)/theld2));
        A(:,ind_remove)=[];
        B(:,ind_remove)=[];
        AA(:,ind_remove)=[]; AA(ind_remove,:)=[];
        delta(ind_remove)=[];
        delta_last(ind_remove)=[];
        f_sample(ind_remove)=[];
        N=length(f_sample);
        etc=min(etc,N);
    end
    
    %% stopping criteria
    erro=  max(max(abs(delta - delta_last)));
    if erro < tol || iter >= maxiter
        converged = true;
    end
    iter = iter + 1;
    
end

res_x=1./abs(delta);
res_x=sqrt(res_x)*norm_y;






