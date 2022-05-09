function mu=GSL(Y)
%% The proposed GSL method

%% Normalized
[N,T]=size(Y);
Norm_y=norm(Y,'fro')/sqrt(N*T);
Y=Y/Norm_y;

%%  Initialization
a=1e-10;
b=a;
gamma=ones(N,1);
Z=ones(N,3)/3;
alpha=1;
iter = 0;
maxiter=100;
converged=false;

while ~converged
    %%  Update x
    deltal=[gamma(2:length(gamma));gamma(1)];
    deltar=[gamma(end);gamma(1:length(gamma)-1)];
    gamma=(  Z(:,1).*deltar + Z(:,2).*gamma +  Z(:,3).*deltal   );
    sigma=1./(alpha + gamma );
    mu= alpha*(sigma.*  Y  );
    
    %% Update alpha
    resid=Y-mu;
    alpha=( T*N + 2*a )/( 2*b +  norm(resid, 'fro')^2  + T*  sum(sigma)  );
    
    %%  Update gamma
    sum_mu=sum( mu.*conj(mu), 2);
    mu2=sum_mu + T*real(sigma);
    mu2l=[mu2(2:length(mu2));mu2(1)];
    mu2r=[mu2(end);mu2(1:length(mu2)-1)];
    Z_2=Z(:,2);
    Z_1l=[Z(2:length(mu2),1);Z(1,1)];
    Z_3r=[Z(end,3);Z(1:length(mu2)-1,3)];
    mu2_combine= Z_1l.*mu2l + Z_2.*mu2 + Z_3r.*mu2r;
    Z_shape=(Z_1l + Z_2 + Z_3r);
    c_k=(Z_shape*T+  2*a    );   %%%% a*1000;
    d_k=(  mu2_combine  +  2*b   );
    gamma=c_k./d_k;   
     
    if iter>=30
        [~,ind_sort]=sort(gamma,'descend');
        gamma(ind_sort(1:round(end*0.9)))=1e3;
        if iter>=maxiter-5
            gamma(ind_sort(1:round(end*0.8)))=1e5;
        end
    end
    if iter >= maxiter
        converged = true;
    end
    iter = iter + 1;
    
end

mu=mu*Norm_y;