function Fn=search_P(Y,Fs)


 
N=length(Y);
F = ([1:N]-1)*Fs/N; 
  
y_diff_hibert = abs(hilbert([0;diff(Y)]));
y_diff_corr_fft_abs=abs(fft(y_diff_hibert));    %%%% method 1: diff  Hilbert 
y_diff_corr_fft_abs=y_diff_corr_fft_abs(1:round(N/2));

set_remove=[1];
for ii=2:N
   if y_diff_corr_fft_abs(  ii  )< y_diff_corr_fft_abs(  ii-1  )
       set_remove=[set_remove,ii];
   else
       break;
   end
 end
 y_diff_corr_fft_abs(set_remove)=0;
 [~,ind]=sort(y_diff_corr_fft_abs,'descend');
 
 Num=10;
 F_active= F(ind(1:Num));

 %% remove closed values;
  for ii=1:Num
       F_choose = F_active(ii);
       for jj=ii+1:Num
          if abs(F_active(jj)-F_choose)< 10 %F_choose*0.05
              F_active(jj)=0;
          end
       end      
  end
 F_active(F_active==0)=[];
 Num=length(F_active);
       
      
%% chose the optimal ind 
Z=zeros(Num,Num); 
C=zeros(Num,Num); 
 for ii=1:Num
    F_choose = F_active(ii);
    for jj=1: Num
       modmod=  mod(F_active(jj),F_choose);
       if    min(modmod, abs(F_choose-modmod) )       <  min(5, F_choose*0.05)   &&   F_active(jj)/F_choose<=10
          Z(ii,jj)=1;
          C(ii,jj)= round(F_active(jj)/F_choose);
       end
    end
 end
 
 [~,Z_ind]=max(sum(Z,2));
%%
%  sum_Z=sum(Z,2);
%  Z_ind=find(sum_Z==max(sum_Z));
%  [~,temp_ind]=min(F_active(Z_ind));
%  Z_ind=Z_ind(temp_ind);
%%
 
 
 z_choose=Z(Z_ind,:); 
 c_choose=C(Z_ind,:); 
 F_res=F_active(z_choose==1);   
 C_res=c_choose(z_choose==1);     
 
 Fn=sum(F_res)/sum(C_res);
 
 
 
 end
 
 
 

 
 
 
 
 
 
 
 
 
 
%  ind= min(ind(1:4));
%  Fn=F(ind);




% 
% 
% figure;
% stem(F(1:round(N/2)),  (y_diff_corr_fft_abs(1:round(N/2)))/ max(y_diff_corr_fft_abs))
% axis([0 1500 0 1])
% title(['FFT of direction diff signal'])
% 
%  
 
 
 
 
 