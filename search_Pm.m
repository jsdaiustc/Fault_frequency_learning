function [Fn,unactive]=search_Pm(Pm,f_sample)



unactive=0;
[~,ind]=sort(Pm,'descend');
Num=min(10,length(Pm));
F_active= f_sample(ind(1:Num));

 %% remove closed values;
  for ii=1:Num
       F_choose = F_active(ii);
       for jj=ii+1:Num
          if abs(F_active(jj)-F_choose)< 10  
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
 
 [Z_max,Z_ind]=max(sum(Z,2));

 
 
%%
 z_choose=Z(Z_ind,:); 
 c_choose=C(Z_ind,:); 
 F_res=F_active(z_choose==1);   
 C_res=c_choose(z_choose==1);     
 Fn=sum(F_res)/sum(C_res);
 
 
 if Z_max<=1
     unactive=1;
 elseif Z_max==2
   ind11= find(z_choose==1);
   if abs(ind11(1)-ind11(2))>5
       unactive=1;
   end
 end
   
 end
 
 
 
 
 
 
 