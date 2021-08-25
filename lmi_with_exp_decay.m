function [ K , P ] = lmi_with_exp_decay( A , B , p , alpha, r )
%LMI_QUASI_LPV_DESIGN Summary of this function goes here
%   Detailed explanation goes here

n = size(A,1);
m = size(B,2);

nv = size(A,3);
np = size(p,1);

setlmis([]);

XX = lmivar( 1 , [n  1] );
FF = lmivar( 1 , [n  1] );
GG = lmivar( 2 , [m  n] );

lmi = newlmi;

% F > 0
lmiterm([-lmi 1 1 FF],1,1);

%He{AF + BG + alpha F} < 0
for i=1:nv
    lmi = newlmi;
    lmiterm( [lmi 1 1  FF] ,  A(:,:,i) , 1 , 's' );
    lmiterm( [lmi 1 1  FF] ,  alpha    , 1 , 's' );
    lmiterm( [lmi 1 1  GG] ,  B(:,:,i) , 1 , 's' );
end

for i=1:np
        lmi = newlmi;
        lmiterm([-lmi 1 1  0] , 1          );
        lmiterm([-lmi 1 2 FF] , p(i,:) , 1 );
        lmiterm([-lmi 2 2 FF] , 1      , 1 );
end

% [ -rF   A(d)F + B(d)G ]  < 0  
% [  *       -rF        ]  

for i=1:nv
   lmi = newlmi; 
   lmiterm( [lmi 1 1  FF] , -r , 1 );
   lmiterm( [lmi 1 2  FF] , A(:,:,i) , 1 );
   lmiterm( [lmi 1 2  GG] , B(:,:,i) , 1 ); 
   lmiterm( [lmi 2 2  FF] , -r , 1 );
end

lmi = newlmi;
lmiterm([-lmi 1 1 XX] , 1 , 1  );
lmiterm([-lmi 1 2  0] , eye(n) );
lmiterm([-lmi 2 2 FF] , 1 , 1  );

lmis=getlmis;

no=decnbr(lmis);
co=zeros(no,1);

for j=1:no
    Xj = defcx(lmis,j,XX);
    co(j) = trace(Xj);
end

[~,gopt] = mincx( lmis , co , [1e-5 2000 1e9 100 0]);
  
F = dec2mat(lmis,gopt,FF);
G = dec2mat(lmis,gopt,GG);

K = G/F;
P = inv(F);

end