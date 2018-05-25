%SIMPLIFY AND TEST
clc
clear all
syms n1 n2 n1c n1d n2c n2d beta1 beta2 alpha1 alpha2 b c T;

%n2=1-n1;
m=n1;
n2 = 1-n1;
n1d=1-n1c;
n2d=1-n2c;
T=1;
%n1=n1c+n1d;
%n2=n2c+n2d;
%T = 0.1; %separation of scales in replicator dynamics
    b=12;
    c=.5
    stablefull = []
    counter=1;
    for alpha1 = 6.1:10:26.1
        for alpha2=10.1:10:20.1
            for beta1=0.6:.4:1
                for beta2 = 0.6:.4:1
                    
                    p1c(n1c,n2c,n1)= n1*( n1c*(b-c)  +  n1d*(-c) )+(n2)*( n2c*beta1*beta2*alpha1*n1*(b-c)  +  n2c*beta2*(1-beta1)*b  +  n2c*beta1*(1-beta2)*(-c) + n2d*beta1*(-c));
                    p1d(n1c,n2c,n1)= n1*( n1c*b ) + (n2)*( n2c*beta2*b );
                    p2c(n1c,n2c,n1) = (n2)*( n2c*(b-c) + n2d*(-c) ) + (n1)*( n1c*beta1*beta2*alpha2*n2*(b-c) + n1c*(1-beta1)*beta2*(-c) + n1c*(1-beta2)*beta1*(b) + n1d*beta2*(-c));
                    p2d(n1c,n2c,n1) = (n2)*b*n2c + (n1) * n1c * beta1*(b);
                    
                    p1bar(n1c,n2c,n1) = n1c*p1c(n1c,n2c,n1) + n1d*p1d(n1c,n2c,n1)
                    p2bar(n1c,n2c,n1) = n2c*p2c(n1c,n2c,n1) + n2d*p2d(n1c,n2c,n1);
                    pbar(n1c,n2c,n1) = n1*p1bar(n1c,n2c,n1)+n2*p2bar(n1c,n2c,n1);
                    
                    n1cdot(n1c,n2c,n1)=n1c*(p1c-p1bar);
                    %n1ddot(n1c,n2c,n1)=n1d*(p1d-p1bar);
                    n2cdot(n1c,n2c,n1)=n2c*(p2c-p2bar);
                    %n2ddot(n1c,n2c,n1)=n2d*(p2d-p2bar)
                    n1dot(n1c,n2c,n1)=n1*T*(p1bar-pbar);
                    
                    
                    vars = [n1c,n2c,n1];
                    dots = [n1cdot,n2cdot,n1dot];
                    shayni = jacobian(dots,vars);
                     
                    [n1cstar,n2cstar,n1star]=solve(n1cdot==0,n2cdot==0,n1dot==1);
                    for fp = 1:length(n1cstar)
                        fp
                        if real(double(n1cstar(fp)))<1 &real(double(n1cstar(fp)))>0 &real( double(n2cstar(fp)))<1 &real(double(n2cstar(fp)))>0 &real(double(n1star(fp)))<1 &real(double(n1star(fp)))>0
                            lam = eig(shayni(n1cstar(fp),n2cstar(fp),n1star(fp)));
                            internal = fp
                            if real(double(lam(1)))<0 & real(double(lam(2)))<0 & real(double(lam(3)))<0
                                
                                
                                stablefull=[stablefull;[fp,alpha1,alpha2,beta1,beta2,n1star(fp),n1cstar(fp),n2cstar(fp)]]
                                
                            end
                        end
                      
                    end
                end
            end
        end
    end
stablefull
% for fp =2:length(n1cstar)
%     fp
%     if fp==10 | fp ==12
%         continue
%     end
%     location = [n1cstar(fp),n2cstar(fp),n1star(fp)]
%     
%     %rank(shayni(n1cstar(fp),0,n2cstar(fp),0,n1star(fp)))
%     [veceasy,lameasy] = eig(shayni(n1cstar(fp),n2cstar(fp),n1star(fp),beta1,beta1,alpha1,alpha1,b,c))
%     [vec,lam] = eig(shayni(n1cstar(fp),n2cstar(fp),n1star(fp)))
%     vec= simplify(vec)
%     lam = simplify(lam)
% end
% 
% %#one-dimensional
% % jacqui = jacobian(n1dot,n1)
% % n1star1dim11 = simple(solve(n1dot(1,1,n1)==0,n1))
% % n1star1dim10 = simple(solve(n1dot(1,0,n1)==0,n1))
% % n1star1dim01 = simple(solve(n1dot(0,1,n1)==0,n1))
% % n1star1dim00 = simple(solve(n1dot(0,0,n1)==0,n1))
% % 
% % lam110 = simple(eig(jacqui(1,1,0)))
% % lam111 = simple(eig(jacqui(1,1,1)))
% % lam11m = simple(eig(jacqui(1,1,n1star1dim11(3))))
% % 
% % lam100 = simple(eig(jacqui(1,0,0)))
% % lam101 = simple(eig(jacqui(1,0,1)))
% % lam10m = simple(eig(jacqui(1,0,n1star1dim10(3))))
% % 
% % lam010 = simple(eig(jacqui(0,1,0)))
% % lam011 = simple(eig(jacqui(0,1,1)))
% % lam01m = simple(eig(jacqui(0,1,n1star1dim01(3))))
% % 
% % lam000 = eig(jacqui(0,0,n1star1dim00))