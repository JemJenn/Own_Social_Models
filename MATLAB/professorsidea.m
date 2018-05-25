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

p1c(n1c,n2c,n1,beta1,beta2,alpha1,alpha2,b,c)= n1*( n1c*alpha1*(b-c)/n1  +  n1d*(-c) )+(n2)*( n2c*beta1*beta2*(b-c)  +  n2c*beta2*(1-beta1)*b  +  n2c*beta1*(1-beta2)*(-c) + n2d*beta1*(-c));
p1d(n1c,n2c,n1,beta1,beta2,alpha1,alpha2,b,c)= n1*( n1c*b ) + (n2)*( n2c*beta2*b );
p2c(n1c,n2c,n1,beta1,beta2,alpha1,alpha2,b,c) = (n2)*( n2c*alpha2*(b-c)/n2 + n2d*(-c) ) + (n1)*( n1c*beta1*beta2*(b-c) + n1c*(1-beta1)*beta2*(-c) + n1c*(1-beta2)*beta1*(b) + n1d*beta2*(-c));
p2d(n1c,n2c,n1,beta1,beta2,alpha1,alpha2,b,c) = (n2)*b*n2c + (n1) * n1c * beta1*(b);

p1bar(n1c,n2c,n1,beta1,beta2,alpha1,alpha2,b,c) = n1c*p1c(n1c,n2c,n1,beta1,beta2,alpha1,alpha2,b,c) + n1d*p1d(n1c,n2c,n1,beta1,beta2,alpha1,alpha2,b,c)
p2bar(n1c,n2c,n1,beta1,beta2,alpha1,alpha2,b,c) = n2c*p2c(n1c,n2c,n1,beta1,beta2,alpha1,alpha2,b,c) + n2d*p2d(n1c,n2c,n1,beta1,beta2,alpha1,alpha2,b,c);
pbar(n1c,n2c,n1,beta1,beta2,alpha1,alpha2,b,c) = n1*p1bar(n1c,n2c,n1,beta1,beta2,alpha1,alpha2,b,c)+n2*p2bar(n1c,n2c,n1,beta1,beta2,alpha1,alpha2,b,c);

n1cdot(n1c,n2c,n1,beta1,beta2,alpha1,alpha2,b,c)=n1c*(p1c-p1bar)
%n1ddot(n1c,n2c,n1,beta1,beta2,alpha1,alpha2,b,c)=n1d*(p1d-p1bar);
n2cdot(n1c,n2c,n1,beta1,beta2,alpha1,alpha2,b,c)=n2c*(p2c-p2bar)
%n2ddot(n1c,n2c,n1,beta1,beta2,alpha1,alpha2,b,c)=n2d*(p2d-p2bar)
n1dot(n1c,n2c,n1,beta1,beta2,alpha1,alpha2,b,c)=n1*(p1bar-pbar)
%n2dot(n1c,n2c,n1,beta1,beta2,alpha1,alpha2,b,c)=n2*T*(p2bar-pbar)
%assume(n1*(n1c+n1d) + n2*(n2c+n2d) == 1);
% 
% n1cnull=symfun(solve(n1cdot==0,n1c),[n2c,n1,beta1,beta2,alpha1,alpha2,b,c]);
% n2cnull=symfun(solve(n2cdot==0,n2c),[n1c,n1,beta1,beta2,alpha1,alpha2,b,c]);
% %   [grid1,grid2]=meshgrid(0:0.05:1,0:0.05:1);
% %quiver(grid1,grid2,n1cnull(grid2))
% 
% [n1starstar(beta1,beta2,alpha1,alpha2,b,c),n1cstarstar(beta1,beta2,alpha1,alpha2,b,c),n2cstarstar(beta1,beta2,alpha1,alpha2,b,c)]  = solve(n1dot==0,n1cdot==0,n2cdot==0,n1,n1c,n2c);
 [n1star,n1cstar,n2cstar] = solve(n1dot==0,n1cdot==0,n2cdot==0,n1,n1c,n2c);
% % 
% n1cstar=simplify(n1cstar);
% n2cstar=simplify(n2cstar);
% n1star =simplify(n1star);
defectionPoints = [];
extinctionPoints =[];
locationsOfInterest =[];
locationsOfParticularInterest=[];
totalCooperation=[];
partialCooperation=[];
problemPoints=[]
%alpha2=alpha1;
interiorPoints=[];
for fp=1:length(n1cstar)    
    if (n1star(fp)~=0 & n1star(fp)~=1 & n1cstar(fp)~=0 & n1cstar(fp)~=1 & n2cstar(fp)~= 0 & n2cstar(fp)~=1)
        interiorPoints = [interiorPoints;[fp,n1cstar(fp),n2cstar(fp),n1star(fp)]]
    end    
    if (n1star(fp)==0 | n1star(fp)==1)
        extinctionPoints=[extinctionPoints ; [fp,n1cstar(fp),n2cstar(fp),n1star(fp)]];
    elseif (n1cstar(fp) == 0 & n2cstar(fp) ==0 )
        defectionPoints = [defectionPoints ; [fp,n1cstar(fp),n2cstar(fp),n1star(fp)]];
    elseif (n1cstar(fp)~=0 & n2cstar(fp)~=0)
        if (n1cstar(fp)==1 & n2cstar(fp)==1)
            totalCooperation = [totalCooperation;[fp,n1cstar(fp),n2cstar(fp),n1star(fp)]];
        else
            partialCooperation = [partialCooperation;[fp,n1cstar(fp),n2cstar(fp),n1star(fp)]];
        end
        locationsOfParticularInterest = [locationsOfParticularInterest;[fp,n1cstar(fp),n2cstar(fp),n1star(fp)]];
    else
        locationsOfInterest = [locationsOfInterest;[fp,n1cstar(fp),n2cstar(fp),n1star(fp)]];
    end
end
extinctionPoints
defectionPoints
locationsOfInterest
partialCooperation
totalCooperation
interiorPoints
locationsOfParticularInterest;
vars = [n1c,n2c,n1];
dots = [n1cdot,n2cdot,n1dot];
shayni = jacobian(dots,vars);
vecsAndlams=[];

% for fp =1:length(n1cstar)
%     fp
%     location = [n1cstar(fp),n2cstar(fp),n1star(fp)]
%     if location(3)==1 | location(3)==0
%         continue
%     end
%     if fp==1|fp==19|fp==27|fp==28
%         continue
%     end
%     %rank(shayni(n1cstar(fp),0,n2cstar(fp),0,n1star(fp),beta1,beta2,alpha1,alpha2,b,c))
%     eigenvalues = eig(shayni(n1cstar(fp),n2cstar(fp),n1star(fp),beta1,beta2,alpha1,alpha2,b,c));
%     %vec= simplify(vec)
%     lam = simplify(eigenvalues)
% end


%parameterise and numerical
% eigvec=[]
% alpha1=4
% alpha2=4
% beta1=0.5
% beta2=0.5
% b=2
% c=0.5
for alpha1 = 1.05:0.5:4.05
    for alpha2 = 1.05:0.5:4.05
        for beta1 = 0.1:0.2:0.9
            for beta2 = 0.1:0.2:0.9
                for b = 1.1:1:4.1
                    for c=0.51:0.5:3.1
                        for fp =3:6
                            for i = 1:50
                                [n1starnum,n1cstarnum,n2cstarnum]  = vpasolve([n1dot(n1c,n2c,n1,beta1,beta2,alpha1,alpha2,b,c)==0,n1cdot(n1c,n2c,n1,beta1,beta2,alpha1,alpha2,b,c)==0,n2cdot(n1c,n2c,n1,beta1,beta2,alpha1,alpha2,b,c)==0],[n1,n1c,n2c],'random',true);
                                if (n1starnum~=0 & n1starnum~=1 & n1cstarnum~=0 & n1cstarnum~=1 & n2cstarnum~= 0 & n2cstarnum~=1)
                                    
                                    
                                    fp;
                                    %                             n1cstar = n1cstarstar(beta1,beta2,alpha1,alpha2,b,c);
                                    %                             n2cstar = n2cstarstar(beta1,beta2,alpha1,alpha2,b,c);
                                    %                             n1star = n1starstar(beta1,beta2,alpha1,alpha2,b,c);
                                    %
                                    %                             z=0.5;
                                    location = [n1cstarnum,n2cstarnum,n1starnum]
                                    %                             if location(3)==1 | location(3)==0
                                    %                                 continue
                                    %                             end
                                    %                             %rank(shayni(n1cstar(fp),0,n2cstar(fp),0,n1star(fp),beta1,beta2,alpha1,alpha2,b,c))
                                    %eigenvalues = eval(eig(shayni(n1cstar(fp),n2cstar(fp),n1star(fp),beta1,beta2,alpha1,alpha2,b,c)))
                                    eigenvalues = eval(eig(shayni(n1cstarnum,n2cstarnum,n1starnum,beta1,beta2,alpha1,alpha2,b,c)));
                                    if (imag(eigenvalues(1))~=0 | imag(eigenvalues(2))~=0 |imag(eigenvalues(3))~=0)
                                        amazing =location;
                                        if amazing(1)>0 & amazing(1)<1 & amazing(2)>0 & amazing(3)<1 & amazing(3)>0 & amazing(2)<1
                                            vamazingpt = amazing
                                            eigenvalues=eigenvalues
                                            param=[beta1,beta2,alpha1,alpha2,b,c]
                                            
                                        end
                                    end
                                end
                            end
                
                            %vec= simplify(vec)
                            %lam = simplify(lam)
                        end
                    end
                end
            end
        end
    end
end

% %#one-dimensional
% jacqui = jacobian(n1dot,n1)
% n1star1dim11 = simplify(solve(n1dot(1,1,n1,beta1,beta2,alpha1,alpha2,b,c)==0,n1))
% n1star1dim10 = simplify(solve(n1dot(1,0,n1,beta1,beta2,alpha1,alpha2,b,c)==0,n1))
% n1star1dim01 = simplify(solve(n1dot(0,1,n1,beta1,beta2,alpha1,alpha2,b,c)==0,n1))
% n1star1dim00 = simplify(solve(n1dot(0,0,n1,beta1,beta2,alpha1,alpha2,b,c)==0,n1))
% 
% lam110 = simplify(eig(jacqui(1,1,0,beta1,beta2,alpha1,alpha2,b,c)))
% lam111 = simplify(eig(jacqui(1,1,1,beta1,beta2,alpha1,alpha2,b,c)))
% lam11m = simplify(eig(jacqui(1,1,n1star1dim11(3),beta1,beta2,alpha1,alpha2,b,c)))
% 
% lam100 = simplify(eig(jacqui(1,0,0,beta1,beta2,alpha1,alpha2,b,c)))
% lam101 = simplify(eig(jacqui(1,0,1,beta1,beta2,alpha1,alpha2,b,c)))
% lam10m = simplify(eig(jacqui(1,0,n1star1dim10(3),beta1,beta2,alpha1,alpha2,b,c)))
% 
% lam010 = simplify(eig(jacqui(0,1,0,beta1,beta2,alpha1,alpha2,b,c)))
% lam011 = simplify(eig(jacqui(0,1,1,beta1,beta2,alpha1,alpha2,b,c)))
% lam01m = simplify(eig(jacqui(0,1,n1star1dim01(3),beta1,beta2,alpha1,alpha2,b,c)))
% 
% lam000 = eig(jacqui(0,0,n1star1dim00,beta1,beta2,alpha1,alpha2,b,c))