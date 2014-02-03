function [spkStress]=finiteStrainSMA(elast,poiss,beta,R,h,epsi,TempI,Temp,q,nrtol,nrmax)

%function [spkStress]=finiteStrainSMA(elast,poiss,beta,R,h,epsi,TempI,Temp,q,delta,tsteps,nrtol,nrmax)
%
%calculate second Piola-Kirchhoff stress tensors for SMA under shear
%deformation, according to algorithm in arghavani et al, where:
%
%elast : elastic/youngs modulus (MPa)
%poiss : Poisson ratio
%beta  : ? (MPa/C)
%R     : ? (MPa)
%h     : ? (MPa)
%epsi  : maximum transf. strain norm in uniaxial test (%)
%TempI : Temperature reference (?)
%Temp  : Temperature
%q     : small parameter ~ 10^-4 for C^t tensor estimate at nucleation
%nrtol : newton-raphson tolerance (default : 10^-4)
%nrmax : maximum newton-raphson steps (default : 100)
%
%written/tested WJB 10/11 Matlab R2010a/RHEL5

if (nargin<12)
    nrtol=1e-4;
end

if (nargin < 13)
    nrmax=100;
end


%sanity checks

if (nrmax < 1) || (nrmax>1000)
    error('please specify 1<nrmax<1000');
end

if (nrtol < 1e-6) || (nrtol > 1e-2)
    error('please specify 10^-6 < nrtol < 1e-2');
end

if (epsi < 0) || (epsi > 100)
    error('please specify 0 < epsi < 100');
end

if (poiss < 0) || (poiss > 1)
    error('please specify 0 < poiss < 1');
end

if (q < 1e-6) || (q > 1e-2)
    error('please specify 10^-6 < q < 1e-2');
end


%to fraction
epsi = epsi/100;

%tau
tempDelta = Temp-TempI;
tau = beta * ((tempDelta + abs(tempDelta)) / 2);

%lame's first const
lambda	= (elast*poiss) / ((1+poiss)*(1-2*poiss));

%lame's second const
mu	= (elast) / (2*(1+poiss));

%initialize shear fraction
kappa=0;

%create deformation gradient "F" and Green-Lagrange strain "E" tensors
defGrad = getFTensor(kappa);
strTens = getETensor(kappa);

%calculate right Cauchy-Green (RCG=="C") deformation tensor
rightC	= getCTensor(defGrad);

%transformative right CG @ nucleation
transRCG = ones(3,3);
transRCGprior = ones(3,3);

%calculate initial Cauchy stress compts
alpha = getCauchyStressCompts(rightC,ones(3,3),lambda,mu);

%init multipliers
gamma = 0;
deltaZeta=0;

%transformative strain @ nucleation
transE = zeros(3,3);
transEprior = zeros(3,3);

%test condition, algorithm start (nucleation)
test = (norm(getStressDirnTensor(alpha,rightC)) < (tau + R)) && (norm(transE)==0);

%main loops

for i=1:1
    
    %main algorithm
    
    if test
        
        %set transRCG to unity
        %we are done
        transRCG=ones(3,3);
        
    else
        
        %transRCG ==prior
        transRCG = transRCGprior;
        
        %normal tensor
        N = getNTensor(getStressDirnTensor(alpha,rightC),transE);
        
        %calculate Cauchy stress compts
        alpha = getCauchyStressCompts(rightC,transRCG,lambda,mu);
        
        %calculate SPK, limit function
        spkStress = getSPKTensor(alpha,rightC,transRCG);
        
        X = getXTensor(h,tau,transE,gamma,N);
        Y = getYTensor(rightC,spkStress,transRCG,X);
        
        f = getLimitFunction(Y,R);
        
        if (f < 0)
            
            %set solution
            %ie., transRCG remains unchanged from previous code block
            
        elseif getCompletionCondition(alpha,rightC,transRCG,tau,R)
            
            %set solution to unity
            transRCG = ones(3,3);
            
        else
            
            %solve
            
            %initialize transRCG for this step to prior; done
            
            %init deltaZeta
            deltaZeta=0;
            
            %init transRCG for nucleation case
            if (norm(transEprior)==0)
                transRCG = getInitTransRCG(alpha,rightC,q);
            end
            
            [transRCG,deltaZeta] = solveFirstSystem(transRCG,deltaZeta,rightC,nrtol,nrmax);
            
            if (norm(getETensorFromC(transRCG))<epsi)
                %then leave transRCG alone; solution for this step
            else
                %solve second system, using first system results as init
                transRCG = solveSecondSystem(transRCG,deltaZeta,rightC,nrtol,nrmax);
            end
        end
    end
    
    %calculate Cauchy stress compts
    alpha = getCauchyStressCompts(rightC,transRCG,lambda,mu);
    %calculate SPK stress
    spkStress = getSPKTensor(alpha,rightC,transRCG);
    %
    
end

end %main function

function alpha = getCauchyStressCompts(rightC,transRCG,lambda,mu)

%return cauchy stress compts using Saint-Venant Kirchhoff strain energy
%function

alpha(1) = (lambda/4)*(sum(rightC.*inv(transRCG))-3) - mu/2;
alpha(2) = mu/2;
alpha(3) = 0;

end

function output = getCompletionCondition(alpha,rightC,transRCG,tau,R)

%return the completion condition test

constA = getStressDirnTensor(alpha,rightC);
constB = log(transRCG);

output = norm(constA + R.*constB./norm(constB)) < tau;

end

function C = getCTensor(F)
%return right Cauchy-Green tensor

C	= F'*F;

end

function output = getDeviator(input)

%return the deviatoric

[eignvc eignvl]=eig(input);

output = input - (1/3).*eignvl;


end

function E = getETensor(kappa)

%return strain in shear deformation
E = [0 kappa 0; kappa kappa^2 0; 0 0 0];

end

function F = getFTensor(kappa)

%return deformation gradiant under shear
F=[1 kappa 0; 0 1 0; 0 0 1];

end

function output = getInitTransRCG(alpha,rightC,q)

%return initial transf RCG for nucleation case

const = getStressDirnTensor(alpha,rightC);
N = const./norm(const);

output = ones(3,3) + 2*q.*N;


end

function output = getLimitFunction(Y,R)

%return the limit function
output = norm(getDeviator(Y))-R;

end

function N = getNTensor(devProdCS,transE)

%return N

if (norm(transE)==0)
    N = devProdCS ./ norm(devProdCS);
else
    N = transE./ norm(transE);
end

end

function output = getSPKTensorElastic(alpha,rightC)

%return second Piola-Kirchhoff (SPK) stress tensor w/ transRCG==1
output = 2*(alpha(1)*rightC + alpha(2)*rightC + alpha(3)*rightC*rightC);


end

function output =getSPKTensor(alpha,rightC,transRCG)

%return SPK w/ transRCG != 1
constA = inv(transRCG);
constB = rightC / constA;
output = 2*(alpha(1) * constA + alpha(2)*constA*constB + alpha(3)*constA * (constB^2));

end

function output = getStressDirnTensor(alpha,rightC)

%return stress tensor direction
spkStress = getSPKTensorElastic(alpha,rightC);

prodCS = rightC*spkStress;
output = getDeviator(prodCS);


end

function output = getETensorFromC(C)

%return E based on C

output = (C - ones(3,3))./2;

end

function X = getXTensor(h,tau,transE,gamma,N)

%return the X tensor
X = h*transE + (tau + gamma)*N;

end

function Y = getYTensor(rightC,spkStress,transRCG,X)

%return the Y Tensor
Y = rightC*spkStress - transRCG*X;


end

function getElements()

%write tensor element string
%i,j,k,l runs 1->3

for i=1:3
    for j=1:3
        
        for k=1:3
            
            
            for l=1:3
                
                base = ['output(',num2str((i-1)*3+k),',',num2str((j-1)*3+l),')='];
                
                str=base;
                if delta(i,k)
                    str=[str,'0.5*x(',num2str(l),',',num2str(j),')'];
                end
                
                if delta(i,l)
                    
                    
                    if ~strcmp(str,base)
                        str=[str,'+'];
                    end
                    
                    str=[str,'0.5*x(',num2str(k),',',num2str(j),')'];
                end
                
                if delta(j,l)
                    
                    if ~strcmp(str,base)
                        str=[str,'+'];
                    end
                    
                    str=[str,'0.5*x(',num2str(i),',',num2str(k),')'];
                end
                
                if delta(k,j)
                    if ~strcmp(str,base)
                        str=[str,'+'];
                    end
                    
                    str=[str,'0.5*x(',num2str(i),',',num2str(l),')'];
                end
                
                if strcmp(str,base)
                    str=[str,'0;'];
                else
                    str = [str,';'];
                end
                
                disp(str);
                
            end
        end
    end
    
end

end

function output=delta(i,j)

output = 0;

if (i==j)
    output=1;
end


end

function output = getIsoScalarsDerivative(eigv,flag,ind1,ind2)

%return the scalar constants for tensor derivative
%function specified in flag

test 	= sum(eigv<0);

if (test>0)
    error(['eigv not > 0 for input matrix to getIsoScalarsDerivative']);
end

output=zeros(1,6);

switch flag
    case 1
        constA = (sqrt(eigv(ind1))-sqrt(eigv(ind2)));
    case 2
        constA = (log(eigv(ind1))-log(eigv(ind2)));
    case 3
        constA = (exp(eigv(ind1))-exp(eigv(ind2)));
end

constB = eigv(ind1)-eigv(ind2);
constC = constB^2;
constD = constC*constB;

switch flag
    case 1
        constE = 1/sqrt(eigv(ind2));
        constF = 1/sqrt(eigv(ind1));
    case 2
        constE = 1/eigv(ind2);
        constF = 1/eigv(ind1);
    case 3
        constE = exp(eigv(ind2));
        constF = exp(eigv(ind1));
end


output(1) =  constA/constC  - constE / constB;
output(2) = 2*eigv(ind2)*constA / constC - (eigv(ind1)+eigv(ind2))/constB * constE;
output(3) = 2*constA / constD - (constE+constF)/constC;
output(4) = eigv(ind2)*output(3);
output(5) = output(4);
output(6) = eigv(ind2)*output(4);


end

function output = getScalarsProductDerivative(eigvx,eigvy,ind1,ind2)

%return the scalar constants for tensor product derivative

output=zeros(1,6);

constA = eigvx(ind1)+eigvx(ind2);
constB = eigvx(ind1)-eigvx(ind2);
constC = constB^2;
constD = constC*constB;
constE = eigvy(ind1)-eigvy(ind2);
constF = (eigvy(ind1)/eigvx(ind2) + eigvy(ind2)/eigvx(ind1) - eigvy(ind1)/eigvx(ind1) + eigvy(ind2)/eigvx(ind2));
constG = (eigvy(ind2)/eigvx(ind1) - eigvy(ind2)/eigvx(ind2));
constH = (eigvy(ind1)/eigvx(ind2) - eigvy(ind2)/eigvx(ind2));


%CHECK ME; possible typo in paper Int. J. Numer. Engng 2004 : 61:880-895

output(1) =  constE/constC  - (1 / constB) * constG;
output(2) = 2*eigvx(ind2)*constE / constC + (constA / constB) * constG;
output(3) = 2*constE / constD -(1/constC)*constF ;
output(4) = 2*eigvx(ind2) * constE / constD +(eigvx(ind2)/constC)*constF + (1/constB)*constH;
output(5) = 2*eigvx(ind2) * constE / constD +(eigvx(ind2)/constC)*constF + (1/constB)*constG;
output(6) = 2*eigvx(ind2) * eigvx(ind2) * constE / constD + (eigvx(ind1)*eigvx(ind2))/constC *(eigvy(ind1)/eigvx(ind1) + eigvy(ind2)/eigvy(ind1)) - eigvx(ind2)*eigvx(ind2) / constC *(eigvy(ind1)/eigvx(ind1) + eigvy(ind2)/eigvy(ind1)) - constA/constB * (eigvy(ind2)/eigvx(ind2));

end



function output = getTensorSqDerivative(x)

%return the derivative dX^2/dX
%nb: this bit needs to be optimized

output = zeros(9,9);
output(1,1)=0.5*x(1,1)+0.5*x(1,1)+0.5*x(1,1)+0.5*x(1,1);
output(1,2)=0.5*x(2,1)+0.5*x(1,2);
output(1,3)=0.5*x(3,1)+0.5*x(1,3);
output(2,1)=0.5*x(2,1)+0.5*x(1,2);
output(2,2)=0;
output(2,3)=0;
output(3,1)=0.5*x(3,1)+0.5*x(1,3);
output(3,2)=0;
output(3,3)=0;
output(1,4)=0.5*x(1,2)+0.5*x(1,2);
output(1,5)=0.5*x(2,2)+0.5*x(1,1);
output(1,6)=0.5*x(3,2);
output(2,4)=0.5*x(2,2)+0.5*x(1,1);
output(2,5)=0.5*x(1,2)+0.5*x(1,2);
output(2,6)=0.5*x(1,3);
output(3,4)=0.5*x(3,2);
output(3,5)=0.5*x(1,3);
output(3,6)=0;
output(1,7)=0.5*x(1,3)+0.5*x(1,3);
output(1,8)=0.5*x(2,3);
output(1,9)=0.5*x(3,3)+0.5*x(1,1);
output(2,7)=0.5*x(2,3);
output(2,8)=0;
output(2,9)=0.5*x(1,2);
output(3,7)=0.5*x(3,3)+0.5*x(1,1);
output(3,8)=0.5*x(1,2);
output(3,9)=0.5*x(1,3)+0.5*x(1,3);
output(4,1)=0.5*x(2,1)+0.5*x(2,1);
output(4,2)=0.5*x(1,1)+0.5*x(2,2);
output(4,3)=0.5*x(2,3);
output(5,1)=0.5*x(1,1)+0.5*x(2,2);
output(5,2)=0.5*x(2,1)+0.5*x(2,1);
output(5,3)=0.5*x(3,1);
output(6,1)=0.5*x(2,3);
output(6,2)=0.5*x(3,1);
output(6,3)=0;
output(4,4)=0;
output(4,5)=0.5*x(1,2)+0.5*x(2,1);
output(4,6)=0;
output(5,4)=0.5*x(1,2)+0.5*x(2,1);
output(5,5)=0.5*x(2,2)+0.5*x(2,2)+0.5*x(2,2)+0.5*x(2,2);
output(5,6)=0.5*x(3,2)+0.5*x(2,3);
output(6,4)=0;
output(6,5)=0.5*x(3,2)+0.5*x(2,3);
output(6,6)=0;
output(4,7)=0;
output(4,8)=0.5*x(1,3);
output(4,9)=0.5*x(2,1);
output(5,7)=0.5*x(1,3);
output(5,8)=0.5*x(2,3)+0.5*x(2,3);
output(5,9)=0.5*x(3,3)+0.5*x(2,2);
output(6,7)=0.5*x(2,1);
output(6,8)=0.5*x(3,3)+0.5*x(2,2);
output(6,9)=0.5*x(2,3)+0.5*x(2,3);
output(7,1)=0.5*x(3,1)+0.5*x(3,1);
output(7,2)=0.5*x(3,2);
output(7,3)=0.5*x(1,1)+0.5*x(3,3);
output(8,1)=0.5*x(3,2);
output(8,2)=0;
output(8,3)=0.5*x(2,1);
output(9,1)=0.5*x(1,1)+0.5*x(3,3);
output(9,2)=0.5*x(2,1);
output(9,3)=0.5*x(3,1)+0.5*x(3,1);
output(7,4)=0;
output(7,5)=0.5*x(3,1);
output(7,6)=0.5*x(1,2);
output(8,4)=0.5*x(3,1);
output(8,5)=0.5*x(3,2)+0.5*x(3,2);
output(8,6)=0.5*x(2,2)+0.5*x(3,3);
output(9,4)=0.5*x(1,2);
output(9,5)=0.5*x(2,2)+0.5*x(3,3);
output(9,6)=0.5*x(3,2)+0.5*x(3,2);
output(7,7)=0;
output(7,8)=0;
output(7,9)=0.5*x(1,3)+0.5*x(3,1);
output(8,7)=0;
output(8,8)=0;
output(8,9)=0.5*x(2,3)+0.5*x(3,2);
output(9,7)=0.5*x(1,3)+0.5*x(3,1);
output(9,8)=0.5*x(2,3)+0.5*x(3,2);
output(9,9)=0.5*x(3,3)+0.5*x(3,3)+0.5*x(3,3)+0.5*x(3,3);




end

function [eigv,eigp] = getEigenprojections(input)

%return eigenprojections & values for
%constructing spectral expansions
%we could use eig but it's nice to see the math

%invariants

first 	= trace(input);
secon 	= 0.5 * ((trace(input)^2 - trace(input*input)));
third	= det(input);

%eigenvalues

numer 	= (-2*first^3 + 9*first*secon - 27*third) / 54;
denom 	= (first^2 - 3*secon) / 9;
theta	= acos(numer / sqrt(denom^3));


constA	= -2 * sqrt(denom);
constB	= first/3;

eigv(1) = constA * cos(theta/3) + constB;
eigv(2)	= constA * cos((theta+2*pi)/3) + constB;
eigv(3)	= constA * cos((theta-2*pi)/3) + constB;

%eigenprojections

eigp=zeros(3,3,3);

for i=1:3
    
    index 	= complement(i,[1,2,3]);
    
    if (eigv(i) ~= eigv(index(1))) && (eigv(i) ~= eigv(index(2)))
        
        eigp(i,:,:) = eigv(i) / (2*eigv(i)^3 - first * eigv(i)^2 - third) * (input^2 - (first-eigv(i)).*input + (third / eigv(i)));
        
    else if (eigv(i) ~= eigv(index(1))) && (eigv(index(1)) == eigv(index(2)))
            
            eigp(index(1),:,:) = ones(3,3) - eigv(i) / (2*eigv(i)^3 - first * eigv(i)^2 - third) * (input^2 - (first-eigv(i)).*input + (third / eigv(i)));
            
        else
            eigp(i,:,:) = ones(3,3);
            
        end
        
    end
    
end

end


function [output] = getLogFunction(input)


%return the log function of input using spectral decomp


[eigv,eigp] = getEigenprojections(input);

output 	= zeros(3,3);
output 	= output + log(eigv(1)).*eip(1,:,:);

if (eigv(1) ~= eigv(2))
    output = output + log(eigv(2)).*eip(2,:,:);
end


if (eigv(3) ~= eigv(2)) && (eigv(3) ~= eigv(1))
    output = output + log(eigv(3)).*eip(3,:,:);
end


end

function [output] = getExpFunction(input)


%return the exp function of input using spectral decomp


[eigv,eigp] = getEigenprojections(input);

output 	= zeros(3,3);
output 	= output + exp(eigv(1)).*eip(1,:,:);

if (eigv(1) ~= eigv(2))
    output = output + exp(eigv(2)).*eip(2,:,:);
end


if (eigv(3) ~= eigv(2)) && (eigv(3) ~= eigv(1))
    output = output + exp(eigv(3)).*eip(3,:,:);
end


end

function [output] = getSqrtFunction(input)


%return the sqrt function of input using spectral decomp


[eigv,eigp] = getEigenprojections(input);

test 	= sum(eigv<0);

if (test>0)
    error(['eigv not > 0 for input matrix :',num2str(input)]);
end

output 	= zeros(3,3);
output 	= output + sqrt(eigv(1)).*eip(1,:,:);

if (eigv(1) ~= eigv(2))
    output = output + sqrt(eigv(2)).*eip(2,:,:);
end


if (eigv(3) ~= eigv(2)) && (eigv(3) ~= eigv(1))
    output = output + sqrt(eigv(3)).*eip(3,:,:);
end


end

function [output] = getIsoTensorFunctionDeriv(input,flag)


%return the isotropic tensor function deriv of input using spectral decomp
%flag=1 -> sqrt
%flag=2 -> log
%flag=3 -> exp


%construct dx^2/dx

dx2dx = getTensorSqDerivative(input)

%get eigenprojections

[eigv,eigp] = getEigenprojections(input);

%sqrt & ln
if (flag < 3)
    test 	= sum(eigv<0);
    
    if (test>0)
        error(['eigv not > 0 for input matrix :',num2str(input)]);
    end
end

%init
output = zeros(9,9);


%build up products for construction

eigProjProduct(1,:,:) = kron(eigp(1,:,:),eigp(1,:,:));
eigProjProduct(2,:,:) = kron(eigp(2,:,:),eigp(2,:,:));
eigProjProduct(3,:,:) = kron(eigp(3,:,:),eigp(3,:,:));

%CHECK THIS
ident = ones(9,9);

%cartesian functions & derivatives
switch flag
    case 1
        constA = sqrt(eigv);
        constB = 1./constA;
    case 2
        constA = log(eigv);
        constB = 1./eigv
    case 3
        constA = exp(eigv);
        constB = constA;
end


if (eigv(1) ~= eigv(2)) && (eigv(2) ~= eigv(3))
    
    
    
    for i=1:3
        
        index 	= complement(i,[1,2,3]);
        
        output = output + constA(i) / ((eigv(i)-eigv(index(1))) * eigv(i)-eigv(index(2))) *(dx2dx - (eigv(index(1))+eigv(index(2))).*ident -(2*eigv(i)-eigv(index(1))-eigv(index(2))).*eigProjProduct(i,:,:) - (eigv(index(1))-eigv(index(2)))*(eigProjProduct(index(1),:,:)-eigProjProduct(index(2),:,:)) + constB(i).*eigProjProduct(i,:,:));
        
    end
    
    
    
elseif (eigv(1) == eigv(2)) && (eigv(2) == eigv(3))
    
    output = ident*constA;
    
    
else
    
    if (eigv(1) == eigv(2)) && (eigv(2) ~= eigv(3))
        
        ind1=3; ind2=1;
        
    elseif (eigv(2) == eigv(3)) && (eigv(1) ~= eigv(2))
        
        ind1=1; ind2=2;
        
    elseif (eigv(1) ==eig(3)) && (eigv(2) ~= eigv(3))
        
        ind1=2; ind2=3;
    end
    
    scalars = getIsoScalarsDerivative(eigv,flag,ind1,ind2)
    
    
    output = scalars(1).*dx2dx - scalars(2).*ident - scalars(3).*kron(input,input)+scalars(4).*kron(input,ones(3,3)) + scalars(5).*kron(ones(3,3),input) - scalars(6).*ident;
    
    
end

end

function [output] = getProductTensorFunctionDeriv(X,Y)


%return the partial derivative dy/dx of a tensor function product y=a.x


%construct dx^2/dx

dx2dx = getTensorSqDerivative(X)

%get x eigenprojections

[eigv,eigp] = getEigenprojections(X);

%and for y

[eigvy,eigpy]=getEigenProjections(Y);


%init
output = zeros(9,9);


%build up permuted products for construction

eigProjProduct = zeros(9,9,9);

eigProjProduct(1,:,:) = kron(eigp(1,:,:),eigp(1,:,:));
eigProjProduct(5,:,:) = kron(eigp(2,:,:),eigp(2,:,:));
eigProjProduct(9,:,:) = kron(eigp(3,:,:),eigp(3,:,:));
eigProjProduct(2,:,:) = kron(eigp(1,:,:),eigp(2,:,:));
eigProjProduct(3,:,:) = kron(eigp(1,:,:),eigp(3,:,:));
eigProjProduct(4,:,:) = kron(eigp(2,:,:),eigp(1,:,:));
eigProjProduct(6,:,:) = kron(eigp(2,:,:),eigp(3,:,:));
eigProjProduct(7,:,:) = kron(eigp(3,:,:),eigp(1,:,:));
eigProjProduct(8,:,:) = kron(eigp(3,:,:),eigp(2,:,:));

%CHECK THIS
ident = ones(9,9);


if (eigv(1) ~= eigv(2)) && (eigv(2) ~= eigv(3))
    
    
    
    for i=1:3
        
        index 	= complement(i,[1,2,3]);
        ind = (i-1)*3 + i;
        index2 = complement(ind,[1,5,9]);
        
        output = output + eigvy(i) / ((eigvx(i)-eigvx(index(1))) * eigvx(i)-eigv(index(2))) *(dx2dx - (eigvx(index(1))+eigvx(index(2))).*ident -(2*eigvx(i)-eigvx(index(1))-eigvx(index(2))).*eigProjProduct(ind,:,:) - (eigv(index(1))-eigv(index(2)))*(eigProjProduct(index2(1),:,:)-eigProjProduct(index2(2),:,:)));
        
        
        
        
    end
    for j=1:3
        for k=1:3
            output = output + eigvy(i) / eigvx(j) .* eigProjProduct((i-1)*3+j,:,:);
        end
    end
    
    
elseif (eigv(1) == eigv(2)) && (eigv(2) == eigv(3))
    
    output = ident*eigvy(1)/eigvx(2); 
    
else
    
    if (eigv(1) == eigv(2)) && (eigv(2) ~= eigv(3))
        
        ind1=3; ind2=1;
        
    elseif (eigv(2) == eigv(3)) && (eigv(1) ~= eigv(2))
        
        ind1=1; ind2=2;
        
    elseif (eigv(1) ==eig(3)) && (eigv(2) ~= eigv(3))
        
        ind1=2; ind2=3;
    end
    
    scalars = getScalarsProductDerivative(eigvx,eigvy,ind1,ind2)
    
    
    output = scalars(1).*dx2dx - scalars(2).*ident - scalars(3).*kron(X,X)
    +scalars(4).*kron(X,ones(3,3)) + scalars(5).*kron(ones(3,3),X) - scalars(6).*ident;
    
    
end

end


function [output] = getSqFunction(input)


%return the square of input using spectral decomp
%this is just a test of the logic


[eigv,eigp] = getEigenprojections(input);

test 	= sum(eigv<0);

if (test>0)
    error(['eigv not > 0 for input matrix :',num2str(input)]);
end

output 	= zeros(3,3);
output 	= output + eigv(1)*eigv(1).*eip(1,:);

if (eigv(1) ~= eigv(2))
    output = output + eigv(2)*eigv(2).*eip(2,:);
end


if (eigv(3) ~= eigv(2)) && (eigv(3) ~= eigv(1))
    output = output + eigv(3)*eigv(3).*eip(3,:);
end


end

function [transRCG,deltaZeta] = solveFirstSystem(transRCG,deltaZeta,rightC,nrtol,nrmax)

%perform newton raphson, unsaturated transformation case
%quit when average delta is < nrtol, abort if loops > nrmax
%FIXME : better convergence eg., adaptive step

for i=1:nrmax
    
    %calculate
    
    
    %update
    
    
    %test
    
    
    if 0
        break
    end
    
end


end

function transRCG = solveSecondSystem(transRCG,deltaZeta,rightC,nrtol,nrmax)

%perform newton raphson, the saturated transformation case
%quit when average delta is < nrtol, abort if loops > nrmax
%FIXME : as above

for i=1:nrmax
    
    %calculate
    
    
    %update
    
    
    %test
    
    if 0
        break
    end
    
end

end


function output = updateStepNR(jacobian,func)

%return parameter deltas
output = - jacobian \ func;

end



%
% TENSOR ELEMENTS FOR NR STEPS
%

function output =tensorComptU4(transRCG,rightC)

%return U_ijkl, partial derivative U wrt transRCG; fourth order
%chain rule; dU/dC^t = dU/dC * dC/dC^t

%isotropic

dUdC = getIsoTensorFunctionDeriv(rightC,1);

%product; rightC  = C^e * C^t

dCdCt = getProductTensorFunctionDeriv(transRCG,rightC);

output = dUdC * dCdCt;


end

function output = tensorComptU2(C)

%return U_ij

output = getSqrtFunction(C);

end

function output = tensorComptH4(W)
%return H_ijkl, partial derivative H wrt W; fourth order

output = getIsoTensorFunctionDeriv(W,2);

end


function output = tensorComptH2(transRCGPrior,U)

%return H_ij = log(W); W= U*inv(transRCGprior)*U

constB = U * inv(transRCGPrior) * U;

output = getLogFunction(constB);

end

function output = tensorComptG2(Q)

%return G_ij = exp(Q); Q=delta_gamma * inv(U^t) * A * inv(U^t)

output = getExpFunction(Q);

end

function output = tensorComptG4(Q)

%return G_ijkl = dG/dQ

output = getIsoTensorFunctionDeriv(G,3);

end

function output = tensorComptA2(Y,transRCG)

%return A_ij

const = getDeviator(Y);

output = 2.*const ./norm(const) * transRCG;

end


function output = tensorComptA4(transRCG,A)

%return A_ijkl = dA/dC^t

output = getProductTensorFunctionDeriv(transRCG,A);

%nb; compare with/use explicit expression

end

function output = tensorComptA4tilde(transRCG,A,rightC)

%return A_ijkl = dA/dC = dA/dC^t * dC^t/dC

dAdCt = getProductTensorFunctionDeriv(transRCG,A);
dCtdC = getProductTensorFunctionDeriv(rightC,transRCG);

%nb; an explicit expression exists, should use this

end

function output = tensorComptZ2(Y)

%return Z_ij

const = getDeviator(Y);
output = const ./norm(const);

end

function output = tensorComptZ4(Y,Z)

%return Z_ijkl
output = 1./norm(getDeviator(Y)).*(ones(9,9)-kron(Z,Z'));

end

function output = tensorComptN4(transE,N)

%return N_ijkl
output = 1./norm(transE).*(ones(9,9) - kron(N,N));

end

function output = tensorComptR12(U4,H2,U2,dGamma,A4)

%return R_12 = U4 * H2 * U2 + U2 * H2 * U4 + U2 * H4 * U2 + dGamma * A4

output = U4 * kron(H2*U2,ones(3,3)) + kron(U2*H2,ones(3,3)) * U4 + kron(U2,ones(3,3))*H4 * kron(U2,ones(3,3)) + dGamma .* A4;

end

function output = tensorComptR14(dGamma,Z4,transRCG,N)

%return R_14 = -2*dGamma*Z4 * C2 * (C*N-1/3 * ||C*N|| * I)

constA = -2*dGamma.*Z4*kron(transRCG,ones(3,3));
constB = (C*N - 1/3.*sum(sum(C*N)).*ones(3,3));

%contract

constC = constA .* repmat(constB,3,3);

output(1,1) = sum(sum(constC(1:3,1:3)));
output(1,2) = sum(sum(constC(1:3,4:6)));
output(1,3) = sum(sum(constC(1:3,7:9)));

output(2,1) = sum(sum(constC(4:6,1:3)));
output(2,2) = sum(sum(constC(4:6,4:6)));
output(2,3) = sum(sum(constC(4:6,7:9)));


output(3,1) = sum(sum(constC(7:9,1:3)));
output(3,2) = sum(sum(constC(7:9,4:6)));
output(3,3) = sum(sum(constC(7:9,7:9)));
end

function output = tensorComptR21(Z2,Y4tilde)

%return R_21 = Z2 * Y4tilde


constA = kron(Z2,ones(3,3))*Y4tilde;

output = constA(1:3,1:3)+constA(4:6,4:6)+constA(7:9,7:9);

end

function output = tensorComptR22(Z2,Y4)

%return R_22 = Z2 * Y4


constA = kron(Z2,ones(3,3))*Y4;

output = constA(1:3,1:3)+constA(4:6,4:6)+constA(7:9,7:9);

end

function output = tensorComptR24(Z,transRCG,N)

%return R_24 = -Z*C*N

output = -sum(sum(Z.*(transRCG*N)));

end

function output = tensorComptR12PT1(U4tilde,G4,G2,U2tilde)

%return R_12 = U4tilde * G2 * U2tilde + U2tilde * G4 * U2tilde + U2tilde *
%G * U4tilde for PT1 equation system

output = U4 * kron(G2*U2tilde,eye(3,3)) + kron(U2tilde,eye(3,3))*G4*kron(U2tilde,eye(3,3)) + kron(U2tilde*G2,eye(3,3)) * U4tilde;

end

function output = tensorComptR13PT1(U2tilde,A2,G2)

%return R_13 for PT1 system

output = U2tilde * U2tilde * A2 * U2tilde * G2 * U2tilde;

end

function output = tensorComptR14PT1(U2tilde,G4,Z4,N2,transRCG,dGamma)

%return R_14 = -2*dGamma*U2*U2*G4*U2*U2*Z4*C*(C*N-1/3 * ||C*N|| * I)

constA = (C*N - 1/3.*sum(sum(C*N)).*eye(3,3));

%tic

a=U2tilde;
b=G4;
c=Z4;
d=transRCG;

% expand 2nd  rank tensor U_im
%put m tiles in 3rd d
%i index has period 27 w/ matrix row index

foo(:,:,1)=kron(a(:,1),ones(9,27)); 
foo(:,:,2)=kron(a(:,2),ones(9,27)); 
foo(:,:,3)=kron(a(:,3),ones(9,27)); 

% expand 4th rank tensor G_mnpq
%put m tiles in 3rd d
%n index has period 27 w/ matrix column index
%p index has period 3 w/ block matrix row index
%q index has period 3 w/ block matrix column index

foo2(:,:,1)=[repmat(b(1:3,1:3),9,3), repmat(b(1:3,4:6),9,3), repmat(b(1:3,7:9),9,3)]; 
foo2(:,:,2)=[repmat(b(4:6,1:3),9,3), repmat(b(4:6,4:6),9,3), repmat(b(4:6,7:9),9,3)]; 
foo2(:,:,3)=[repmat(b(7:9,1:3),9,3), repmat(b(7:9,4:6),9,3), repmat(b(7:9,7:9),9,3)]; 

%expand 4th rank tensor Z_rmpq
%put m tiles in 3rd d
%r index has period 9 w/ block matrix row index
%p index has period 3 w/ block matrix row index
%q index has period 3 w/ block matrix column index

foo3(:,:,1)=repmat([repmat(c(1:3,1:3),1,9); repmat(c(4:6,1:3),1,9); repmat(c(7:9,1:3),1,9)],3,1); 
foo3(:,:,2)=repmat([repmat(c(1:3,4:6),1,9); repmat(c(4:6,4:6),1,9); repmat(c(7:9,4:6),1,9)],3,1); 
foo3(:,:,3)=repmat([repmat(c(1:3,7:9),1,9); repmat(c(4:6,7:9),1,9); repmat(c(7:9,7:9),1,9)],3,1); 

%expand 2nd rank tensor C_ms
%put m tiles in 3rd d
%s index has period 3 w/ block matrix column index

foo4(:,:,1)=repmat(kron(d(1,:),ones(27,3)),1,3); 
foo4(:,:,2)=repmat(kron(d(2,:),ones(27,3)),1,3); 
foo4(:,:,3)=repmat(kron(d(3,:),ones(27,3)),1,3);

%CONTRACT m; U_im * G_mnpq * Z_rmpq * C_ms ->  F_inrspq

F = sum(foo.*foo2.*foo3.*foo4,3);

%clear foo foo2 foo3 foo4

%reshape F_inrspq to suit contraction w/ n
%put n tiles in 3rd d

foo5(:,:,1) = F(:,1:9);
foo5(:,:,2) = F(:,10:18);
foo5(:,:,3) = F(:,19:27); 


%CONTRACT n; U_nj * F_inrspq = G_ijrspq

G = [(foo5(:,:,1).*a(1,1) + foo5(:,:,2).*a(2,1) + foo5(:,:,3).*a(3,1)), 	%j==1 tile
	(foo5(:,:,1).*a(1,2) + foo5(:,:,2).*a(2,2) + foo5(:,:,3).*a(3,2)), 	%j==2 tile
	(foo5(:,:,1).*a(1,3) + foo5(:,:,2).*a(2,3) + foo5(:,:,3).*a(3,3)) ];	%j==3 tile
%reshape G_ijrspq to suit contraction w/ r
%put r tiles in 3rd d; r has period 9 w/ row index

foo6(:,:,1) = G(1:3,:);			%r==1 tiles
foo6(:,:,2) = G(4:6,:);			%r==2 tiles
foo6(:,:,3) = G(7:9,:);			%r==3 tiles
foo6(:,:,4) = G(10:12,:);		%r==1 tiles
foo6(:,:,5) = G(13:15,:);		%r==2 tiles
foo6(:,:,6) = G(16:18,:);		%r==3 tiles
foo6(:,:,7) = G(19:21,:);		%r==1 tiles
foo6(:,:,8) = G(22:24,:);		%r==2 tiles
foo6(:,:,9) = G(25:27,:);		%r==3 tiles

%CONTRACT r; U_pr * G_ijrspq = H_ijpspq

H = [(foo6(:,:,1).*a(1,1) + foo6(:,:,2).*a(1,2) + foo6(:,:,3).*a(1,3));		%p==1 tiles
	(foo6(:,:,1).*a(2,1) + foo6(:,:,2).*a(2,2) + foo6(:,:,3).*a(2,3)); 		%p==2 tiles
	(foo6(:,:,1).*a(3,1) + foo6(:,:,2).*a(3,2) + foo6(:,:,3).*a(3,3)); 		%p==3 tiles
	(foo6(:,:,4).*a(1,1) + foo6(:,:,5).*a(1,2) + foo6(:,:,6).*a(1,3)); 		%p==1 tiles
	(foo6(:,:,4).*a(2,1) + foo6(:,:,5).*a(2,2) + foo6(:,:,6).*a(2,3)); 		%p==2 tiles
	(foo6(:,:,4).*a(3,1) + foo6(:,:,5).*a(3,2) + foo6(:,:,6).*a(3,3)); 		%p==3 tiles
	(foo6(:,:,7).*a(1,1) + foo6(:,:,8).*a(1,2) + foo6(:,:,9).*a(1,3)); 		%p==1 tiles
	(foo6(:,:,7).*a(2,1) + foo6(:,:,8).*a(2,2) + foo6(:,:,9).*a(2,3)); 		%p==2 tiles
	(foo6(:,:,7).*a(3,1) + foo6(:,:,8).*a(3,2) + foo6(:,:,9).*a(3,3))];		%p==3 tiles

%reshape H_ijpspq to suit contraction w/ s
%put s tiles in 3rd d; s has period 9 w/ row index

foo7(:,:,1) = H(:,1:3);			%s==1 tiles
foo7(:,:,2) = H(:,4:6);			%s==2 tiles
foo7(:,:,3) = H(:,7:9);			%s==3 tiles
foo7(:,:,4) = H(:,10:12);		%s==1 tiles
foo7(:,:,5) = H(:,13:15);		%s==2 tiles
foo7(:,:,6) = H(:,16:18);		%s==3 tiles
foo7(:,:,7) = H(:,19:21);		%s==1 tiles
foo7(:,:,8) = H(:,22:24);		%s==2 tiles
foo7(:,:,9) = H(:,25:27);		%s==3 tiles

%CONTRACT s; U_sq * H_ijpspq = I_ijpqpq

I = [(foo7(:,:,1).*a(1,1) + foo7(:,:,2).*a(2,1) + foo7(:,:,3).*a(3,1)), 	%q==1 tiles
	(foo7(:,:,1).*a(1,2) + foo7(:,:,2).*a(2,2) + foo7(:,:,3).*a(3,2)), 		%q==2 tiles
	(foo7(:,:,1).*a(1,3) + foo7(:,:,2).*a(2,3) + foo7(:,:,3).*a(3,3)), 		%q==3 tiles
	(foo7(:,:,4).*a(1,1) + foo7(:,:,5).*a(2,1) + foo7(:,:,6).*a(3,1)), 		%q==1 tiles
	(foo7(:,:,4).*a(1,2) + foo7(:,:,5).*a(2,2) + foo7(:,:,6).*a(3,2)), 		%q==2 tiles
	(foo7(:,:,4).*a(1,3) + foo7(:,:,5).*a(2,3) + foo7(:,:,6).*a(3,3)), 		%q==3 tiles
	(foo7(:,:,7).*a(1,1) + foo7(:,:,8).*a(2,1) + foo7(:,:,9).*a(3,1)), 		%q==1 tiles
	(foo7(:,:,7).*a(1,2) + foo7(:,:,8).*a(2,2) + foo7(:,:,9).*a(3,2)), 		%q==2 tiles
	(foo7(:,:,7).*a(1,3) + foo7(:,:,8).*a(2,3) + foo7(:,:,9).*a(3,3))]; 	%q==3 tiles

%CONTRACT p; I_ijpqpq = J_ijqq
%CONTRACT q; J_ijqq = output_ij

foo8 = kron(ones(9,9),constA).*I;
J = [sum(foo8(1:9,:)); sum(foo8(10:18,:)); sum(foo8(19:27,:))];

output = -2.*dGamma.*[sum(J(:,1:9),2), sum(J(:,10:18),2), sum(J(:,19:27),2)];


end







