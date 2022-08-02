function dlambdadt=adjoint_Eqn(t,lambda,x,tvec,u,paras,epsilon)

x=interp1(tvec,x,t);
Abeta=x(1);
ptau=x(2);
tauo=x(3);
Nd=x(4);
C=x(5);
u=pchip(tvec,u,t);
epsilon=pchip(tvec,epsilon,t);


Gx=[paras.lambda_beta*(1-2*Abeta/paras.K_Abeta)-u 0 0 0 0;...
    paras.lambda_ptau*(1-ptau/paras.K_ptau)       -paras.lambda_ptau*Abeta/paras.K_ptau 0 0 0;...
    0 0 0 0 0;...
    0 paras.lambda_Ndptau*(1-Nd/paras.K_Nd) paras.lambda_NdTO*(1-Nd/paras.K_Nd) -(paras.lambda_NdTO*tauo+paras.lambda_Ndptau*ptau)/paras.K_Nd 0;...
    0 paras.lambda_Ct*(1-C/paras.K_C) 0 paras.lambda_CN*(1-C/paras.K_C) -(paras.lambda_CN*Nd+paras.lambda_Ct*ptau)/paras.K_C];
% 
 dlambdadt=-[epsilon*u.^2;0;0;0;1]-Gx'*lambda;
%dlambdadt=zeros(5,1);

% %dlambdadt(1)=-(paras.lambda_beta*(1-2*Abeta/paras.K_Abeta)-u)*lambda(1)-paras.alpha*u;
% dlambdadt(1)=-(paras.lambda_beta*(1-2*Abeta/paras.K_Abeta)-u)*lambda(1)-paras.alpha;
% dlambdadt(2)=-(paras.lambda_ptau*(1-ptau/paras.K_ptau))*lambda(1)+paras.lambda_ptau*Abeta/paras.K_ptau*lambda(2);
% dlambdadt(3)=0;
% dlambdadt(4)=-paras.lambda_Ndptau*(1-Nd/paras.K_Nd)*lambda(2)-paras.lambda_NdTO*(1-Nd/paras.K_Nd)*lambda(3)+(paras.lambda_NdTO*tauo+paras.lambda_Ndptau*ptau)/paras.K_Nd*lambda(4);
% dlambdadt(5)=-paras.alpha1-paras.lambda_C*(1-C/paras.K_C)*lambda(4)+paras.lambda_C*Nd/paras.K_C*lambda(5);
% 
