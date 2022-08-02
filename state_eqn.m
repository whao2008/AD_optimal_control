function dxdt=state_eqn(t,x,paras,tvec,u)
Abeta=x(1);
ptau=x(2);
tauo=x(3);
Nd=x(4);
C=x(5);
%CM=u(6);
 
u=pchip(tvec,u,t);

dxdt=[paras.lambda_beta*Abeta*(1-Abeta/paras.K_Abeta)-u*Abeta;
      paras.lambda_ptau*Abeta*(1-ptau/paras.K_ptau);
      paras.lambdatauO;
      (paras.lambda_NdTO*tauo+paras.lambda_Ndptau*ptau)*(1-Nd/paras.K_Nd);
      (paras.lambda_CN*Nd+paras.lambda_Ct*ptau)*(1-C/paras.K_C)];
      %paras.lambda_CMNd*Nd*CI*(1-CM/paras.K_CI)];
      
      