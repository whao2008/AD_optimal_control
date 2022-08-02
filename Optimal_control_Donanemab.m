paras.initial=[107.6,0.47,0,0.5,27.6];
paras.lambda_beta=.18;
paras.K_Abeta=200;

paras.lambda_ptau=1e-2;
paras.K_ptau=10;

paras.lambdatauO=1.73;

paras.lambda_NdTO=5e-5;
paras.lambda_Ndptau=.1;
paras.K_Nd=1;

paras.lambda_CN=10;
paras.K_C=100;


paras.lambda_Ct=10;
delta = 0.001;
treat_age=70;
T=treat_age+1.5;
alpha=1;

epsilon0=.1;
umax=1.04*2;

x0=paras.initial; %reset initial conditions
opts=odeset('NonNegative',1,'AbsTol',1e-13,'RelTol',1e-13);% Impose nonnegative constraint on parameters
dt=1e-3;
%u=w;
tvec=[treat_age:dt:T];
M=(T-treat_age)/dt;
u=ones(M+1,1)*0;
solx=ode45(@(t,x) state_eqn(t,x,paras,tvec,u),tvec,x0,opts);
x = deval(solx,tvec)';
str={'A_\beta','\tau_p','\tau_o','N','C','u'};
for i=1:5
    if i<3
        subplot(3,2,i);
    elseif i>3
        subplot(3,2,i-1);
    else
        continue
    end
    hold on
    plot((tvec-treat_age)*365/7,(x(:,i)),'k','linewidth',3)
    set(gca, 'YScale', 'log')

    ylabel(str{i})
end
x_no_treatment=deval(solx,[treat_age:dt:T])';


x0=paras.initial; %reset initial conditions
u=ones(M+1,1);
x=zeros(M+1,5);
lambda=zeros(M+1,5);
w=zeros(size(w));
eta=1;
test=-1;
cur_it=0;
epsilon= epsilon0*(exp(-2*(tvec-T0)))';

%epsilon= epsilon0*ones(size(epsilon));
while(test < 0)

    oldu = u;
    oldx = x;
    oldlambda = lambda;

    solx=ode45(@(t,x) state_eqn(t,x,paras,tvec,u),tvec,x0,opts);
    x = deval(solx,tvec)';

    sollambda= ode45(@(t,lambda) adjoint_Eqn(t,lambda,x,tvec,u,paras,epsilon),[T T0],[alpha 0 0 0 alpha]);
    lambda = deval(sollambda,tvec)';

    Abeta=x(:,1);
    lambda1=lambda(:,1);

    temp=(Abeta.*lambda1)./(epsilon.*Abeta)/2;
    u1 = max(min(umax,temp),0);% max(0,temp);
    u=0.9*oldu+0.1*u1;
    test=min([delta*norm(u,1)-norm(oldu-u,1) delta*norm(x,1)-norm(oldx-x,1) delta*norm(lambda,1)-norm(oldlambda-lambda,1)]);
    if test>0
        J=alpha*Abeta(end)+alpha*x(end,5)+trapz(tvec,x(:,5))+trapz(tvec,epsilon.*u.^2.*Abeta);
    end
    cur_it=cur_it+1;
    if cur_it>300
        break
    end
end
hold on
str={'A_\beta','\tau_p','\tau_o','N','C','u'};
for i=1:5
    if i<3
        subplot(3,2,i);
    elseif i>3
        subplot(3,2,i-1);
    else
        continue
    end
    hold on
    plot((tvec-T0)*365/7,(x(:,i)),colors{1},'linewidth',3)
    xlabel('Time (week)')
    xlim([0 78])
    ylabel(str{i})
    set(gca, 'YScale', 'log')
    set(gca,'fontsize',20)
end


subplot(3,2,5);
hold on
plot((tvec-T0)*365/7,u,colors{1},'linewidth',3)
xlabel('Time (week)')
%    set(gca, 'YScale', 'log')

xlim([0 78])
ylabel('u')
set(gca,'fontsize',20)

subplot(3,2,6)
i=5;
hold on
plot((tvec-T0)*365/7,((x(:,i)-x_no_treatment(:,i))./x_no_treatment(:,i))*100,colors{1},'linewidth',3)
plot((tvec-T0)*365/7,zeros(size(tvec)),'k','linewidth',3)

xlabel('Time (week)')
xlim([0 78])
ylabel([str{i} ' change %'])
%set(gca, 'YScale', 'log')
set(gca,'fontsize',20)

