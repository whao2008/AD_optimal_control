ParasValues;
delta = 0.001;
T=100;
umax=.12;
epsilon0=10;
year=0;
colors={'g','r','k','m'};
for T0=[10:20:70 100]
    x0=paras.initial; %reset initial conditions
    opts=odeset('NonNegative',1);% Impose nonnegative constraint on parameters
    dt=1e-1;
    %u=w;
    tvec=[0:dt:T0];
    M=T0/dt;
    u=zeros(M+1,1);
    solx=ode45(@(t,x) state_eqn(t,x,paras,tvec,u),tvec,x0,opts);
    x = deval(solx,tvec)';
    str={'A_\beta','\tau_p','\tau_o','N','C','u'};
    if year
        for i=1:5
            subplot(3,2,i);
            hold on
            plot(tvec,x(:,i),'b','linewidth',3)
            xlabel('Time (year)')
            ylabel(str{i})
            set(gca,'fontsize',20)
        end
    end

    x0=x(end,:);
    M=(T-T0)/dt;
    u=ones(M+1,1);
    x=zeros(M+1,5);
    lambda=zeros(M+1,5);
    tvec=[T0:dt:T];
    w=1e-3*exp(1e-2*tvec)';
    w=zeros(size(w));
    eta=1;
    test=-1;
    cur_it=0;
    epsilon= epsilon0*(1./(tvec-T0+1).^2)';

    while(test < 0)

        oldu = u;
        oldx = x;
        oldlambda = lambda;

        solx=ode45(@(t,x) state_eqn(t,x,paras,tvec,u),tvec,x0,opts);
        x = deval(solx,tvec)';

        sollambda= ode45(@(t,lambda) adjoint_Eqn(t,lambda,x,tvec,u,paras,epsilon),[T T0],[0 0 0 0 0]);
        lambda = deval(sollambda,tvec)';

        Abeta=x(:,1);
        lambda1=lambda(:,1);

        %u=u+eta*(Abeta.*lambda1-(w+2*epsilon*u));
        %temp=(Abeta.*lambda1+paras.alpha*Abeta)./epsilon/2;
        temp=(Abeta.*lambda1)./(epsilon.*(Abeta+paras.K_Abeta))/2;
        u1 = max(min(umax,temp),0);% max(0,temp);
        u = 0.5*(u1 + oldu);
        u=0.9*oldu+0.1*u1;
        test=min([delta*norm(u,1)-norm(oldu-u,1) delta*norm(x,1)-norm(oldx-x,1) delta*norm(lambda,1)-norm(oldlambda-lambda,1)])
        cur_it=cur_it+1
        if cur_it>300
            break
        end
    end
    hold on
    str={'A_\beta','\tau_p','\tau_o','N','C','u'};
    for i=1:5
        subplot(3,2,i);
        hold on
        if year
            plot(tvec,x(:,i),colors{(T0-10)/20+1},'linewidth',3)
            xlabel('Time (year)')
        else
            plot((tvec-T0)*365/7,x(:,i),colors{(T0-10)/20+1},'linewidth',3)
            xlabel('Time (week)')
        end
        xlim([0 100])
        ylabel(str{i})
        set(gca, 'YScale', 'log')
        set(gca,'fontsize',20)
    end
    subplot(3,2,6);
    hold on
    if year
        plot(tvec,u,colors{(T0-10)/20+1},'linewidth',3)

        xlabel('Time (year)')
    else
        plot((tvec-T0)*365/7,u,colors{(T0-10)/20+1},'linewidth',3)
        xlabel('Time (week)')
    end
    set(gca, 'YScale', 'log')

    xlim([0 100])
    ylabel('u')
    set(gca,'fontsize',20)
end
