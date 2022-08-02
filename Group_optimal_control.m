addpath('../Old_code/')
load paras_values
diseae_type=1;%AD
diseae_type=3;%MCI

%best fit
disease_id=find(paras_list(end,:)==diseae_type);
[~,id]=min(paras_list(end-2,disease_id));
best_ind=disease_id(id);


C_change=[];

for id=disease_id
    if id==best_ind
        draw_pic=1;
    else
        draw_pic=0;
    end
    x=paras_list(:,id);
    paras.initial=x([3 6 8 12 16])';

    paras.lambda_beta=x(1);
    paras.K_Abeta=x(2);

    paras.lambda_ptau=x(4);
    paras.K_ptau=x(5);

    paras.lambdatauO=x(7);

    paras.lambda_NdTO=x(9);
    paras.lambda_Ndptau=x(10);
    paras.K_Nd=x(11);

    paras.lambda_CN=x(13);
    paras.lambda_Ct=x(14);
    paras.K_C=x(15);

    tmp=[];
    figure(1)
    epsilon0=5;
    treat_age=60;
    umax=0.12;%low dosage
    colors={'b'};
    Optimal_control_week
    tmp=[tmp ((x(end,end)-x_no_treatment(end,end))./x_no_treatment(end,end))*100 J];
    umax=0.21;%high dosage
    colors={'r'};
    Optimal_control_week
    tmp=[tmp ((x(end,end)-x_no_treatment(end,end))./x_no_treatment(end,end))*100 J];
    umax=1.04*2;%high dosage
    epsilon0=0.1;
    colors={'g'};
    Optimal_control_week
    tmp=[tmp ((x(end,end)-x_no_treatment(end,end))./x_no_treatment(end,end))*100 J];

    figure(2)
    epsilon0=5;
    treat_age=70;
    umax=0.12;%low dosage
    colors={'b'};
    Optimal_control_week
    tmp=[tmp ((x(end,end)-x_no_treatment(end,end))./x_no_treatment(end,end))*100 J];
    umax=0.21;%high dosage
    colors={'r'};
    Optimal_control_week
    tmp=[tmp ((x(end,end)-x_no_treatment(end,end))./x_no_treatment(end,end))*100 J];
    umax=1.04*2;%high dosage
    epsilon0=0.1;
    colors={'g'};
    Optimal_control_week
    tmp=[tmp ((x(end,end)-x_no_treatment(end,end))./x_no_treatment(end,end))*100 J];

    figure(3)
    epsilon0=5;
    treat_age=60;
    umax=0.12;%low dosage
    colors={'b'};
    Optimal_control_year
    tmp=[tmp ((x(end,end)-x_no_treatment(end,end))./x_no_treatment(end,end))*100 J];
    umax=0.21;%high dosage
    colors={'r'};
    Optimal_control_year
    tmp=[tmp ((x(end,end)-x_no_treatment(end,end))./x_no_treatment(end,end))*100 J];
    umax=1.04*2;%high dosage
    epsilon0=0.1;
    colors={'g'};
    Optimal_control_year
    tmp=[tmp ((x(end,end)-x_no_treatment(end,end))./x_no_treatment(end,end))*100 J];

    figure(4)
    epsilon0=5;
    treat_age=70;
    umax=0.12;%low dosage
    colors={'b'};
    Optimal_control_year
    tmp=[tmp ((x(end,end)-x_no_treatment(end,end))./x_no_treatment(end,end))*100 J];
    umax=0.21;%high dosage
    colors={'r'};
    Optimal_control_year
    tmp=[tmp ((x(end,end)-x_no_treatment(end,end))./x_no_treatment(end,end))*100 J];
    umax=1.04*2;%high dosage
    epsilon0=0.1;
    colors={'g'};
    Optimal_control_year
    tmp=[tmp ((x(end,end)-x_no_treatment(end,end))./x_no_treatment(end,end))*100 J];

    C_change=[C_change;tmp];
    id
end

for i=1:size(C_change,1)
    %tmp=['&' num2str(patient_id(paras_list(end-1,disease_id(i))))];
    if i==1
        tmp=['\multirow{' num2str(size(C_change,1)+1) '}{*}{78 weeks }'];
    else
        tmp=[];
    end

    tmp=[tmp '&' num2str(Age((paras_list(end-1,disease_id(i)))))];
    tmp=[tmp '&' char(Gender((paras_list(end-1,disease_id(i)))))];
    for j=1:2:12
        tmp=[tmp '&$' num2str(abs(C_change(i,j)*1e2),2) '\times10^{-2}$'];
    end
    tmp=[tmp '\\\cline{2-9}'];
    disp(tmp)
end

tmp=['& \multicolumn{2}{|c|}{Avg}'];
for j=1:2:12
    tmp=[tmp '&$' num2str(abs(mean(C_change(:,j))*1e2),2) '\times10^{-2}$'];
end
tmp=[tmp '\\\hline'];
disp(tmp)
for i=1:size(C_change,1)
    if i==1
        tmp=['\multirow{' num2str(size(C_change,1)+1) '}{*}{10 years }'];
    else
        tmp=[];
    end

    tmp=[tmp '&' num2str(Age((paras_list(end-1,disease_id(i)))))];
    tmp=[tmp '&' char(Gender((paras_list(end-1,disease_id(i)))))];
    for j=13:2:size(C_change,2)
        tmp=[tmp '&$' num2str(abs(C_change(i,j)),2) '$'];
    end
    tmp=[tmp '\\\cline{2-9}'];
    disp(tmp)
end

tmp=['& \multicolumn{2}{|c|}{Avg}'];
for j=13:2:size(C_change,2)
    tmp=[tmp '&$' num2str(abs(mean(C_change(:,j))),2) '$'];
end
tmp=[tmp '\\\hline'];
disp(tmp)