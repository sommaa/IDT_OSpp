%                                          ___  ________  _________                                                   %
%                                         |\  \|\   ___ \|\___   ___\                                                 %
%                                         \ \  \ \  \_|\ \|___ \  \_|                                                 %
%                                          \ \  \ \  \ \\ \   \ \  \                                                  %
%                                           \ \  \ \  \_\\ \   \ \  \                                                 %
%                                            \ \__\ \_______\   \ \__\                                                %
%                                             \|__|\|_______|    \|__|                                                %
%                                                                                                                     %
%                       Author: Andrea Somma;                                                                         %
%                       Ignition delay time; Combustion course (Prof: Marco Mehl)(2021-2022);                         %
%                       Politecnico of Milan.                                                                         %

fclose('all'); clc; clear all; close all;

%% inputs
MoleOx = "O2 0.21 N2 0.79";
MoleFuel = "IC8H16 1";
species_youwant = "OH";
T_range = ["600 K" "610 K" "625 K" "635 K" "645 K" "650 K" "665 K" "675 K" "685 K" "700 K" "715 K" "725 K" "750 K" "765 K" ...
    "775 K" "785 K" "800 K" "810 K" "825 K" "835 K" "845 K" "850 K" "865 K" "875 K" "885 K" "900 K" "925 K" "950 K" ...
    "1000 K" "1050 K" "1075 K" "1100 K" "1125 K" "1150 K" "1175 K" "1200 K" "1225 K" "1250 K"  "1275 K" "1300 K"  ...
    "1325 K" "1350 K" "1375 K" "1450 K" "1500 K" "1525 K" "1550 K" "1600 K"];%;
T_range_num = [600 610 625 635 645 650 665 675 685 700 715 725 750 765 775 785 800 810 825 835 845 850 865 875 885 ... 
    900 925 950 1000 1050 1075 1100 1125 1150 1175 1200 1225 1250 1275 1300 1325 1350 1375 1450 1500 1525 1550 1600];
P_range = ["20 atm" "25 atm" "30 atm" "32.5 atm" "35 atm" "37.5 atm" "40 atm" "45 atm" "50 atm" "65 atm" "70 atm"];
phi = ["0.75" "0.8" "0.9" "1" "1.1"];
mech_file = 'C7.txt';
trans_file = 'C7_T.txt';
thermo_file = 'C7_TH.txt'; lines_to_skip = 2;

%% test full mech

mkdir("complete_model_test")

RunFileID = fopen('./complete_model_test/Run.bat', 'wt');


for i = 1:length(T_range)
    for j = 1:length(P_range)
        for k = 1:length(phi)
            counter = regexprep(strcat(T_range(i)," ",P_range(j)," ",phi(k)), ' ', '_');
            OUT_file = "1";
            testname = strcat(counter,".dic");
            
            % Open the file for reading in text mode.
            fileID = fopen('./complete_model_test/input2.dic', 'rt');
            % Open the file for reading in text mode.
            outputFileID = fopen(strcat('./complete_model_test/',testname), 'wt');
            % Read the first line of the file.
            textLine = fgetl(fileID);
            lineCounter = 1;
            
            while ischar(textLine)
                if contains(textLine, '@Kinetics ')
                    textLineOut = sprintf('%s', textLine, strcat('../kin/',mech_file,';'));
                elseif contains(textLine, '@Thermodynamics ')
                    textLineOut = sprintf('%s', textLine, strcat('../kin/',thermo_file,';'));
                elseif contains(textLine, '@Output ')
                    textLineOut = sprintf('%s', textLine, strcat('../',OUT_file,';'));
                elseif contains(textLine, '@EndTime')
                    textLineOut = sprintf('%s', textLine, "2 s",';');
                elseif contains(textLine, '@Temperature') && lineCounter < 70
                    textLineOut = sprintf('%s', textLine,T_range(i),';');
                elseif contains(textLine, '@Pressure') && lineCounter < 70
                    textLineOut = sprintf('%s', textLine, P_range(j),';');
                elseif contains(textLine, '@FuelMoles ')
                    textLineOut = sprintf('%s', textLine, strcat(" ",MoleFuel),';');
                elseif contains(textLine, '@OxidizerMoles ')
                    textLineOut = sprintf('%s', textLine, strcat(" ",MoleOx),';');
                elseif contains(textLine, '@Species ')
                    textLineOut = sprintf('%s', textLine, species_youwant,';');
                elseif contains(textLine, '@OutputFolder')
                    textLineOut = sprintf('%s', textLine, counter,';');
                elseif contains(textLine, '@EquivalenceRatio')
                    textLineOut = sprintf('%s', textLine, strcat(" ",phi(k)),';');
                else
                    textLineOut = sprintf('%s %f', textLine);
                end
                % Output line of text to output file
                fprintf(outputFileID, '%s\n', textLineOut);
                
                % Read the next line.
                textLine = fgetl(fileID);
                lineCounter = lineCounter + 1;
            end
            % All done reading all lines, so close the files.
            fclose(fileID);
            fclose(outputFileID);
            
            textLine = strcat('"%OPENSMOKEPP_EXE_FOLDER%\OpenSMOKEpp_BatchReactor.exe" --input'," ",testname);
            textLineOut = sprintf('%s', textLine);
            fprintf(RunFileID, '%s\n', textLineOut);
            
            
        end
    end
end

fclose(RunFileID);
! cd complete_model_test\ && Run.bat & exit/b

clc

%% range high
T_range_h = ["1000 K" "1050 K" "1075 K" "1100 K" "1125 K" "1150 K" "1175 K" "1200 K" "1225 K" "1250 K"  "1275 K" "1300 K" ...
    "1325 K" "1350 K" "1375 K" "1450 K" "1500 K" "1525 K" "1550 K" "1600 K"];
T_range_h_num = [1000 1050 1075 1100 1125 1150 1175 1200 1225 1250 1275 1300 1325 1350 1375 1450 1500 1525 1550 1600];

%%
cycle_counter = 0;

for i = 1:length(T_range_h)
    for j = 1:length(P_range)
        for k = 1:length(phi)
            cycle_counter = cycle_counter + 1;
            cartel = regexprep(strcat(T_range_h(i)," ",P_range(j)," ",phi(k)), ' ', '_');
            complete_dir_IDT = strcat('./complete_model_test/',cartel,'/IDT.out');
            
            completeID_IDT = fopen(complete_dir_IDT,'rt');
            
            textLine_complete_IDT = fgetl(completeID_IDT);
            lineCounter = 1;
            
            while ischar(textLine_complete_IDT)
                K = split(textLine_complete_IDT);
                
                if strcmp(K(1),'OH(max)')
                    tau_compl = K(2);
                end
                textLine_complete_IDT = fgetl(completeID_IDT);
            end
            
            tau_c_1_h(cycle_counter) = str2num(cell2mat(tau_compl));
            tau_c_3d_1_h(i,j,k) =  str2num(cell2mat(tau_compl));
            fclose('all');
        end
    end
end

%% ga range high
lb = [1e-8,-1.4,13000,-1.1];
ub = [1e-4,-0.7,28000,-0.1];

options= optimoptions('ga','ConstraintTolerance',1e-11,'FunctionTolerance', 1e-11,...
    'MaxGeneration',100,'UseParallel',false,'PopulationSize',100,...
    'PlotFcn', {@gaplotbestf,@gaplotstopping,@gaplotrange,@gaplotselection,@gaplotscorediversity});

coeff1h = ga(@(K)min1(K,T_range_h,P_range,phi,tau_c_1_h),4,[],[],[],[],lb,ub,[],options);

%% plot range high

for i = 1:length(T_range_h)
    for j = 1:length(P_range)
        for k = 1:length(phi)
            
            tmp1 = split(T_range_h(i));
            T = str2num(tmp1(1));
            tmp1 = split(P_range(j));
            P = str2num(tmp1(1))*101325;
            phi_num = str2num(phi(k));
            tau_eq_1(i,j,k) = coeff1h(1)*P^coeff1h(2)*exp(coeff1h(3)/T)*phi_num^coeff1h(4);
            
        end
    end
end

cycle_counter = 0;
clf

for i = 1:length(P_range)
    for k = 1:length(phi)
        cycle_counter = cycle_counter+1;
        figure(cycle_counter)
        semilogy(1000./T_range_h_num,tau_eq_1(:,i,k))
        hold on
        semilogy(1000./T_range_h_num,tau_c_3d_1_h(:,i,k))
        title(strcat(P_range(i),phi(k)))
        legend("model", "OSpp")
    end
end

%% range low

cycle_counter = 0;

T_range_l = ["685 K" "700 K" "715 K" "725 K" "750 K" "765 K" "775 K"];
T_range_l_num = [685 700 715 725 750 765 775];
P_range_l = ["14 atm" "16 atm" "18 atm" "20 atm" "22 atm" "24 atm" "26 atm" "28 atm"];

for i = 1:length(T_range_l)
    for j = 1:length(P_range_l)
        for k = 1:length(phi)
            cycle_counter = cycle_counter + 1;
            cartel = regexprep(strcat(T_range_l(i)," ",P_range_l(j)," ",phi(k)), ' ', '_');
            complete_dir_T = strcat('./complete_model_test/',cartel,'/Output.out');
            
            Out = readmatrix(complete_dir_T,'FileType','text');
            
            tau_c_1_l(cycle_counter) = nan;
            if Out(end,5) <2400
                tau_c_1_l(cycle_counter) = nan;
            else
                for kkk = 1:length(Out(:,5))-1
                    slope_1(kkk) = (Out(kkk+1,5)-Out(kkk,5))/(Out(kkk+1,1)-Out(kkk,1));
                    if slope_1(kkk) == inf
                        slope_1(kkk) = slope_1(kkk-1)+0.01;
                    end
                end
                
                t_c_1_l(cycle_counter) = nan;
                written = 0;
                
                for jjj = 4:length(slope_1)
                    if slope_1(jjj)<slope_1(jjj-1) && slope_1(jjj-1)<= slope_1(jjj-2) && slope_1(jjj-2)<= slope_1(jjj-3) ...
                            && slope_1(jjj-1) ~= 0 && slope_1(jjj-2) ~= 0 && slope_1(jjj-3) ~= 0 && slope_1(jjj) ~= 0 && ...
                            Out(jjj,5) < 1000 && written == 0 && Out(jjj,5)>T_range_l_num(i)+1 
                        t_1st = Out(jjj,1);
                        tau_c_1_l(cycle_counter) = t_1st;
                        written = 1;
                    end
                end
            end
            clear slope_1
            tau_c_3d_l(i,j,k) =  tau_c_1_l(cycle_counter);
            fclose('all');
        end
    end
end

%% ga range low

lb = [1e-9,-0.6,12000,-2,-0.8];
ub = [8e-6,-0.4,17000,5,-0.3];

options= optimoptions('ga','ConstraintTolerance',1e-10,'FunctionTolerance', 1e-10,...
    'MaxGeneration',100,'UseParallel',false,'PopulationSize',100,...
    'PlotFcn', {@gaplotbestf,@gaplotstopping,@gaplotrange,@gaplotselection,@gaplotscorediversity});

coeffl = ga(@(K)min2(K,T_range_l,P_range_l,phi,tau_c_1_l),5,[],[],[],[],lb,ub,[],options);

%% plot range low

for i = 1:length(T_range_l)
    for j = 1:length(P_range_l)
        for k = 1:length(phi)
            
            tmp1 = split(T_range_l(i));
            T = str2num(tmp1(1));
            tmp1 = split(P_range_l(j));
            P = str2num(tmp1(1))*101325;
            phi_num = str2num(phi(k));
            tau_eq_l(i,j,k) = coeffl(1)*P^coeffl(2)*exp(coeffl(3)/T)*T^coeffl(4)*phi_num^coeffl(5);
            
        end
    end
end

cycle_counter = 0;

for i = 1:length(P_range)
    for k = 1:length(phi)
        cycle_counter = cycle_counter+1;
        figure(cycle_counter)
        semilogy(1000./T_range_l_num,tau_eq_l(:,i,k))
        hold on
        semilogy(1000./T_range_l_num,tau_c_3d_l(:,i,k))
        title(strcat(P_range(i),phi(k)))
        legend("model", "OSpp")
    end
end

%% NTC
T_range_NTC = ["600 K" "610 K" "625 K" "635 K" "645 K" "650 K" "665 K" "675 K" "685 K" "700 K" "715 K" "725 K" "750 K" ...
    "765 K" "775 K" "785 K" "800 K" "810 K" "825 K" "835 K" "845 K" "850 K" "865 K" "875 K" "885 K" "900 K" "925 K" ...
    "950 K" "1000 K" "1050 K" "1075 K" "1100 K" "1125 K" "1150 K" "1175 K" "1200 K" "1225 K" "1250 K"  "1275 K" "1300 K"...
    "1325 K" "1350 K" "1375 K" "1450 K" "1500 K" "1525 K" "1550 K" "1600 K"];
T_range_NTC_num = [600 610 625 635 645 650 665 675 685 700 715 725 750 765 775 785 800 810 825 835 845 850 865 875 885 ...
    900 925 950 1000 1050 1075 1100 1125 1150 1175 1200 1225 1250 1275 1300 1325 1350 1375 1450 1500 1525 1550 1600];
cycle_counter = 0;

for i = 1:length(T_range_NTC)
    for j = 1:length(P_range)
        for k = 1:length(phi)
            cycle_counter = cycle_counter + 1;
            cartel = regexprep(strcat(T_range_NTC(i)," ",P_range(j)," ",phi(k)), ' ', '_');
            complete_dir_IDT = strcat('./complete_model_test/',cartel,'/IDT.out');
            
            completeID_IDT = fopen(complete_dir_IDT,'rt');
            
            textLine_complete_IDT = fgetl(completeID_IDT);
            lineCounter = 1;
            
            while ischar(textLine_complete_IDT)
                K = split(textLine_complete_IDT);
                
                if strcmp(K(1),'OH(max)')
                    tau_compl = K(2);
                end
                textLine_complete_IDT = fgetl(completeID_IDT);
            end
            
            tau_c_1(cycle_counter) = str2num(cell2mat(tau_compl));
            tau_c_3d_1(i,j,k) =  str2num(cell2mat(tau_compl));
            fclose('all');
        end
    end
end

%%
lb = [900,-1.38,0.0001,300,0.001,-6];
ub = [960,-1.33,0.8,18800,0.008,-4];

options= optimoptions('ga','ConstraintTolerance',1e-10,'FunctionTolerance', 1e-10,...
    'MaxGeneration',120,'UseParallel',false,'PopulationSize',120,...
    'PlotFcn', {@gaplotbestf,@gaplotstopping,@gaplotrange,@gaplotselection,@gaplotscorediversity});

coeff_whole = ga(@(K)min3(K,T_range_NTC,P_range,phi,tau_c_1,coeffl,coeff1h),6,[],[],[],[],lb,ub,[],options);

%% plot whole
clf

for i = 1:length(T_range_NTC)
    for j = 1:length(P_range)
        for k = 1:length(phi)
            Teq= coeff_whole(1); omega= coeff_whole(2); key= coeff_whole(3); C0= coeff_whole(4); mu= coeff_whole(5); ...
                sigma=coeff_whole(6);
            A1= coeffl(1); n1= coeffl(2); B1= coeffl(3); m1= coeffl(4); beta1= coeffl(5);
            Ah= coeff1h(1); nh= coeff1h(2); Bh= coeff1h(3); betah= coeff1h(4);
            
            
            tmp1 = split(P_range(j));
            P = str2num(tmp1(1))*101325;
            tmp1 = split(T_range_NTC(i));
            T = str2num(tmp1(1));
            phi_num = str2num(phi(k));
            
            DTcf = omega*(T-Teq*P^key*phi_num^mu*(100/(99+phi_num))^sigma);
            Tcf = T+0.5*(DTcf +sqrt((DTcf)^2+C0));
            DPcf = DTcf/Tcf*P;
            Pcf = P+DPcf;
            tau1 = A1*P^n1*exp(B1/T)*T^m1*phi_num^beta1;
            tauh_cf = Ah*Pcf^nh*exp(Bh/Tcf)*phi_num^betah;
            tauh = Ah*P^nh*exp(Bh/T)*phi_num^betah;
            tau_eq(i,j,k) = tau1+tauh_cf*(1-tau1/tauh);
            tau_eq_l(i,j,k) = coeffl(1)*P^coeffl(2)*exp(coeffl(3)/T)*T^coeffl(4)*phi_num^coeffl(5);
            tau_eq_1(i,j,k) = coeff1h(1)*P^coeff1h(2)*exp(coeff1h(3)/T)*phi_num^coeff1h(4);
        end
    end
end

cycle_counter = 0;

for i = 1:length(P_range)
    for k = 1:length(phi)
        cycle_counter = cycle_counter + 1;
        
        figure(cycle_counter)
        semilogy(10000./T_range_NTC_num,tau_c_3d_1(:,i,k),"LineWidth",2.2)
        hold on
        semilogy(10000./T_range_NTC_num,tau_eq(:,i,k),"LineWidth",1.4)
        semilogy(10000./T_range_NTC_num,tau_eq_l(:,i,k),"LineWidth",1.5,'LineStyle','-.')
        semilogy(10000./T_range_NTC_num,tau_eq_1(:,i,k),"LineWidth",1.5,'LineStyle','-.')
        title(strcat(P_range(i),phi(k)))
        legend("OSpp","model tot", "model low", "model high",'Location','southeast')
        xlabel('10000/T [K]')
        ylabel('Time [s]')
        grid on
    end
end

%% errore percentuale
cycle_counter = 0;
close all

for i = 1:length(T_range)
    for j = 1:length(P_range)
        for k = 1:length(phi)
            cycle_counter = cycle_counter + 1;
            tau_line(cycle_counter) = tau_eq(i,j,k);
            if tau_c_1(cycle_counter)/tau_line(cycle_counter) < 0.3 || tau_c_1(cycle_counter)/tau_line(cycle_counter) > 2
                tau_c_1(cycle_counter) = 0;
                tau_c_1(cycle_counter) = 0;
            end
        end
    end
end

loglog(tau_c_1,tau_line, '.', 'markersize', 5)
hold on
x = linspace(1e-6,0.4);
loglog(x,x,"LineWidth",1.2)
loglog(x,x*1.3,"LineWidth",1.2)
loglog(x,x*0.7,"LineWidth",1.2)
xlabel("correlated IDT",'FontSize',17)
ylabel("OSpp IDT",'FontSize',17)
legend({"","0% error","+30% error","-30% error"},'FontSize',12)
grid on