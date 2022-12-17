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

function Err = min3(K,T_range,P_range,phi,tau_c,coeffl,coeff1h)
Teq= K(1); omega= K(2); key= K(3); C0= K(4); mu= K(5); sigma= K(6);
A1= coeffl(1); n1= coeffl(2); B1= coeffl(3); m1= coeffl(4);
Ah= coeff1h(1); nh= coeff1h(2); Bh= coeff1h(3);

cycle_counter = 0;

for i = 1:length(T_range)
    for j = 1:length(P_range)
        for k = 1:length(phi)
            cycle_counter = cycle_counter + 1;
            
            tmp1 = split(P_range(j));
            P = str2num(tmp1(1))*101325;
            tmp1 = split(T_range(i));
            T = str2num(tmp1(1));
            phi_num = str2num(phi(k));
            
            DTcf = omega*(T-Teq*P^key*phi_num^mu*(100/(99+phi_num))^sigma);
            Tcf = T + 0.5*(DTcf + sqrt(DTcf^2+C0));
            DPcf = DTcf/Tcf*P;
            Pcf = P+DPcf;
            tau1 = A1*P^n1*exp(B1/T)*T^m1;
            tauh_cf = Ah*Pcf^nh*exp(Bh/(Tcf));
            tauh = Ah*P^nh*exp(Bh/T);
            tau_eq(cycle_counter) = tau1+tauh_cf*(1-tau1/tauh);
            weights(cycle_counter) = 1;
            
        end
    end
end

res = zeros(1,cycle_counter);

for jj = 1:cycle_counter
    if tau_c(jj)<0.8 && tau_c(jj) ~= 0
        res(jj) = (abs(tau_eq(jj)-tau_c(jj)))/tau_c(jj)*weights(jj);
    elseif tau_eq(jj) < 0
        res(jj) = res(jj)*10000;
    else
        res(jj) = 0;
    end
end

Err = sum(res);
