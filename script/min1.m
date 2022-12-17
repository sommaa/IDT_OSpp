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

function Err = min1(K,T_range,P_range,phi,tau_c)
Ah= K(1); nh= K(2); Bh= K(3); betah = K(4);

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
            
            tau_h(cycle_counter) = Ah*P^nh*exp(Bh/T)*phi_num^betah;
            weights(cycle_counter) = 1;
            
            if tau_h(cycle_counter) > 0.2*1e3 || tau_h(cycle_counter) < 1e-70
                tau_h(cycle_counter) = 100000;
            end
            if tau_c(cycle_counter) > tau_h(cycle_counter)
                weights(cycle_counter) = 1.2;
            end

        end
    end
end

res = zeros(1,cycle_counter);
for jj = 1:cycle_counter
    res(jj) = abs(tau_h(jj)-tau_c(jj))/abs(tau_c(jj))*weights(cycle_counter);
    %if priority_index(cycle_counter) == 1
    %    res(jj)= res(jj)*10;
    %end
end

Err = sum(res);
