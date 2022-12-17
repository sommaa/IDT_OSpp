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

function Err = min2(K,T_range,P_range,phi,tau_c_l)
A1= K(1); n1= K(2); B1= K(3); m1= K(4); beta1= K(5); 

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
            
            tau_l_eq(cycle_counter) = A1*P^n1*exp(B1/T)*T^m1*phi_num^beta1;
            weights(cycle_counter) = 1;
            if tau_c_l(cycle_counter) > tau_l_eq(cycle_counter)
                weights(cycle_counter) = 1.8;
            end
        end
    end
end

res = zeros(1,cycle_counter);

for jj = 1:cycle_counter
    if ~isnan(tau_c_l(jj))
        res(jj) = (abs(tau_l_eq(jj)-tau_c_l(jj)))/tau_c_l(jj)*weights(jj);
    else
        res(jj) = 0;
    end
end

Err = sum(res);
