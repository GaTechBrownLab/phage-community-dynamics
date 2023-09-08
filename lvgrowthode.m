function dNdt = lvgrowthode(~,N,p,numsp,scale)
% Sarah Sundius, February 9, 2023
%
% Function to define the Lotka-Volterra growth model
% Inputs:   N = species densities, oriented (r,c) = (time,species)
%           p = parameter set (to optimize/fit) 
%           numsp = number of species
%           scale = boolean to use log scale or not (1 = log scaled)
%
% Outputs:  dNdt = Lotka-Volterra model

if numsp == 1   
    % Define the model, p = [r k]
    if scale
        dNdt = p(1,1) - (p(1,1)/p(1,2))*exp(N);
        %dNdt = p(1,1) + p(1,2)*exp(N);
    else
        dNdt = p(1,1)*N*(1 - N/p(1,2));
    end

elseif numsp == 2 
    % Define the model, p = [r1 b11 b12; r2 b21 b22]
    if scale
        dNdt(1,1) = p(1,:)*[1; N(1); N(2)];
        dNdt(2,1) = p(2,:)*[1; N(1); N(2)];
    else
        dNdt(1,1) = p(1,:)*[N(1); N(1).^2; N(1).*N(2)];
        dNdt(2,1) = p(2,:)*[N(2); N(1).*N(2); N(2).^2];
    end 

elseif numsp == 3
    % Define the model, p = [r1 b11 b12 b13; r2 b21 b22 b23; r3 b31 b32 b33]
    if scale
        dNdt(1,1) = p(1,:)*[1; N(1); N(2); N(3)];
        dNdt(2,1) = p(2,:)*[1; N(1); N(2); N(3)];
        dNdt(3,1) = p(3,:)*[1; N(1); N(2); N(3)];
    else
        dNdt(1,1) = p(1,:)*[N(1); N(1).^2; N(1).*N(2); N(1).*N(2)];
        dNdt(2,1) = p(2,:)*[N(2); N(1).*N(2); N(2).^2; N(2).*N(3)];
        dNdt(3,1) = p(3,:)*[N(3); N(1).*N(3); N(2).*N(3); N(3).^2];
    end 

elseif numsp == 4
    % Define the model, p = [r1 b11 b12 b13 b14; r2 b21 b22 b23 b24; r3 b31 b32 b33 b34; r4 b41 b42 b43 b44]
    if scale
        dNdt(1,1) = p(1,:)*[1; N(1); N(2); N(3); N(4)];
        dNdt(2,1) = p(2,:)*[1; N(1); N(2); N(3); N(4)];
        dNdt(3,1) = p(3,:)*[1; N(1); N(2); N(3); N(4)];
        dNdt(4,1) = p(4,:)*[1; N(1); N(2); N(3); N(4)];
    else
        dNdt(1,1) = p(1,:)*[N(1); N(1).^2; N(1).*N(2); N(1).*N(3); N(1).*N(4)];
        dNdt(2,1) = p(2,:)*[N(2); N(1).*N(2); N(2).^2; N(2).*N(3); N(2).*N(4)];
        dNdt(3,1) = p(3,:)*[N(3); N(1).*N(3); N(2).*N(3); N(3).^2; N(3).*N(4)];
        dNdt(4,1) = p(4,:)*[N(4); N(1).*N(4); N(2).*N(4); N(3).*N(4); N(4).*2];
    end 
    
end

end