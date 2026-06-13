function out = reward_v2(P,Nap,Nue,H,connection)
out = 0;
for i=1:Nap
    for j =1:Nue
        if connection(i, j) == 1
            received_power = P(i)*H(i,j);
            interference_power =sum(H(i,connection(i, :) == 1)) - received_power;
            noise_power = randn() + 1i * randn();  % Complex Gaussian noise
            SINR=abs(received_power / (interference_power + noise_power^2));
            out = out+log(1+SINR);
        end
    end
end
end

