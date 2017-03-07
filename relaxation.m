function [ A, b ] = relaxation( t, T1, T2 )
    A = diag([exp(-t/T2) exp(-t/T2) exp(-t/T1)]);
    b = [0 0 1-exp(-t/T1)]'; % Assuming Mo = 1
end

