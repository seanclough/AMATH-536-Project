function theo_matrix = plot_theo()
    epsilon = .85;
    lambda = 0;
    alpha = 2*epsilon+2*lambda-1;
    gamma = 1.26;
    N_0 = 400;
    r_0 = 1/365;
    r = zeros(1,6);
    for i=1:6
        r(i) = 1/365*(gamma)^i;
    end
    theo_matrix = zeros(201,6);
    for i=1:6
        %\prod_{l=1,l\neq j}^i(\gamme^i-gamma^l)
        result1 = ones(1,i);
        % u is used in place of j
        for u=1:i
            % o is used in place of l
            for o=1:i
                if o ~= u
                    result1(u) = result1(u)*(gamma^u-gamma^o);
                end
            end
        end
        for j=1:201
            result2 = 0;
            % m is used in place of j
            for m=1:i
                result2 = result2+(-1)^i/(gamma^m*result1(m))*(exp(-alpha*gamma^m/N_0*20*(j-1))-exp(-alpha*gamma^i/N_0*20*(j-1)));
            end
            theo_matrix(j,i) = N_0/gamma^i/alpha^i*(2*epsilon)^(i-1)*(1-exp(-alpha*gamma^i/N_0*20*(j-1)))+N_0*sqrt(gamma^(i*(i-1)))*(2*epsilon)^(i-1)/alpha^i*result2;
        end
    end
    for i=1:6
        plot(0:20:4000,theo_matrix(:,i))
    end
end