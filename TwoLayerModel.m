function [T,uz] = TwoLayerModel (t,Lambda,p,S)
   
   % Standardize vector orientations 
        t = t(:);
        p = p(:).';

    % Extract trial parameters and construct material shorthand
        alpha_f = p(1); % Film thermal diffusivity [m^2/s]
        alpha_s = p(2); % Substrate thermal diffusivity [m^2/s]
        R = p(3); % Interfacial thermal resistance [m^2 K/W]
        L = S.L; % Film thickness [m]
        q = 2*pi/Lambda; % grating wavevector [1/m]

        % Conductivities must be recalculated whenever either trial diffusivity changes
            kf = S.Cf*alpha_f; % Film thermal conductivity [W/(m-K)]
            ks = S.Cs*alpha_s; % Substrate thermal conductivity [W/(m-K)]

        nu_f = S.nuf;
        nu_s = S.nus;
        mu_f = S.Ef/(2*(1 + nu_f)); % Film shear modulus [Pa]
        mu_s = S.Es/(2*(1 + nu_s)); % Substrate shear modulus [Pa]

        % Effective coupling
            thermo_f = S.alphathf*(1 + nu_f)/(1 - nu_f);
            thermo_s = S.alphaths*(1 + nu_s)/(1 - nu_s);

    % Build the symmetric angular-frequency grid
        nPositive = S.Nw/2;
        omegaPositive = logspace(log10(S.wmin), log10(S.wmax), nPositive);
        omega = [-fliplr(omegaPositive), omegaPositive];
        temperaturePositive = complex(zeros(size(omegaPositive)));
        displacementPositive = complex(zeros(size(omegaPositive)));

    % Assemble the elastic boundary-condition matrix
        exp2qL = exp(2*q*L);
        shearRatio = mu_s/mu_f;

        elasticMatrix = [1, 2*nu_f-2,     1,       2-2*nu_f,              0,           0; ...
                         1, 2*nu_f-1,     -1,      2*nu_f-1,              0,           0; ...
                         1, q*L-3+4*nu_f, -exp2qL, exp2qL*(4*nu_f-q*L-3), -1,          3-4*nu_s-q*L; ...
                         1, q*L,          exp2qL,  q*L*exp2qL,            -1,          -q*L; ...
                         1, q*L-2+2*nu_f, exp2qL,  exp2qL*(q*L+2-2*nu_f), -shearRatio, (2-2*nu_s-q*L)*shearRatio; ...
                        1, q*L-1+2*nu_f, -exp2qL, exp2qL*(2*nu_f-q*L-1), -shearRatio, (1-2*nu_s-q*L)*shearRatio];

       elasticSolver = decomposition(elasticMatrix,"lu");    % Solve the coupled thermal and quasi-static elastic problems
      
        for frequencyIndex = 1:nPositive
            currentOmega = omegaPositive(frequencyIndex);
            beta_f = sqrt(q^2 + 1i*currentOmega/alpha_f);
            beta_s = sqrt(q^2 + 1i*currentOmega/alpha_s);
            r_f = 1i*currentOmega/alpha_f;
            r_s = 1i*currentOmega/alpha_s;

            conductivityRatio = ks*beta_s/(kf*beta_f);
            interfaceFactor = 1 - conductivityRatio + ks*beta_s*R;
            thermalDenominator = 2*conductivityRatio + interfaceFactor*(1 - exp(-2*beta_f*L));

            A_f = S.Q0/(kf*beta_f) * (2*conductivityRatio + interfaceFactor)/thermalDenominator;
            B_f = S.Q0/(kf*beta_f) * interfaceFactor*exp(-2*beta_f*L)/thermalDenominator;
            A_f_interface = A_f*exp((q - beta_f)*L);
            B_f_interface = B_f*exp((q + beta_f)*L);
            A_s_interface = S.Q0/(kf*beta_f) *2*exp((q - beta_f)*L)/thermalDenominator;

            % Thermal forcing of the elastic system
                elasticForcing = [ thermo_f*beta_f*(A_f - B_f)/r_f; ...
                                   thermo_f*q*(A_f + B_f)/r_f; ...
                                   thermo_f*q*(A_f_interface + B_f_interface)/r_f - thermo_s*q*A_s_interface/r_s; ...
                                   thermo_f*beta_f*(A_f_interface - B_f_interface)/r_f - thermo_s*beta_s*A_s_interface/r_s; ...
                                   thermo_f*beta_f*(A_f_interface - B_f_interface)/r_f - shearRatio*thermo_s*beta_s*A_s_interface/r_s; ...
                                   thermo_f*q*(A_f_interface + B_f_interface)/r_f - shearRatio*thermo_s*q*A_s_interface/r_s ];

            % Solve for the six elastic coefficients 
                elasticCoefficients = elasticSolver\elasticForcing;

            % Surface values at z = 0
                temperaturePositive(frequencyIndex) = A_f + B_f;
                displacementPositive(frequencyIndex) = thermo_f*beta_f*(B_f - A_f)/r_f + elasticCoefficients(1) + elasticCoefficients(3);
        end

    % Reconstruct the full spectrum and transform back to time
        temperatureSpectrum = [conj(fliplr(temperaturePositive)), temperaturePositive];
        displacementSpectrum = [conj(fliplr(displacementPositive)), displacementPositive];

      T = real(inverseFourierPiecewiseLinear(omega,temperatureSpectrum,t))/(2*pi);
      uz = real(inverseFourierPiecewiseLinear(omega,displacementSpectrum,t))/(2*pi);
end

function x = inverseFourierPiecewiseLinear(omega, spectrum, t)
   % Performs an inverse Fourier transform using piecewise-linear interpolation.    
    omega = omega(:).';
    spectrum = spectrum(:).';
    t = t(:);
    
    deltaOmega = diff(omega);
    spectralSlope = diff(spectrum)./deltaOmega;
    
    x = complex(zeros(size(t)));
    
    for timeIndex = 1:numel(t)
        currentTime = t(timeIndex);
        scaledInterval = currentTime*deltaOmega;
        useSeries = abs(scaledInterval) < 1e-4;
        integralConstant = complex(zeros(size(deltaOmega)));
        integralLinear = complex(zeros(size(deltaOmega)));
        interval = deltaOmega(~useSeries);
        phase = exp(1i*currentTime*interval);
        integralConstant(~useSeries) = (phase - 1)/(1i*currentTime);
        integralLinear(~useSeries) = phase.*(interval/(1i*currentTime) + 1/currentTime^2) - 1/currentTime^2;

        interval = deltaOmega(useSeries);
        integralConstant(useSeries) = interval + 1i*currentTime*interval.^2/2 - currentTime^2*interval.^3/6;
        integralLinear(useSeries) = interval.^2/2 + 1i*currentTime*interval.^3/3 - currentTime^2*interval.^4/8;

        x(timeIndex) = sum(exp(1i*currentTime*omega(1:end-1)).* (spectrum(1:end-1).*integralConstant  + spectralSlope.*integralLinear));
    end
end
