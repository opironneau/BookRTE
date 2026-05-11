function solveIQ()
    % solve Milne problem and multi-group radiative transfer
    % Translated from C++ to MATLAB

    clc;
    
    %% Constants
    pi_val = atan(1)*4; 
    stefan = (pi_val^2)^2/15;
    
    Ce = 2.0; 
    Cs = 2.0e-6;
    rho0 = 1.0; 
    drho = -0.7; % density gradient
    TOA = 1.2; 
    Z = TOA*rho0*(1+drho*TOA/2); % optical thickness
    B0 = 1.4744e-8; 
    T00_val = 4798.;
    
    Te = (273+18)/T00_val; 
    Ts = 5778/T00_val;
    q0 = -0.3;
    mus = 0.5;
    beta_val = 0.5;
    
    verbose = false;
    Nz = 60; % nb points in tau
    dz = Z/(Nz-1);
    kmax = 9;  % nb fixed point iterations
    jmaxmax = 600; % max of max nb of points for integration in nu range
    newton_iters = 50; % to compute the temperature
    epsdycho = 1e-4; 
    epsnewton = 1.e-10;  % precision for dychotomy & Newton
    kappamin = 0.001;  % if kappa read is too small max it with kappamin
    nt = 5; 
    z1 = Z*0.5; z2 = Z*0.7; z3 = 0.8*Z;  % altitudes for cloud and Rayleigh scattering
    nu1 = 0.5; nu2 = 1.0; % range of Rayleigh scattering
    cloud = 1.0;

    %% Preallocated Arrays
    kappanu = zeros(1, jmaxmax);
    nu = zeros(1, jmaxmax);
    J0 = zeros(jmaxmax, Nz); J0old = zeros(jmaxmax, Nz);
    J2 = zeros(jmaxmax, Nz); J2old = zeros(jmaxmax, Nz);
    K0 = zeros(jmaxmax, Nz); K0old = zeros(jmaxmax, Nz);
    K2 = zeros(jmaxmax, Nz); K2old = zeros(jmaxmax, Nz);
    T = zeros(1., Nz);
    kappaz = zeros(jmaxmax, Nz);
    kappasum = zeros(jmaxmax, Nz);
    
    % Paths - Adjust as needed
    basedir = '/Users/pironneau/Dropbox/aranger/TeX2026/BookVRTE/prog4/';
    mykappafile = fullfile(basedir, '_kappa.txt'); % 1 - kappa
    myresulttemperature = fullfile(basedir, 'temperatureyyx');
    
    jmax = 0; % Will be set by readkappa

    %% Main Execution Block
    jmax = readkappa(mykappafile);
    
    for K = 0:1  % C++ used for(int K=0;K<2;K++)
        
        for j = 1:jmax
            if K == 1
                if (nu(j) > 3/18) && (nu(j) < 3/14)
                    kappanu(j) = max(0.5, 1.2*kappanu(j)); % Visible R
                end
            elseif K == 2 % If you expand loop to K=2
                kappanu(j) = 0.5; % constant kappa
            end
        end
        
        resultfilename = [myresulttemperature, num2str(K), '.txt'];
        fprintf('\n iterations \t [T] near earth [T] far ||G|| and ||S||\n');
        
        t0 = tic;
        
        % Run the multi-block fixed point iterations
        multiBlock(Te/2);
        
        cpu_time = toc(t0);
        fprintf(' Time CPU = %10.6f\n', cpu_time);
        fprintf('\n tau\t [T]:\n');
        
        % Write to file and console
        fileID = fopen(resultfilename, 'w');
        for i = 1:Nz-1
            z_alt = backtoz((i-1)*dz);
            T_celcius = T(i)*T00_val - 273;
            fprintf('%f %f\n', z_alt, T_celcius);
            if fileID > 0
                fprintf(fileID, '%f %f\n', z_alt, T_celcius);
            end
        end
        if fileID > 0
            fclose(fileID);
        end
    end

    %% Nested Helper Functions 
    % Note: These have access to the main variables initialized above.
    
    function val = kapparhoz(nu_val, z_val)
        val = 1 + cloud*((z_val>z1)-(z_val>z2)) * ...
            ((nu_val>3/10)*(nu_val<3/1) + (nu_val<3/20));
    end

    function val = as(z_val, nu_val)
        val = 0.3 ...
            + 0.3*4*(z2-z_val)*(z_val-z1)*(z_val>z1)*(z_val<z2)/((z1-z2)^2) ...
            + 0.3*16*(nu_val<nu2)*(nu_val>nu1)*(((nu_val-nu1)*(nu_val-nu2)/((nu1-nu2)^2))^2)*(z_val>z3);
    end

    function res = expint_E1(t)
        if nargin < 1, t = 1; end
        t1 = abs(t);
        Kexpint = 14 + (t1-1)*4;
        gaNtaua = 0.577215664901533;
        
        if t1 < 1e-5
            res = 0; return;
        end
        if t1 > 4
            if verbose, disp('argument of E_1>2.5'); end
            tx = 1/t1;
            res = exp(-t1)*tx * (1 +(-1+(2+(-6+(24+(-120+720*tx)*tx)*tx)*tx)*tx)*tx);
            return;
        end
        
        ak = t1;
        soNtaue = -gaNtaua - log(t1) + ak;
        for k = 2:(floor(Kexpint)-1)
            ak = ak * (-t1*(k-1)/(k^2));
            soNtaue = soNtaue + ak;
        end
        res = soNtaue;
    end

    function res = expint_E2(t)
        if nargin < 1, t = 1; end
        t1 = abs(t); res = exp(-t1) - t1*expint_E1(t1);
    end
    function res = expint_E3(t)
        if nargin < 1, t = 1; end
        t1 = abs(t); res = (exp(-t1) - t1*expint_E2(t1))/2;
    end
    function res = expint_E4(t)
        if nargin < 1, t = 1; end
        t1 = abs(t); res = (exp(-t1) - t*expint_E3(t1))/3;
    end
    function res = expint_E5(t)
        if nargin < 1, t = 1; end
        t1 = abs(t); res = (exp(-t1) - t*expint_E4(t1))/4;
    end

    function val = Bsun(nu_val)
        val = (nu_val^2)*nu_val/(exp(nu_val/Ts) - 1);
    end
    function val = Bearth(nu_val)
        val = (nu_val^2)*nu_val/(exp(nu_val/Te) - 1);
    end
    function val = BB(nu_val, T_val)
        val = (nu_val^2)*nu_val/(exp(min(nu_val/T_val, 100)) - 1);
    end
    function val = dBB(nu_val, T_val)
        a = exp(nu_val/T_val);
        if a > 1e100
            val = (nu_val^2/T_val) * 1e-30;
        else
            val = a * (nu_val^2/max(1.e-30, a - 1)^2) / T_val;
        end
    end

    function num_freqs = readkappa(filename)
        fprintf('reading kappa(nu) from file %s\n', filename);
        if ~isfile(filename)
            % create a dummy array for the sake of not crashing if file isn't present
            warning(['Cannot open file: ', filename, '. Using dummy data.']);
            data = [linspace(1, 10, 50)', rand(50,1)*0.5+0.1];
        else
            data = load(filename);
        end
        
        num_freqs = min(size(data, 1), jmaxmax);
        
        for j_idx = 1:num_freqs
            waveno = data(j_idx, 1);
            kappaux = data(j_idx, 2);
            kappanu(j_idx) = max(kappaux, kappamin);
            nu(j_idx) = 3 / waveno;
            kappasum(j_idx, 1) = 0; % MATLAB 1-based index (0 altitude)
            
            for it = 1:Nz
                z_val = (it-1)*dz;
                kappaz(j_idx, it) = kapparhoz(nu(j_idx), z_val);
                if it > 1
                    kappasum(j_idx, it) = kappasum(j_idx, it-1) + ...
                        (kappaz(j_idx, it) + kappaz(j_idx, it-1))*dz/2;
                end
            end
        end
        fprintf('Number of frequencies %d\n', num_freqs);
    end

    function updateJK()
        for jnu = 1:jmax
            H = zeros(1, Nz);
            S = zeros(1, Nz);
            kappanuj = kappanu(jnu); 
            nuj = nu(jnu);
            c0 = -q0 * mus * Cs * Bsun(nuj) * exp(-kappanuj * kappasum(jnu, Nz) / mus);
            
            for i = 1:Nz
                z_val = (i-1)*dz;
                kappazj = kappaz(jnu, i) * kappanuj;
                asz = as(z_val, nuj);
                sigs = kappazj * asz;
                siga = kappazj * (1 - asz);
                
                S(i) = sigs * J0old(jnu, i) + siga * BB(nuj, T(i));
                H(i) = 9*beta_val*sigs/8 * (J2old(jnu, i) - J0old(jnu, i)/3 - K0old(jnu, i) + K2old(jnu, i));
                c0 = c0 - q0 * expint_E2(kappasum(jnu, i)*kappanuj) * S(i) * dz;
            end
            
            for i = 1:Nz
                ksz = kappasum(jnu, i) * kappanuj;
                J0z = Ce * Bearth(nuj) * expint_E3(ksz)/2 ...
                    + Cs * Bsun(nuj) * exp(-(kappasum(jnu, Nz) - kappasum(jnu, i))*kappanuj/mus)/2 ...
                    + c0 * expint_E2(ksz)/2;
                J2z = Ce * Bearth(nuj) * expint_E5(ksz)/2 ...
                    + Cs * (mus^2) * Bsun(nuj) * exp(-(kappasum(jnu, Nz) - kappasum(jnu, i))*kappanuj/mus)/2 ...
                    + c0 * expint_E4(ksz)/2;
                K0z = 0; K2z = 0;
                
                for j_idx = 2:Nz
                    Hj = (H(j_idx) + H(j_idx-1))/2;
                    Sj = (S(j_idx) + S(j_idx-1))/2 - Hj/3;
                    aux = kappanuj * (kappasum(jnu, i) - (kappasum(jnu, j_idx) + kappasum(jnu, j_idx-1))/2);
                    
                    expE1ij = expint_E1(aux);
                    expE3ij = expint_E3(aux);
                    expE5ij = expint_E5(aux);
                    
                    J0z = J0z + (expE1ij*Sj + expE3ij*Hj)*dz/2;
                    J2z = J2z + (expE3ij*Sj + expE5ij*Hj)*dz/2;
                    K0z = K0z - (expE1ij - expE3ij)*Hj*dz/2;
                    K2z = K2z - (expE3ij - expE5ij)*Hj*dz/2;
                end
                
                J0(jnu, i) = J0z;
                J2(jnu, i) = J2z;
                K0(jnu, i) = K0z;
                K2(jnu, i) = K2z;
            end
        end
    end

    function myeq = root(rhs, T0_root, i)
        myeq = -rhs;
        for j_idx = 2:jmax
            myeq = myeq + (1 - as((i-1)*dz, nu(j_idx))) * kappanu(j_idx) * ...
                BB((nu(j_idx) + nu(j_idx-1))/2, T0_root) * (nu(j_idx) - nu(j_idx-1));
        end
    end

    function rhs = rhsTeq(i)
        rhs = 0;
        for j_idx = 2:jmax
            rhs = rhs + (1 - as((i-1)*dz, nu(j_idx))) * kappanu(j_idx) * J0(j_idx, i) * (nu(j_idx) - nu(j_idx-1));
        end
    end

    function T_ret = getTbydycho(rhs, i, Tstart)
        T0_root = Tstart; T1_root = 2*T0_root;
        if rhs == 0
            T_ret = 0; return;
        end
        counter = 0;
        myeq0 = root(rhs, T0_root, i);
        myeq1 = root(rhs, T1_root, i);
        
        while (myeq0 > 0 && counter < 100)
            T0_root = T0_root / 2; counter = counter + 1;
            myeq0 = root(rhs, T0_root, i);
            if verbose, fprintf('%f down %f\n', T0_root, myeq0); end
        end
        while (myeq1 < 0 && counter < 100)
            T1_root = 2 * T1_root; counter = counter + 1;
            myeq1 = root(rhs, T1_root, i);
            if verbose, fprintf('%f up %f\n', T1_root, myeq1); end
        end
        
        if verbose
            myeq0 = root(rhs, T0_root, i);
            myeq1 = root(rhs, T1_root, i);
            if myeq0 > myeq1, disp('BUG in dychotomy'); end
        end
        
        while (T1_root - T0_root > epsdycho && counter < 100)
            Taux = (T1_root + T0_root) / 2; counter = counter + 1;
            myeq0 = root(rhs, Taux, i);
            if myeq0 > 0
                T1_root = Taux; 
            else
                T0_root = Taux;
            end
            if verbose
                myeq0_disp = root(rhs, T0_root, i);
                myeq1_disp = root(rhs, T1_root, i);
                fprintf('%f middle %f %f %f\n', T0_root, T1_root, myeq0_disp, myeq1_disp);
            end
        end
        if counter > 99, disp('Divergence in dichotomy'); end
        T_ret = (T1_root + T0_root) / 2;
    end

    function genT()
        for i = 1:Nz
            rhs = rhsTeq(i);
            T(i) = getTbydycho(rhs, i, T(i));
            presfunc = 1.0;
            inewton = 0;
            
            while (inewton < newton_iters && abs(presfunc) > epsnewton)
                inewton = inewton + 1;
                T0_val = T(i);
                left = 0; deriv = 0;
                
                for j_idx = 2:jmax
                    dnu = nu(j_idx) - nu(j_idx-1);
                    nu11 = (nu(j_idx) + nu(j_idx-1)) / 2;
                    kappaux = (1 - as((i-1)*dz, nu(j_idx))) * (kappanu(j_idx) + kappanu(j_idx-1)) / 2;
                    left = left + kappaux * BB(nu11, T0_val) * dnu;
                    deriv = deriv + kappaux * dBB(nu11, T0_val) * dnu;
                end
                
                presfunc = rhs - left;
                if abs(deriv) > 1e-10
                    T(i) = T0_val + presfunc/deriv;
                end
                
                if verbose
                    fprintf('%d T= %f residue= %f deriv= %f Tstefan= %f\n', ...
                        i, T(i), presfunc, deriv, sqrt(sqrt(rhs*15*2))/pi_val);
                end
            end
            if inewton >= newton_iters
                disp('Newton precision doubtful');
            end
        end
    end

    function multiBlock(initT)
        T(:) = initT; % initialize
        J0old(:) = 0; J2old(:) = 0; K0old(:) = 0; K2old(:) = 0;
        
        for k = 1:kmax
            updateJK();
            genT();
            
            % Copy values
            J0old = J0; J2old = J2; K0old = K0; K2old = K2;
            
            normG = 0;
            for j_idx = 2:jmax
                for i = 1:Nz
                    normG = normG + abs(J0(j_idx, i)) * (nu(j_idx) - nu(j_idx-1)) * dz;
                end
            end
            % MATLAB maps T[2] to T(3)
            fprintf('k= %d %f  %f\n', k-1, T(3)*4798-273, normG);
        end
    end

    function z_out = backtoz(z_val)
        if drho == 0
            z_out = z_val;
        else
            z_out = (sqrt(1 + 2*drho*z_val/rho0) - 1) / drho;
        end
    end

end