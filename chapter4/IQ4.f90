!=======================================================================
! MAIN PROGRAM use gfortran file
!=======================================================================
        program main
        implicit none
        
        ! Define double precision kind
        integer, parameter :: wp = selected_real_kind(15, 307)
        
        ! Constants
        real(wp), parameter :: pi = 4.0_wp * atan(1.0_wp)
        real(wp), parameter :: stefan = ((pi**2)**2) / 15.0_wp
        
        real(wp), parameter :: Ce = 2.0_wp
        real(wp), parameter :: Cs = 2.0e-6_wp
        real(wp), parameter :: rho0 = 1.0_wp
        real(wp), parameter :: drho = -0.7_wp
        real(wp), parameter :: TOA = 1.2_wp
        real(wp), parameter :: Z = TOA * rho0*(1.0_wp+drho*TOA/2.0_wp)
        real(wp), parameter :: B0 = 1.4744e-8_wp
        real(wp), parameter :: T0_val = 4798.0_wp
        
        real(wp), parameter :: Te = (273.0_wp + 18.0_wp) / T0_val
        real(wp), parameter :: Ts = 5778.0_wp / T0_val
        real(wp), parameter :: q0 = -0.3_wp
        real(wp), parameter :: mus = 0.5_wp
        real(wp), parameter :: beta_val = 0.5_wp
        
        logical, parameter  :: verbose = .false.
        integer, parameter  :: Nz = 100
        real(wp), parameter :: dz = Z / real(Nz - 1, wp)
        integer, parameter  :: kmax = 14
        integer, parameter  :: jmaxmax = 600
        integer, parameter  :: newton_iters = 50
        real(wp), parameter :: epsdycho = 1.0e-4_wp
        real(wp), parameter :: epsnewton = 1.0e-10_wp
        real(wp), parameter :: kappamin = 0.001_wp
        
        real(wp), parameter :: z1 = Z * 0.5_wp
        real(wp), parameter :: z2 = Z * 0.7_wp
        real(wp), parameter :: z3 = 0.8_wp * Z
        real(wp), parameter :: nu1 = 0.5_wp
        real(wp), parameter :: nu2 = 1.0_wp
        real(wp), parameter :: cloud = 1.0_wp

        ! Global Arrays (1-based indexing in Fortran)
        real(wp), dimension(jmaxmax) :: kappanu, nu
        real(wp), dimension(jmaxmax,Nz)::J0,J0old,J2,J2old,K0,K0old,K2,K2old
        real(wp), dimension(Nz) :: T
        real(wp), dimension(jmaxmax, Nz) :: kappaz, kappasum
        
        integer :: jmax
        
        integer :: K_idx, j, i
        character(len=256) :: basedir,mykappafile
        character(len=256) :: myresulttemp,resultfilename
        character(len=10) :: K_str
        real(wp) :: t_start, t_end
!	/Users/pironneau/Dropbox/aranger/TeX2026/BookVRTE        
        basedir='prog4/'
        mykappafile = trim(basedir) // '_kappa.txt'
        myresulttemp = trim(basedir) // 'temperatureyyx'
        
        call readkappa(mykappafile)
        
        do K_idx = 0, 1
          do j = 1, jmax
            if (K_idx == 1) then
              if (nu(j)>3.0_wp/18.0_wp.and.nu(j)<3.0_wp/14.0_wp) then
                        kappanu(j) = max(0.5_wp, 1.2_wp * kappanu(j))
                end if
                else if (K_idx == 2) then 
                            kappanu(j) = 0.5_wp
                end if
            end do
                
            write(K_str, '(I0)') K_idx
            resultfilename = trim(myresulttemp) // trim(K_str) // '.txt'
                
            print*,char(10),"iter [T] near [T] far ||G|| and ||S||"
                
            call cpu_time(t_start)
                
            call multiBlock(Te / 2.0)
                
            call cpu_time(t_end)
            print '(" Time CPU = ", F10.6)', t_end - t_start
                
            print *, char(10), " tau          [T]:"
            open(unit=20, file=trim(resultfilename), status='unknown')
            do i = 1, Nz - 1
              print *, backtoz(real(i - 1, wp) * dz), T(i) * T0_val - 273.0_wp
         write(20, *) backtoz(real(i - 1, wp) * dz), T(i) * T0_val - 273.0_wp
            end do
                close(20)          
        end do
		end program main



        !------------------------------------------------------------
        ! Helper Functions for Logical -> Numeric Conversion & Math
        !------------------------------------------------------------
        function kapparhoz(nu_val, z_val) result(res)
                real(wp), intent(in) :: nu_val, z_val
                real(wp) :: res, cloud_term, freq_term
                
                cloud_term = 0.0_wp
                if (z_val > z1) cloud_term = cloud_term + 1.0_wp
                if (z_val > z2) cloud_term = cloud_term - 1.0_wp
                
                freq_term = 0.0_wp
        if (nu_val>0.3_wp.and.nu_val<3.0_wp)freq_term=freq_term+1.0_wp
                if (nu_val < 0.15_wp) freq_term = freq_term + 1.0_wp
                
                res = 1.0_wp + cloud * cloud_term * freq_term
        end function kapparhoz

        function as_func(z_val, nu_val) result(res)
                real(wp), intent(in) :: z_val, nu_val
                real(wp) :: res, term1, term2, term3
                
                term1 = 0.3_wp
                term2 = 0.0_wp
        if (z_val > z1 .and. z_val < z2) then
        term2 =0.3_wp*4.0_wp*(z2-z_val)*(z_val-z1)/((z1-z2)**2)
                end if
                term3 = 0.0_wp
        if (nu_val < nu2 .and. nu_val > nu1 .and. z_val > z3) then
            term3=0.3_wp*16.0_wp*(((nu_val-nu1)*(nu_val-nu2)/((nu1-nu2)**2))**2)
             end if
                res = term1 + term2 + term3
        end function as_func

        !------------------------------------------------------------
        ! Exponential Integral Functions
        !------------------------------------------------------------
        function expint_E1(t_in) result(res)
                real(wp), intent(in), optional :: t_in
                real(wp) :: t_val, t1, tx, ak, soNtaue
                integer :: Kexpint, k
                real(wp), parameter :: gaNtaua = 0.577215664901533_wp
                
                if (present(t_in)) then
                        t_val = t_in
                else
                        t_val = 1.0_wp
                end if
                
                t1 = abs(t_val)
                if (t1 < 1.0e-5_wp) then
                        res = 0.0_wp
                        return
                end if
                
                if (t1 > 4.0_wp) then
                        if (verbose) print *, "argument of E_1 > 2.5"
                        tx = 1.0_wp / t1
                        res = exp(-t1)*tx * (1.0_wp + (-1.0_wp + (2.0_wp + (-6.0_wp + &
                                 (24.0_wp + (-120.0_wp + 720.0_wp*tx)*tx)*tx)*tx)*tx)*tx)
                        return
                end if
                
                Kexpint = 14 + int((t1 - 1.0_wp) * 4.0_wp)
                ak = t1
                soNtaue = -gaNtaua - log(t1) + ak
                
                do k = 2, Kexpint - 1
                        ak = ak * (-t1 * real(k - 1, wp) / real(k**2, wp))
                        soNtaue = soNtaue + ak
                end do
                res = soNtaue
        end function expint_E1

        function expint_E2(t_in) result(res)
                real(wp), intent(in), optional :: t_in
                real(wp) :: t_val, t1
                if (present(t_in)) then; t_val = t_in; else; t_val = 1.0_wp; end if
                t1 = abs(t_val)
                res = exp(-t1) - t1 * expint_E1(t1)
        end function expint_E2

        function expint_E3(t_in) result(res)
                real(wp), intent(in), optional :: t_in
                real(wp) :: t_val, t1
                if (present(t_in)) then; t_val = t_in; else; t_val = 1.0_wp; end if
                t1 = abs(t_val)
                res = (exp(-t1) - t1 * expint_E2(t1)) / 2.0_wp
        end function expint_E3

        function expint_E4(t_in) result(res)
                real(wp), intent(in), optional :: t_in
                real(wp) :: t_val, t1
                if (present(t_in)) then; t_val = t_in; else; t_val = 1.0_wp; end if
                t1 = abs(t_val)
                res = (exp(-t1) - t_val * expint_E3(t1)) / 3.0_wp
        end function expint_E4

        function expint_E5(t_in) result(res)
                real(wp), intent(in), optional :: t_in
                real(wp) :: t_val, t1
                if (present(t_in)) then; t_val = t_in; else; t_val = 1.0_wp; end if
                t1 = abs(t_val)
                res = (exp(-t1) - t_val * expint_E4(t1)) / 4.0_wp
        end function expint_E5

        !------------------------------------------------------------
        ! Blackbody Functions
        !------------------------------------------------------------
        function Bsun(nu_val) result(res)
                real(wp), intent(in) :: nu_val
                real(wp) :: res
                res = (nu_val**2) * nu_val / (exp(nu_val / Ts) - 1.0_wp)
        end function Bsun

        function Bearth(nu_val) result(res)
                real(wp), intent(in) :: nu_val
                real(wp) :: res
                res = (nu_val**2) * nu_val / (exp(nu_val / Te) - 1.0_wp)
        end function Bearth

        function BB(nu_val, T_val) result(res)
                real(wp), intent(in) :: nu_val, T_val
                real(wp) :: res, exp_val
                exp_val = exp(min(nu_val / T_val, 100.0_wp))
                res = (nu_val**2) * nu_val / (exp_val - 1.0_wp)
        end function BB

        function dBB(nu_val, T_val) result(res)
                real(wp), intent(in) :: nu_val, T_val
                real(wp) :: res, a
                a = exp(nu_val / T_val)
                if (a > 1.0e100_wp) then
                        res = (nu_val**2 / T_val) * 1.0e-30_wp
                else
                        res = a * (nu_val**2 / max(1.0e-30_wp, a - 1.0_wp)**2) / T_val
                end if
        end function dBB

        function backtoz(z_val) result(res)
                real(wp), intent(in) :: z_val
                real(wp) :: res
                if (abs(drho) < 1.0e-12_wp) then
                        res = z_val
                else
                        res = (sqrt(1.0_wp + 2.0_wp * drho * z_val / rho0) - 1.0_wp) / drho
                end if
        end function backtoz

        !------------------------------------------------------------
        ! Main Physics Routines
        !------------------------------------------------------------
        subroutine readkappa(filename)
                character(len=*), intent(in) :: filename
                integer :: io_status, j_idx, it
                real(wp) :: waveno, kappaux, z_val
                
                print *, "Reading kappa(nu) from file ", trim(filename)
                
                open(unit=10, file=trim(filename), status='old', action='read', iostat=io_status)
                if (io_status /= 0) then
                        print *, "Cannot open file: ", trim(filename)
                        stop
                end if
                
                j_idx = 0
                do while (j_idx < jmaxmax)
                        read(10, *, iostat=io_status) waveno, kappaux
                        if (io_status < 0) exit ! EOF
                        
                        j_idx = j_idx + 1
                        kappanu(j_idx) = max(kappaux, kappamin)
                        nu(j_idx) = 3.0_wp / waveno
                        kappasum(j_idx, 1) = 0.0_wp
                        
                        do it = 1, Nz
                                z_val = real(it - 1, wp) * dz
                                kappaz(j_idx, it) = kapparhoz(nu(j_idx), z_val)
                                if (it > 1) then
                                        kappasum(j_idx, it) = kappasum(j_idx, it - 1) + &
                                                (kappaz(j_idx, it) + kappaz(j_idx, it - 1)) * dz / 2.0_wp
                                end if
                        end do
                end do
                
                close(10)
                jmax = j_idx
                print *, "Number of frequencies ", jmax
        end subroutine readkappa

        subroutine updateJK()
                integer :: jnu, i, j_idx
                real(wp) :: H_arr(Nz), S_arr(Nz)
                real(wp) :: kappanuj, nuj, c0, z_val, kappazj, asz, sigs, siga
                real(wp) :: ksz, J0z, J2z, K0z, K2z, Hj, Sj, aux
                real(wp) :: expE1ij, expE3ij, expE5ij
                
                do jnu = 1, jmax
                        kappanuj = kappanu(jnu)
                        nuj = nu(jnu)
                        c0 = -q0 * mus * Cs * Bsun(nuj) * exp(-kappanuj * kappasum(jnu, Nz) / mus)
                        
                        do i = 1, Nz
                                z_val = real(i - 1, wp) * dz
                                kappazj = kappaz(jnu, i) * kappanuj
                                asz = as_func(z_val, nuj)
                                sigs = kappazj * asz
                                siga = kappazj * (1.0_wp - asz)
                                
                                H_arr(i) = 9.0_wp * beta_val * sigs / 8.0_wp * &
                                                   (J2old(jnu, i) - J0old(jnu, i)/3.0_wp - K0old(jnu, i) + K2old(jnu, i))
                                S_arr(i) = sigs * J0old(jnu, i) + siga * BB(nuj, T(i)) - H_arr(i) /3.
                                c0 = c0 - q0 * expint_E2(kappasum(jnu, i) * kappanuj) * S_arr(i) * dz
                        end do
                        
                        do i = 1, Nz
                                ksz = kappasum(jnu, i) * kappanuj
                                J0z = Ce * Bearth(nuj) * expint_E3(ksz)/2.0_wp &
                                        + Cs * Bsun(nuj) * exp(-(kappasum(jnu, Nz) - kappasum(jnu, i))*kappanuj/mus)/2.0_wp &
                                        + c0 * expint_E2(ksz)/2.0_wp
                                
                                J2z = Ce * Bearth(nuj) * expint_E5(ksz)/2.0_wp &
                                        + Cs * (mus**2) * Bsun(nuj) * exp(-(kappasum(jnu, Nz) - kappasum(jnu, i))*kappanuj/mus)/2.0_wp &
                                        + c0 * expint_E4(ksz)/2.0_wp
                                
                                K0z = 0.0_wp
                                K2z = 0.0_wp
                                
                                do j_idx = 2, Nz
                                        Hj = (H_arr(j_idx) + H_arr(j_idx - 1)) / 2.0_wp
                                        Sj = (S_arr(j_idx) + S_arr(j_idx - 1)) / 2.0_wp
                                        aux = kappanuj * (kappasum(jnu, i) - (kappasum(jnu, j_idx) + kappasum(jnu, j_idx - 1)) / 2.0_wp)
                                        
                                        expE1ij = expint_E1(aux)
                                        expE3ij = expint_E3(aux)
                                        expE5ij = expint_E5(aux)
                                        
                                        J0z = J0z + (expE1ij * Sj + expE3ij * Hj) * dz / 2.0_wp
                                        J2z = J2z + (expE3ij * Sj + expE5ij * Hj) * dz / 2.0_wp
                                        K0z = K0z - (expE1ij - expE3ij) * Hj * dz / 2.0_wp
                                        K2z = K2z - (expE3ij - expE5ij) * Hj * dz / 2.0_wp
                                end do
                                
                                J0(jnu, i) = J0z
                                J2(jnu, i) = J2z
                                K0(jnu, i) = K0z
                                K2(jnu, i) = K2z
                        end do
                end do
        end subroutine updateJK

        function root_func(rhs, T0_root, i) result(myeq)
                real(wp), intent(in) :: rhs, T0_root
                integer, intent(in) :: i
                real(wp) :: myeq, z_val
                integer :: j_idx
                
                myeq = -rhs
                z_val = real(i - 1, wp) * dz
                do j_idx = 2, jmax
                        myeq = myeq + (1.0_wp - as_func(z_val, nu(j_idx))) * kappanu(j_idx) * &
                                   BB((nu(j_idx) + nu(j_idx - 1)) / 2.0_wp, T0_root) * (nu(j_idx) - nu(j_idx - 1))
                end do
        end function root_func

        function rhsTeq(i) result(rhs)
                integer, intent(in) :: i
                real(wp) :: rhs, z_val
                integer :: j_idx
                
                rhs = 0.0_wp
                z_val = real(i - 1, wp) * dz
                do j_idx = 2, jmax
                        rhs = rhs + (1.0_wp - as_func(z_val, nu(j_idx))) * kappanu(j_idx) * &
                                  J0(j_idx, i) * (nu(j_idx) - nu(j_idx - 1))
                end do
        end function rhsTeq

        function getTbydycho(rhs, i, Tstart) result(T_ret)
                real(wp), intent(in) :: rhs, Tstart
                integer, intent(in) :: i
                real(wp) :: T_ret, T0_root, T1_root, Taux, myeq0, myeq1
                integer :: counter
                
                T0_root = Tstart
                T1_root = 2.0_wp * T0_root
                if (abs(rhs) < 1.0e-30_wp) then
                        T_ret = 0.0_wp
                        return
                end if
                
                counter = 0
                myeq0 = root_func(rhs, T0_root, i)
                myeq1 = root_func(rhs, T1_root, i)
                
                do while (myeq0 > 0.0_wp .and. counter < 100)
                        T0_root = T0_root / 2.0_wp
                        counter = counter + 1
                        myeq0 = root_func(rhs, T0_root, i)
                        if (verbose) print *, T0_root, " down ", myeq0
                end do
                
                do while (myeq1 < 0.0_wp .and. counter < 100)
                        T1_root = 2.0_wp * T1_root
                        counter = counter + 1
                        myeq1 = root_func(rhs, T1_root, i)
                        if (verbose) print *, T1_root, " up ", myeq1
                end do
                
                if (verbose .and. myeq0 > myeq1) print *, "BUG in dychotomy"
                
                do while (T1_root - T0_root > epsdycho .and. counter < 100)
                        Taux = (T1_root + T0_root) / 2.0_wp
                        counter = counter + 1
                        myeq0 = root_func(rhs, Taux, i)
                        if (myeq0 > 0.0_wp) then
                                T1_root = Taux
                        else
                                T0_root = Taux
                        end if
                        if (verbose) print *, T0_root, " middle ", T1_root
                end do
                
                if (counter > 99) print *, "Divergence in dichotomy"
                T_ret = (T1_root + T0_root) / 2.0_wp
        end function getTbydycho

        subroutine genT()
                integer :: i, j_idx, inewton
                real(wp) :: rhs, presfunc, T0_val, left_val, deriv, dnu, nu11, kappaux, z_val
                
                do i = 1, Nz
                        rhs = rhsTeq(i)
                        T(i) = getTbydycho(rhs, i, T(i))
                        presfunc = 1.0_wp
                        inewton = 0
                        
                        do while (inewton < newton_iters .and. abs(presfunc) > epsnewton)
                                inewton = inewton + 1
                                T0_val = T(i)
                                left_val = 0.0_wp
                                deriv = 0.0_wp
                                z_val = real(i - 1, wp) * dz
                                
                                do j_idx = 2, jmax
                                        dnu = nu(j_idx) - nu(j_idx - 1)
                                        nu11 = (nu(j_idx) + nu(j_idx - 1)) / 2.0_wp
                                        kappaux = (1.0_wp - as_func(z_val, nu(j_idx))) * &
                                                          (kappanu(j_idx) + kappanu(j_idx - 1)) / 2.0_wp
                                        left_val = left_val + kappaux * BB(nu11, T0_val) * dnu
                                        deriv = deriv + kappaux * dBB(nu11, T0_val) * dnu
                                end do
                                
                                presfunc = rhs - left_val
                                if (abs(deriv) > 1.0e-10_wp) T(i) = T0_val + presfunc / deriv
                                
                                if (verbose) print *, i, " T= ", T(i), " residue= ", presfunc
                        end do
                        if (inewton >= newton_iters) print *, "Newton precision doubtful"
                end do
        end subroutine genT

        subroutine multiBlock(initT)
                real(wp), intent(in) :: initT
                integer :: k, i, j_idx
                real(wp) :: normG
                
                T = initT
                J0old = 0.0_wp
                J2old = 0.0_wp
                K0old = 0.0_wp
                K2old = 0.0_wp
                
                do k = 1, kmax
                        call updateJK()
                        call genT()
                        
                        J0old = J0
                        J2old = J2
                        K0old = K0
                        K2old = K2
                        
                        normG = 0.0_wp
                        do j_idx = 2, jmax
                                do i = 1, Nz
                                        normG = normG + abs(J0(j_idx, i)) * (nu(j_idx) - nu(j_idx - 1)) * dz
                                end do
                        end do
                        ! T(3) corresponds to T[2] in 0-based C++
                        print *, "k= ", k - 1, T(3) * 4798.0_wp - 273.0_wp, normG
                end do
        end subroutine multiBlock