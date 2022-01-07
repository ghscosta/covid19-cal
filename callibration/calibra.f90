program corona
!LEGENDA: 0 = SUSCET�VEL ; 1 = EXPOSTO ; 2 = ASSINTOM�TICO ; 3 = INFECTADO ; 4 = RECUPERADO ; 5 = CONFIRMADO
implicit none

real*8, parameter :: kcont = 15.0	 !N�MERO M�DIO DE CONTATOS FORA DE CASA
integer, parameter :: tmax = 14     !DIAS DE SIMULA��O
integer, parameter :: tmax2 = 28
integer :: iseed,amostra,cidade      !SEED, AMOSTRA E CIDADE
integer, parameter :: ninte = 100     !N�MERO DE ITERA��ES PARA CALIBRAGEM
integer, parameter :: nsampt = 1     !N�MERO DE ITERA��ES PARA CALIBRAGEM
logical :: flag

character*20 :: day_ini
integer :: cik,dias,diaInicial,diaFinal
character*20 :: datas(30)

!TAXAS DE TRANSI��O ENTRE OS DIFERENTES COMPARTIMENTOS (CRUAS)
real*8 :: gu, ga, teste_U, teste_A,undNot,b,c,d,f,alpha,aA,aU

!PAR�METROS DE MODULA��O DO N�MERO DE CONTATOS E ISOLAMENTO SOCIAL
real*8 :: pmuni
integer :: conf0,conf1,conf2

!VARI�VEIS PARA CONTROLAR INDIV�DUOS NAS CIDADES
real*8 :: Ns,Ne,Na,Ni,Nr,Nag,Nc,Nrc

!TAXAS DE TRANSI��O DO MODELO
real*8 :: rdif,rsea,rsei,rea,rai,rir,pq,rtot,rpmax,ramax,rar,rsec,rcr,ric,ale,rac
real*8 :: r0,r0up,r0do,rt,rtup,rtdo,und,undup,undo
real*8 :: betaA,betaU,betaAdo,betaUdo,betaAup,betaUup

!VARI�VEIS DE MINIMIZA��O PARA CALIBRAGEM
real*8 :: ref(0:tmax2),Ncc
real*8 :: dif(0:tmax2)
real*8 :: dift,dift0
real*8 :: difer,difer0,dp
real*8 :: minimo(0:tmax),beta_min
real*8 ::dif_min
logical :: varia, escrev
integer :: Nam
real*8 :: Nsmin,Ncmin,Nemin,Namin,Nrmin,Nimin

!VARI�VEIS AUXILIARES
integer :: patin,patout,contif,nsamp
real*8 :: dt,t
integer :: Neff,tgt,Nefft,c0,i,j,inte
integer :: ab1,ab3,tser1,tser,tser3,pat,tser12
character*20 :: ab2,str,amt,day,city,dayz
real*8 :: rnd,p,q,rnd2,z,cont,Nout,tser2,tser4,tser5,ctf,taxA,taxU,taxC
integer :: tper1,tper2,t_sun,t_sat,tx
real*8 ::ax1,ax2,ax3,ncal,Nconf(0:tmax2),zc,bcd,incub,zero,lethNom,lethReal,fih,incub2,incub3,incub1
real*8 :: testing,lee,recuperados
character*20 :: mn,ab,label,amos
integer :: a,iErr,ct,cini,cifi,ncid,nga,nsa

!Inputs
character*199 :: timeSeries,fatSeries,dataIni
integer :: pop

    amostra = 1
    iseed = 1239686
        
    read*, timeSeries,fatSeries,dataIni,pop

    day_ini = dataIni    
    ct = 0

    open(133,file=timeSeries)
    open(266,file=fatSeries)

             do
                read(266,*,IOSTAT = iErr)ab,lee
                if(ab == dataIni)exit
                ct=ct+1
            end do

            lethNom = 1d0*lee

            ct = 0
            ref = 0
            do
                read(133,*,IOSTAT = iErr)ab,a
                if(ab == day_ini)exit
                ct=ct+1
            end do

            do j = 0,tmax2
                read(133,*,IOSTAT = iErr)ab,a
                if(iErr /= 0 )exit
                ref(j) = a
            end do

            conf0 = ref(0)
            conf1 = ref(7)

            close(133)
            close(266)

            call calcula_rec()

            flag = .false.

            call taxas()

            pmuni = 0.01
            dp = 0.01

            dift0 = 10000000000000d0
            dif_min = dift0
            difer0 = 0
            minimo  = 0

            do nsa = 1,ninte

                dif = 0
                Nconf = 0

                call cond_ini()
                if(Ns == 0) go to 188
                call calcula_taxas()
                call rk4(tmax)
                call minimiza_quadrados()
                call atualiza_profilaxia()

            end do

            dp = 0.001

            do nsa = 1,ninte

                dif = 0
                Nconf = 0

                call cond_ini()
                if(Ns == 0) go to 188
                call calcula_taxas()
                call rk4(tmax)

                call minimiza_quadrados()
                call atualiza_profilaxia()

            end do

            flag = .True.

            do nga = 1,1000

                    write(amos,*)nga

                    call taxas()
                    call cond_ini()
                    call calcula_taxas()

                    open(277,file = 'epiQuantities.dat')
                    open(278,file = 'hiddenCompart.dat')
                    zero = 0d0
                    write(277,'(7F20.10)')zero,Nc,ref(0),1d0*r0*(Ns/pop),1d0*Ns/pop,1d0*Nr/Nc,1d0*nga


                    call rk4(tmax2)

            end do

        188 close(133)
        close(266)

	contains

    subroutine calcula_rec()

       real*8 :: lethPart,casos0,casos1
       integer :: linha
       real*8 :: a,b2,alpha2,c_bas,delta_bas
       character*20 :: date

       open(133,file=timeSeries)
       open(266,file=fatSeries)

       incub = 6.4
       lethReal = 0.75

       b = 2d0/incub
       c = 2d0/incub
       d = 1d0/3.2
       f = 1d0/(1d0/d + 1d0/c)
       alpha = 0.06
       aA = alpha
       aU = alpha
       fih = 1d0*b/(b+f)

       read(133,*)date,casos0
       read(266,*)date,lethPart

       recuperados = 0
       do
           do linha = 1,7
                read(133,*)date,casos1
                read(266,*)date,lethPart
           end do

           fih = 1d0*b/(b+f)

            alpha2 = 1.0
            a = -1.0 * alpha2 * alpha2
            b2 =  (1.0 + alpha2)

            c_bas = lethreal/(fih*lethNom)
            delta_bas = (b2*b2 - 4.0*a*c_bas)**(1d0/2d0)
            teste_U = -1.0*b2 + delta_bas
            teste_U = -teste_U / (2.0*a)
            teste_A = teste_U * alpha2


            if(teste_U .gt. 0.9) teste_U = 0.9
            if(teste_A .gt. 0.9) teste_A = 0.9

            undNot = (1.0-teste_A)*(1.0-fih*teste_U)
            undNot = 1.0*undNot/(teste_A + (1.0 - teste_A)*fih*teste_U)

           recuperados = recuperados + undNot*(casos1-casos0)
           casos0 = casos1

           if(date == day_ini) exit

           if(day_ini == "2021-01-01") then
               recuperados = 0d0
               exit
           endif

       end do

    end subroutine

    subroutine taxas()

        real*8 :: a,b2,alpha2,c_bas,delta_bas

        if(flag .eqv. .False.) then

            incub = 6.4
            lethReal = 0.75

            b = 2d0/incub
            c = 2d0/incub
            d = 1d0/3.2
            f = 1d0/(1d0/d + 1d0/c)
            alpha = 0.06
            aA = alpha
            aU = alpha

            fih = 1d0*b/(b+f)

            alpha2 = 1.0
            a = -1.0 * alpha2 * alpha2
            b2 =  (1.0 + alpha2)

            c_bas = lethreal/(fih*lethNom)
            delta_bas = (b2*b2 - 4.0*a*c_bas)**(1d0/2d0)
            teste_U = -1.0*b2 + delta_bas
            teste_U = -teste_U / (2.0*a)
            teste_A = teste_U * alpha2

            if(teste_U .gt. 0.9) teste_U = 0.9
            if(teste_A .gt. 0.9) teste_A = 0.9

            undNot = (1.0-teste_A)*(1.0-fih*teste_U)
            undNot = 1.0*undNot/(teste_A + (1.0 - teste_A)*fih*teste_U)

            gu = abs(d*teste_U/(1d0-teste_U))
            ga = abs((f+c)*teste_A/(1d0-teste_A))

        else

            pmuni = beta_min*(0.9d0 + 0.2d0*ran2(iseed))
            incub1 = gammadist(12.8d0,4d0)
            incub2 = gammadist(12.8d0,4d0)
            incub3 = gammadist(12.8d0,4d0)

            lethReal = 0.5 + 0.5*ran2(iseed)

            b = 1d0/incub1
            c = 1d0/incub2
            d = 1d0/incub3

            f = 1d0/(1d0/d + 1d0/c)
            alpha = 0.06
            aA = alpha
            aU = alpha

            fih = 1d0*b/(b+f)

            alpha2 = 1.0
            a = -1.0 * alpha2 * alpha2
            b2 =  (1.0 + alpha2)

            c_bas = lethreal/(fih*lethNom)
            delta_bas = (b2*b2 - 4.0*a*c_bas)**(1d0/2d0)

            teste_U = -1.0*b2 + delta_bas
            teste_U = -teste_U / (2.0*a)
            teste_A = teste_U * alpha2
            
            if(teste_U .gt. 0.9) teste_U = 0.9
            if(teste_A .gt. 0.9) teste_A = 0.9

            undNot = (1.0-teste_A)*(1.0-fih*teste_U)
            undNot = 1.0*undNot/(teste_A + (1.0 - teste_A)*fih*teste_U)

            gu = abs(d*teste_U/(1d0-teste_U))
            ga = abs((f+c)*teste_A/(1d0-teste_A))

        endif


    end subroutine

    subroutine cond_ini() !CONDI��O INICIAL

        character*20 :: str
        integer :: ab1,ab3,ab4,ab5

        t = 0

        Ns = 0 ; Ne = 0 ; Ni = 0 ; Na = 0 ; Nr = 0 ; Nc = 0 ; Nrc = 0

        Nag = pop
        Ns = pop
        cont = 1
        z = 1d0

        Nr = pop
        Ns = Ns-Na-Ne-Nr-Ni-Nc
        call eigen_CI()

    end subroutine

    subroutine eigen_CI()

        real*8 :: M(3,3),vec(3),vec0(3),difvec(3),ndvec,profi
        integer :: i
        real*8 :: mu_A,beta_U,beta_R,alpha_C,alpha_R,beta_T,alpha_T,beta_C
        real*8 :: lambda_A,lambda_U,norm,norm0,maxi,dconf,tdel

        mu_A = 1d0*b
        beta_U = 1d0*c
        beta_R = 1d0*f
        beta_C = 1d0*ga
        alpha_C = 1d0*gu
        alpha_R = 1d0*d
        beta_T = beta_U+beta_R+beta_C
        alpha_T = alpha_C+alpha_R

        vec(1) = 0.4 ; vec(2) = 0.77 ; vec(3) = 0.991
        profi = pmuni

        lambda_A = aA*kcont*profi
        lambda_U = aU*kcont*profi

        M = 0

        M(1,1) = -1d0*mu_A  ; M(2,1) =  1d0*lambda_A  ; M(3,1) = 1d0*lambda_U
        M(1,2) =  1d0*mu_A  ; M(2,2) =  -1d0*beta_T   ; M(3,2) = 0
        M(1,3) =  0         ; M(2,3) =  1d0*beta_U    ; M(3,3) = -1d0*alpha_T

        maxi = min(M(1,1),M(2,2),M(3,3))

        M(1,1) = M(1,1)-1d0*maxi
        M(2,2) = M(2,2)-1d0*maxi
        M(3,3) = M(3,3)-1d0*maxi

        do i = 1,1000

            vec0 = vec
            vec = matmul(vec0,M)

            norm = sum(vec*vec)
            norm0= sum(vec0*vec0)

            norm = sqrt(norm)
            norm0 = sqrt(norm0)

            vec = 1d0*vec/norm
            vec0 = 1d0*vec0/norm0

            difvec = abs(vec-vec0)
            ndvec = sum(difvec*difvec)

        end do

        dconf = (conf1-conf0)
        tdel = 7d0

        Nc = int(1d0*conf0)
        Ni = int(1d0*dconf/(tdel)/(beta_C*vec(2)/vec(3)+alpha_C))
        Na = int(1d0*Ni*vec(2)/vec(3))
        Ne = int(1d0*Ni*vec(1)/vec(3))
        Nr = int(recuperados)

        Ns = pop-Nr-Ni-Na-Nc-Ne
        if(Ns .lt. 0) then
         
            Ns = 0
            Ni = 0 ; Na = 0 ; Ne = 0
            Nr = Ns-Nc
        end if
    end subroutine

    subroutine calcula_taxas()

        real*8 :: ale

        betaA = 1d0*aA*pmuni*kcont*pop/Ns
        betaU = 1d0*aU*pmuni*kcont*pop/Ns

        r0 =1d0*betaU/(c+f+ga)*(betaA/betaU+c/(d+gu))

    end subroutine

    subroutine minimiza_quadrados()
        Nconf(0) = conf0
        do j = 0,tmax
            dif(j) = (ref(j) - Nconf(j))*(ref(j) - Nconf(j))
          
        end do

        dift = sum(dif)

        if(dift .lt. dif_min) then
            minimo = Nconf(0:tmax)
            beta_min = pmuni
            dif_min = dift
        end if

    end subroutine

    subroutine atualiza_profilaxia()

        difer = dift0 - dift

        dift0 = dift

        if(difer .gt. -0.001)pmuni = pmuni+1d0*dp

        if(difer .lt. 0) then

            dp = dp*(-1d0)
            pmuni = pmuni+1d0*dp

        end if

        if(pmuni .lt. 0) pmuni=0

        difer0 = difer

    end subroutine

    subroutine rk4(par)

        real*8 :: FS1,FS2,FS3,FS4,FE1,FE2,FE3,FE4,FA1,FA2,FA3,FA4
        real*8 :: FU1,FU2,FU3,FU4,FC1,FC2,FC3,FC4,FR1,FR2,FR3,FR4
        integer :: par
        real*8, parameter :: dt = 0.1

        tser1 = 1d0

        do while (t .lt. par)

            FS1 = edo_S(1d0*Ns,1d0*Na,1d0*Ni,1d0*pop)
            FE1 = edo_E(1d0*Ns,1d0*Ne,1d0*Na,1d0*Ni,1d0*pop)
            FA1 = edo_A(1d0*Ne,1d0*Na)
            FU1 = edo_U(1d0*Na,1d0*Ni)
            FC1 = edo_C(1d0*Na,1d0*Ni)
            FR1 = edo_R(1d0*Na,1d0*Ni)

            FS2 = edo_S(1d0*Ns+FS1*dt/2d0,1d0*Na+FA1*dt/2d0,1d0*Ni+FU1*dt/2d0,1d0*pop)
            FE2 = edo_E(1d0*Ns+FS1*dt/2d0,1d0*Ne+FE1*dt/2d0,1d0*Na+FA1*dt/2d0,1d0*Ni+FU1*dt/2d0,1d0*pop)
            FA2 = edo_A(1d0*Ne+FE1*dt/2d0,1d0*Na+FA1*dt/2d0)
            FU2 = edo_U(1d0*Na+FA1*dt/2d0,1d0*Ni+FU1*dt/2d0)
            FC2 = edo_C(1d0*Na+FA1*dt/2d0,1d0*Ni+FU1*dt/2d0)
            FR2 = edo_R(1d0*Na+FA1*dt/2d0,1d0*Ni+FU1*dt/2d0)

            FS3 = edo_S(1d0*Ns+FS2*dt/2d0,1d0*Na+FA2*dt/2d0,1d0*Ni+FU2*dt/2d0,1d0*pop)
            FE3 = edo_E(1d0*Ns+FS2*dt/2d0,1d0*Ne+FE2*dt/2d0,1d0*Na+FA2*dt/2d0,1d0*Ni+FU2*dt/2d0,1d0*pop)
            FA3 = edo_A(1d0*Ne+FE2*dt/2d0,1d0*Na+FA2*dt/2d0)
            FU3 = edo_U(1d0*Na+FA2*dt/2d0,1d0*Ni+FU2*dt/2d0)
            FC3 = edo_C(1d0*Na+FA2*dt/2d0,1d0*Ni+FU2*dt/2d0)
            FR3 = edo_R(1d0*Na+FA2*dt/2d0,1d0*Ni+FU2*dt/2d0)

            FS4 = edo_S(1d0*Ns+FS3*dt,1d0*Na+FA3*dt,1d0*Ni+FU3*dt,1d0*pop)
            FE4 = edo_E(1d0*Ns+FS3*dt,1d0*Ne+FE3*dt,1d0*Na+FA3*dt,1d0*Ni+FU3*dt,1d0*pop)
            FA4 = edo_A(1d0*Ne+FE3*dt,1d0*Na+FA3*dt)
            FU4 = edo_U(1d0*Na+FA3*dt,1d0*Ni+FU3*dt)
            FC4 = edo_C(1d0*Na+FA3*dt,1d0*Ni+FU3*dt)
            FR4 = edo_R(1d0*Na+FA3*dt,1d0*Ni+FU3*dt)

            Ns = Ns+1d0*dt*(FS1+2d0*FS2+2d0*FS3+FS4)/6d0
            Ne = Ne+1d0*dt*(FE1+2d0*FE2+2d0*FE3+FE4)/6d0
            Na = Na+1d0*dt*(FA1+2d0*FA2+2d0*FA3+FA4)/6d0
            Ni = Ni+1d0*dt*(FU1+2d0*FU2+2d0*FU3+FU4)/6d0
            Nc = Nc+1d0*dt*(FC1+2d0*FC2+2d0*FC3+FC4)/6d0
            Nr = Nr+1d0*dt*(FR1+2d0*FR2+2d0*FR3+FR4)/6d0

            t = t+dt

            if(t .gt. tser1) then
                Nconf(tser1) = Nconf(tser1)+Nc

                if(flag .eqv. .true.) then
                    write(277,'(7F20.10)')real(tser1),Nc,ref(tser1),1d0*r0*(Ns/pop),1d0*Ns/pop,1d0*Nr/Nc,1d0*nga
                    write(278,'(8F20.10)')real(tser1),Ns,Ne,Na,Ni,Nc,Nr,1d0*nga
                end if

                tser1=tser1+1

            endif

        end do

    end subroutine

    function edo_S(dS,dA,dU,dS0)
        real*8:: dS,dA,dU,dS0
        real*8::edo_S
        edo_S = -1d0*betaU*dU*dS/dS0 -1d0*betaA*dA*dS/dS0
        return

    end function

    function edo_E(dS,dE,dA,dU,dS0)
        real*8:: dS,dE,dA,dU,dR,dC,dS0
        real*8:: edo_E
        edo_E = +1d0*betaU*dU*dS/dS0 +1d0*betaA*dA*dS/dS0-b*dE

        return

    end function

    function edo_A(dE,dA)
        real*8:: dS,dE,dA,dU,dR,dC,dS0
        real*8:: edo_A
        edo_A = 1d0*b*dE-1d0*ga*dA-1d0*c*dA-1d0*f*dA

        return

    end function

    function edo_U(dA,dU)
        real*8:: dS,dE,dA,dU,dR,dC,dS0
        real*8:: edo_U
        edo_U = 1d0*c*dA-1d0*gu*dU-1d0*d*dU

        return

    end function

    function edo_C(dA,dU)
        real*8:: dS,dE,dA,dU,dR,dC,dS0
        real*8:: edo_C
        edo_C = 1d0*gu*dU+1d0*ga*dA

        return

    end function

    function edo_R(dA,dU)

        real*8:: dS,dE,dA,dU,dR,dC,dS0
        real*8:: edo_R

        edo_R = 1d0*d*dU+1d0*f*dA

        return

    end function

    function ran2(idum)
        INTEGER :: idum
        REAL*8 ::ran2
        Integer, PARAMETER ::IM1=2147483563,IM2=2147483399,IMM1=IM1-1,&
             IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,&
             NTAB=32,NDIV=1+IMM1/NTAB
        real*8,parameter   ::EPS=1.2e-7,RNMX=1.-EPS,AM=1./IM1
        INTEGER            ::idum2,j,k,iv(NTAB),iy
        SAVE iv,iy,idum2
        DATA idum2/123456789/, iv/NTAB*0/, iy/0/
        if (idum.le.0) then
           idum=max(-idum,1)
           idum2=idum
           do j=NTAB+8,1,-1
              k=idum/IQ1
              idum=IA1*(idum-k*IQ1)-k*IR1
              if (idum.lt.0) idum=idum+IM1
              if (j.le.NTAB) iv(j)=idum
           enddo
           iy=iv(1)
        endif
        k=idum/IQ1
        idum=IA1*(idum-k*IQ1)-k*IR1
        if (idum.lt.0) idum=idum+IM1
        k=idum2/IQ2
        idum2=IA2*(idum2-k*IQ2)-k*IR2
        if (idum2.lt.0) idum2=idum2+IM2
        j=1+iy/NDIV
        iy=iv(j)-idum2
        iv(j)=idum
        if(iy.lt.1)iy=iy+IMM1
        ran2=min(AM*iy,RNMX)
        return
    end function ran2

    function gammadist(alpha,beta)

        real*8  gammadist
        real*8  x,alpha,beta,xm,varx
        real*8  xmax,beta_to_alp,alpha_m_1,g_alpha,fmax,fgamma

        xm=alpha/beta
        varx=alpha/beta**2
        xmax=10*(xm+sqrt(varx))

        !alguns parametros
        beta_to_alp=beta**alpha
        alpha_m_1=alpha-1
        g_alpha=gamma(alpha)
        x=(alpha-1)/beta !moda
        fmax=beta_to_alp*x**alpha_m_1*exp(-beta*x)/g_alpha

        11 x=xmax*ran2(iseed)
        fgamma=beta_to_alp*x**alpha_m_1*exp(-beta*x)/g_alpha
        if (fmax*ran2(iseed)>fgamma) goto 11
        gammadist = x

    end function gammadist


    end program
