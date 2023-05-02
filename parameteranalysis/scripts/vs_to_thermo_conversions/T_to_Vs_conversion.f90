! Convert S-wave velocity to mantle temperature using parameterisation of Yamauchi & Takei 2016.
! Read in parameters (z in km, vs in km/s).
! Input file "input.zVs", contains 'readmax' entries.


!--------------------------------------------------
! Anelasticity model parameters and switches
!--------------------------------------------------
module commondat
	implicit none
	real, public :: mu0,dmudT,dmudP,eta0,E,Va,solgrad,sol50
	integer, public :: switchdens,switchsol,switchgs,ns
end module commondat
!--------------------------------------------------
! Main functions
!--------------------------------------------------
module transdat
	implicit none
	real, public :: zT,VsT
end module transdat
!--------------------------------------------------
! Yamauchi & Takei 2016 intrinsic parameters
!--------------------------------------------------
module extradat
	implicit none
	real, public :: Ab,alpha,tauP,Teta,beta,delphi,gamma,lambdaphi
    data Ab,alpha,tauP,Teta,beta,delphi,gamma,lambdaphi / 0.664,0.38,6.e-5,0.94,0.,0.,5.,0./
end module extradat
!--------------------------------------------------
! Grose & Afonso 2013 expansivity
!--------------------------------------------------
module afonsodat
	implicit none
	real, public :: K0,KT,grun,p0,a0,a1
    data K0,KT,grun,p0,a0,a1 / 130.e9,4.8,6.,3330.,2.832e-5,0.758e-8 /
end module afonsodat
!--------------------------------------------------
! Reference mantle values
!--------------------------------------------------
module refdat
	implicit none
	real, public :: R,Pr,TKr,TK0,d,dr,rhor,alphaT,bmod,freq,pi,mn,T50
    data R,Pr,TKr,TK0,d,dr,rhor,alphaT,bmod,freq,pi,mn,T50 &
	/ 8.3145,1.5e9,1473.,273.,1.e-3,1.e-3,3291.,3.59e-5,115.2,0.01,3.1415927,3.,1326./
end module refdat
!--------------------------------------------------
! McKenzie & Bickle 1988 solidus calculation parameters
!--------------------------------------------------
module meltdat
	implicit none
	real, public :: melt,ifail,Vs0,dVsdT
end module meltdat

!--------------------------------------------------
! Main program
!--------------------------------------------------

program VsT
	use commondat
	integer, parameter :: readmax = 10000000 
	integer, parameter :: npar=11
	real :: depsol
	real, dimension(npar) :: pin
	real, dimension(readmax) :: z,T,Vs,rhoo,Qo,etao
	real, dimension(3,readmax) :: gin,Ts
	integer :: o, g, h, j, i, k
	integer :: narg, cptArg !#of arg & counter of arg
	character(30)  ::  name
	character(100)  ::  gsf = "average_profile_dannberg_2017.txt"


	narg=command_argument_count()
	!Loop over the arguments
	if(narg.gt.0)then
! 	!loop across options
	  do cptArg=1,narg
	    call get_command_argument(cptArg,name)
	    read(name,*)pin(cptArg)
	  end do
	else
	  write(*,*)"no input parameters given"
	  stop
	end if 

	mu0=pin(1)
	dmudT=pin(2)
	dmudP=pin(3)
	eta0=pin(4)
	E=pin(5)
	Va=pin(6)
	solgrad=pin(7)
    sol50=pin(8)
	switchdens=int(pin(9))
	switchsol=int(pin(10))
	switchgs=int(pin(11))

	
	!Load variable grain size profile if required
	if (switchgs.eq.1) then
	  h=1
	  do while(h.le.readmax)
	    open(54, file=gsf)
	    read(54,*,end=777) gin(1,h),gin(2,h),gin(3,h)
	    h=h+1
	  enddo
777	  continue
	  g=h-1
	  close(54)
	endif
	
	T0=0.
	
	open(1,file='input.zT')
	j=1
	do while(j.le.readmax)
	  read (1,*,end=888) z(j), T(j)
	  call Vs_calc(z(j),Vs(j),T(j),Qo(j),etao(j),rhoo(j))
	  j=j+1
	end do
888	continue
	k=j-1
	
	open(2,file='output_T_to_Vs.zTVsQvd')	
	i=1
	do while (i.le.k)
		write (2,*) z(i), T(i), Vs(i), Qo(i), etao(i), rhoo(i)
		i=i+1
	end do
	close(1)
	close(2)
	
! -------------------------------------------------
contains
!--------------------------------------------------

!--------------------------------------------------
! Brent to iteratively minimise misfit between calculated and predicted dV
!--------------------------------------------------
  function brent(xmin,xin)

	use afonsodat
	  
    real :: brent,xin,xmin
    integer, parameter :: itmax=500
    real, parameter :: CGOLD=.3819660
    real, parameter :: ZEPS=1.0e-10
    real, parameter :: AX=1.
    real, parameter :: BX=1.5
    real, parameter :: CX=2.0
    real, parameter :: tol=1.48e-8

!     Given a function f, and given a bracketing triplet of abscissas AX, bx, cx (such that bx is
!     between AX and cx, and f(bx) is less than both f(AX) and f(CX)), this routine isolates
!     the minimum to a fractional precision of about tol using Brent’s method. The abscissa of
!     the minimum is returned as xmin, and the minimum function value is returned as brent,
!     the returned function value.
!     Parameters: MAXimum allowed number of iterations; golden ratio; and a small number that
!     protects against trying to achieve fractional accuracy for a minimum that happens to be
!     exactly zero.
    integer :: iter
    real :: a,b,dd,e,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm
    a=min(AX,CX)        !a and b must be in ascending order, though the input
    b=max(AX,CX)    !abscissas need not be.
    v=BX          !Initializations...
    w=v
    x=v
    e=0.    !This will be the distance moved on the step before last.
    fx=abs((K0*(3./2.)*(x**(7./3.)-x**(5./3.))*(1.+(((3./4.)*(KT-4.))*(x**(2./3.)-1.))))-xin)
    fv=fx
    fw=fx
    do iter=1,itmax   !Main program loop.
!       write(*,*) 'iteration', iter
      xm=0.5*(a+b)
      tol1=(tol*abs(x))+ZEPS
      tol2=2.*tol1
      if(abs(x-xm).le.(tol2-0.5*(b-a))) goto 3    !Test for done here.
      if(abs(e).gt.tol1) then        !Construct a trial parabolic fit.
        r=(x-w)*(fx-fv)
        q=(x-v)*(fx-fw)
        p=(x-v)*q-(x-w)*r
        q=2.*(q-r)
        if(q.gt.0.) p=-p
        q=abs(q)
        etemp=e
        e=dd
	if(abs(p).ge.abs(0.5*q*etemp).or.p.le.q*(a-x).or.p.ge.q*(b-x)) goto 1
	!The above conditions determine the acceptability of the parabolic fit. Here it is o.k.:
	dd=p/q !Take the parabolic step.
	u=x+dd
	if(u-a.lt.tol2 .or. b-u.lt.tol2) dd=sign(tol1,xm-x)
	goto 2     !Skip over the golden section step.
      endif
1     if(x.ge.xm) then    !We arrive here for a golden section step, which we take
        e=a-x     !into the larger of the two segments.
      else
        e=b-x
      endif
      dd=CGOLD*e           !Take the golden section step.
2     if(abs(dd).ge.tol1) then     !Arrive here with dd computed either from parabolic fit, or
        u=x+dd              !else from golden section.
      else
        u=x+sign(tol1,dd)
      endif
      fu=abs((K0*(3./2.)*(u**(7./3.)-u**(5./3.))*(1.+(((3./4.)*(KT-4.))*(u**(2./3.)-1.))))-xin) !This is the one function evaluation per iteration,
      if(fu.le.fx) then    !and now we have to decide what to do with our function
        if(u.ge.x) then     !evaluation. Housekeeping follows:
          a=x
        else
          b=x
        endif
        v=w
        fv=fw
        w=x
        fw=fx
        x=u
        fx=fu
      else
        if(u.lt.x) then
          a=u
        else
          b=u
        endif
        if(fu.le.fw .or. w.eq.x) then
	  v=w
	  fv=fw
	  w=u
	  fw=fu
	else if(fu.le.fv .or. v.eq.x .or. v.eq.w) then
	  v=u
	  fv=fu
	endif
      endif      !Done with housekeeping. Back for another iteration.
    enddo 
    write(*,*) 'brent exceed maximum iterations'
3   xmin=x    !Arrive here ready to exit with best values.
    brent=fx
    return
    end function brent  
      
!--------------------------------------------------
! For Vs calculation
!--------------------------------------------------
subroutine Vs_calc(z,Vs,T,Q,eta,rho)
	!     use melt, only: vs0,dvsodt
	use commondat
	use extradat
	use refdat
	use afonsodat
	
	real :: TK, Pf, Pg, Tsol, sigmap, Ap, test, val, exponent
	real :: rho, eta, Q, Aeta, Tn, tauM, tau,tauS, Ju, J1, J2
	real :: intalphaT,rhoP0,alphaP0,dV0,misfit

	Pg=(z/30.)
	Pf=Pg*1.e9
	TK=T+273.
	
	if (switchsol.eq.0) then
	  Tsol=sol50+(solgrad*(z-50.))
	else if (switchsol.eq.1) then
	  call tsolid(z,Tsol)
	end if 
	
	Tn=TK/(Tsol+273.)
	if (Tn.lt.Teta) then
	  Aeta=1.
	else if (Tn.ge.Teta.and.Tn.lt.1.) then
	  Aeta=exp((-1.*((Tn-Teta)/(Tn-(Tn*Teta))))*log(gamma))
	else
	  Aeta=(1./gamma)*exp(-delphi)
	endif
	
	if (switchgs.eq.1) then
	  call gsize(gs, z)
	  d=gs
	endif
	
	eta=((((d/dr)**mn)*(eta0*exp((E/R)*(1./TK-1./TKr))*&
	exp((Va/R)*(Pf/TK-Pr/TKr))))*Aeta)
	Ju=1./(1.e9*(mu0+(dmudP*Pg)+(dmudT*T)))

	if (T.le.0..or.eta.gt.1.e40) then
		eta=1.e40
	end if
	tauM=eta*Ju
	tau=(3.*z*1000.)/4200.
	tauS=tau/(2.*pi*tauM)

	if (Tn.lt.0.91) then
	  Ap=0.01
	else if (Tn.ge.0.91.and.Tn.lt.0.96) then
	  Ap=0.01+(0.4*(Tn-0.91))
	else if (Tn.ge.0.96.and.Tn.lt.1.) then
	  Ap=0.03
	else
	  Ap=0.03+beta
	endif 
	if (Tn.lt.0.92) then
	  sigmap=4.
	else if (Tn.ge.0.92.and.Tn.lt.1.) then
	  sigmap=4.+(37.5*(Tn-0.92))
	else
	  sigmap=7.
	endif  
	J1=Ju*(1.+((Ab*(tauS**alpha))/alpha)+((sqrt(2.*pi)/2.)*Ap*sigmap&
	*(1.-erf((log(tauP/tauS))/(sqrt(2.)*sigmap)))))
	
	if (switchdens.eq.0) then
		rho=rhor*((1.-alphaT*(T-600.))+(Pg/bmod))
	else if (switchdens.eq.1) then
	    dV0=(-7.334963115431564676e-23*(Pf**2))+(7.510867653681621105e-12*Pf)+1.000184023114681908e+00  
	    alphaP0=dV0*exp((grun+1)*((dV0**(-1.))-1.))
	    rhoP0=p0*dV0
	    intalphaT=(a0*(TK-273.))+((a1/2.)*((TK**2.)-(273.**2.)))
	    rho=rhoP0*(1.-(alphaP0*intalphaT))
	endif
    
	Vs=1./(sqrt(rho*J1)*1000.)
	J2=(Ju*(pi/2.)*((Ab*(tauS**alpha))+(Ap*exp(-1.*(((log(tauP/tauS))**2.)&
	/(2.*sigmap**2.))))))+(Ju*tauS)
	Q=J2/J1 
	return
end subroutine Vs_calc

!--------------------------------------------------
! For variable grain size calculations
!--------------------------------------------------

  subroutine gsize(g,depg)
      implicit none
      real,intent(out) :: g
      real, intent(in)  :: depg
      integer :: depint
      depint=nint(depg)
      g=gin(2,depint+1)

    return
  end subroutine gsize
  
!--------------------------------------------------
! Calculating McKenzie & Bickle 1988 profiles
!--------------------------------------------------

  subroutine tsolid(z,Tsol)
    implicit none
    real, intent(in) :: z 
    real, intent(out) :: Tsol 
    real :: P, T1, Pd, dPdT, dT, err, dTs, Ts
    integer :: N
    real, parameter :: T0=1373.1
    real, parameter :: Bs=1.2e-2! Convert S-wave velocity to mantle temperature using parameterisation of Yamauchi & Takei 2016.
    real, parameter :: As=4.968e-4
    real, parameter :: Cs=136.05
    P=z/30.
    T1=(Cs*P)+T0
    N=0
35  Pd=(T1-T0)/Cs+As*exp(Bs*(T1-T0))-P
    N=N+1
    dPdT=1./Cs+As*Bs*exp(Bs*(T1-T0))
    dT=-Pd/dPdT
    T1=T1+dT
    err=abs(dT/T1)
    if (abs(dT).gt.0.1.and.N.lt.10) then
	goto 35
    end if
    Ts=T1
    dPdT=1./Cs+As*Bs*exp(Bs*(Ts-T0))
    dTs=1./dPdT
    Tsol=Ts-273.
    return
  end subroutine tsolid
  

end program VsT


  