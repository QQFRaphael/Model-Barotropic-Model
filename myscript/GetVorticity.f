cc this program is used to calculate velocity potential
cc  or stream function from wind field 
cc  using spherical function expansion
cc procedure: 1.wind>>divergence and vorticty
cc            2.prepare base spherical function (Legendre polynominals
cc              and sine and cosine functions)
cc            3.expansion of divergence and vorticity 
cc            4.expansion coefficients for vel potential or stream func
cc            5.vel potential and stream function by composition
cc resolution:2.5x2.5
cc
                  parameter (ix=145,iy=73,it=480,nytt=40)
	parameter (lp=21,mp=21)
	dimension u(ix,iy,it),v(ix,iy,it)
	dimension um(ix,iy),vm(ix,iy)
	dimension div(ix,iy),vor(ix,iy)
	dimension pov(ix,iy),str(ix,iy)
cc
	dimension alm(0:lp,0:mp),blm(0:lp,0:mp)
	dimension Alf(0:lp,0:mp),Blf(0:lp,0:mp)
	dimension spx(ix,0:mp),cpx(ix,0:mp)
	dimension plm(iy,0:lp,0:mp),spy(iy),cpy(iy),pnlm(0:lp,0:mp)
cc
	r0=6.37E+6
	cc0=3.1415926/180.0
cc
cc prepare sine and cosine base functions
	do i=1,ix
	do m=0,mp
	spx(i,m)=sin(float(m*(i-1))*2.5*cc0)
	cpx(i,m)=cos(float(m*(i-1))*2.5*cc0)
	enddo
	enddo
cc prepare lengedre polynominals 
	do j=1,iy
	cpy(j)=cos(float(j-iy/2-1)*2.5*cc0)
	spy(j)=sin(float(j-iy/2-1)*2.5*cc0)
	enddo
	do l=1,lp
	do m=0,l
	pm2=1.0
	if(l.eq.m) then
	lww=1
	else
	lww=l-m
	endif
	do l2=1,lww
	pm2=pm2/sqrt(float(l2))
	enddo
	pm2=pm2/sqrt(float(2*l+1))
	do l1=1,l+m
	pm2=pm2*sqrt(float(l1))
	enddo
	pnlm(l,m)=pm2*sqrt(2.0)
	enddo
	enddo
cc
	do j=1,iy
c++++
ccc m=0
	plm(j,0,0)=1.0
	plm(j,1,0)=spy(j)
	do l=1,lp-1
	plm(j,l+1,0)=(2*l+1)*spy(j)*plm(j,l,0)-(l+0)*plm(j,l-1,0)
	plm(j,l+1,0)=plm(j,l+1,0)/float(l+1-0)
	enddo
ccc m=1
	plm(j,1,1)=cpy(j)
	do l=1,lp-1
	plm(j,l+1,1)=plm(j,l-1,1)+(2*l+1)*cpy(j)*plm(j,l,1-1)
	enddo
ccc
	do l=1,lp-1
	do m=2,l+1
	plm(j,l+1,m)=plm(j,l-1,m)+(2*l+1)*cpy(j)*plm(j,l,m-1)
	enddo
	enddo
c++++
	enddo
	print *,'finish polynominals'
cc
	open(11,file='GetVorticity.clim',
     .    form='unformatted',access='direct',recl=144*73)
	irecl=0
cc
	open(31,file='umon200.dat',
     .    form='unformatted',access='direct',recl=144*73)
cc
	open(32,file='vmon200.dat',
     .    form='unformatted',access='direct',recl=144*73)
cc
	irecp=0
	do k=1,it
cc+++++
	irecp=irecp+1
	read(31,rec=irecp) ((u(i,j,k),i=1,ix-1),j=1,iy)
        read(32,rec=irecp) ((v(i,j,k),i=1,ix-1),j=1,iy)
cc add one point
	do j=1,iy
	u(ix,j,k)=u(1,j,k)
	v(ix,j,k)=v(1,j,k)
	enddo
cc+++++
	enddo
cc
	do 100 mon=1,12
	print *,mon
cc+++++
cc climatology
	do i=1,ix
	do j=1,iy
	um(i,j)=0.0
	vm(i,j)=0.0
	do k=mon,it,12  !1948-2001
	um(i,j)=um(i,j)+u(i,j,k)
	vm(i,j)=vm(i,j)+v(i,j,k)
	enddo
	um(i,j)=um(i,j)/nytt !54.0
	vm(i,j)=vm(i,j)/nytt !54.0
	enddo
	enddo
cc calculate divergence and vorticity
cc div*r0,vor*r0
	do i=2,ix-1
	do j=2,iy-1
        div(i,j)=(um(i+1,j)-um(i-1,j))/(cpy(j)*5.0*cc0)
     .   +(vm(i,j+1)*cpy(j+1)-vm(i,j-1)*cpy(j-1))/(cpy(j)*5.0*cc0)
        vor(i,j)=(vm(i+1,j)-vm(i-1,j))/(cpy(j)*5.0*cc0)
     .   -(um(i,j+1)*cpy(j+1)-um(i,j-1)*cpy(j-1))/(cpy(j)*5.0*cc0)
	enddo
	enddo
	do j=2,iy-1
        div(1,j)=(um(2,j)-um(ix-1,j))/(cpy(j)*5.0*cc0)
     .   +(vm(1,j+1)*cpy(j+1)-vm(1,j-1)*cpy(j-1))/(cpy(j)*5.0*cc0)
        vor(1,j)=(vm(2,j)-vm(ix-1,j))/(cpy(j)*5.0*cc0)
     .   -(um(1,j+1)*cpy(j+1)-um(1,j-1)*cpy(j-1))/(cpy(j)*5.0*cc0)
        div(ix,j)=div(1,j)
        vor(ix,j)=vor(1,j)
	enddo
	do i=1,ix
        div(i,1)=0.0
        div(i,iy)=0.0
        vor(i,1)=0.0
        vor(i,iy)=0.0
	enddo
        print *,'finish vorticity/divergence'
cc expansion of div and vor according to spherical functions
cc to get expansion coef alm and blm
cc11 vorticity
	do l=1,lp
	do m=0,l
	alm(l,m)=0.0
	blm(l,m)=0.0
	do i=1,ix-1
	do j=1,iy
        alm(l,m)=alm(l,m)+r0*vor(i,j)*plm(j,l,m)*cpx(i,m)*cpy(j)
        blm(l,m)=blm(l,m)+r0*vor(i,j)*plm(j,l,m)*spx(i,m)*cpy(j)
	enddo
	enddo
	alm(l,m)=alm(l,m)/pnlm(l,m)/3.1415926
     .                   *2.5*2.5*cc0**2
	blm(l,m)=blm(l,m)/pnlm(l,m)/3.1415926
     .                   *2.5*2.5*cc0**2
	enddo
	alm(l,0)=alm(l,0)/2.0
	enddo
        print *,'finish alm,blm'
cc get Alf and Blf
	do l=1,lp
	do m=0,l
	Alf(l,m)=r0*alm(l,m)/pnlm(l,m)/float(l*(l+1))
	Blf(l,m)=r0*blm(l,m)/pnlm(l,m)/float(l*(l+1))
	enddo
	enddo
        print *,'finish Alf,Blf'
cc composition of Alf and Blf to get vel potential and stream func
	do i=1,ix-1
	do j=1,iy
        str(i,j)=0.0
	do l=1,lp
	do m=0,l
        str(i,j)=str(i,j)+(Alf(l,m)*cpx(i,m)+Blf(l,m)*spx(i,m))
     .                   *plm(j,l,m)
	enddo
	enddo
	enddo
	enddo
cc expansion of div and vor according to spherical functions
cc to get expansion coef alm and blm
cc22 divergence
	do l=1,lp
	do m=0,l
	alm(l,m)=0.0
	blm(l,m)=0.0
	do i=1,ix-1
	do j=1,iy
        alm(l,m)=alm(l,m)+r0*div(i,j)*plm(j,l,m)*cpx(i,m)*cpy(j)
        blm(l,m)=blm(l,m)+r0*div(i,j)*plm(j,l,m)*spx(i,m)*cpy(j)
	enddo
	enddo
	alm(l,m)=alm(l,m)/pnlm(l,m)/3.1415926
     .                   *2.5*2.5*cc0**2
	blm(l,m)=blm(l,m)/pnlm(l,m)/3.1415926
     .                   *2.5*2.5*cc0**2
	enddo
	alm(l,0)=alm(l,0)/2.0
	enddo
        print *,'finish alm,blm'
cc get Alf and Blf
	do l=1,lp
	do m=0,l
	Alf(l,m)=r0*alm(l,m)/pnlm(l,m)/float(l*(l+1))
	Blf(l,m)=r0*blm(l,m)/pnlm(l,m)/float(l*(l+1))
	enddo
	enddo
        print *,'finish Alf,Blf'
cc composition of Alf and Blf to get vel potential and stream func
	do i=1,ix-1
	do j=1,iy
        pov(i,j)=0.0
	do l=1,lp
	do m=0,l
        pov(i,j)=pov(i,j)+(Alf(l,m)*cpx(i,m)+Blf(l,m)*spx(i,m))
     .                   *plm(j,l,m)
	enddo
	enddo
	enddo
	enddo
cc pov/r0,str/r0
	do i=1,ix-1
	do j=1,iy
        pov(i,j)=pov(i,j)/r0
        str(i,j)=str(i,j)/r0
	enddo
	enddo
cc output
        irecl=irecl+1
        write(11,rec=irecl) ((div(i,j)/r0,i=1,ix-1),j=1,iy)
        irecl=irecl+1
        write(11,rec=irecl) ((vor(i,j)/r0,i=1,ix-1),j=1,iy)
        irecl=irecl+1
        write(11,rec=irecl) ((pov(i,j)*(-1.0),i=1,ix-1),j=1,iy)
        irecl=irecl+1
        write(11,rec=irecl) ((str(i,j)*(-1.0),i=1,ix-1),j=1,iy)
cc+++++
 100    continue
cc
	stop
	end


