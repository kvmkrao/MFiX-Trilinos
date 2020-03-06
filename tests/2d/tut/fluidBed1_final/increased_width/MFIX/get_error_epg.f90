
implicit real*8(a-h,o-z)
parameter(n=11016)
real*8 a(n),b(n),er(n)
real*8 c(n),d(n),ern(n),tm(n),tt(n)
real*8 t1,t2,tmp,sum1,avg
real*8 l2norm,l2nort
integer i1,i2

open(1,file='fort.800')
open(2,file='../TRI/fort.800')
open(3,file='errorn.txt')

do j=1,10000
   l2norm = 0.0d0
   l2nort = 0.0d0
        do i=1,918
	read(1,*) tm(j), i1, ii, a(i)
	read(2,*) tt(j), i2, jj, b(i)
	er(i) = 0.0d0
        if(i.gt.205) then 
	  er(i) = abs(a(i)-b(i))
	  if(a(i).lt.1.0001) l2norm = l2norm + a(i)*a(i)
          if(b(i).lt.1.0001) l2nort = l2nort + b(i)*b(i)
	end if

	end do 

!	if((tm(j).gt.0.7213).and.(tm(j).lt.0.72140)) then
	if((tm(j).gt.0.2014).and.(tm(j).lt.0.2030)) then
!	   if(i1.eq.31) then 
	   if(i1.eq.26) then 
	     do k=1,918
	   	c(k) = a(k)
	     write(101,*) tm(j),i1,k,c(k)
	     end do 
	   write(*,*) tm(j),i1
!	   exit
	   goto 999
	   end if  
	end if 

!	if(tt(j).gt.0.72148.and.tt(j).lt.0.72150) then
	if(tt(j).gt.0.2020.and.tt(j).lt.0.2040) then
!	   if(i2.eq.14) then 
	   if(i2.eq.16) then 
	   do k=1,918
	       d(k) = b(k)
	       write(102,*) tt(j),i1,k,d(k)
	   end do 
	   write(*,*) tt(j), i2
!	   exit
	   end if 
	end if
     write(3,*) tm(j), tt(j), i1,i2, maxval(er),minval(er),l2norm,l2nort

enddo

999	continue

l2norm = 0.0d0
l2nort = 0.0d0

tmp  = 0.0d0
sum1 = 0.0d0
do i=1,918
        er(i)=0.0d0
	if(i.gt.205) then 
	tmp   = tmp + (d(i) - c(i)) * (d(i) -c(i))
	er(i) = abs(d(i)) - abs(c(i))
	sum1  = sum1 + abs(er(i))
	write(11,*) i, d(i),c(i)
	end if
end do 

avg = sum1/float(918)

sum1 = 0.0
do i=215,918
    sum1  = sum1 + (avg-er(i))**2
end do 

!write(*,*) "l2 norm", tmp,maxval(er),dsqrt(sum1/918.00)
write(*,*) "l2 norm", tmp,maxval(er),sqrt(tmp/918.00)

end
