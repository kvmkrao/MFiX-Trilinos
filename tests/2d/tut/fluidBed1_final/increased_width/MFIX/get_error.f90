
implicit real*8(a-h,o-z)
parameter(n=11016)
real*8 a(n),b(n),er(n)

open(1,file='fort.700')
open(2,file='../TRI/fort.700')
open(3,file='error1.txt')

do j=1,50000
   do i=1,980
	read(1,*) ii,a(i)
	read(2,*) jj,b(i)
	er(i) = abs(a(i)-b(i))
   end do
   write(3,*) j, maxval(er),minval(er),maxval(a),maxval(b)
enddo

end
