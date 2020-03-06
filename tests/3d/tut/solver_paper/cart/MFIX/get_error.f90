
implicit real*8(a-h,o-z)
parameter(n=11016)
real*8 a(n),b(n),er(n)

open(1,file='fort.111')
open(2,file='../TRI/fort.211')
open(3,file='error1.txt')

do j=1,50000
   do i=1,11016
	read(1,*)a(i)
	read(2,*)b(i)
	er(i) = abs(a(i)-b(i))
   end do
   write(3,*) j, maxval(er),minval(er),maxval(a),maxval(b)
enddo

end
