        program main
        parameter(is=41,js=41)
        implicit double precision (a-h,o-z)
        complex*16 z1(is,js)
        complex*16 d1(is,js)
c       z1(1,2)=1.190797402D-4
c       z1(1,2)=4.0D-6+z1(1,2)
 	read(88) z1
 	print *,z1(1,2),z1(1,12)
 	print *,z1(41,2),z1(41,12)
 	z1(1,2)=z1(1,2)+1.190797402D-4
        write(2)z1
        write(2)d1

        write(2)z1
        write(2)d1
        stop
        end
