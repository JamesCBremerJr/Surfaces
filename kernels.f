cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       Main evaluation routine.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


        subroutine eval_kernel(ikernel,zk,x,y,val)
        implicit double precision (a-h,o-z)
        dimension x(13),y(13)
        double complex zk,val
c
c       Evaluate one of the kernels associated with Laplace's equation
c       or the Helmholtz equation at wavenumer zk.
c
        if (ikernel .eq. 1) then
        call ksingle0(x,y,val)
        return
        endif
c
        if (ikernel .eq. 2) then
        call kdouble0(x,y,val)
        return
        endif
c
        if (ikernel .eq. 3) then
        call ksingle0prime(x,y,val)
        return
        endif
c
        if (ikernel .eq. 4) then
        call ksingle(x,y,zk,val)
        return
        endif
c
        if (ikernel .eq. 5) then
        call kdouble(x,y,zk,val)
        return
        endif
c
        if (ikernel .eq. 6) then
        call ksingleprime(x,y,zk,val)
        return
        endif
c
        end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       Kernels for Laplace's equation.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


        subroutine ksingle0(x,y,val)
        implicit double precision (a-h,o-z)
        dimension x(13),y(13)
c
        data pi      / 3.14159265358979323846264338327950288d0 /
        data over4pi / 0.0795774715459476678844418816862571810d0 /
c
        double complex val
c
        xwht  =  x(1)
        x1    =  x(2)
        x2    =  x(3)
        x3    =  x(4)
c        dx11  =  x(5)
c        dx21  =  x(6)
c        dx31  =  x(7)
c        dx12  =  x(8)
c        dx22  =  x(9)
c        dx32  =  x(10)
c        dnx1  =  x(11)
c        dnx2  =  x(12)
c        dnx3  =  x(13)
c
        ywht  =  y(1)
        y1    =  y(2)
        y2    =  y(3)
        y3    =  y(4)
c        dy11  =  y(5)
c        dy21  =  y(6)
c        dy31  =  y(7)
c        dy12  =  y(8)
c        dy22  =  y(9)
c        dy32  =  y(10)
c        dny1  =  y(11)
c        dny2  =  y(12)
c        dny3  =  y(13)
c
        dd    = sqrt((x1-y1)**2+(x2-y2)**2+(x3-y3)**2)
c
        val   = -xwht*ywht*over4pi/dd
c
        end


        subroutine kdouble0(x,y,val)
        implicit double precision (a-h,o-z)
        dimension x(13),y(13)
c
        data pi      / 3.14159265358979323846264338327950288d0 /
        data over4pi / 0.0795774715459476678844418816862571810d0 /
c
        double complex val
c
        xwht  =  x(1)
        x1    =  x(2)
        x2    =  x(3)
        x3    =  x(4)
c        dx11  =  x(5)
c        dx21  =  x(6)
c        dx31  =  x(7)
c        dx12  =  x(8)
c        dx22  =  x(9)
c        dx32  =  x(10)
c        dnx1  =  x(11)
c        dnx2  =  x(12)
c        dnx3  =  x(13)
c
        ywht  =  y(1)
        y1    =  y(2)
        y2    =  y(3)
        y3    =  y(4)
c        dy11  =  y(5)
c        dy21  =  y(6)
c        dy31  =  y(7)
c        dy12  =  y(8)
c        dy22  =  y(9)
c        dy32  =  y(10)
        dny1  =  y(11)
        dny2  =  y(12)
        dny3  =  y(13)
c
        dd    = sqrt((x1-y1)**2+(x2-y2)**2+(x3-y3)**2)
        val   = ((x1-y1)*dny1 + (x2-y2)*dny2 + (x3-y3)*dny3)/dd**3
        val   = -val*xwht*ywht*over4pi
c
        end



        subroutine ksingle0prime(x,y,val)
        implicit double precision (a-h,o-z)
        dimension x(13),y(13)
c
        data pi      / 3.14159265358979323846264338327950288d0 /
        data over4pi / 0.0795774715459476678844418816862571810d0 /
c
        double complex val
c
        xwht  =  x(1)
        x1    =  x(2)
        x2    =  x(3)
        x3    =  x(4)
c        dx11  =  x(5)
c        dx21  =  x(6)
c        dx31  =  x(7)
c        dx12  =  x(8)
c        dx22  =  x(9)
c        dx32  =  x(10)
        dnx1  =  x(11)
        dnx2  =  x(12)
        dnx3  =  x(13)
c
        ywht  =  y(1)
        y1    =  y(2)
        y2    =  y(3)
        y3    =  y(4)
c        dy11  =  y(5)
c        dy21  =  y(6)
c        dy31  =  y(7)
c        dy12  =  y(8)
c        dy22  =  y(9)
c        dy32  =  y(10)
c        dny1  =  y(11)
c        dny2  =  y(12)
c        dny3  =  y(13)
c
        dd    = sqrt((x1-y1)**2+(x2-y2)**2+(x3-y3)**2)
        val   = ((y1-x1)*dnx1 + (y2-x2)*dnx2 + (y3-x3)*dnx3)/dd**3
        val   = -val*xwht*ywht*over4pi
c
        end




cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       Kernels for Helmholtz equation
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


        
        subroutine ksingle(x,y,zk,val)
        implicit double precision (a-h,o-z)
        double complex zk,val,ima
c
        data ima     / (0.0d0,1.0d0) /
        data pi      / 3.14159265358979323846264338327950288d0 /
        data over4pi / 0.0795774715459476678844418816862571810d0 /
c
        dimension x(13),y(13)
c
        xwht  =  x(1)
        x1    =  x(2)
        x2    =  x(3)
        x3    =  x(4)
c        dx11  =  x(5)
c        dx21  =  x(6)
c        dx31  =  x(7)
c        dx12  =  x(8)
c        dx22  =  x(9)
c        dx32  =  x(10)
c        dnx1  =  x(11)
c        dnx2  =  x(12)
c        dnx3  =  x(13)
c
        ywht  =  y(1)
        y1    =  y(2)
        y2    =  y(3)
        y3    =  y(4)
c        dy11  =  y(5)
c        dy21  =  y(6)
c        dy31  =  y(7)
c        dy12  =  y(8)
c        dy22  =  y(9)
c        dy32  =  y(10)
        dny1  =  y(11)
        dny2  =  y(12)
        dny3  =  y(13)
c
        dd    = sqrt((x1-y1)**2+(x2-y2)**2+(x3-y3)**2)
c
        val   = exp(ima*zk*dd)/dd
        val   = -val*xwht*ywht*over4pi
c                
        end


        subroutine kdouble(x,y,zk,val)
        implicit double precision (a-h,o-z)
        dimension x(13),y(13)
        double complex val,ima,zk
c
        data ima     / (0.0d0,1.0d0) /
        data pi      / 3.14159265358979323846264338327950288d0 /
        data over4pi / 0.0795774715459476678844418816862571810d0 /
c
        xwht  =  x(1)
        x1    =  x(2)
        x2    =  x(3)
        x3    =  x(4)
c        dx11  =  x(5)
c        dx21  =  x(6)
c        dx31  =  x(7)
c        dx12  =  x(8)
c        dx22  =  x(9)
c        dx32  =  x(10)
c        dnx1  =  x(11)
c        dnx2  =  x(12)
c        dnx3  =  x(13)
c
        ywht  =  y(1)
        y1    =  y(2)
        y2    =  y(3)
        y3    =  y(4)
c        dy11  =  y(5)
c        dy21  =  y(6)
c        dy31  =  y(7)
c        dy12  =  y(8)
c        dy22  =  y(9)
c        dy32  =  y(10)
        dny1  =  y(11)
        dny2  =  y(12)
        dny3  =  y(13)
c
        dd    = sqrt((x1-y1)**2+(x2-y2)**2+(x3-y3)**2)
        val   = ((x1-y1)*dny1 + (x2-y2)*dny2 + (x3-y3)*dny3)/dd**3
        val   = val * exp(ima*zk*dd)*(1-ima*zk*dd)
        val   = -val*xwht*ywht*over4pi
c
        end




        subroutine ksingleprime(x,y,zk,val)
        implicit double precision (a-h,o-z)
        dimension x(13),y(13)
        double complex val,ima,zk
c
        data ima     / (0.0d0,1.0d0) /
        data pi      / 3.14159265358979323846264338327950288d0 /
        data over4pi / 0.0795774715459476678844418816862571810d0 /
c
        xwht  =  x(1)
        x1    =  x(2)
        x2    =  x(3)
        x3    =  x(4)
c        dx11  =  x(5)
c        dx21  =  x(6)
c        dx31  =  x(7)
c        dx12  =  x(8)
c        dx22  =  x(9)
c        dx32  =  x(10)
        dnx1  =  x(11)
        dnx2  =  x(12)
        dnx3  =  x(13)
c
        ywht  =  y(1)
        y1    =  y(2)
        y2    =  y(3)
        y3    =  y(4)
c        dy11  =  y(5)
c        dy21  =  y(6)
c        dy31  =  y(7)
c        dy12  =  y(8)
c        dy22  =  y(9)
c        dy32  =  y(10)
c        dny1  =  y(11)
c        dny2  =  y(12)
c        dny3  =  y(13)
c
        dd    = sqrt((x1-y1)**2+(x2-y2)**2+(x3-y3)**2)
        val   = ((y1-x1)*dnx1 + (y2-x2)*dnx2 + (y3-x3)*dnx3)/dd**3
        val   = val * exp(ima*zk*dd)*(1-ima*zk*dd)
        val   = -val*xwht*ywht*over4pi
c
        end




        subroutine ksingletang1(x,y,zk,val)
        implicit double precision (a-h,o-z)
        dimension x(13),y(13)
        double complex val,ima,zk
c
        data ima     / (0.0d0,1.0d0) /
        data pi      / 3.14159265358979323846264338327950288d0 /
        data over4pi / 0.0795774715459476678844418816862571810d0 /
c
        xwht  =  x(1)
        x1    =  x(2)
        x2    =  x(3)
        x3    =  x(4)
        dx11  =  x(5)
        dx21  =  x(6)
        dx31  =  x(7)
c        dx12  =  x(8)
c        dx22  =  x(9)
c        dx32  =  x(10)
        dnx1  =  x(11)
        dnx2  =  x(12)
        dnx3  =  x(13)
c
        ywht  =  y(1)
        y1    =  y(2)
        y2    =  y(3)
        y3    =  y(4)
c        dy11  =  y(5)
c        dy21  =  y(6)
c        dy31  =  y(7)
c        dy12  =  y(8)
c        dy22  =  y(9)
c        dy32  =  y(10)
c        dny1  =  y(11)
c        dny2  =  y(12)
c        dny3  =  y(13)
c
        dd    = sqrt((x1-y1)**2+(x2-y2)**2+(x3-y3)**2)
        val   = ((y1-x1)*dx11 + (y2-x2)*dx21 + (y3-x3)*dx31)/dd**3
        val   = val * exp(ima*zk*dd)*(1-ima*zk*dd)
        val   = -val*xwht*ywht*over4pi
c
        end


        subroutine ksingletang2(x,y,zk,val)
        implicit double precision (a-h,o-z)
        dimension x(13),y(13)
        double complex val,ima,zk
c
        data ima     / (0.0d0,1.0d0) /
        data pi      / 3.14159265358979323846264338327950288d0 /
        data over4pi / 0.0795774715459476678844418816862571810d0 /
c
        xwht  =  x(1)
        x1    =  x(2)
        x2    =  x(3)
        x3    =  x(4)
c        dx11  =  x(5)
c        dx21  =  x(6)
c        dx31  =  x(7)
        dx12  =  x(8)
        dx22  =  x(9)
        dx32  =  x(10)
        dnx1  =  x(11)
        dnx2  =  x(12)
        dnx3  =  x(13)
c
        ywht  =  y(1)
        y1    =  y(2)
        y2    =  y(3)
        y3    =  y(4)
c        dy11  =  y(5)
c        dy21  =  y(6)
c        dy31  =  y(7)
c        dy12  =  y(8)
c        dy22  =  y(9)
c        dy32  =  y(10)
c        dny1  =  y(11)
c        dny2  =  y(12)
c        dny3  =  y(13)
c
        dd    = sqrt((x1-y1)**2+(x2-y2)**2+(x3-y3)**2)
        val   = ((y1-x1)*dx12 + (y2-x2)*dx22 + (y3-x3)*dx32)/dd**3
        val   = val * exp(ima*zk*dd)*(1-ima*zk*dd)
        val   = -val*xwht*ywht*over4pi
c
        end
