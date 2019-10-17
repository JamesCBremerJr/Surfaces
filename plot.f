        implicit double precision (a-h,o-z)
        dimension zs1(2,10 000)
        dimension zs2(2,10 000)
        dimension zs3(2,10 000)
        dimension zs4(2,10 000)
c
        n1  = 120
c
        pi = acos(-1.0d0)
        dd = n1
        dd = 2*pi/dd
        do 1000 j=1,n1
        t = dd*(j-1)
        x = 2*cos(t)
        y = sin(t)
        zs1(1,j)=x
        zs1(2,j)=y
 1000 continue
c
        n2 = 200
        dd = n2
        dd = 2*pi/dd
        do 1100 j=1,n2
        t = dd*(j-1)
        x = 4+cos(t)
        y = 5*sin(t)
        zs2(1,j)=x
        zs2(2,j)=y
 1100 continue
c
        n3 = 77
        dd = n3
        dd = 2*pi/dd
        do 1200 j=1,n3
        t = dd*(j-1)
        x =-1+cos(t)
        y =-5*sin(t)
        zs3(1,j)=x
        zs3(2,j)=y
 1200 continue
c
        n4 = 177
        dd = n4
        dd = 2*pi/dd
        do 1300 j=1,n4
        t = dd*(j-1)
        x = 17*cos(t)
        y = -7+sin(t)
        zs4(1,j)=x
        zs4(2,j)=y
 1300 continue

        iplot-point=1
        call plot_points("Points on an ellipse.*",iplot,n1,zs1)
c
        iplot=2
        call plot_points2("Points on two ellipses*",iplot,n1,zs1,n2,zs2)
c
        iplot=3
        call plot_points3("Points on three ellipses*",iplot,n1,zs1,
     -    n2,zs2,n3,zs3)
c
        iplot=4
        call plot_points4("Points on four ellipses*",iplot,n1,zs1,
     -    n2,zs2,n3,zs3,n4,zs4)

c
c$$$        iout=3
c$$$        n = 5
c$$$c
c$$$        do 3000 i=1,n
c$$$c
c$$$        zs3(1,i) = i
c$$$        zs3(2,i) = i
c$$$c
c$$$        zs3(3,i) = i+1
c$$$        zs3(4,i) = i+1
c$$$ 3000 continue
c$$$c
c$$$        call plot_tikz_rectangles(iout,n,zs3)

        end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       This file contains code for generating primitive plots.  The
c       following subroutines produce GNUPLOT files:
c
c   plot_points - plot a collection of points in the plane
c
c   plot_points2 - plot two separate collections of points in the plane
c
c   plot_points3 - plot three collections of points in the plane
c
c   plot_points4 - plot four collections of points in the plane
c
c   plot_points_on_curve - plot a curve with a collection of points
c       lying on it
c
c       The following subroutines produce TIKZ output suitable for
c       producing publications:
c
c
c   plot_tikz_rectangles - produce a plot displying a disjoint 
c       collection of rectangles in the plane
c
c
c   plot3d_points - plot a collection of points in space
c
c   plot_triangles - construct a postscript file which draws a 
c       collection of triangles in the *plane*
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


        subroutine plot_triangles(iout,n,verts)
        implicit double precision (a-h,o-z)
        dimension verts(6,1)
        character*11 filename
c
 1000 format ("tris",I3.3,".eps")
 1100 format (".01 setlinewidth")
 1200 format ("newpath")
 1300 format ("stroke")
c 1400 format (E22.15," ",E22.15," scale")
c 1500 format (E22.15," ",E22.15," translate")
c
 2000 format (E22.15," ",E22.15," moveto")
 2100 format (E22.15," ",E22.15," lineto")
c
c       Determine the scale of the plot.
c
        xmin = 1d50
        xmax =-1d50
        ymin = 1d50
        ymax =-1d50
c
        do 4000 i=1,n
        x1 = verts(1,i)
        y1 = verts(2,i)
        x2 = verts(3,i)
        y2 = verts(4,i)
        x3 = verts(5,i)
        y3 = verts(6,i)
        ymax = max(ymax,y1)
        ymax = max(ymax,y2)
        ymax = max(ymax,y3)
c
        xmax = max(xmax,x1)
        xmax = max(xmax,x2)
        xmax = max(xmax,x3)
c
        ymin = min(ymin,y1)
        ymin = min(ymin,y2)
        ymin = min(ymin,y3)
c
        xmin = min(xmin,y1)
        xmin = min(xmin,y2)
        xmin = min(xmin,y3)

 4000 continue
c
        ddx = 7.5*72.0d0/(xmax-xmin)
        ddy = 7.5*72.0d0/(ymax-ymin)
c
        write (filename,1000) iout
c
        iw = 20
        open(iw,FILE=filename)

c
        do 3000 i=1,n
        x1 = (verts(1,i)-xmin)*ddx+36
        y1 = (verts(2,i)-ymin)*ddy+36*3
        x2 = (verts(3,i)-xmin)*ddx+36
        y2 = (verts(4,i)-ymin)*ddy+36*3
        x3 = (verts(5,i)-xmin)*ddx+36
        y3 = (verts(6,i)-ymin)*ddy+36*3
c
        write (iw,1200)
        write(iw,1100)
        write (iw,2000) x1,y1
        write (iw,2100) x2,y2
        write (iw,2100) x3,y3
        write (iw,2100) x1,y1
        write (iw,1300)
 3000 continue
c
        close(iw)
c
        end


        subroutine plot_tikz_rectangles(iout,n,zs)
        implicit double precision (a-h,o-z)
        dimension zs(4,n)
        character*11 filename
c
c       Produce a TIKZ file which generates a picture of a collection
c       of rectanges in the plane.
c
c                             Input Parameters:
c
c   iout - the name of the output file will be tikz${iout}.tex
c   n - the number of rectangles
c   zs - a (4,n) array each column of which gives the coordinates of
c       a triangle.
c
c
 0100 format("tikz",I3.3,".tex")
 0200 format("\draw [fill=gray!50] (",E12.5,",",E12.5,") rectangle (",
     -       E12.5,",",E12.5, ");")
 0300 format("\draw (",E12.5,",",E12.5,") rectangle (",
     -       E12.5,",",E12.5, ");")

c
        iw = 20
        write(filename,0100) iout
        open(iw,FILE=filename)
c
        xmin = 1d50
        xmax =-1d50
        ymin = 1d50
        ymax =-1d50
c
        do 1000 i=1,n
        x1 = zs(1,i)
        y1 = zs(2,i)
        x2 = zs(3,i)
        y2 = zs(4,i)
c
        xmin = min(xmin,x1)
        xmin = min(xmin,x2)
        xmax = max(xmax,x1)
        xmax = max(xmax,x2)
        ymin = min(ymin,y1)
        ymin = min(ymin,y2)
        ymax = max(ymax,y1)
        ymax = max(ymax,y2)
c
 1000 continue
c
        write(iw,*) "\begin{tikzpicture}[scale=1]"

        write(iw,0300) xmin,ymin,xmax,ymax

        do 1100 i=1,n
        write(iw,0200) zs(1,i),zs(2,i),zs(3,i),zs(4,i)
 1100 continue
c
        write(iw,*) "\end{tikzpicture}"

        end


        subroutine plot_points(title,iplot,n,zs)
        implicit double precision (a-h,o-z)
        dimension zs(2,1)
        character*1 title(1),AST
        character*5 filename1
        character*9 filename2
        character*80 command
        character*1 quote
c
        data quote /'"'/, AST / '*' /
c
c       Produce a plot showing a user-specified collection of points
c       in the plane. 
c
c                         Input Parameters:
c
c   title - an asterisk terminated string to use as a title
c   iplot - the filename of the GNU output file will be gn${iplot}
c   n - number of points
c   zs - a (2,n) array whose columns given the coordinates of the points
c       to plot
c
c
 1000 format("gn",I3.3)
 1100 format("gn",I3.3,".dat")
 1200 format((E12.5,2X,E12.5))
c
 1300 format(" plot [",E12.5,":",E12.5,"] [",E12.5,":",E12.5," ] ",
     1   '"',A9,'"'," with points pointtype 7 pointsize .25")
 1400 format("gnuplot -persist ",A9)
 1500 format(" set title ",'"',80A1,'"')
c
c       Find the length of the title.
c
        len = 0
        do 2000 j=1,10 000
        if (title(j) .eq. AST) goto 2100
        len=len+1
 2000 continue
 2100 continue
c
c       Open the data file, gn???.dat, and the GNUplot script file,
c       gn???.
c
        write(filename1,1000) iplot
        write(filename2,1100) iplot
c
        iw1=20
        iw2=21
c
        open(iw1,FILE=filename1)
        open(iw2,FILE=filename2)
c
c       Write the coordinates of the points to the data file and 
c       close it.
c
        write(iw2,1200) (zs(1,j),zs(2,j),j=1,n)
        close(iw2)
c
c       Decide on the length of the axes.
c
        xmin = 1d50
        xmax =-1d50
        ymin = 1d50
        ymax =-1d50
c
        do 2200 j=1,n
        xmax = max(zs(1,j),xmax)
        xmin = min(zs(1,j),xmin)        
        ymax = max(zs(2,j),ymax)
        ymin = min(zs(2,j),ymin)
 2200 continue
c
        dlen = max((xmax-xmin),(ymax-ymin))
        dlen = dlen*1.1d0
c
        xmin = (xmax+xmin)/2-dlen/2
        xmax = xmin+dlen
        ymin = (ymax+ymin)/2-dlen/2
        ymax = ymin+dlen
c
        write(iw1,*)    "set terminal wxt"
        write(iw1,*)    "unset key"
c
        if (len .gt. 0) then
        write (iw1,1500) (title(j),j=1,len)
        endif
c
        write(iw1,1300) xmin,xmax,ymin,ymax,filename2
c
        close(iw1)
c
        end


        subroutine plot_points2(title,iplot,n1,zs1,n2,zs2)
        implicit double precision (a-h,o-z)
        dimension zs1(2,1),zs2(2,1)
c
        character*1 title(1),AST
        character*5 filename1
c
        character*9 filename2
        character*10 filename3
c
        character*80 command
        character*1 quote
c
        data quote /'"'/, AST / '*' /
c
c       Produce a plot showing  two user-specified collection of points
c       in the plane. 
c
c                         Input Parameters:
c
c   title - an asterisk terminated string to use as a title
c   iplot - the filename of the GNU output file will be gn${iplot}
c   n? - number of points in the ?th collection
c   zs? - a (2,n) array whose columns given the coordinates of the points
c       in the ?th collection
c
 1600 format("set style line 1 pointtype 7 ps .25 lc rgb ",'"',
     -     "red0",'"')
 1700 format("set style line 2 pointtype 7 ps .25 lc rgb ",'"',
     -    "blue0",'"')
c
 1000 format("gn",I3.3)
 1100 format("gn",I3.3,".dat")
 1150 format("gn",I3.3,".dat2")
c
 1200 format((E12.5,2X,E12.5))
c
 1300 format(" plot [",E12.5,":",E12.5,"] [",E12.5,":",E12.5," ] ",
     -   '"',A9,'"'," with points ls 1,",
     -   '"',A10,'"'," with points ls 2")

 1400 format("gnuplot -persist ",A9)
 1500 format(" set title ",'"',80A1,'"')
c
c       Find the length of the title.
c
        len = 0
        do 2000 j=1,10 000
        if (title(j) .eq. AST) goto 2100
        len=len+1
 2000 continue
 2100 continue
c
c       Open the data files, gn???.dat, and the GNUplot script file,
c       gn???.
c
        write(filename1,1000) iplot
        write(filename2,1100) iplot
        write(filename3,1150) iplot
c
        iw1=20
        iw2=21
        iw3=22
c
        open(iw1,FILE=filename1)
        open(iw2,FILE=filename2)
        open(iw3,FILE=filename3)
c
c       Write the coordinates of the points to the data files and 
c       close them.
c
        write(iw2,1200) (zs1(1,j),zs1(2,j),j=1,n1)
        write(iw3,1200) (zs2(1,j),zs2(2,j),j=1,n2)
        close(iw2)
        close(iw3)
c
c       Decide on the length of the axes.
c
        xmin = 1d50
        xmax =-1d50
        ymin = 1d50
        ymax =-1d50
c
        do 2200 j=1,n1
        xmax = max(zs1(1,j),xmax)
        xmin = min(zs1(1,j),xmin)        
        ymax = max(zs1(2,j),ymax)
        ymin = min(zs1(2,j),ymin)
 2200 continue
c
        do 2250 j=1,n2
        xmax = max(zs2(1,j),xmax)
        xmin = min(zs2(1,j),xmin)        
        ymax = max(zs2(2,j),ymax)
        ymin = min(zs2(2,j),ymin)
 2250 continue
c
        dlen = max((xmax-xmin),(ymax-ymin))
        dlen = dlen*1.1d0
c
        xmin = (xmax+xmin)/2-dlen/2
        xmax = xmin+dlen
        ymin = (ymax+ymin)/2-dlen/2
        ymax = ymin+dlen
c
        write(iw1,*)    "set terminal wxt"
        write(iw1,*)    "unset key"
c
        if (len .gt. 0) then
        write (iw1,1500) (title(j),j=1,len)
        endif
c
        write(iw1,1600)
        write(iw1,1700)
        write(iw1,1300) xmin,xmax,ymin,ymax,filename2,filename3
        close(iw1)
c
        end



        subroutine plot_points3(title,iplot,n1,zs1,n2,zs2,n3,zs3)
        implicit double precision (a-h,o-z)
        dimension zs1(2,1),zs2(2,1),zs3(2,1)
c
        character*1 title(1),AST
        character*5 filename1
c
        character*9 filename2
        character*10 filename3
        character*10 filename4
c
        character*80 command
        character*1 quote
c
        data quote /'"'/, AST / '*' /
c
c       Produce a plot showing  two user-specified collection of points
c       in the plane. 
c
c                         Input Parameters:
c
c   title - an asterisk terminated string to use as a title
c   iplot - the filename of the GNU output file will be gn${iplot}
c   n? - number of points in the ?th collection
c   zs? - a (2,n) array whose columns given the coordinates of the points
c       in the ?th collection
c
 1600 format("set style line 1 pointtype 7 ps .25 lc rgb ",'"',
     -    "red0",'"')
 1700 format("set style line 2 pointtype 7 ps .25 lc rgb ",'"',
     -    "blue0",'"')
 1800 format("set style line 3 pointtype 7 ps .25 lc rgb ",'"',
     -    "gray0",'"')
c
 1000 format("gn",I3.3)
 1100 format("gn",I3.3,".dat")
 1150 format("gn",I3.3,".dat2")
 1175 format("gn",I3.3,".dat3")
c
 1200 format((E12.5,2X,E12.5))
c
 1300 format(" plot [",E12.5,":",E12.5,"] [",E12.5,":",E12.5," ] ",
     -   '"',A9,'"',"  with points ls 1, ",
     -   '"',A10,'"'," with points ls 2, ",
     -   '"',A10,'"'," with points ls 3")

 1400 format("gnuplot -persist ",A9)
 1500 format(" set title ",'"',80A1,'"')
c
c       Find the length of the title.
c
        len = 0
        do 2000 j=1,10 000
        if (title(j) .eq. AST) goto 2100
        len=len+1
 2000 continue
 2100 continue
c
c       Open the data files, gn???.dat, and the GNUplot script file,
c       gn???.
c
        write(filename1,1000) iplot
        write(filename2,1100) iplot
        write(filename3,1150) iplot
        write(filename4,1175) iplot
c
        iw1=20
        iw2=21
        iw3=22
        iw4=23
c
        open(iw1,FILE=filename1)
        open(iw2,FILE=filename2)
        open(iw3,FILE=filename3)
        open(iw4,FILE=filename4)
c
c       Write the coordinates of the points to the data files and 
c       close them.
c
        write(iw2,1200) (zs1(1,j),zs1(2,j),j=1,n1)
        write(iw3,1200) (zs2(1,j),zs2(2,j),j=1,n2)
        write(iw4,1200) (zs3(1,j),zs3(2,j),j=1,n3)
c
        close(iw2)
        close(iw3)
        close(iw4)
c
c       Decide on the length of the axes.
c
        xmin = 1d50
        xmax =-1d50
        ymin = 1d50
        ymax =-1d50
c
        do 2200 j=1,n1
        xmax = max(zs1(1,j),xmax)
        xmin = min(zs1(1,j),xmin)        
        ymax = max(zs1(2,j),ymax)
        ymin = min(zs1(2,j),ymin)
 2200 continue
c
        do 2250 j=1,n2
        xmax = max(zs2(1,j),xmax)
        xmin = min(zs2(1,j),xmin)        
        ymax = max(zs2(2,j),ymax)
        ymin = min(zs2(2,j),ymin)
 2250 continue
c
        do 2275 j=1,n3
        xmax = max(zs3(1,j),xmax)
        xmin = min(zs3(1,j),xmin)        
        ymax = max(zs3(2,j),ymax)
        ymin = min(zs3(2,j),ymin)
 2275 continue
c
c
        dlen = max((xmax-xmin),(ymax-ymin))
        dlen = dlen*1.1d0
c
        xmin = (xmax+xmin)/2-dlen/2
        xmax = xmin+dlen
        ymin = (ymax+ymin)/2-dlen/2
        ymax = ymin+dlen
c
        write(iw1,*)    "set terminal wxt"
        write(iw1,*)    "unset key"
c
        if (len .gt. 0) then
        write (iw1,1500) (title(j),j=1,len)
        endif
c
        write(iw1,1600)
        write(iw1,1700)        
        write(iw1,1800)

        write(iw1,1300) xmin,xmax,ymin,ymax,filename2,filename3,
     -    filename4
        close(iw1)
c
        end



        subroutine plot_points4(title,iplot,n1,zs1,n2,zs2,n3,zs3,
     -    n4,zs4)
        implicit double precision (a-h,o-z)
        dimension zs1(2,1),zs2(2,1),zs3(2,1),zs4(2,1)
c
        character*1 title(1),AST
        character*5 filename1
c
        character*9 filename2
        character*10 filename3
        character*10 filename4
        character*10 filename5

c
        character*80 command
        character*1 quote
c
        data quote /'"'/, AST / '*' /
c
c       Produce a plot showing  two user-specified collection of points
c       in the plane. 
c
c                         Input Parameters:
c
c   title - an asterisk terminated string to use as a title
c   iplot - the filename of the GNU output file will be gn${iplot}
c   n? - number of points in the ?th collection
c   zs? - a (2,n) array whose columns given the coordinates of the points
c       in the ?th collection
c
 1600 format("set style line 1 pointtype 7 ps .5 lc rgb ",'"',
     -    "red0",'"')
 1700 format("set style line 2 pointtype 7 ps .5 lc rgb ",'"',
     -    "blue0",'"')
 1800 format("set style line 3 pointtype 7 ps .5 lc rgb ",'"',
     -    "gray0",'"')
 1900 format("set style line 4 pointtype 7 ps .5 lc rgb ",'"',
     -    "green0",'"')

c
 1000 format("gn",I3.3)
 1100 format("gn",I3.3,".dat")
 1150 format("gn",I3.3,".dat2")
 1175 format("gn",I3.3,".dat3")
 1190 format("gn",I3.3,".dat4")
c
 1200 format((E12.5,2X,E12.5))
c
 1300 format(" plot [",E12.5,":",E12.5,"] [",E12.5,":",E12.5," ] ",
     -   '"',A9,'"',"  with points ls 1,",
     -   '"',A10,'"'," with points ls 2,",
     -   '"',A10,'"'," with points ls 3,",
     -   '"',A10,'"'," with points ls 4")

 1400 format("gnuplot -persist ",A9)
 1500 format(" set title ",'"',80A1,'"')
c
c       Find the length of the title.
c
        len = 0
        do 2000 j=1,10 000
        if (title(j) .eq. AST) goto 2100
        len=len+1
 2000 continue
 2100 continue
c
c       Open the data files, gn???.dat, and the GNUplot script file,
c       gn???.
c
        write(filename1,1000) iplot
        write(filename2,1100) iplot
        write(filename3,1150) iplot
        write(filename4,1175) iplot
        write(filename5,1190) iplot

c
        iw1=20
        iw2=21
        iw3=22
        iw4=23
        iw5=24
c
        open(iw1,FILE=filename1)
        open(iw2,FILE=filename2)
        open(iw3,FILE=filename3)
        open(iw4,FILE=filename4)
        open(iw5,FILE=filename5)
c
c       Write the coordinates of the points to the data files and 
c       close them.
c
        write(iw2,1200) (zs1(1,j),zs1(2,j),j=1,n1)
        write(iw3,1200) (zs2(1,j),zs2(2,j),j=1,n2)
        write(iw4,1200) (zs3(1,j),zs3(2,j),j=1,n3)
        write(iw5,1200) (zs4(1,j),zs4(2,j),j=1,n4)
c
        close(iw2)
        close(iw3)
        close(iw4)
        close(iw5)

c
c       Decide on the length of the axes.
c
        xmin = 1d50
        xmax =-1d50
        ymin = 1d50
        ymax =-1d50
c
        do 2200 j=1,n1
        xmax = max(zs1(1,j),xmax)
        xmin = min(zs1(1,j),xmin)        
        ymax = max(zs1(2,j),ymax)
        ymin = min(zs1(2,j),ymin)
 2200 continue
c
        do 2250 j=1,n2
        xmax = max(zs2(1,j),xmax)
        xmin = min(zs2(1,j),xmin)        
        ymax = max(zs2(2,j),ymax)
        ymin = min(zs2(2,j),ymin)
 2250 continue
c
        do 2275 j=1,n3
        xmax = max(zs3(1,j),xmax)
        xmin = min(zs3(1,j),xmin)        
        ymax = max(zs3(2,j),ymax)
        ymin = min(zs3(2,j),ymin)
 2275 continue
c
        do 2290 j=1,n4
        xmax = max(zs4(1,j),xmax)
        xmin = min(zs4(1,j),xmin)        
        ymax = max(zs4(2,j),ymax)
        ymin = min(zs4(2,j),ymin)
 2290   continue
c
        dlen = max((xmax-xmin),(ymax-ymin))
        dlen = dlen*1.1d0
c
        xmin = (xmax+xmin)/2-dlen/2
        xmax = xmin+dlen
        ymin = (ymax+ymin)/2-dlen/2
        ymax = ymin+dlen
c
        write(iw1,*)    "set terminal wxt"
        write(iw1,*)    "unset key"
c
        if (len .gt. 0) then
        write (iw1,1500) (title(j),j=1,len)
        endif
c
        write(iw1,1600)
        write(iw1,1700)        
        write(iw1,1800)
        write(iw1,1900)

        write(iw1,1300) xmin,xmax,ymin,ymax,filename2,filename3,
     -    filename4,filename5
        close(iw1)
c
        end



        subroutine plot_points_on_curve(title,iplot,n,zs,m,zscurve)
        implicit double precision (a-h,o-z)
        dimension zs(2,1),zscurve(2,1)
        character*1 title(1),AST
        character*5 filename1
        character*9 filename2
        character*10 filename3
        character*80 command
        character*1 quote
c
c       Plot a collection of points lying on a curve.
c
        data quote /'"'/, AST / '*' /
c
c
c
c                         Input Parameters:
c
c                         Output Parameters:
c
c
 1000 format("gn",I3.3)
 1100 format("gn",I3.3,".dat")
 1150 format("gn",I3.3,".dat2")
 1200 format((E12.5,2X,E12.5))
c
 1300 format(" plot [",E12.5,":",E12.5,"] [",E12.5,":",E12.5," ] ",
     -   '"',A10,'"'," with lines ls 1,",
     -   '"',A9,'"'," with points ls 2")
 1400 format("gnuplot -persist ",A9)
 1500 format(" set title ",'"',80A1,'"')
c
 1600 format("set style line 1 lt 1 lw 1 lc rgb ",'"',"gray90",'"')
 1700 format("set style line 2 pointtype 7 ps .5 lc rgb ",'"',
     -    "gray0",'"')
c
c
c       Find the length of the title.
c
        len = 0
        do 2000 j=1,10 000
        if (title(j) .eq. AST) goto 2100
        len=len+1
 2000 continue
 2100 continue
c
c       Open the data file, gn???.dat, and the GNUplot script file,
c       gn???.
c
        write(filename1,1000) iplot
        write(filename2,1100) iplot
        write(filename3,1150) iplot
c     
        iw1=20
        iw2=21
        iw3=22
c
        open(iw1,FILE=filename1)
        open(iw2,FILE=filename2)
        open(iw3,FILE=filename3)
c
c       Write the coordinates of the points to the data files.
c
        write(iw2,1200) (zs(1,j),zs(2,j),j=1,n)
        close(iw2)
c
        write(iw3,1200) (zscurve(1,j),zscurve(2,j),j=1,m)
        close(iw3)
c
c       Decide on the length of the axes.
c
        xmin = 1d50
        xmax =-1d50
        ymin = 1d50
        ymax =-1d50
c
        do 2200 j=1,n
        xmax = max(zs(1,j),xmax)
        xmin = min(zs(1,j),xmin)
        ymax = max(zs(2,j),ymax)
        ymin = min(zs(2,j),ymin)
 2200 continue
c
        do 2250 j=1,m
        xmax = max(zscurve(1,j),xmax)
        xmin = min(zscurve(1,j),xmin)
        ymax = max(zscurve(2,j),ymax)
        ymin = min(zscurve(2,j),ymin)
 2250 continue
c
        dlen = max((xmax-xmin),(ymax-ymin))
        dlen = dlen*1.1d0
c
        xmin = (xmax+xmin)/2-dlen/2
        xmax = xmin+dlen
        ymin = (ymax+ymin)/2-dlen/2
        ymax = ymin+dlen
c
        write(iw1,*)    "set terminal wxt"
        write(iw1,*)    "unset key"
c
        if (len .gt. 0) then
        write (iw1,1500) (title(j),j=1,len)
        endif
c
        write(iw1,1600) 
        write(iw1,1700) 
        write(iw1,1300) xmin,xmax,ymin,ymax,filename3,filename2
c
        close(iw1)
c
        end



        subroutine plot3d_points(title,iplot,n,zs)
        implicit double precision (a-h,o-z)
        dimension zs(3,n)
c
        character*1 title(1),AST
        character*5 filename1
        character*9 filename2
c
        character*1 quote
        data quote /'"'/, AST / '*' /
c
 1000 format("gn",I3.3)
 1100 format("gn",I3.3,".dat")
 1200 format(E12.5,2X,E12.5,2X,E12.5)
c 1300 format(" splot [",E12.5,":",E12.5,"] [",E12.5,":",E12.5," ] ",
c     1   '"',A9,'"'," with points pointtype 7 pointsize .25")
 1300 format(" splot ",
     1   '"',A9,'"'," with points pointtype 7 pointsize .25")

 1400 format("gnuplot -persist ",A9)
 1500 format(" set title ",'"',80A1,'"')
c
c       Find the length of the title.
c
        len = 0
        do 2000 j=1,10 000
        if (title(j) .eq. AST) goto 2100
        len=len+1
 2000 continue
 2100 continue
c
c       Open the data file, gn???.dat, and the GNUplot script file,
c       gn???.
c
        write(filename1,1000) iplot
        write(filename2,1100) iplot
c
        iw1=20
        iw2=21
c
        open(iw1,FILE=filename1)
        open(iw2,FILE=filename2)
c
c       Write the coordinates of the points to the data file and 
c       close it.
c
        do 0100 j=1,n
        write(iw2,1200) zs(1,j),zs(2,j),zs(3,j)
 0100 continue
c
        close(iw2)
c
        write(iw1,*)    "set terminal wxt"
        write(iw1,*)    "unset key"
c
        if (len .gt. 0) then
        write (iw1,1500) (title(j),j=1,len)
        endif
c
        write(iw1,1300) filename2
c
        close(iw1)
        
        end


        



        subroutine plot3d_points2(title,iplot,n1,zs1,n2,zs2)
        implicit double precision (a-h,o-z)
        dimension zs1(3,1),zs2(3,1)
c
        character*1  title(1),AST
        character*5  filename1
        character*9  filename2
        character*10 filename3
c
        character*1 quote
        data quote /'"'/, AST / '*' /
c
 1000 format("gn",I3.3)
 1100 format("gn",I3.3,".dat")
 1150 format("gn",I3.3,".dat2")

 1200 format(E12.5,2X,E12.5,2X,E12.5)
c 1300 format(" splot [",E12.5,":",E12.5,"] [",E12.5,":",E12.5," ] ",
c     1   '"',A9,'"'," with points pointtype 7 pointsize .25")
c 1300 format(" splot ",
c     1   '"',A9,'"'," with points pointtype 7 pointsize .25")
c
 1300 format(" splot ",
     -   '"',A9,'"'," with points pointtype 7 pointsize .25,",
     -   '"',A10,'"'," with points pointtype 7 pointsize .25")

 1400 format("gnuplot -persist ",A9)
 1500 format(" set title ",'"',80A1,'"')
c
c       Find the length of the title.
c
        len = 0
        do 2000 j=1,10 000
        if (title(j) .eq. AST) goto 2100
        len=len+1
 2000 continue
 2100 continue
c
c       Open the data file, gn???.dat, and the GNUplot script file,
c       gn???.
c
        write(filename1,1000) iplot
        write(filename2,1100) iplot
        write(filename3,1150) iplot
c
        iw1=20
        iw2=21
        iw3=22
c
        open(iw1,FILE=filename1)
        open(iw2,FILE=filename2)
        open(iw3,FILE=filename3)
c
c       Write the coordinates of the points to the data file and 
c       close it.
c
        do 0100 j=1,n1
        write(iw2,1200) zs1(1,j),zs1(2,j),zs1(3,j)
 0100 continue
c
        do 0200 j=1,n2
        write(iw3,1200) zs2(1,j),zs2(2,j),zs2(3,j)
 0200 continue
c
        close(iw2)
        close(iw3)
c
        write(iw1,*)    "set terminal wxt"
        write(iw1,*)    "unset key"
c
        if (len .gt. 0) then
        write (iw1,1500) (title(j),j=1,len)
        endif
c
        write(iw1,1300) filename2,filename3
c
        close(iw1)
c
        end


        


        subroutine plot3d_points3(title,iplot,n1,zs1,n2,zs2,n3,zs3)
        implicit double precision (a-h,o-z)
        dimension zs1(3,1),zs2(3,1),zs3(3,1)
c
        character*1  title(1),AST
        character*5  filename1
        character*9  filename2
        character*10 filename3
        character*10 filename4
c
        character*1 quote
        data quote /'"'/, AST / '*' /
c
 1000 format("gn",I3.3)
 1100 format("gn",I3.3,".dat")
 1150 format("gn",I3.3,".dat2")
 1170 format("gn",I3.3,".dat3")


 1200 format(E12.5,2X,E12.5,2X,E12.5)
c 1300 format(" splot [",E12.5,":",E12.5,"] [",E12.5,":",E12.5," ] ",
c     1   '"',A9,'"'," with points pointtype 7 pointsize .25")
c 1300 format(" splot ",
c     1   '"',A9,'"'," with points pointtype 7 pointsize .25")
c
 1300 format(" splot ",
     -   '"',A9,'"',"  with points pointtype 7 pointsize .25,",
     -   '"',A10,'"'," with points pointtype 7 pointsize .25,",
     -   '"',A10,'"'," with points pointtype 7 pointsize .25")
c
 1400 format("gnuplot -persist ",A9)
 1500 format(" set title ",'"',80A1,'"')
c
c       Find the length of the title.
c
        len = 0
        do 2000 j=1,10 000
        if (title(j) .eq. AST) goto 2100
        len=len+1
 2000 continue
 2100 continue
c
c       Open the data file, gn???.dat, and the GNUplot script file,
c       gn???.
c
        write(filename1,1000) iplot
        write(filename2,1100) iplot
        write(filename3,1150) iplot
        write(filename4,1170) iplot

c
        iw1=20
        iw2=21
        iw3=22
        iw4=23
c
        open(iw1,FILE=filename1)
        open(iw2,FILE=filename2)
        open(iw3,FILE=filename3)
        open(iw4,FILE=filename4)

c
c       Write the coordinates of the points to the data file and 
c       close it.
c
        do 0100 j=1,n1
        write(iw2,1200) zs1(1,j),zs1(2,j),zs1(3,j)
 0100 continue
c
        do 0200 j=1,n2
        write(iw3,1200) zs2(1,j),zs2(2,j),zs2(3,j)
 0200 continue
c
        do 0300 j=1,n3
        write(iw4,1200) zs3(1,j),zs3(2,j),zs3(3,j)
 0300 continue
c
        close(iw2)
        close(iw3)
c
        write(iw1,*)    "set terminal wxt"
        write(iw1,*)    "unset key"
c
        if (len .gt. 0) then
        write (iw1,1500) (title(j),j=1,len)
        endif
c
        write(iw1,1300) filename2,filename3,filename4
c
        close(iw1)
c
        end


        
