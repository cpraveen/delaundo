      program main
!
!     advancing front method for simply connected convex domain
!
!     this code is written with an emphasis on clarity, rather than
!     performance.  it is not an efficient method for triangulation
!     by advancing front.
!
!     known bugs:
!
!     1.  since the code is written in single precision arithmetic,
!         the triangulation algorithm may fail when two points are
!         too close.
!
!     definition of variables
!
!       cell(n,i)      cell(1,i), cell(2,i) and cell(3,i) are the node numbers
!                        of the three nodes comprising cell number i
!
!       icell          number of current cells
!
!       imax           maximum number of cells
!
!       free(j)        list of nodes currently on face or interior, j=1,jfree
!
!       jfree          number of nodes currently in front plus interior
!
!       jmax           maximum number of nodes
!
!       kface          number of faces currently in front
!
!       kmax           maximum number of faces
!
!       hash(j,k)      hash table for sides on front
!                        hash(j,0) is number of entries for key "j"
!                        hash(j,k) is kth entry and equals index "k" for
!                                  this side in the heap table heap(.,k)      
!
!       heap(1,k)      node number of beginning of face number k in heap
!     
!       heap(2,k)      node number of end of face number k in heap
!
!       hmax           maximum number of entries per key in hash table
!
!       mmod           modulus used for key is mmod+1
!
!       np(j)          conversion table of node numbers used for plotting
!                        intermediate results
!
!       x(j), y(j)     cartesian components of node number j
!
!     assumptions
!
!       all triangles are ordered CCW
!
      integer    hmax
      parameter (imax=25000,hmax=10000,jmax=12500,kmax=50000,mmod=99)
      integer    cell(3,imax), free(jmax), hash(0:mmod,0:hmax), 
     >           heap(2,kmax), np(jmax) 
      real       x(jmax), y(jmax)

      ! initialize
      call initial(free,hash,hmax,heap,icell,imax,jfree,jmax,
     >             kface,mmod,x,y)
      print *,'How many iterations do you want ?'
      read  *, itermax
      iter = 0
      do while ( kface.gt.0 )
        iter = iter + 1 
        print *,'Starting iteration number ',iter 
        ! always select the first entry in heap (see subroutine distance)   
        j1 = heap(1,1)  ! select "from" node on face
        j2 = heap(2,1)  ! select "to" node on face
        ! find third point
        call third (cell,free,heap,icell,j1,j2,j3,jfree,jmax,
     >              kface,np,x,y)
        !  modify front
        call front(cell,free,hash,hmax,heap,icell,imax,j1,j2,j3,
     >             jfree,kface,kmax,mmod) 
        if ( iter.eq.itermax ) then
          print *,'Max iterations in advancing front. Stop.'
          call output(free,icell,jfree,jmax,kface,cell,heap,np,x,y)
          stop
        end if  
      end do
      call output(free,icell,jfree,jmax,kface,cell,heap,np,x,y)  ! write output
      stop
      end

      SUBROUTINE ADD(heap,j1,j2,kface,kmax)
!
!     add face j1-j2 to heap
!
      integer heap(2,*) 

      if ( kface.lt.kmax ) then
        kface = kface + 1
        heap(1,kface) = j1
        heap(2,kface) = j2
      else
        print *,'Maximum number of faces exceeded.  Stop.'
        stop
      end if
      return
      end

      SUBROUTINE CIRCLE(x1,x2,y1,y2,x3,y3,xc,yc,r2)
!
!     computes the center of the circumscribing circle and square of radius
!
      xc = ((y3-y2)*(x1**2-x2**2+y1**2-y2**2)-
     >      (y1-y2)*(x3**2-x2**2+y3**2-y2**2))/ 
     >     (2.0*((x1-x2)*(y3-y2)-(x3-x2)*(y1-y2)))
      if ( y3.ne.y2 ) then
        yc = (0.5*(x3**2-x2**2+y3**2-y2**2)-(x3-x2)*xc)/(y3-y2)
      else
        yc = (0.5*(x1**2-x2**2+y1**2-y2**2)-(x1-x2)*xc)/(y1-y2) 
      end if   
      r2 = (x1-xc)**2+(y1-yc)**2  
      return
      end
!
      SUBROUTINE DISTANCE(dist,hash,hmax,heap,jn,j3,kn,mmod) 
!
!     computes distance of node jn from node j3
!    
      integer  dist, hash(0:mmod,0:hmax), heap(2,*), hmax, kn(2)

      mkey = mmod + 1
      ! search for distance one assuming order jn-j3
      key = mod(jn,mkey)
      j   = 1
      do while ( j.le.hash(key,0) )
        k = hash(key,j)
        if ( (heap(2,k).eq.j3).and.(heap(1,k).eq.jn) ) then
          dist  = 1
          kn(1) = k
          return 
        else
          j = j+1
        end if
      end do
      ! search for distance one assuming order j3-jn
      key = mod(j3,mkey)
      j   = 1
      do while ( j.le.hash(key,0) )
        k = hash(key,j)
        if ( (heap(2,k).eq.jn).and.(heap(1,k).eq.j3) ) then
          dist  = 1
          kn(1) = k
          return 
        else
          j = j+1
        end if
      end do
      ! search for distance two assuming order jn-jm-j3
      key = mod(jn,mkey)
      do i = 1,hash(key,0)
        k  = hash(key,i)     ! face number
        if ( k.ne.1 )  then  ! the first face cannot be j1-j2
          jm = heap(2,k)       ! "to" node on face
          if ( heap(1,k).eq.jn ) then 
            nkey = mod(jm,mkey)  ! new key 
            j = 1
            do while ( j.le.hash(nkey,0) )
              m = hash(nkey,j)
              if ( (heap(2,m).eq.j3).and.(heap(1,m).eq.jm) ) then
                dist  = 2
                kn(1) = k
                kn(2) = m
                return 
              else
                j = j+1
              end if
            end do
          end if
        end if 
      end do
      ! search for distance two assuming order j3-jm-jn
      key = mod(j3,mkey)
      do i = 1,hash(key,0)
        k  = hash(key,i)     ! face number
        jm = heap(2,k)       ! "to" node on face
        if ( heap(1,k).eq.j3 ) then
          nkey = mod(jm,mkey)  ! new key 
          j = 1
          do while ( j.le.hash(nkey,0) )
            m = hash(nkey,j)
            if ( m.ne.1 ) then  ! the second face cannot be j1-j2
              if ( (heap(2,m).eq.jn).and.(heap(1,m).eq.jm) ) then
                dist  = 2
                kn(1) = m
                kn(2) = k
                return 
              else
                j = j+1
              end if
            else
              j = j + 1
            end if
          end do
        end if
      end do
      ! distance is greater than 2, thus set dist = 0
      dist = 0  
      return
      end

      SUBROUTINE FDELETE(d1,d2,k1,k2,hash,heap,hmax,kface,mmod)
!
!     deletes faces k = 1, k1(i) (i = 1,d1), and k2(i) (i=1,d2) from heap
!
      integer d1,d2,hash(0:mmod,0:hmax),heap(2,*),hmax,kd(6),k1(2),k2(2)

      ! add faces to be deleted from hash and heap to array kd()
      kd(1) = 1
      do k = 1,d1  ! this loop not executed if d1 = 0
        kd(k+1) = k1(k)
      end do
      do k = 1,d2  ! this loop not executed if d2 = 0
        kd(k+d1+1) = k2(k)
      end do
      n = 1 + d1 + d2  ! n is number of faces to be deleted
      ! reorder the face numbers in ascending order
      do j = 2,n
        k = kd(j)
        do i = j-1,1,-1
          if ( kd(i).le.k ) goto 10
          kd(i+1) = kd(i)
        end do
        i = 0  
10      kd(i+1) = k
      end do
      kd(n+1) = kface + 1  ! last face to "delete" is beyond end of list
      ! delete faces and compress heap
      k = 0
      do i = 1,n
        do j = kd(i)+1,kd(i+1)-1
          k = k + 1
          heap(1,k) = heap(1,j)
          heap(2,k) = heap(2,j)
        end do
      end do
      kface = k
      return
      end
      SUBROUTINE FRONT(cell,free,hash,hmax,heap,icell,imax,j1,j2,j3,
     >                 jfree,kface,kmax,mmod) 
!
!     modifies front and adds cells
!
      integer cell(3,*), d1, d2, free(*), hash(0:mmod,0:hmax), 
     >        heap(2,*), hmax, k1(2), k2(2)
      logical bndry, done


      mkey = mmod + 1
      ! determine if j3 is on boundary or interior
      bndry = .false.
      done  = .false. 
      k     = 1
      do while ( .not.done )
        if ( heap(1,k).eq.j3 ) then
          bndry = .true.
          done  = .true.
        else if ( k.lt.kface ) then 
          k = k + 1
        else
          done = .true.
        end if
      end do
      !
      ! point j3 is in interior
      !
      if ( .not.bndry ) then
        d1 = 0
        d2 = 0
        call fdelete(d1,d2,k1,k2,hash,heap,hmax,kface,mmod)  
        call add(heap,j1,j3,kface,kmax)
        call add(heap,j3,j2,kface,kmax)
        call rehash(hash,heap,hmax,kface,mmod)
        ! jfree and free(j) do not change
        cell(1,icell+1) = j1  ! add new cell
        cell(2,icell+1) = j2
        cell(3,icell+1) = j3
        icell = icell + 1
        return
      end if
      !  
      ! point j3 is on front
      !
      ! determine distance from j1 to j3 and j2 to j3
      call distance(d1,hash,hmax,heap,j1,j3,k1,mmod) 
      call distance(d2,hash,hmax,heap,j2,j3,k2,mmod) 
      if ( (d1.eq.1).and.(d2.eq.1) ) then
        ! case 1
        call fdelete(d1,d2,k1,k2,hash,heap,hmax,kface,mmod)
        call rehash(hash,heap,hmax,kface,mmod)
        call ndelete(free,hash,hmax,heap,jfree,j1,mmod)
        call ndelete(free,hash,hmax,heap,jfree,j2,mmod)
        call ndelete(free,hash,hmax,heap,jfree,j3,mmod)
        cell(1,icell+1) = j1  
        cell(2,icell+1) = j2  
        cell(3,icell+1) = j3  
        icell = icell + 1
      else if ( (d1.eq.1).and.(d2.eq.2) ) then
        ! case 2a
        j4 = heap(1,k2(2))
        call fdelete(d1,d2,k1,k2,hash,heap,hmax,kface,mmod)
        call rehash(hash,heap,hmax,kface,mmod)
        call ndelete(free,hash,hmax,heap,jfree,j1,mmod)
        call ndelete(free,hash,hmax,heap,jfree,j2,mmod)
        call ndelete(free,hash,hmax,heap,jfree,j3,mmod)
        call ndelete(free,hash,hmax,heap,jfree,j4,mmod)
        cell(1,icell+1) = j1
        cell(2,icell+1) = j2
        cell(3,icell+1) = j3
        cell(1,icell+2) = j2
        cell(2,icell+2) = j4
        cell(3,icell+2) = j3
        icell = icell + 2 
      else if ( (d1.eq.2).and.(d2.eq.1) ) then 
        ! case 2b
        j4 = heap(1,k1(1))
        call fdelete(d1,d2,k1,k2,hash,heap,hmax,kface,mmod)
        call rehash(hash,heap,hmax,kface,mmod)
        call ndelete(free,hash,hmax,heap,jfree,j1,mmod)
        call ndelete(free,hash,hmax,heap,jfree,j2,mmod)
        call ndelete(free,hash,hmax,heap,jfree,j3,mmod)
        call ndelete(free,hash,hmax,heap,jfree,j4,mmod)
        cell(1,icell+1) = j1
        cell(2,icell+1) = j2
        cell(3,icell+1) = j3
        cell(1,icell+2) = j1
        cell(2,icell+2) = j3
        cell(3,icell+2) = j4
        icell = icell + 2 
      else if ( (d1.eq.1).and.(d2.eq.0) ) then
        ! case 3a
        call fdelete(d1,d2,k1,k2,hash,heap,hmax,kface,mmod)
        call add(heap,j3,j2,kface,kmax)
        call rehash(hash,heap,hmax,kface,mmod)
        call ndelete(free,hash,hmax,heap,jfree,j1,mmod)
        cell(1,icell+1) = j1
        cell(2,icell+1) = j2
        cell(3,icell+1) = j3
        icell = icell + 1
      else if ( (d1.eq.0).and.(d2.eq.1) ) then
        ! case 3b
        call fdelete(d1,d2,k1,k2,hash,heap,hmax,kface,mmod)
        call add(heap,j1,j3,kface,kmax)
        call rehash(hash,heap,hmax,kface,mmod)
        call ndelete(free,hash,hmax,heap,jfree,j2,mmod)
        cell(1,icell+1) = j1
        cell(2,icell+1) = j2
        cell(3,icell+1) = j3
        icell = icell + 1
      else if ( (d1.eq.2).and.(d2.eq.2) ) then
        ! case 4
        j4 = heap(1,k1(1))
        j5 = heap(1,k2(2))
        call fdelete(d1,d2,k1,k2,hash,heap,hmax,kface,mmod)
        call rehash(hash,heap,hmax,kface,mmod)
        call ndelete(free,hash,hmax,heap,jfree,j1,mmod)
        call ndelete(free,hash,hmax,heap,jfree,j2,mmod)
        call ndelete(free,hash,hmax,heap,jfree,j3,mmod)
        call ndelete(free,hash,hmax,heap,jfree,j4,mmod)
        call ndelete(free,hash,hmax,heap,jfree,j5,mmod)
        cell(1,icell+1) = j1
        cell(2,icell+1) = j2
        cell(3,icell+1) = j3
        cell(1,icell+2) = j1
        cell(2,icell+2) = j3
        cell(3,icell+2) = j4
        cell(1,icell+3) = j2
        cell(2,icell+3) = j5
        cell(3,icell+3) = j3
        icell = icell + 3
      else if ( (d1.eq.0).and.(d2.eq.2) ) then
        ! case 5a
        j4 = heap(2,k2(1))
        call fdelete(d1,d2,k1,k2,hash,heap,hmax,kface,mmod)
        call add(heap,j1,j3,kface,kmax)
        call rehash(hash,heap,hmax,kface,mmod)
        call ndelete(free,hash,hmax,heap,jfree,j2,mmod)
        call ndelete(free,hash,hmax,heap,jfree,j4,mmod)
        cell(1,icell+1) = j1
        cell(2,icell+1) = j2
        cell(3,icell+1) = j3
        cell(1,icell+2) = j2
        cell(2,icell+2) = j4
        cell(3,icell+2) = j3
        icell = icell + 2
      else if ( (d1.eq.2).and.(d2.eq.0) ) then   
        ! case 5b
        j4 = heap(1,k1(1))
        call fdelete(d1,d2,k1,k2,hash,heap,hmax,kface,mmod)
        call add(heap,j3,j2,kface,kmax)
        call rehash(hash,heap,hmax,kface,mmod)
        call ndelete(free,hash,hmax,heap,jfree,j1,mmod)
        call ndelete(free,hash,hmax,heap,jfree,j4,mmod)
        cell(1,icell+1) = j1
        cell(2,icell+1) = j2
        cell(3,icell+1) = j3
        cell(1,icell+2) = j1
        cell(2,icell+2) = j3
        cell(3,icell+2) = j4
        icell = icell + 2
      else if ( (d1.eq.0).and.(d2.eq.0) ) then 
        ! case 6
        call fdelete(d1,d2,k1,k2,hash,heap,hmax,kface,mmod)
        call add(heap,j1,j3,kface,kmax)
        call add(heap,j3,j2,kface,kmax)
        call rehash(hash,heap,hmax,kface,mmod)
        cell(1,icell+1) = j1
        cell(2,icell+1) = j2
        cell(3,icell+1) = j3
        icell = icell + 1
      end if
      if ( icell.gt.imax ) then
        print *,'Maximum number of cells exceeded.  Stop.'
        stop
      end if
      return
      end

      SUBROUTINE INITIAL(free,hash,hmax,heap,icell,imax,jfree,jmax,
     >                   kface,mmod,x,y)
!
!     initializes arrays free, hash, heap, x and y
!      
      integer free(*), hash(0:mmod,0:hmax), heap(2,*), hmax 
      logical done, found
      real    x(*), y(*)

      icell = 0
      mkey  = mmod + 1  
      print *,'**********************************'
      print *,'*      Advancing Front Code      *'
      print *,'*                                *'
      print *,'*   Generation of Unstructured   *'
      print *,'*     Grid Inside Rectangle      *'
      print *,'**********************************'
      print *,' '
      print *,'The number of boundary points should be an even integer.'
      print *,'What is the number of points on the boundary ?'
      read  *, ip
      ! convert to even number
      if ( 2*(ip/2).ne.ip ) then
        ip = 2*(ip/2)
        print *,'Changed number of boundary points even value of ',ip 
      end if
      il = ip/4 + 1
      jl = ip/2 + 2 - il
      print *,'What is the total number of points ?'
      read  *, jfree
      ! checks
      if ( jfree.le.ip ) then
        print *,'Total points must exceed number of boundary points.'
        stop
      else if ( jfree.gt.jmax ) then
        print *,'Increase jmax and recompile.'
        stop
      else if ( 2*jfree.gt.imax ) then     
        ! the estimated number of cells is 2*jfree
        print *,'Increase imax and recompile.' 
        stop
      end if
      dx    = 1.0
      dy    = 1.0
      ! boundary points
      jn = 1  
      kf = 1
      do i = 1,il-1
        x(jn) = i*dx
        y(jn) = dy
        heap(1,kf) = jn
        heap(2,kf) = jn + 1
        jn    = jn + 1
        kf    = kf + 1
      end do
      do j = 1,jl-1
        x(jn) = il*dx
        y(jn) = j*dy
        heap(1,kf) = jn
        heap(2,kf) = jn + 1
        jn    = jn + 1
        kf    = kf + 1
      end do  
      do i = il,2,-1
        x(jn) = i*dx
        y(jn) = jl*dy
        heap(1,kf) = jn
        heap(2,kf) = jn + 1
        jn    = jn + 1
        kf    = kf + 1
      end do 
      do j = jl,2,-1
        x(jn) = dx
        y(jn) = j*dy
        heap(1,kf) = jn
        heap(2,kf) = jn + 1
        jn    = jn + 1
        kf    = kf + 1
      end do
      heap(2,kf-1) = 1
      kface = kf - 1
      ! interior points
      print *,'What is the initial seed for random number generator ?'
      read  *, iflag
      ran1     = rand(iflag)
      iflag    = 0
      xlen2    = ((il-1)*dx)**2
      xmin2    = 1.e-8
      itermax  = 10000
      do j = jn,jfree
        found  = .false.
        iter   = 1
        do while ( .not.found )
          if ( iter.gt.itermax ) then
            print *,'Cannot find random point.  Stop.'
            stop
          else 
            ran1  = rand(iflag)
            x(j)  = dx*(1+ (il-1)*ran1)
            ran2  = rand(iflag)
            y(j)  = dy*(1+ (il-1)*ran2)
            done = .false.
            jj = 1
            do while ( .not.done )
              if ( jj.eq.j ) then
                done  = .true.
                found = .true.
              else   
                ratio = ((x(j)-x(jj))**2 + (y(j)-y(jj))**2)/xlen2
                if ( ratio.lt.xmin2 ) then
                  done  = .true.   
                  found = .false.
                else
                  jj = jj + 1 
                end if              
              end if
            end do
            iter = iter + 1
          end if
        end do
      end do
      ! initialize array free(j)
      do j = 1,jfree
        free(j) = j
      end do
      call rehash(hash,heap,hmax,kface,mmod)  ! create hash table
      return
      end
      SUBROUTINE NDELETE(free,hash,hmax,heap,jfree,jn,mmod)
!
!     removes node jn from array free() if it is not connected
!     to any faces in the front
!
      integer free(*), hash(0:mmod,0:hmax), heap(2,*), hmax 
      logical delete, done

      mkey = mmod + 1
      key  = mod(jn,mkey)
      delete = .true.
      done = .false. 
      j    = 0
      do while ( .not.done )
        j = j + 1
        if ( j.gt.hash(key,0) ) then 
          done = .true.
        else 
          k = hash(key,j)
          if ( heap(1,k).eq.jn ) then
            delete = .false.
            done   = .true.
          end if 
        end if 
      end do
      if ( delete ) then
        jr = 0
        done = .false.
        do while ( .not.done )
          jr = jr + 1
          if ( jn.eq.free(jr) ) done = .true.
          if ( jr.gt.jfree ) then
            print *,'Error in subroutine NDELETE.  Stop.'
            stop
          end if
        end do 
        do j = jr,jfree-1
          free(j) = free(j+1)  
        end do
        jfree = jfree-1
      end if
      return
      end
      SUBROUTINE OUTPUT(free,icell,jfree,jmax,kface,cell,heap,np,x,y)
!
!     write output files in TECPLOT(TM) format
!     TECPLOT is a trademark of Amtec Engineering, Bellevue, WA
!
!     plot1 includes cells (zone 1), front (zone 2) and the 
!     free points (zone 3).  if all points have been triangulated,
!     then only zone 1 is created.
!
!     plot2 creates geometry file with each triangle shaded.  
!     this permits confirmation of the region covered by the
!     triangulation.
!
      character*80 filename
      integer      cell(3,*), free(*), heap(2,*), np(*)
      real         x(*), y(*)
      logical      found

      ! np(j) takes the absolute node number and gives the new node number
      do j = 1,jmax
        np(j) = 0 
      end do
      jnode = 0
      do i = 1,icell
        do k = 1,3 
          if ( np(cell(k,i)).eq.0 ) then
            jnode       = jnode + 1
            np(cell(k,i)) = jnode 
          end if
        end do
      end do
      filename = 'plot1'
      open(1,file=filename,status='unknown',access='sequential',
     >     form='formatted',err=1000)
      ! cells
      write(1,100) jnode, icell
100   format('TITLE = "Advancing Front "',/,'VARIABLES = X,Y',/,
     >       'ZONE T="Cells", F=FEPOINT, ET=TRIANGLE, N=',i5,
     >       ', E=',i5)
      do j = 1, jnode
        ! find newly numbered nodes in consecutive order 
        jj    = 0  
        found = .false. 
        do while ( .not.found )
          jj = jj + 1
          if ( jj.gt.jmax ) then
            print *,'Error in plot subroutine.  Stop.'
            stop
          else if ( np(jj).eq.j ) then
            found = .true.
          end if
        end do
        write(1,200)  x(jj), y(jj)
      end do
200   format(2(e11.4e3,1x))
      do i = 1, icell
        write(1,300)  np(cell(1,i)),np(cell(2,i)),np(cell(3,i))
      end do
300   format(4i10)
      ! front
      if ( kface.gt.0 ) then
        write(1,400) kface+1
400     format('ZONE T="Front", F=POINT, I=',i5)
        do k = 1,kface
          write(1,200) x(heap(1,k)), y(heap(1,k))   
        end do
        write(1,200) x(heap(2,kface)),y(heap(2,kface))
      end if 
      ! front and interior points for scatter plot
      if ( kface.gt.0 ) then
        write(1,500) jfree
500     format('ZONE T="Free Points", F=POINT, I=',i5)
        do j = 1,jfree
          write(1,200) x(free(j)), y(free(j))
        end do
      end if
      close(1) 
      ! geometry record
      ! all triangles should be colored blue and any gaps would be evident
      filename = 'plot2'
      open(1,file=filename,status='unknown',access='sequential',
     >     form='formatted',err=1000)
      igeo = 16           ! number of geometries per record
      irec = icell/igeo   ! number of records with igeo geometries 
      npts = 3
      do index = 1,irec
        write(1,600) igeo
600     format('GEOMETRY T=LINE,M=GRID,C=YELLOW,F=POINT,FC=BLUE',/,i5)
        do i = (index-1)*igeo+1,index*igeo
          write(1,300) npts
          write(1,200) x(cell(1,i)), y(cell(1,i))
          write(1,200) x(cell(2,i)), y(cell(2,i))
          write(1,200) x(cell(3,i)), y(cell(3,i))
        end do
      end do
      irem = icell-irec*igeo  ! number of geometries in last record
      if ( irem.gt.0 ) then
        write(1,600) irem
        do i = irec*igeo+1,icell
          write(1,300) npts
          write(1,200) x(cell(1,i)), y(cell(1,i))
          write(1,200) x(cell(2,i)), y(cell(2,i))
          write(1,200) x(cell(3,i)), y(cell(3,i))
        end do
      end if
      print *,'Wrote geometry file for ',icell,' cells'
      close(1) 
      return
1000  print *,'Error opening file.  Stop.'
      stop
      end

     






      
      LOGICAL FUNCTION OVERLAP(cell,icell,j1,j2,j3)
!
!     determines if cell j1-j2-j3 overlaps any existing cells
!
      integer cell(3,*)
      logical done

      overlap = .false.
      done    = .false.
      i       = 1
      do while ( .not.done )
        if ( i.gt.icell ) then
          done = .true.
        else
          i1 = cell(1,i)
          i2 = cell(2,i)
          i3 = cell(3,i)
          if ( ((i1.eq.j1).and.(i2.eq.j2)) .or.
     >         ((i1.eq.j2).and.(i2.eq.j3)) .or.
     >         ((i1.eq.j3).and.(i2.eq.j1)) .or.
     >         ((i2.eq.j1).and.(i3.eq.j2)) .or.
     >         ((i2.eq.j2).and.(i3.eq.j3)) .or.
     >         ((i2.eq.j3).and.(i3.eq.j1)) .or.
     >         ((i3.eq.j1).and.(i1.eq.j2)) .or.
     >         ((i3.eq.j2).and.(i1.eq.j3)) .or.
     >         ((i3.eq.j3).and.(i1.eq.j1)) ) then
            overlap = .true. 
            done    = .true.
          else 
            i = i + 1 
          end if
        end if
      end do
      return
      end


      SUBROUTINE REHASH(hash,heap,hmax,kface,mmod)
!
!     reconstructs hash table 
!
      integer hash(0:mmod,0:hmax), heap(2,*), hmax

      mkey = mmod + 1
      do i = 0,mmod
        hash(i,0) = 0
      end do
      do k = 1,kface
        key         = mod(heap(1,k),mkey)  ! key for hash table
        hash(key,0) = hash(key,0) + 1      ! increment no. entries
        if ( hash(key,0).le.hmax ) then
          hash(key,hash(key,0)) = k        ! add entry 
        else
          print *,'Maximum number of hash entries exceeded. Stop.'
          stop
        end if 
      end do
      return 
      end
      SUBROUTINE THIRD(cell,free,heap,icell,j1,j2,j3,jfree,jmax,
     >                 kface,np,x,y)
!
!     find third point to form triangle
!
!     definition of variables
!
!       inside   = .true. if there could be a point inside the circumscribing
!                         circle formed by (x1,y1), (x2,y2) and (x3,y3)  
!
!       done     = .true. if all points on list have been checked to
!                         determine if there is a point inside the circle
!
      parameter(toler=1.e-7)
      integer cell(3,*), free(*), heap(2,*), np(*) 
      logical done, inside, overlap 
      real    x(*), y(*)

      x1 = x(j1)
      x2 = x(j2)
      y1 = y(j1)
      y2 = y(j2)
      ! select a trial point on front or interior 
      call trial(cell,free,heap,icell,jfree,j1,j2,j3,jmax,kface,np,x,y)
      x3 = x(j3)
      y3 = y(j3)
      ! examine all available points other than j1, j2 and j3 to see 
      ! if there is a point within the circumscribing circle, and if so,
      ! use this point as the next choice for j3
      inside = .true.
      do while ( inside ) 
        ! center of circumscribing circle and square of radius
        call circle(x1,x2,y1,y2,x3,y3,xc,yc,r2)
        ! test all points in the available list
        done = .false.
        j    = 1
        do while ( .not.done )
          if ( j.gt.jfree ) then
            done   = .true.      !  all points have been checked, and
            inside = .false.     !  there are no points inside circle  
          else
            jn = free(j)
            ! jn is checked as possible third node iff it is not j1, j2,
            ! or j3, and if the cell j1-j2-j3 would not overlap other cells
            if ( (jn.ne.j1).and.(jn.ne.j2).and.(jn.ne.j3) ) then
              d2 = (x(jn)-xc)**2+(y(jn)-yc)**2 
              if ( d2 .lt. r2 ) then
                cr = (x2-x1)*(y3-y1)-(y2-y1)*(x3-x1)
                if ( (cr.gt.toler).and.
     >               (.not.overlap(cell,icell,j1,j2,jn)) ) then 
                  done   = .true.  ! no need to check remainder of points.
                  inside = .true.  ! there is a point within the circle.
                  j3   = jn        ! use this point as the next choice.
                  x3   = x(j3)
                  y3   = y(j3)   
                else
                  j = j + 1
                end if
              else
                j = j + 1 
              end if  
            else
              j = j + 1         
            end if
          end if 
        end do
      end do  
      return
      end
!
      SUBROUTINE TRIAL(cell,free,heap,icell,jfree,j1,j2,j3,
     >                 jmax,kface,np,x,y)
!
!     selects a trial node j3
!
      parameter(toler=1.e-7)
      integer cell(3,*), free(*), heap(2,*), np(*)            
      real    x(*), y(*)
      logical done, overlap

      x1 = x(j1)
      x2 = x(j2)
      y1 = y(j1)
      y2 = y(j2)
      done = .false.
      j  = 1 
      do while ( .not.done )
        if ( j.gt.jfree )  then
          print *,'Cannot find eligible trial point.  Stop.'
          call output(free,icell,jfree,jmax,kface,cell,heap,np,x,y)
          stop
        else 
          j3 = free(j)
          x3 = x(j3)     
          y3 = y(j3) 
          !  check to insure that
          !  1) j3 is not the same as j1 or j2
          !  2) j1-j2-j3 forms a triangle in CCW sense
          !  3) j1-j2-j3 does not overlap an existing triangle
          if ( (j3.ne.j1).and.(j3.ne.j2).and.
     >         (.not.overlap(cell,icell,j1,j2,j3)) ) then
            if ( ((x2-x1)*(y3-y1)-(y2-y1)*(x3-x1)).gt.toler ) then
              done = .true.  
            else 
              j = j + 1
            end if
          else
            j = j + 1
          end if
        end if 
      end do
      return
      end
!
