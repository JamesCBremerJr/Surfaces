        implicit double precision (a-h,o-z)

        end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       This is the end of the debugging code and beginning of the
c       thread safe heap allocation/deallocation code.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       This file contains a primitive code for managing a heap.
c   
c      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        subroutine heap_init(ier,maxchunks,heap,lheap)
        implicit double precision (a-h,o-z)
        dimension heap(1)
c
        ier = 0
c
c       Allocate memory for the contorl structures.
c      

c
c       Initialize the control structures.
c
        end


        subroutine heap_allocate(ier,heap,lchunk,ichunk)
        implicit double precision (a-h,o-z)
c
c       Allocate a chunk from the heap.
c
        end


        subroutine heap_deallocate(ier,heap,ichunk)
        implicit double precision (a-h,o-z)
c
c       Deallocate a chunk specified by the user via its index.
c
        end


        subroutine heap_pointer(ier,heap,ichunk,iptr)
        implicit double precision (a-h,o-z)
c
c       Fetch 
c
        end


        subroutine heap_gc(ier,heap)
        implicit double precision (a-h,o-z)
c
c       Perform garbage collection on the heap.
c
        end
