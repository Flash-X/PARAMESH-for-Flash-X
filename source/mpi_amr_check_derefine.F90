!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

      subroutine amr_check_derefine (mype)

! By K. Olson (NASA/GSFC and GMU), 4/97
!! MODIFICATIONS
!!  2022-12-01 Klaus Weide  Include Flashx_mpi_implicitNone.fh

      Use paramesh_dimensions
      Use physicaldata
      Use tree
      Use timings
      Use Paramesh_comm_data, ONLY : amr_mpi_meshComm

#include "Flashx_mpi_implicitNone.fh"

      Integer, intent(in) :: mype

! local variables ----------------------------------------------------------

      Integer :: i,j,k
      Integer :: nodetype2(maxblocks_tr)
      Integer :: nodetype_recv(maxblocks_tr)
      Integer :: tnodetype_recv(maxblocks_tr)
      Integer :: nodetype_send(maxblocks_tr)
      Integer :: ipar,ipar_proc
      Integer :: reqr(maxblocks_tr)
      Integer :: statr(MPI_STATUS_SIZE,maxblocks_tr)
      Integer :: isg,ierr,neighs,neighr,isend,jsend,ksend
      Integer :: nderefines, loop_count
      Logical :: derefine_chi(nchild,maxblocks_tr)
      Logical :: refine_par(maxblocks_tr)
      Logical :: trefine_par(maxblocks_tr)
      logical :: unmodified(maxblocks_tr)
      logical :: derefine0(maxblocks_tr)

! --------------------------------------------------------------------------

      no_of_Calls_check_derefine = no_of_Calls_check_derefine + 1

!-----If the block is marked for derefinement and it is not
!-----a leaf block Then Do not derefine it
      Do i = 1,lnblocks
         If (derefine(i).and.nodetype(i).ne.1) derefine(i) = .FALSE.
      End Do 

!-----Allow the child blocks to only derefine if their parent is NOT
!-----marked for refinement.
      neighr = 0
      Do i = 1,lnblocks
         refine_par(i) = .FALSE.
         If (parent(1,i) > 0) Then
            If (parent(2,i).ne.mype) Then
               neighr = neighr + 1
               Call MPI_IRECV(refine_par(i),                           & 
                              1,                                       & 
                              MPI_LOGICAL,                             & 
                              parent(2,i),                             & 
                              i,                                       & 
                              amr_mpi_meshComm,                          & 
                              reqr(neighr),                            & 
                              ierr)
            Else
               refine_par(i) = refine(parent(1,i))
            End If
         End If
      End Do 

      neighs = 0
      Do i = 1,lnblocks
         Do j = 1,nchild
          If (child(1,j,i) > 0) Then
             If (child(2,j,i).ne.mype) Then
                neighs = neighs + 1
                Call MPI_SSEND(refine(i),                              & 
                               1,                                      & 
                               MPI_LOGICAL,                            & 
                               child(2,j,i),                           & 
                               child(1,j,i),                           & 
                               amr_mpi_meshComm,                         & 
                               ierr)
            End If
          End If
         End Do 
      End Do 

      If (neighr > 0) Then
         Call MPI_WAITALL (neighr, reqr, statr, ierr)
      End If

      Do i = 1,lnblocks
         If (nodetype(i) == 1.and.derefine(i)) Then
            If (refine_par(i)) derefine(i)=.False.
         End If
      End Do 

!----Turn off refine flags of non-leaf blocks
      Do i = 1,lnblocks
         If (nodetype(i) > 1.and.refine(i)) refine(i) = .False.
      End Do 

!-----Check neighbors to check if OK to derefine
!-----set nodetype2 = 2 if it either has children or it is marked for
!-----refinement
!-----First initialize the record of parents of leaf blocks whose children 
!-----have all had their derefine flags left unmodified.
      unmodified = .True.
      loop_count = 0
      derefine0 = derefine
      nderefines = 1

      Do While (loop_count < 2.and.nderefines > 0)

      If (loop_count == 1) derefine = derefine0

      Do i = 1,lnblocks
        nodetype2(i) = 1
        If ( (child(1,1,i) >= 1.or.refine(i)) .and.                    & 
              unmodified(i)                ) Then ! this node has children 
!        If ( (child(1,1,i) >= 1.and.unmodified(i)) .or.                    & 
!              refine(i)                ) Then ! this node has children 
                                                  ! or it is marked for 
                                                  ! refinement Then its
                                                  ! type is 2
          nodetype2(i) = 2
        End If
      End Do 

!----Check for neighboring blocks which are more than one level of refinement
!----different
!----cycle through block faces
      Do k = 1,1+2*k3d
      Do j = 1,1+2*k2d
      Do i = 1,3

         If (i == 2 .and. j == 1+k2d .and. k == 1+k3d) then
            ! Don't do anything !!!
         Else

         isend = 3-i+1
         jsend = (1+k2d*2)-j+1
         ksend = (1+k3d*2)-k+1

         nodetype_recv(1:lnblocks) = 0

         neighr = 0
         Do isg = 1,lnblocks
            nodetype_recv(isg) = 0
            If (surr_blks(1,i,j,k,isg) > -1) Then
               If (surr_blks(2,i,j,k,isg).ne.mype) Then
                  neighr = neighr + 1
                  Call MPI_IRECV(nodetype_recv(isg),                   & 
                                 1,                                    & 
                                 MPI_INTEGER,                          & 
                                 surr_blks(2,i,j,k,isg),               & 
                                 surr_blks(1,i,j,k,isg),               & 
                                 amr_mpi_meshComm,                       & 
                                 reqr(neighr),                         & 
                                 ierr)
                  mess_counter_chk_deref = mess_counter_chk_deref + 1
               Else
                  nodetype_recv(isg) = nodetype2(surr_blks(1,i,j,k,isg))
               End If
            End If
         End Do 

!--------send nodetype2 to neigh if neighbor is off processor and nodetype2 = 2
         neighs = 0
         Do isg = 1,lnblocks
            If (surr_blks(1,isend,jsend,ksend,isg) > -1) Then
               If (surr_blks(2,isend,jsend,ksend,isg).ne.mype) Then
                  neighs = neighs + 1
                  Call MPI_SSEND(nodetype2(isg),                       & 
                                 1,                                    & 
                                 MPI_INTEGER,                          & 
                                 surr_blks(2,isend,jsend,ksend,isg),   & 
                                 isg,                                  & 
                                 amr_mpi_meshComm,                       & 
                                 ierr)
               End If
            End If
         End Do 
         
         If (neighr > 0) Then
            Call MPI_WAITALL (neighr, reqr, statr, ierr)
         End If
      
!         Do isg = 1,lnblocks
!            If (nodetype_recv(isg) == 2) nodetype2(isg) = 2
!         End Do 

         End If ! End If (i == 1 .and. j == 2 .and. ...)

!--------Now reset derefine flags based on value of nodetype2
         Do isg = 1,lnblocks
            If (nodetype_recv(isg) == 2 .and. derefine(isg)) Then
               derefine(isg) = .FALSE.
            End If
         End Do 

      End Do  ! End Do i
      End Do  ! End Do j
      End Do  ! End Do k

!-----If a block Do es not have a parent (i.e. = -1) Then you can't derefine
!-----it further so if it is marked for derefinement turn derefine off
      Do i = 1,lnblocks
         If (derefine(i).and.parent(1,i) < 0) derefine(i) = .FALSE.
      End Do 

!-----check if all siblings are also marked for derefinement, if not Then
!-----Don't derefine this block
!-----parents collect messages from children and count the number of children
!-----marked for derefinement (stored in nodetype_recv).

      neighr = 0
      Do isg = 1,lnblocks
         Do j = 1,nchild
            derefine_chi(j,isg) = .FALSE.
            If (child(1,j,isg) > -1) Then
            If (child(2,j,isg).ne.mype) Then
               neighr = neighr + 1
!--------------derefine_chi(j,isg) - this is just junk
               Call MPI_IRECV(derefine_chi(j,isg),                     & 
                              1,                                       & 
                              MPI_LOGICAL,                             & 
                              child(2,j,isg),                          & 
                              child(1,j,isg),                          & 
                              amr_mpi_meshComm,                          & 
                              reqr(neighr),                            & 
                              ierr)
            Else
               derefine_chi(j,isg) = derefine(child(1,j,isg))
            End If
            End If
         End Do 
      End Do 

!-----Children send a message to parent if they are marked for derefinement
      neighs = 0
      nodetype_recv(:) = 0    ! using this variable as a counter here

      Do i = 1,lnblocks
            ipar = parent(1,i) ! parent of i
            ipar_proc = parent(2,i) ! processor parent is stored on
         If (ipar > -1) Then
            If (ipar_proc.ne.mype) Then
               neighs = neighs + 1
               Call MPI_SSEND(derefine(i),                             & 
                              1,                                       & 
                              MPI_LOGICAL,                             & 
                              ipar_proc,                               & 
                              i,                                       & 
                              amr_mpi_meshComm,                          & 
                              ierr)
            End If
         End If
      End Do 

      If (neighr > 0) Then
         Call MPI_WAITALL (neighr, reqr, statr, ierr)
      End If
      
      Do i = 1,lnblocks
         Do j = 1,nchild
            If (derefine_chi(j,i)) Then
               nodetype_recv(i) = nodetype_recv(i) + 1
            End If
         End Do 
      End Do 
      nodetype_send(1:lnblocks) = nodetype_recv(1:lnblocks)

!-----Now parent sends nodetype_recv to its children if nodetype_recv = nchild
!-----record modifications to the derefine state of parents, so this
!-----info can be used in any further loop_count iterations
      Do isg = 1,lnblocks
         If (nodetype(isg) == 2.and.nodetype_send(isg) == nchild       & 
            .and. unmodified(isg)) Then
                         unmodified(isg) = .False.
                         work_block(isg) = gr_btWorkDefaultLeaf
         End If
      End Do 

      neighr = 0
      Do isg = 1,lnblocks
         If (parent(1,isg) > -1) Then
            If (parent(2,isg).ne.mype) Then
               neighr = neighr + 1
               Call MPI_IRECV(nodetype_recv(isg),                      & 
                              1,                                       &   
                              MPI_INTEGER,                             & 
                              parent(2,isg),                           & 
                              isg,                                     & 
                              amr_mpi_meshComm,                          & 
                              reqr(neighr),                            & 
                              ierr)
            End If
         End If
      End Do 

      neighs = 0
      Do isg = 1,lnblocks
         Do j = 1,nchild
            If (child(1,j,isg) >= 1) Then
               If (child(2,j,isg).ne.mype) Then
                  neighs = neighs + 1
                  Call MPI_SSEND(nodetype_send(isg),                   & 
                                 1,                                    & 
                                 MPI_INTEGER,                          & 
                                 child(2,j,isg),                       & 
                                 child(1,j,isg),                       & 
                                 amr_mpi_meshComm,                       & 
                                 ierr)
               Else
                  nodetype_recv(child(1,j,isg)) = nodetype_send(isg)
               End If
            End If
         End Do 
      End Do 
      
      If (neighr > 0) Then
         Call MPI_WAITALL (neighr, reqr, statr, ierr)
      End If

!-----Now loop though the blocks one final time and if nodetype_recv .ne. nchild
!-----and if derefine = .TRUE., then Don't derefine
      Do isg = 1,lnblocks
         If (derefine(isg).and.nodetype_recv(isg).ne.nchild) Then
            derefine(isg) = .FALSE.
         End If
      End Do 
      
      nderefines = 0
      Do isg = 1,lnblocks
        If ( nodetype(isg) == 1 .and.                                  & 
           (derefine0(isg).and.(.not.derefine(isg))) )                 & 
           nderefines = nderefines + 1
      End Do 

      Call comm_int_sum_to_all(nderefines,nderefines)

      loop_count = loop_count+1

      End Do  ! End Do While (loop_count < 2.and.nderefines > 0)

      Return
      End Subroutine amr_check_derefine


