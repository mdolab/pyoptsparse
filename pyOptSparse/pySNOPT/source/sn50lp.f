!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
*     File  sn50lp.f
*
*     s5LP     s5BSx    s5dgen   s5FixX   s5FixS   s5finE   s5getB
*     s5grdE   s5hs     s5Inf    s5InfE   s5IniE   s5LG     s5LPit
*     s5pric   s5rc     s5rcE    s5setE   s5setp   s5setx   s5step
*
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5LP
     &   ( iExit, Prob, ProbTag, Elastc, subopt,
     &     LPlog, needLU, needx,
     &     m, n, nb, nDegen, itLP, itLPmax, itn,
     &     lEmode, lvlInf, PrtLvl,
     &     minimz, iObj, sclObj, ObjAdd, tolFP, tolLP, tolx,
     &     nInf, sInf, wtInf, piNorm, rgNorm,
     &     ne, nlocA, locA, indA, Acol,
     &     hElast, hEstat, hfeas, hs, kBS,
     &     Ascale, bl, bu, blBS, buBS,
     &     gBS, pi, rc, nrhs0, nrhs, rhs, x, xBS, xFreez,
     &     iy, iy1, y, y1, y2,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      external
     &     LPlog
      logical
     &     Elastc, needLU, needx
      integer
     &     Prob, iExit, iObj, m, mBS, minimz, nrhs0, nrhs, n, nb, ne,
     &     nlocA, nInf, nS, itLP, itLPmax, itn, lEmode, lvlInf, PrtLvl,
     &     nDegen, lencw, leniw, lenrw, subopt, hElast(nb),
     &     hEstat(nb), hs(nb), locA(nlocA), indA(ne), kBS(m+1),
     &     hfeas(m+1), iy(nb), iy1(nb), iw(leniw)
      double precision
     &     sclObj, ObjAdd, tolFP, tolLP, tolx, sInf, wtInf, piNorm,
     &     rgNorm, Acol(ne), Ascale(nb), bl(nb), bu(nb), blBS(m+1),
     &     buBS(m+1), rc(nb), rhs(nrhs0), x(nb), xBS(m+1), xFreez(nb),
     &     gBS(m+1), pi(m), y(nb), y1(nb), y2(nb), rw(lenrw)
      character
     &     ProbTag*20, cw(lencw)*8

!     ==================================================================
!     s5LP   solves a linear program.
!
!     The problem type can be:
!       type = 'FP '  feasible point only
!       type = 'FPE'  feasible point for equalities only
!       type = 'FPS'  feasible point for QP subproblem
!       type = 'LP '  LP problem
!
!     The optimization can pass through the following phases:
!
!       Phase 1               find a feasible point for all variables
!
!       Elastic Phase 1       make the non-elastic variables feasible
!                             while allowing infeasible elastics
!
!       Phase 2               minimize the objective
!
!       Elastic Phase 2       minimize a composite objective while
!                             keeping the non-elastics feasible
!
!                             In this phase, lvlInf means the following:
!
!                 lvlInf = 0  zero     weight on the infeasibilities
!                                      (infeasibillities are ignored)
!                          1  finite   weight on the infeasibilities
!                          2  infinite weight on the infeasibilities
!                                      (the objective is ignored)
!
!     On exit:
!
!      iExit         Status
!      -----         ------
!        -3          Too many iterations
!        -2          LP is unbounded
!        -1          Non-elastic variables are infeasible
!         0          LP solution found
!        >0          Fatal error
!
!     The array kBS is a permutation on the column indices.
!     kBS(1  :m )  holds the column indices of the basic variables.
!     Superbasics have been temporarily fixed at their current value.
!
!     30 Sep 1991: First version based on Qpsol routine lpcore.
!     20 Jul 1996: Slacks changed to be the row value.
!     06 Aug 1996: First Min Sum version.
!     14 Jul 1997: Thread-safe version.
!     24 Dec 1999: Suboptimization option added.
!     01 Aug 2003: snEXIT and snPRNT adopted.
!     24 Dec 2003: pi checked for NaN and Inf entries.
!     07 May 2006: s4ksav now handles negative values of hs.
!     ==================================================================
      character
     &     str*115
      logical
     &     checkx, chkFea, chkpi, gotE, feasbl, gotg, gotR, incres,
     &     jstFea, jstPh1, LUok, needf, needv, needpi, newB, newLU,
     &     optiml, prt10, prtLog, prtSum
      integer
     &     inform, itnfix, itnlim, jq, jBq, jBr,
     &     jSq, jSr, kchk, kDegen, kfac, klog, kObj, kp, kPrc, kPrPrt,
     &     ksav, kSumm, leng, lenL0, lenL, lenU0, lenU, LUitn, LUmod,
     &     LUsiz0, LUmax, LUreq, linesL, linesS, MnrHdP, MnrHdS, nBS,
     &     nElast, nFac, nFreez, nInfE, nonOpt, nSwap, nnH, neg,
     &     nfix(2), PrintL, PrintS, toldj1, toldj2, toldj3, typeLU
      double precision
     &     Anorm, Bgrwth, Bold, featol, condHz, djq, djqPrt, infBnd,
     &     ObjLP, ObjPrt, pivot, rowerr, sgnObj, sInfE, step, tolx0,
     &     tolinc, weight, dummy(1)
*     ------------------------------------------------------------------
      double precision   zero
      parameter         (zero   = 0.0d+0)
      integer            Intern,     Extern
      parameter         (Intern = 0, Extern = 1)
      integer            Check
      parameter         (Check  = 1)
      integer            FP,         FPE,        FPS
      parameter         (FP     = 0, FPE    = 3, FPS    = 4)
      integer            B
      parameter         (B      = 0)
      integer            WithB
      parameter         (WithB  = 1)
      integer            Init,       Optml,      Cycle
      parameter         (Init   = 0, Optml  = 1, Cycle = 2)

      parameter         (toldj1 = 184) ! phase 1 dj tol for p.p.
      parameter         (toldj2 = 185) ! phase 2 dj tol for p.p.
      parameter         (toldj3 = 186) ! current optimality tol
      parameter         (lenL0  = 171) ! size of L0
      parameter         (lenU0  = 172) ! size of initial  U
      parameter         (lenL   = 173) ! size of current  L
      parameter         (lenU   = 174) ! size of current  U
      parameter         (kObj   = 205) ! xBS(kObj) is the obj. slack
      parameter         (LUitn  = 215) ! itns since last factorize
      parameter         (LUmod  = 216) ! number of LU mods
      parameter         (PrintL = 218) ! (on/off) log     status
      parameter         (PrintS = 219) ! (on/off) summary status
      parameter         (linesL = 220) ! # lines in log     file
      parameter         (linesS = 221) ! # lines in summary file
      parameter         (MnrHdP = 223) ! >0 => Minor heading for iPrint
      parameter         (MnrHdS = 225) ! >0 => Minor heading for iSumm
*     ------------------------------------------------------------------
      infBnd    = rw( 70) ! definition of an infinite bound

      kchk      = iw( 58) ! check (row) frequency
      kfac      = iw( 59) ! factorization frequency
      ksav      = iw( 60) ! save basis map
      klog      = iw( 61) ! log/print frequency
      kSumm     = iw( 62) ! Summary print frequency
      kDegen    = iw( 63) ! max. expansions of featol
      itnlim    = iw( 89) ! limit on total iterations
      nFac      = iw(210) ! # of LU factorizations

      if (nFac .gt. 0) then
         LUsiz0    = iw(lenL0) + iw(lenU0)
         LUmax     = 2*LUsiz0
      end if

      nnH    = 0                ! local value of nnH
      neg    = 0
      leng   = 1
      nS     = 0                ! local value of nS
      mBS    = m + 1

      prt10  =              PrtLvl .ge. 10
      prtLog =              PrtLvl .ge.  1  .and.
     &         (mod( itLP,  klog ) .eq.  0  .and.  itLP   .ne. 0  .or.
     &                      klog   .eq.  1                      )
      prtSum =              PrtLvl .ge.  1  .and.
     &         (mod( itLP,  kSumm) .eq.  0  .and.  itLP   .ne. 0  .or.
     &                      kSumm  .eq.  1                      )
      iw(PrintL) = 0
      iw(PrintS) = 0
      if (prtLog) iw(PrintL) = 1
      if (prtSum) iw(PrintS) = 1

      kPrc   = 0              ! last section scanned in part. pricing
      iExit  = 0
      LUreq  = 0

      chkFea = .true.
      chkpi  = .true.
      feasbl = .false.
      gotE   = .false.
      gotR   = .false.
      jstFea = .false.

*     ------------------------------------------------------------------
*     s5LP operates in either ``Normal'' or ``Elastic'' mode.
*     Everything is normal unless a weighted sum is being minimized or
*     the constraints are infeasible.
*     The logical feasbl refers to the non-elastic variables.
*     wtInf  is the optional parameter Infeasibility Weight.
*     ------------------------------------------------------------------
!     jstPh1 = stop at the end of phase 1 (either regular or elastic)

      jstPh1 = Prob .eq. FP .or.  Prob   .eq. FPE  .or.  Prob  .eq. FPS

!     The phase 2 objective is F1 + wtInf*F2.

      if (Elastc) then
         needf = lvlInf .ne. 2 ! F1 required in phase 2
         needv = lvlInf .ne. 0 ! F2 required in phase 2
      else
         needf = .true.
         needv = .false.
      end if

      needpi = .true.
      newLU  = .true.
      optiml = .false.

      condHz = zero
      ObjLP  = zero
      pivot  = zero
      rgnorm = zero
      step   = zero
      nInfE  = 0
      jq     = 0
      djq    = zero
      jBq    = 0                ! x(jBq) is the incoming   BS
      jBr    = 0                ! x(jBr) is the outgoing   BS
      jSq    = 0                ! x(jSq) is the incoming SBS
      jSr    = 0                ! x(jSr) is the outgoing SBS
      kPrPrt = 0
      sgnObj = minimz
      typeLU = B

      nFreez = 0

      dummy(1)   = zero
      rw(toldj1) = 100.0d+0     ! Used only for LP partial pricing
      rw(toldj2) = 100.0d+0     !

      call s5hs
     &   ( Intern, nb, bl, bu, hs, x )
      call s5dgen
     &   ( inform, Init, PrtLvl, nb, nInf, itn,
     &     featol, tolx, tolinc, hs, bl, bu, x,
     &     itnfix, nfix, tolx0,
     &     iw, leniw, rw, lenrw )

**    ======================Start of main loop==========================
*+    do while (.not. optiml  .and.  iExit .eq. 0)
  100 if       (.not. optiml  .and.  iExit .eq. 0) then
*        ===============================================================
*        Check the initial  x  and move it onto  ( A -I )*x = b.
*        If needLU is true, this will require a basis factorization.
*        ===============================================================
*        If necessary,  factorize the basis  ( B = LU ) and set x.
*        If needLU is false on entry to s5LP, the first call to s2Bfac
*        will try to use existing factors.
*        If needLU is true on entry to s5LP, an LU factorization of
*        type typeLU is computed.
*
*        LUreq =  1  Frequency
*        LUreq =  2  LU nonzeros increased
*        LUreq =  3
*        LUreq =  4
*        LUreq =  5  Singular after LU mod
*        LUreq =  6  Unstable LU mod (growth in new column of U)
*        LUreq =  7  Not enough memory
*        LUreq =  8
*        LUreq =  9
*        LUreq = 10  Row error in setx
*        LUreq = 11  Big  dx   in setx or setpi
!        LUreq = 23  Infeasibility after refactorization
*        ---------------------------------------------------------------
         if (LUreq .gt. 0) needLU = .true.

         if (needx  .or.  needLU) then
            call s2Bfac
     &         ( iExit, typeLU, needLU, newLU, newB,
     &           iObj, itn, PrtLvl, LUreq,
     &           m, mBS, n, nb, nnH, nS, nSwap,
     &           ne, nlocA, locA, indA, Acol,
     &           kBS, hs, bl, bu, blBS, buBS,
     &           nrhs0, nrhs, rhs, x, xBS,
     &           iy, iy1, y, y2,
     &           iw, leniw, rw, lenrw )
            LUsiz0    = iw(lenL0) + iw(lenU0)
            LUmax     = 2*LUsiz0

            if (iExit .gt. 0) go to 100

            gotE   = .false.    ! Check hEstat in elastic mode.
            needpi = .true.     ! Recalculate the pi's.
            needx  = .false.
            chkpi  = .true.
            chkFea = .true.

            if (prt10) iw(MnrHdP) = 1
         end if

         nBS    = m + nS
         nInf   = 0
         sInf   = zero
         call dload ( nBS, zero, gBS, 1 )

         if (Elastc  .and.  .not. gotE) then
*           ------------------------------------------------------------
*           Reset blBS and buBS for any violated elastics.
*           These values are used in s5step.
*           Strictly feasible elastics are returned to normality.
*           ------------------------------------------------------------
            call s5setE
     &         ( nBS, nb, nElast, featol, infBnd,
     &           hElast, hEstat, kBS,
     &           bl, bu, blBS, buBS, xBS )
            gotE  = .true.
         end if

         if ( chkFea ) then

*           In Phase 1 or just after a factorize, check the feasibility
*           of the basic and superbasic non-elastics.
*           jstFea  indicates that we have just become feasible.
*           jstFea is turned off once a step is taken.

            call s5Inf
     &         ( nBS, featol, nInf, sInf, hfeas, blBS, buBS, gBS, xBS )

            if (nInf .gt. 0) then

*              Non-elastics are infeasible.
!              If necessary, switch back to the feasibility phase, after
!              refactorization.
*              Print something if the basis has just been refactorized.

               if ( feasbl ) then
                  call s2tryLU
     &               ( itn, 23, nS, LUreq, LUok, typeLU,
     &                 iw, leniw, rw, lenrw )
                  if (.not. LUok) iExit = 11
                  feasbl = .false.
                  go to 100
               end if

               if (prt10  .and.  iw(LUitn) .eq. 0) then
                  write(str, 1030) itn, nInf, sInf
                  call snPRNT( 21, str, iw, leniw )
               end if
            end if

*           Feasbl = true means that the nonelastics are feasible.
*                    This defines Phase 2.

            if (.not. feasbl) then
               jstFea = nInf .eq. 0
            end if
            feasbl = nInf .eq. 0
            chkFea = nInf .gt. 0
         end if ! chkFea

         if ( Elastc ) then
*           ------------------------------------------------------------
*           Compute the sum of infeasibilities of the elastic variables.
*           ------------------------------------------------------------
            call s5InfE
     &         ( nb, nBS, hEstat, kBS, nInfE, sInfE, bl, bu, x )
            nInf = nInf + nInfE
            sInf = sInf + sInfE
         end if

         ObjLP  = zero
         if (iObj .gt. 0) then
            ObjLP  = xBS(iw(kObj))*sclObj
         end if

         if (feasbl  .and.  jstPh1) then
*           ------------------------------------------------------------
*           The non-elastic variables just became feasible.
*           That is all that is needed.
*           ------------------------------------------------------------
            djqPrt = zero
            pinorm = zero
            call dload ( m, zero, pi, 1 )
            optiml = .true.

         else

            if ( feasbl ) then
!              ---------------------------------------------------------
!              Feasible for the nonelastics.
!              (Elastc = false means no elastics.)
!              ---------------------------------------------------------
               if ( needf ) then
                  if (iObj .ne. 0) then
                     gBS(iw(kObj)) = sgnObj*sclObj
                  end if
               end if

               if (Elastc  .and.  nInfE .gt. 0  .and.  needv) then
                  call s5grdE
     &               ( nb, nBS, wtInf, hEstat, kBS, gBS )
               end if
            end if

            if ( needpi ) then
               call dcopy
     &            ( m, gBS, 1, y, 1 )
               call s5setp
     &            ( inform, m, chkpi, pinorm, y, pi,
     &              iw, leniw, rw, lenrw )
               if (inform .ne. 0) then
                  if (inform .gt. 0) then
                     iExit = inform
                  else          ! pi is infinite or contains a NaN/Inf.
                     call s2tryLU
     &                  ( itn, 11, nS, LUreq, LUok, typeLU,
     &                    iw, leniw, rw, lenrw )
                     if (.not. LUok) iExit = 43
                  end if
                  go to 100
               end if
               needpi = .false.
            end if

*           ============================================================
*           Check for optimality.
*           Find the reduced costs.
*           ============================================================
            if ( feasbl ) then
               rw(toldj3) = tolLP
            else
               rw(toldj3) = tolFP
            end if

            kPrPrt = kPrc
            jq     = 0
            djqPrt = djq
            djq    = zero
            gotg   = .false.
            weight = zero
            if (Elastc  .and.  feasbl) then
               weight = wtInf
            end if

            call s5pric
     &         ( Elastc, feasbl, incres, gotg, subopt,
     &           itn, m, n, nb, leng, neg, nnH,
     &           nS, nFreez, nonOpt, weight, sgnObj, piNorm,
     &           jq, djq, kPrc, rw(toldj1),
     &           ne, nlocA, locA, indA, Acol,
     &           hElast, hs, bl, bu, dummy, pi, rc, x, xFreez,
     &           iw, leniw, rw, lenrw )
            optiml = nonOpt .eq. 0
         end if ! feasbl

         if ( optiml ) then
*           ------------------------------------------------------------
*           Apparently we are optimal.
*           See if any nonbasics have to be set back on their bounds.
*           ------------------------------------------------------------
            call s5dgen
     &         ( inform, Optml, PrtLvl,
     &           nb, nInf, itn,
     &           featol, tolx, tolinc, hs, bl, bu, x,
     &           itnfix, nfix, tolx0,
     &           iw, leniw, rw, lenrw )
            optiml = inform .eq. 0

            if ( optiml ) then

*              So far so good.   Check the row residuals.

               if (iw(LUitn) .gt. 0) then
                  call s5setx
     &               ( inform, Check, itn,
     &                 m, n, nb, nBS, rowerr,
     &                 ne, nlocA, locA, indA, Acol,
     &                 kBS, xBS, nrhs0, nrhs, rhs, x, y, y2,
     &                 iw, leniw, rw, lenrw )
                  optiml = inform .eq. 0
                  LUreq  = inform
               end if
            end if

*           If x is not optimal, set  x  so that ( A  -I )*x = b
*           and check feasibility.

            if (.not. optiml) then
               needx  = .true.
               go to 100
            end if
         end if ! optiml

!        ============================================================
!        Print the details of this iteration.
!        ============================================================
         if (jstPh1  .and.  optiml) then
!           Relax, we are about to exit without printing
         else
            ObjPrt = zero
            if ( feasbl ) then
               if ( needf ) then
                  ObjPrt = ObjAdd + ObjLP
               end if
               if ( needv ) then
                  ObjPrt = ObjPrt + sgnObj*wtInf*sInf
               end if
            end if

            call LPlog
     &         ( Prob, ProbTag,
     &           Elastc, gotR, jstFea, feasbl,
     &           m, mBS, nnH, nS, jSq, jBr, jSr,
     &           iw(linesL), iw(linesS), itn, itLP, kPrPrt, lvlInf,
     &           pivot, step, nInf, sInf, wtInf,
     &           ObjPrt, condHz, djqPrt, rgNorm, kBS, xBS,
     &           iw, leniw )
         end if

         jBq    = 0
         jBr    = 0
         jSq    = 0
         jSr    = 0
         kPrPrt = 0

         if ( optiml ) then
!           ------------------------------------------------------------
!           Convergence.
!           ------------------------------------------------------------
            if (nInf .gt. 0) then

!              No feasible point.
!              Stop or continue in elastic mode, depending on the
!              specified level of infeasibility.

               if (lEmode .gt. 0) then
                  if (Elastc) then
                     if (feasbl) then

*                       Find the final sum of infeasibilities.

                        call s5finE
     &                     ( nBS, nb, nInf, sInf, featol,
     &                       hEstat, kBS, bl, bu, xBS )
                     else

*                       The nonelastic bounds cannot be satisfied
*                       by relaxing the elastic variables. Exit.

                        iExit = -1 ! Infeasible nonelastics
                     end if
                  else
                                ! infeasible constraints in Normal mode.
                                ! Start Elastic Phase 1.
                     if (prt10) then
                        write(str, 8050) itn, ProbTag
                        call snPRNT( 23, str, iw, leniw )
                        write(str, 8060) itn
                        call snPRNT( 23, str, iw, leniw )
                        iw(MnrHdP) = 1
                        iw(MnrHdS) = 1
                     end if

                     Elastc = .true.
                     needf  = lvlInf .ne. 2 ! Need F1 in elastic phase 2
                     needv  = lvlInf .ne. 0 ! Need F2 in elastic phase 2
                     needpi = .true.
                     optiml = .false.
                     djq    = zero
                     step   = zero
                     call s5IniE
     &                  ( nb, nBS, nElast, featol, infBnd,
     &                    hElast, hEstat, kBS,
     &                    bl, bu, blBS, buBS, xBS )
                  end if
                  go to 100
               end if
            end if

            if (prt10 .and. .not. jstPh1) then
               If (jq .ne. 0) then
                  djq    = sgnObj*djq
                  if (klog .eq. 1) then
                     write(str, 1010) djq, jq, rgnorm, pinorm
                     call snPRNT( 31, str, iw, leniw )
                  end if
               else
                  if (klog .eq. 1) then
                     write(str, 1020)          rgnorm, pinorm
                     call snPRNT( 31, str, iw, leniw )
                  end if
               end if
            end if
         else
*           ------------------------------------------------------------
*           A nonbasic has been selected to become superbasic.
*           Compute the vector y such that B y = column jq.
*           ------------------------------------------------------------
*           Unpack column jq into  y1  and solve  B*y = y1.
*           The altered  y1  satisfies  L*y1 = ajq.
*           It is used later to modify L and U.

            call s2unpk
     &         ( jq, m, n, ne, Anorm, nlocA, locA, indA, Acol, y1 )
            call s2Bsol
     &         ( iExit, WithB, m, y1, y, iw, leniw, rw, lenrw )
            if (iExit .gt. 0) return

*           ============================================================
*           Take a simplex step.  A variable will become nonbasic
*           at the new x.
*           ============================================================
            if (itn  .ge. itnlim  .or.  itLP .ge. itLPmax) then
               iExit = -3       ! Excess iterations
               go to 100
            end if

            itLP   = itLP   + 1
            itn    = itn    + 1
            iw(LUitn) = iw(LUitn)  + 1
            newLU  = .false.
            chkpi  = .false.
            jstFea = .false.

*           Decide if we want to print something this iteration.

            prtLog = PrtLvl .ge. 1  .and.  mod( itLP, klog  ) .eq. 0
            prtSum = PrtLvl .ge. 1  .and.  mod( itLP, kSumm ) .eq. 0

            iw(PrintL) = 0
            iw(PrintS) = 0
            if (prtLog) iw(PrintL) = 1
            if (prtSum) iw(PrintS) = 1

*           ------------------------------------------------------------
*           Take a simplex step.
*           The new  x  will still be at a vertex (possibly artificial)
*           Check for unboundedness (inform = 1).
*           ------------------------------------------------------------
            call s5LPit
     &         ( inform, feasbl, incres, needpi, Elastc,
     &           m+1, m, nb, nDegen, LUreq,
     &           kp, jBq, jSq, jBr, jSr, jq,
     &           featol, pivot, step, tolinc,
     &           hElast, hEstat, hfeas, hs, kBS,
     &           bl, bu, blBS, buBS,
     &           x, xBS, y, y1,
     &           iw, leniw, rw, lenrw )
            if (inform .eq. 1) then
               iExit = -2       ! Unbounded direction
               go to 100
            end if

*           Increment featol every iteration.

            featol = featol + tolinc

*           ============================================================
*           Test for error condition and/or frequency interrupts.
*           ============================================================
*           (1) Save a basis map (frequency controlled).
*           (2) Every kdegen iterations, reset featol and move nonbasic
*               variables onto their bounds if they are very close.
*           (3) Refactorize the basis if it has been modified
*               too many times.
*           (4) Update the LU factors of the basis if requested.
*           (5) Check row error (frequency controlled).

            if (mod(itn,ksav) .eq. 0) then
               call s4ksav
     &            ( minimz, m, n, nb, nS, mBS,
     &              itn, nInf, sInf, ObjLP, kBS, hs,
     &              Ascale, bl, bu, x, xBS, cw, lencw, iw, leniw )
            end if

            if (mod( itn, kdegen ) .eq. 0) then
               call s5dgen
     &            ( inform, Cycle, PrtLvl,
     &              nb, nInf, itn,
     &              featol, tolx, tolinc, hs, bl, bu, x,
     &              itnfix, nfix, tolx0,
     &              iw, leniw, rw, lenrw )
               needx  = inform .gt. 0
            end if

            if (LUreq .eq. 0) then
               if (     iw(LUmod) .ge. kfac-1) then
                  LUreq  = 1
               else if (iw(LUmod) .ge. 20  .and.
     &                                iw(lenL)+iw(lenU) .gt. LUmax) then
                  Bgrwth = iw(lenL) + iw(lenU)
                  Bold   = LUsiz0
                  Bgrwth = Bgrwth/Bold
                  if ( prt10 ) then
                     write(str, 1000) Bgrwth
                     call snPRNT( 21, str, iw, leniw )
                  end if
                  LUreq  = 2
               else
                  checkx = mod(iw(LUitn),kchk) .eq. 0
                  if (checkx  .and.  .not. needx) then
                     call s5setx
     &                  ( inform, Check, itn,
     &                    m, n, nb, nBS, rowerr,
     &                    ne, nlocA, locA, indA, Acol,
     &                    kBS, xBS, nrhs0, nrhs, rhs, x, y, y2,
     &                    iw, leniw, rw, lenrw )
                     LUreq  = inform
                  end if
               end if
               if (LUreq .gt. 0) typeLU = B
            end if
         end if ! not optiml

         go to 100
*+    end while
      end if
*     ======================end of main loop============================
*
      call s5hs
     &   ( Extern, nb, bl, bu, hs, x )

      return

 1000 format(' ==> LU file has increased by a factor of', f6.1)
 1010 format(' Biggest dj =', 1p, e11.3, ' (variable', i7, ')',
     &       '    norm rg =',     e11.3, '   norm pi =', e11.3)
 1020 format(   ' Norm rg =', 1p, e11.3, '   norm pi =', e11.3)
 1030 Format(' Itn', i7, ': Infeasible nonelastics.  Num =', i5, 1p,
     &       '   Sum of Infeasibilities =', e8.1 )

 8050 format(' Itn', i7, ': Infeasible ', a)
 8060 format(' Itn', i7, ': Elastic Phase 1 -- making',
     &                   ' non-elastic variables feasible')

      end ! subroutine s5LP

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5BSx
     &   ( Task, nBS, nb, kBS, x, xBS )

      implicit
     &     none
      integer
     &     Task, nBS, nb, kBS(nBS)
      double precision
     &     xBS(nBS), x(nb)

*     =================================================================
*     s5BSx   copies free variables from  xBS  into  x  or vice versa,
*             depending on whether  Task is 'xBS to x' or 'x to xBS'.
*
*     07 Nov 1991: First version based on Minos routine m5bsx.
*     21 Aug 1999: Current version of s5BSx.
*     =================================================================
      integer
     &     j, k
*     ------------------------------------------------------------------
      integer            xtoxBS,     xBStox
      parameter         (xtoxBS = 0, xBStox = 1)
*     ------------------------------------------------------------------

      if (Task .eq. xBStox) then
         do k = 1, nBS
            j     = kBS(k)
            x(j)  = xBS(k)
         end do

      else if (Task .eq. xtoxBS) then
         do k = 1, nBS
            j      = kBS(k)
            xBS(k) = x(j)
         end do
      end if

      end ! subroutine s5BSx

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5dgen
     &   ( inform, Task, PrtLvl, nb, nInf, itn,
     &     featol, tolx, tolinc, hs, bl, bu, x,
     &     itnfix, nfix, tolx0,
     &     iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     Task, inform, itnfix, PrtLvl, nb, nInf, itn, leniw, lenrw,
     &     hs(nb), nfix(2), iw(leniw)
      double precision
     &     featol, tolx, tolinc, tolx0, bl(nb), bu(nb), x(nb), rw(lenrw)

*     ==================================================================
*     s5dgen performs most of the manoeuvres associated with degeneracy.
*     The degeneracy-resolving strategy operates in the following way.
*
*     Over a cycle of iterations, the feasibility tolerance featol
*     increases slightly (from tolx0 to tolx1 in steps of tolinc).
*     This ensures that all steps taken will be positive.
*
*     After kdegen consecutive iterations, nonbasic variables within
*     featol of their bounds are set exactly on their bounds and the
*     basic variables are recomputed to satisfy ( A  -I )*x = b.
*     featol is then reduced to tolx0 for the next cycle of iterations.
*
*
*     If Task = Init, s5dgen initializes the parameters:
*
*     featol  is the current feasibility tolerance.
*     tolx0   is the minimum feasibility tolerance.
*     tolx1   is the maximum feasibility tolerance.
*     tolinc  is the increment to featol.
*     kdegen  is the expand frequency (specified by the user).
*             it is the frequency of resetting featol to tolx0.
*     nDegen  counts the number of degenerate steps (not used here, but
*             incremented by s5step).
*     itnfix  is the last iteration at which an 'Optimal' or 'Cycle'
*             entry set nonbasics onto their bound.
*     nfix(j) counts the number of times an 'Optimal' entry has
*             set nonbasics onto their bound,
*             where j=1 if infeasible, j=2 if feasible.
*
*     tolx0 and tolx1 are both close to the feasibility tolerance tolx
*     specified by the user.  (They must both be less than tolx.)
*
*
*     If Task = Cycle,  s5dgen has been called after a cycle of
*     kdegen iterations.  Nonbasic x(j)s are examined to see if any are
*     off their bounds by an amount approaching featol.  inform returns
*     how many.  Deviations as small as tolz (e.g. 1.0d-11) are not
*     counted. If inform is positive, the basic variables are
*     recomputed.  It is assumed that s5LP or s5QP will then continue
*     iterations.
*
*     itnfix, nfix, tolx0 could be treated as SAVED variables.
*     They are included as arguments to prevent threading problems in a
*     multiprocessing environment.
*
*     If Task = Optml,  s5dgen is being called after a subproblem
*     has been judged optimal, infeasible or unbounded.
*     Nonbasic x(j)s are examined as above.
*
*     07 Nov 1991: First version based on Minos routine m5dgen.
*     12 Jul 1997: Thread-safe version.
*     11 May 2001: Removed summary file printing.
*     01 Aug 2003: snPRNT adopted.
*     01 Aug 2003: Current version of s5dgen.
*     ==================================================================
      character
     &     str*80
      integer
     &     j, kDegen, maxfix
      double precision
     &     b1, b2, d1, d2, eps1, tolx1, tolz
*     ------------------------------------------------------------------
      integer            Init,     Optml,     Cycle
      parameter         (Init = 0, Optml = 1, Cycle = 2)
      double precision   zero
      parameter         (zero = 0.0d+0)
*     ------------------------------------------------------------------
      eps1      = rw(  3) ! eps**(2/3)
      kDegen    = iw( 63) ! max. expansions of featol

      inform    = 0
      if (Task .eq. Init) then

*        Task = Initialize.
*        Initialize at the start of each major iteration.
*        kdegen is the expand frequency      and
*        tolx   is the feasibility tolerance
*        (specified by the user).  They are not changed.
*        nDegen counts the total number of degenerate steps, and is
*        initialized by s5solv or s8solv.

         itnfix  = 0
         nfix(1) = 0
         nfix(2) = 0
         tolx0   = 0.5d+0 *tolx
         tolx1   = 0.99d+0*tolx
         if (kdegen .lt. 99999999) then
            tolinc = (tolx1 - tolx0) / kdegen
         else
            tolinc = zero
         end if
         featol  = tolx0
      else if (Task .eq. Cycle  .or.  Task .eq. Optml) then
*        ---------------------------------------------------------------
*        Task = 'E'nd of cycle or 'O'ptimal.
*        initialize local variables maxfix and tolz.
*        ---------------------------------------------------------------
         maxfix = 2
         tolz   = eps1
         if (Task .eq. Optml) then

*           Task = Optimal.
*           Return with inform = 0 if the last call was at the
*           same itn, or if there have already been maxfix calls
*           with the same state of feasibility.

            if (itnfix .eq. itn   ) return
            if (nInf   .gt.   0   ) then
               j = 1
            else
               j = 2
            end if
            if (nfix(j) .ge. maxfix) return
            nfix(j) = nfix(j) + 1
         end if

*        Set nonbasics on their nearest bound if they are within
*        the current featol of that bound.

         itnfix = itn

         do j = 1, nb
            if (hs(j) .le. 1  .or.  hs(j) .eq. 4) then
               b1    = bl(j)
               b2    = bu(j)
               d1    = abs( x(j) - b1 )
               d2    = abs( x(j) - b2 )
               if (d1 .gt. d2) then
                  b1   = b2
                  d1   = d2
               end if
               if (d1 .le. featol) then
                  if (d1 .gt. tolz) inform = inform + 1
                  x(j) = b1
               end if
            end if
         end do

*        Reset featol to its minimum value.

         featol = tolx0
         if (inform .gt. 0) then

*           The basic variables will be reset.

            if (PrtLvl .ge. 10) then
               write(str, 1000) itn, inform
               call snPRNT( 21, str, iw, leniw )
            end if
         end if
      end if

      return

 1000 format(' Itn', i7, ': Basics recomputed after ', i7,
     &       '  nonbasics set on bound')

      end ! subroutine s5dgen

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5FixX
     &   ( Task, b1, b2, tolx, hs, bl, bu, x )

      implicit
     &     none
      integer
     &     Task, b1, b2, hs(b2)
      double precision
     &     tolx, bl(b2), bu(b2), x(b2)

*     ==================================================================
*     s5FixX  ensures that variables satisfy their simple bounds.
*
*     If Task = xBound, variables x(b1) through x(b2) are made to
*                       satisfy their bounds.
*     If Task = xMove , variables x(b1) through x(b2) are made to
*                       satisfy their bounds. In addition, any nonbasic
*                       variable close to its bound is moved onto it.
*
*     29 Apr 1999: First version of s5FixX.
*     29 Apr 1999: Current version.
*     ==================================================================
      integer
     &     j
      double precision
     &     xj
*     ------------------------------------------------------------------
      integer            xBound,     xMove
      parameter         (xBound = 0, xMove = 1)
*     ------------------------------------------------------------------
      if (Task .eq. xBound) then
         do j = b1, b2
            x(j) = max( x(j), bl(j) )
            x(j) = min( x(j), bu(j) )
         end do
      else if (Task .eq. xMove) then
         do j  = b1, b2
            xj = max( x(j), bl(j) )
            xj = min( xj  , bu(j) )
            if (hs(j) .le. 1) then
               if (xj .le. bl(j) + tolx) xj = bl(j)
               if (xj .ge. bu(j) - tolx) xj = bu(j)
            end if
            x(j)   = xj
         end do
      end if

      end ! subroutine s5FixX

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5fixS
     &   ( Task, m, maxS, mBS, n, nb, nS, hs, kBS,
     &     bl, bu, blBS, buBS, x, xBS )

      implicit
     &     none
      integer
     &     Task, m, maxS, mBS, n, nb, nS, hs(nb), kBS(mBS)
      double precision
     &     bl(nb), bu(nb), x(nb), blBS(mBS), buBS(mBS), xBS(mBS)

!     ==================================================================
!     s5fixS   concerns temporary bounds on superbasic variables.
!     If Task = Fix,  s5fixS sets hs(j) = -1, 0, 1 or 4 for certain
!     superbasic variables.
!
!     If Task = Free, s5fixS changes -1 values to hs(j) = 2.
!
!     30 May 1995: First version of s5fixS.
!     12 Jul 2001: Current version.
!     ==================================================================
      integer
     &     j, k
!     ------------------------------------------------------------------
      integer            Fix,     Free
      parameter         (Fix = 0, Free = 1)
!     ------------------------------------------------------------------

      if (Task .eq. Fix) then
         if (nS .gt. 0) then
!           ------------------------------------------------------------
!           Change superbasic hs(j) to be temporarily fixed.
!           ------------------------------------------------------------
            nS = 0
            do j = 1, nb
               if (hs(j) .eq. 2) then
                  if (bl(j) .eq. bu(j)) then
                     hs(j) =  4
                  else if (x(j) .le. bl(j)) then
                     hs(j) =  0
                  else if (x(j) .ge. bu(j)) then
                     hs(j) =  1
                  else
                     hs(j) = -1
                  end if
               end if
            end do
         end if

      else if (Task .eq. Free) then
!        ---------------------------------------------------------------
!        Free the temporarily fixed structurals.
!        Load the superbasic variables/bounds into xBS, blBS, buBS.
!        ---------------------------------------------------------------
         j = 1
!+       while (j .le. n  .and.  nS .lt. maxS) do
  100    if    (j .le. n  .and.  nS .lt. maxS) then
            if (hs(j) .eq. -1) then
               nS      = nS + 1
               k       = m  + nS
               hs(j)   = 2
               xBS (k) = x(j)
               blBS(k) = bl(j)
               buBS(k) = bu(j)
               kBS (k) = j
            end if
            j  = j + 1
            go to 100
!+       end while
         end if
      end if

      end ! subroutine s5fixS

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5finE
     &   ( nBS, nb, nInfE, sInfE, featol, hEstat, kBS, bl, bu, xBS )

      implicit
     &     none
      integer
     &     nBS, nb, nInfE, hEstat(nb), kBS(nBS)
      double precision
     &     featol, sInfE, bl(nb), bu(nb), xBS(nBS)

*     ==================================================================
*     s5finE  is called at the end of a QP to compute the final number
*     and sum of infeasibilities for the elastic variables. The elastic
*     states of feasible variables are set to zero.
*
*     10 Jun 2001: First version of s5finE.
*     10 Jun 2001: Current version.
*     ==================================================================
      integer
     &     j, jEs, k
      double precision
     &     res
*     ------------------------------------------------------------------
      double precision   zero
      parameter         (zero   = 0.0d+0)
*     ------------------------------------------------------------------
      nInfE = 0
      sInfE = zero

      do k   = 1, nBS
         j   = kBS(k)
         jEs = hEstat(j)

         if (jEs .gt. 0) then
*           ------------------------------------------------------------
*           x(j) is predicted to violate an elastic bound.
*           ------------------------------------------------------------
            if (jEs .eq. 1) then

*              Elastic lower bound.

               res = bl(j)  - xBS(k)
               if (res .gt. featol) then
                  nInfE  =  nInfE + 1
                  sInfE  =  sInfE + res
               else
                  jEs = 0
               end if

            else if (jEs .eq. 2) then

*              Elastic upper bound.

               res = xBS(k) - bu(j)
               if (res .gt. featol) then
                  nInfE  =   nInfE + 1
                  sInfE  =   sInfE + res
               else
                  jEs = 0
               end if
            end if
         end if

         hEstat(j) = jEs

      end do

      end ! subroutine s5finE

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5getB
     &   ( iExit, Start, LPlog, needB, m, maxS, mBS,
     &     n, nb, nnCon, nnJac, nnObj, nName, nS, itQP, itQPmax, itn,
     &     nDegen, numLC, numLIQ, tolFP, tolQP, tolx,
     &     nInf, sInf, wtInf, iObj, sclObj, piNorm, rgNorm,
     &     ne, nlocA, locA, indA, Acol,
     &     hElast, hrtype, hfeas, hs, kBS, Names,
     &     Ascale, bl, bu, blBS, buBS, blSav, buSav,
     &     gBS, pi, rc, nrhs0, nrhs, rhs,
     &     lenx0, nx0, x0, x, xBS,
     &     iy, iy1, y, y1, y2, y3,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      external
     &     LPlog
      logical
     &     needB
      integer
     &     Start, iExit, iObj, itn, itQP, itQPmax, lenx0, nrhs0, nrhs,
     &     nx0, m, maxS, mBS, n, nb, ne, nlocA, nName, nnCon, nnJac,
     &     nnObj, nDegen, nInf, nS, numLC, numLIQ, lencw, leniw, lenrw,
     &     kBS(mBS), hrtype(nb), hfeas(mBS), hElast(nb), hs(nb),
     &     locA(nlocA), indA(ne), iy(nb), iy1(nb), iw(leniw)
      double precision
     &     tolx, tolFP, tolQP, sInf, sclObj, wtInf, piNorm, rgNorm,
     &     Acol(ne), Ascale(nb), bl(nb), bu(nb), blSav(nb),
     &     buSav(nb), rc(nb), x0(lenx0), x(nb), blBS(mBS), buBS(mBS),
     &     gBS(mBS), xBS(mBS), pi(m), rhs(nrhs0), y(nb),
     &     y1(nb), y2(nb), y3(nb), rw(lenrw)
      character
     &     Names(nName)*8, cw(lencw)*8

*     ==================================================================
*     s5getB   finds an initial basis kBS(1:m) for the linear
*     constraints and bounds.
*
*     On entry,
*     Start    is 0,1,2,3: an integer version of the solver's Start
*              (character for sqopt, snopt, snoptc, npopt
*               integer   for snopta).
*
*     In general, the constraints are scaled.
*     First, we attempt to find a  feasible point for the bounds and
*     general linear equalities. This may fail, so there is the chance
*     that the initial basis is not feasible.
*     This difficulty is taken care of later.
*
*      iExit         Status
*      -----         ------
*         -3         Too many iterations
*         -2         variable is unbounded (this should not happen
*         -1         LP is infeasible
*          0         LP solution found
*         >0         Fatal error
*
*     31 Jul 1996: First version of s5getB.
*     12 Jul 1997: Thread-safe version.
*     01 Aug 2003: snEXIT and snPRNT adopted.
*     08 Mar 2004: Implemented gotHes, gotScl for Hot starts.
*     03 Apr 2005: Current version of s5getB.
*     ==================================================================
      logical
     &     Elastc, needLU, needx
      character
     &     ProbTag*20, str*120
      integer
     &     j, iInsrt, iLoadB, iOldB, iStart, lc1,
     &     lCrash, lEmode, lsSave, lvlInf, minimz, MjrPrt, MnrPrt,
     &     nrhsLP0, nrhsLP, nnL, numLEQ, lvlScl,
     &     iCrash, iObjFP, subopt
      double precision
     &     ObjAdd, infBnd, tCrash
*     ------------------------------------------------------------------
      integer            COLD,       BASIS,      WARM,       HOT
      parameter         (COLD   = 0, BASIS  = 1, WARM   = 2, HOT    = 3)
      integer            RowTyp
      parameter         (Rowtyp = 0)
      integer            No
      parameter         (No     =-1)
      integer            Scale
      parameter         (Scale  = 0)
      integer            FPE
      parameter         (FPE    = 3)
      integer            Fix
      parameter         (Fix    = 0)
      integer            Intern,     Extern
      parameter         (Intern = 0, Extern = 1)
      integer            xBound
      parameter         (xBound = 0)

      double precision   zero,            one
      parameter         (zero   = 0.0d+0, one = 1.0d+0)

      parameter         (lvlScl =  75) ! scale option
      parameter         (iCrash =  88) ! Crash option
      integer            gotHes,       gotScl
      parameter         (gotHes = 231, gotScl = 232)
*     ------------------------------------------------------------------
      character          line*4
      data               line /'----'/
*     ------------------------------------------------------------------
      iLoadB = iw(122) ! load file
      iInsrt = iw(125) ! insert file
      iOldB  = iw(126) ! old basis file

      MjrPrt = iw( 92) ! Major print level
      MnrPrt = iw( 93) ! Minor print level

      infBnd = rw( 70) ! definition of an infinite bound
      tCrash = rw( 62) ! crash tolerance.

*     Initialize a few things.

      iExit      = 0
      minimz     = 1
      nnL        = max( nnObj, nnJac )
      nInf       = 0
      ObjAdd     = zero         ! Local to s5getB
      iObjFP     = 0            ! Used for FP calculation
      sclObj     = one

*     Count the linear equality and inequality constraints.

      numLEQ = 0
      do j = n+nnCon+1, nb
         if (bl(j) .eq. bu(j)) numLEQ = numLEQ + 1
      end do
      numLIQ = numLC - numLEQ

*     Initialize the row and column scales.

      if (iw(lvlScl) .gt. 0) then
         call dload
     &      ( nb, one, Ascale, 1 )
      end if

*     ==================================================================
*     Decode Start.
*     ==================================================================
      needB  = .true.

      if (Start .eq. COLD  .or.  Start .eq. BASIS) then
*        --------------------------------
*        Cold start  or  Basis file.
*        --------------------------------
         iStart = 0
         needB  = max( ioldB, iInsrt, iloadB ) .le. 0
         nS     = 0

      else if (Start .eq. WARM) then
*        --------------------------------
*        Warm start.
*        --------------------------------
         iStart = 1
         needB  = .false.

      else if (Start .eq. HOT ) then
*        --------------------------------
*        Hot start.
*        --------------------------------
         iStart = 1
         needB  = .false.
      end if

      call s1page( 1, iw, leniw )

      if (iStart .eq. 0) then
*        ---------------------------------------------------------------
*        Cold start, or Basis file provided.
*        Input a basis file if one exists, thereby defining hs and x.
*        (Otherwise, s2crsh will be called later to define hs.)
*        ---------------------------------------------------------------
*        Initialize x(n+1:nb) and pi(1:m) before the problem is scaled.
*        The basis files initialize all of x.
*        One day they may load pi for nonlinear problems.

         call dload ( m, zero, x(n+1), 1 )
         call dload ( m, zero, pi    , 1 )

         call snPRNT( 21, ' Initial basis', iw, leniw )
         call snPRNT( 21, ' -------------', iw, leniw )

         if ( needB ) then
            call snPRNT( 31, ' No basis file supplied', iw, leniw )

            if (iw(iCrash) .eq. 0) then
               needB  = .false.
               lCrash = 0
               call s2crsh
     &            ( lCrash, MnrPrt, m, n, nb, nnCon,
     &              iw(iCrash), tCrash,
     &              ne, nlocA, locA, indA, Acol,
     &              kBS, hs, hrtype, bl, bu, x,
     &              iw, leniw, rw, lenrw )
            end if
         else
            call s4getB
     &         ( iExit, m, n, nb, nName, nS, iObj,
     &           hs, bl, bu, x, Names, iw, leniw, rw, lenrw )
            if (iExit .gt. 0) then
*              Set blSav and buSav to avoid uninitialized copy.
               call dcopy
     &            ( nb, bl, 1, blSav, 1 )
               call dcopy
     &            ( nb, bu, 1, buSav, 1 )
               go to 900
            end if
         end if
      end if ! iStart = 0

*     ------------------------------------------------------------------
*     Move x inside its bounds.
*     ------------------------------------------------------------------
      call s5FixX
     &   ( xBound, 1, n, tolx, hs, bl, bu, x )

*     ------------------------------------------------------------------
*     Scale the linear part of the constraints.
*     (Any nonlinear elements in A contain fake nonzeros.)
*     ------------------------------------------------------------------
      if (iw(lvlScl) .gt. 0  .and.  numLC .gt. 0) then
         if (iw(gotScl) .eq. 0) then
            iw(gotScl) = 1
            lsSave     = iw(lvlScl)
            iw(lvlScl) = 2
            call s2scal
     &         ( MjrPrt, m, n, nb, nnL, nnCon, nnJac, hrtype,
     &           ne, nlocA, locA, indA, Acol,
     &           Ascale, bl, bu, y, y2,
     &           iw, leniw, rw, lenrw )
            iw(lvlScl) = lsSave
         end if

         call s2scla
     &      ( Scale, m, n, nb, iObj, infBnd, sclObj,
     &        ne, nlocA, locA, indA, Acol,
     &        Ascale, bl, bu, pi, x )
      end if

*     Save the bounds.
*     If required, save the scaled initial point in x0.

      call dcopy
     &   ( nb, bl, 1, blSav, 1 )
      call dcopy
     &   ( nb, bu, 1, buSav, 1 )
      if (nx0 .gt. 0)
     &call dcopy
     &   ( nx0, x, 1, x0   , 1 )

*     ------------------------------------------------------------------
*     Prepare to get feasible for the linear constraints.
*     Relax any nonlinear rows.
*     Allow the linear rows to be elastic.
*     ------------------------------------------------------------------
      if (nnCon .gt. 0  .and.  numLC .gt. 0) then
         call dload
     &      ( nnCon, (-infBnd), bl(n+1), 1 )
         call dload
     &      ( nnCon,   infBnd , bu(n+1), 1 )
      end if

*     ------------------------------------------------------------------
*     Compute a starting basis for the linear constraints.
*     This means attempting to get feasible for the linear equalities.
*     If the E rows are infeasible, s5LP  will take care of it.
*     ------------------------------------------------------------------
      if (MnrPrt .gt. 10) then
         write(str, 1332) (line, j=1,28)
         call snPRNT( 31, str, iw, leniw )
         write(str, 1316) (line, j=1,16)
         call snPRNT( 32, str, iw, leniw )
      end if

      if ( needB ) then
*        ---------------------------------------------------------------
*        Crash is needed to find a basis.
*        ---------------------------------------------------------------
*        Treat Crash 1 the same as Crash 2.

         iw(iCrash)  = max( iw(iCrash), 2 )
         iw(iCrash)  = min( iw(iCrash), 3 )

         if (numLC .gt. 0) then
*           ============================================================
*           Find a feasible point for the linear constraints.
*           ============================================================
*           Crash 2 finds a basis for ALL the LINEAR constraints
*                   (equalities and inequalities).
*           Crash 3 treats linear EQUALITIES separately.
*           ------------------------------------------------------------
            if (iw(iCrash) .eq. 2) then

*              Relax.

            else if (iw(iCrash) .eq. 3) then

               if (numLEQ .gt. 0) then

*                 Find a basis for the linear EQUALITIES.
*                 Linear inequality rows are made to appear to be free.

                  if (MnrPrt .ge. 10) then
                     write(str, 2100) itn
                     call snPRNT( 23, str, iw, leniw )
                  end if
               end if ! numLEQ > 0
            end if ! iCrash = 2 or 3
         end if ! numLC > 0

*        Call Crash even if there are no linear constraints.
*        We haven't done any Solve yet, so we should NOT need
*        call s5hs  ( Extern, nb, bl, bu, hs, x ) .

         lCrash = iw(iCrash)
         call s2crsh
     &      ( lCrash, MnrPrt, m, n, nb, nnCon,
     &        iw(iCrash), tCrash,
     &        ne, nlocA, locA, indA, Acol,
     &        kBS, hs, hrtype, bl, bu, x,
     &        iw, leniw, rw, lenrw )
      end if ! needB

*     ==================================================================
*     1. Set nS to match hs(*).
*     2. Set kBS(m+1) thru kBS(m+nS) to define the initial superbasics.
*     3. Set all nonbasic x to be within bounds, which may change some
*        hs values from 0 to 1.
*     4. Set nonbasic x to be exactly on nearly satisfied bounds.
*        (Some nonbasics may still be between bounds.)
*     ==================================================================
      call s4chek
     &   ( m, maxS, mBS, n, nb, needB, iw(gotHes),
     &     nS, iObj, hs, kBS, bl, bu, x, iw, leniw, rw, lenrw )

      if (needB  .and.  numLC .gt. 0) then
                                ! Fix SBs. This forces simplex steps.
         call s5fixS
     &      ( Fix, m, maxS, mBS, n, nb, nS, hs, kBS,
     &        bl, bu, blBS, buBS, x, xBS )

         if (numLEQ .gt. 0) then
            if (numLIQ .eq. 0) then
               nrhsLP = 0
            else                ! Relax the linear INEQUALITIES.
               call s5LG
     &            ( m, n, nb, nnCon, nrhsLP,
     &              ne, nlocA, locA, indA, Acol,
     &              bl, bu, nrhs0, nrhs, rhs, y3, x, y, rw, lenrw )
            end if

            ProbTag = 'linear equality rows'
            nrhsLP0  =  max(nrhsLP, 1)
            Elastc  = .false.
            lvlInf  =  0
            lEmode  =  0           ! Suspend Elastic mode for equalities
            needLU  = .true.
            needx   =  needLU
            subopt  =  No
            call s5hs
     &         ( Intern, nb, bl, bu, hs, x )
            call s5LP
     &         ( iExit, FPE, ProbTag, Elastc, subopt,
     &           LPlog, needLU, needx,
     &           m, n, nb, nDegen, itQP, itQPmax, itn,
     &           lEmode, lvlInf, MnrPrt,
     &           minimz, iObjFP, sclObj, ObjAdd, tolFP, tolQP, tolx,
     &           nInf, sInf, wtInf, piNorm, rgNorm,
     &           ne, nlocA, locA, indA, Acol,
     &           hElast, hrtype, hfeas, hs, kBS,
     &           Ascale, bl, bu, blBS, buBS,
     &           gBS, pi, rc, nrhsLP0, nrhsLP, y3, x, xBS, x,
     &           iy, iy1, y, y1, y2,
     &           cw, lencw, iw, leniw, rw, lenrw )

*           iExit values are = -3,-2,-1, 0, >0

*           Reset the bounds on the linear inequality rows.

            if (numLIQ .gt. 0) then
               lc1    = n + nnCon + 1
               call dcopy
     &            ( numLC, blSav(lc1), 1, bl(lc1), 1 )
               call dcopy
     &            ( numLC, buSav(lc1), 1, bu(lc1), 1 )
               call s5hs
     &            ( Intern, nb, bl, bu, hs, x )
            end if ! numLIQ > 0

            if (MnrPrt .ge. 10) then
               if (nInf .eq. 0) then

*                 The linear E rows are now satisfied.

                  write(str, 2200) itn
                  call snPRNT( 21, ' ', iw, leniw )
                  call snPRNT( 23, str, iw, leniw )
               else

*                 The linear E rows are infeasible.

                  write(str, 2300) itn
                  call snPRNT( 23, str, iw, leniw )
               end if
            end if
         end if ! numLEQ > 0

         if (iExit .eq. 0) then
*           ------------------------------------------------------------
*           Include the linear INEQUALITIES.
*           ------------------------------------------------------------
            if (numLIQ .gt. 0) then
               if (MnrPrt .ge. 10) then
                  if (numLEQ .gt. 0) then
                     write(str, 2400) itn
                     call snPRNT( 23, str, iw, leniw )
                  else
                     write(str, 2410) itn
                     call snPRNT( 23, str, iw, leniw )
                  end if
               end if

               call s5hs
     &            ( Extern, nb, bl, bu, hs, x )
               call s2Amat
     &            ( RowTyp, MnrPrt, m, n, nb,
     &              nnCon, nnJac, nnObj, iObj,
     &              ne, nlocA, locA, indA, Acol,
     &              bl, bu, hrtype,
     &              iw, leniw, rw, lenrw )
               lCrash = 4
               call s2crsh
     &            ( lCrash, MnrPrt, m, n, nb, nnCon,
     &              iw(iCrash), tCrash,
     &              ne, nlocA, locA, indA, Acol,
     &              kBS, hs, hrtype, bl, bu, x,
     &              iw, leniw, rw, lenrw )
            end if ! numLIQ > 0
         end if
      end if ! needB and numLC > 0

  900 return

 1316 format(1x, 16a4)
 1332 format(1x, 28a4)
 2100 format(' Itn', i7, ': Making linear equality rows feasible')
 2200 format(' Itn', i7, ': Feasible linear equality rows')
 2300 format(' Itn', i7, ': Infeasible linear equality rows')
 2400 format(' Itn', i7, ': Making all linear rows feasible')
 2410 format(' Itn', i7, ': Making the linear rows feasible')

      end ! subroutine s5getB

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5grdE
     &   ( nb, nBS, wtInf, hEstat, kBS, gBS )

      implicit
     &     none
      integer
     &     nb, nBS, hEstat(nb), kBS(nBS)
      double precision
     &     wtInf, gBS(nBS)

*     ==================================================================
*     s5grdE  is called when elastic variables are allowed to violate
*     their true bounds bl and bu.  It updates the gradient gBS
*     to include the gradient of the penalty term.
*
*     On exit,
*     gBS(nBS)    is the rhs for the equations for pi.
*
*     08 Oct 1996: First version of s5grdE.
*     21 Apr 1999: Current version.
*     ==================================================================
      integer
     &     j, jEs, k
*     ------------------------------------------------------------------
      do k = 1, nBS
         j   = kBS(k)
         jEs = hEstat(j)
         if      (jEs .eq. 0) then
*           Relax
         else if (jEs .eq. 1) then
            gBS(k) = gBS(k) - wtInf
         else !  (jEs .eq. 2)
            gBS(k) = gBS(k) + wtInf
         end if
      end do

      end ! subroutine s5grdE

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5hs
     &   ( Mode, nb, bl, bu, hs, x )

      implicit
     &     none
      integer
     &     Mode, nb, hs(nb)
      double precision
     &     bl(nb), bu(nb), x(nb)

*     ==================================================================
*     s5hs   sets the state vector hs.
*     if mode = 'Internal', s5hs sets hs(j) = -1 or 4 for certain
*        nonbasic variables.  This allows s5pric to operate more
*        efficiently.  The internal values of hs are now as follows:
*
*        hs(j) = -1  Nonbasic between bounds (bl     <  x <  bu    )
*        hs(j) =  0  Nonbasic on lower bound (bl-tol <  x <= bl    )
*        hs(j) =  1  Nonbasic on upper bound (bu     <= x <  bu+tol)
*        hs(j) =  2  Superbasic
*        hs(j) =  3  Basic
*        hs(j) =  4  Nonbasic and fixed      (bl     = x  =  bu    )
*
*        where 0 <= tol < the feasibility tolerance.
*
*     if mode = 'External', s5hs changes -1 or 4 values to hs(j) = 0,
*        ready for basis saving and the outside world.
*
*     08 Apr 1992: First version of s5hs.
*     21 Aug 1999: Current version.
*     ==================================================================
      integer
     &     j
*     ------------------------------------------------------------------
      integer            Intern,     Extern
      parameter         (Intern = 0, Extern = 1)
*     ------------------------------------------------------------------
      if (Mode .eq. Intern) then
*        ---------------------------------------------------------------
*        Change nonbasic hs(j) to internal values (including 4 and -1).
*        This may change existing internal values if bl and bu have been
*        changed -- e.g. at the start of each major iteration.
*        ---------------------------------------------------------------
         do j = 1, nb
            if (hs(j) .le. 1) then
               if (bl(j) .eq. bu(j)) then
                  hs(j) =  4
               else if (x(j) .le. bl(j)) then
                  hs(j) =  0
               else if (x(j) .ge. bu(j)) then
                  hs(j) =  1
               else
                  hs(j) = -1
               end if
            end if
         end do

      else if (Mode .eq. Extern) then
*        ---------------------------------------------------------------
*        Change hs to external values.
*        Some nonbasic hs(j) may be 4 or -1.  Change them to 0.
*        ---------------------------------------------------------------
         do j = 1, nb
            if (hs(j) .eq. 4  .or.  hs(j) .eq. -1) hs(j) = 0
         end do
      end if

      end ! subroutine s5hs

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5Inf
     &   ( nBS, featol, nInf, sInf, hfeas, blBS, buBS, gBS, xBS )

      implicit
     &     none
      integer
     &     nBS, nInf, hfeas(nBS)
      double precision
     &     featol, sInf, blBS(nBS), buBS(nBS), gBS(nBS), xBS(nBS)

*     ==================================================================
*     s5Inf  computes the sum and number of infeasibilities.
*
*     hfeas     x(j)                             Meaning
*     -----     -----                             -------
*      -2   infeasible                              x(j) .le. bl(j)-tol
*       0     feasible               bl(j)-tol .le. x(j) .le. bu(j)+tol
*      +2   infeasible                              x(j) .ge. bu(j)+tol
*
*     On exit,
*     nInf        is the number violated non-elastic bounds.
*     gBS(nBS)    is the rhs for the equations for pi.
*
*     29 Oct 1993: First version of s5Inf.
*     21 Aug 1999: Current version.
*     ==================================================================
      integer
     &     k, numInf
      double precision
     &     res, sumInf, xk
*     ------------------------------------------------------------------
      double precision   zero,          one
      parameter         (zero = 0.0d+0, one  = 1.0d+0)
*     ------------------------------------------------------------------
*     Find the number and sum of the non-elastic infeasibilities.
*     ------------------------------------------------------------------
      numInf = 0
      sumInf = zero

*     ------------------------------------------------------------------
*     Set hfeas, used in s5step.
*     ---------------------------------------------------------------
      do 100,    k = 1, nBS
         hfeas(k)  = 0
         xk        = xBS (k)
         res       = blBS(k) - xk
         if (res .le. featol) go to 50
         gBS  (k)  = - one
         hfeas(k)  = - 2
         go to 60

   50    res       = xk - buBS(k)
         if (res .le. featol) go to 100
         gBS  (k)  =   one
         hfeas(k)  =   2

   60    numInf    = numInf + 1
         sumInf    = sumInf + res
  100 continue

      sInf  = sumInf
      nInf  = numInf

      end ! subroutine s5Inf

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5InfE
     &   ( nb, nBS, hEstat, kBS, nInfE, sInfE, bl, bu, x )

      implicit
     &     none
      integer
     &     nb, nBS, nInfE, hEstat(nb), kBS(nBS)
      double precision
     &     sInfE, bl(nb), bu(nb), x(nb)

*     ==================================================================
*     s5InfE  computes the sum and number of elastic infeasibilities.
*
*     On exit,
*     nInfE is the number of elastics allowed to go infeasible.
*
*     20 Aug 1996: First version of s5InfE.
*     21 Aug 1999: Current version.
*     ==================================================================
      integer
     &     j, jEs, k, numInf
      double precision
     &     res, sumInf
*     ------------------------------------------------------------------
      double precision   zero
      parameter         (zero = 0.0d+0)
*     ------------------------------------------------------------------
*     Find the number and sum of the elastic infeasibilities.
*     ------------------------------------------------------------------
      numInf = 0
      sumInf = zero

      do k = 1, nBS
         j   = kBS(k)
         jEs = hEstat(j)

         if      (jEs .eq. 0) then

*           Relax

         else if (jEs .eq. 1) then
            numInf = numInf + 1
            res    = bl(j) - x(j)

            if (res .gt. zero) then
               sumInf = sumInf + res
            end if
         else if (jEs .eq. 2) then
            numInf = numInf + 1
            res    = x(j) - bu(j)

            if (res .gt. zero) then
               sumInf = sumInf + res
            end if
         end if
      end do

      sInfE = sumInf
      nInfE = numInf

      end ! subroutine s5InfE

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5IniE
     &   ( nb, nBS, nElast, featol, infBnd,
     &     hElast, hEstat, kBS,
     &     bl, bu, blBS, buBS, xBS )

      implicit
     &     none
      integer
     &     nb, nBS, nElast, hElast(nb), hEstat(nb), kBS(nBS)
      double precision
     &     featol, infBnd, bl(nb), bu(nb), blBS(nBS), buBS(nBS),
     &     xBS(nBS)

*     ==================================================================
*     s5IniE  initializes hEstat and adjusts the temporary upper and
*     lower bounds for the elastic variables.
*
*     s5IniE is called at the start of elastic mode, i.e., when  Elastc
*     is set  true.  In elastic mode, the elastic variables are allowed
*     to violate their true bounds bl and bu.  The bounds blBS and buBS
*     are redefined for the basic and superbasic elastic variables.
*
*     08 Aug 1996: First version of s5IniE.
*     16 Aug 2000: Current version.
*     ==================================================================
      integer
     &     j, je, k
      double precision
     &     res
*     ------------------------------------------------------------------
      call iload
     &   ( nb, 0, hEstat, 1 )
      nElast = 0

      do   k = 1, nBS
         j   = kBS(k)
         je  = hElast(j)

         if (je .eq. 1  .or.  je .eq. 3) then
            res = bl(j) - xBS(k)
            if (res .gt. featol) then
               nElast    =   nElast + 1
               hEstat(j) =   1
               blBS(k)   = - infBnd
               buBS(k)   =   bl(j)
            end if
         end if

         if (je .eq. 2  .or.  je .eq. 3) then
            res = xBS(k) - bu(j)
            if (res .gt. featol) then
               nElast    =   nElast + 1
               hEstat(j) =   2
               blBS(k)   =   bu(j)
               buBS(k)   =   infBnd
            end if
         end if
      end do

      end ! subroutine s5IniE

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5LG
     &   ( m, n, nb, nnCon, nrhsLP,
     &     ne, nlocA, locA, indA, Acol,
     &     bl, bu, nrhs0, nrhs, rhs, rhsLP, x, y, rw, lenrw )

      implicit
     &     none
      integer
     &     lenrw, m, n, nb, ne, nlocA, nnCon, nrhs0, nrhs, nrhsLP,
     &     locA(nlocA), indA(ne)
      double precision
     &     Acol(ne), bl(nb), bu(nb), x(nb),
     &     rhs(nrhs0), rhsLP(m), y(m), rw(lenrw)

*     ==================================================================
*     s5LG  relaxes the linear inequality constraints ready for a call
*     to s5LP that gets feasible for the linear equality constraints.
*     A rhs is computed that makes the relaxed rows satisfied
*     at the current x.   Then x will not be disturbed more than
*     necessary during a Warm start.
*
*     02 Apr 2005: First version of s5LG.
*     02 Apr 2005: Current version of s5LG.
*     ==================================================================
      integer
     &     i, j
      double precision
     &     eps0, infBnd
*     ------------------------------------------------------------------
      integer            Normal
      parameter         (Normal = 0)
      double precision   zero,            one
      parameter         (zero   = 0.0d+0, one = 1.0d+0)
*     ------------------------------------------------------------------
      eps0   = rw(  2) ! eps**(4/5)
      infBnd = rw( 70) ! definition of an infinite bound

      nrhsLP = m
      call s2Aprd
     &   ( Normal, eps0, ne, nlocA, locA, indA, Acol,
     &     one, x, n, zero, y, m )
      call daxpy
     &   ( nrhsLP, (-one), x(n+1), 1, y, 1 )
      if (nrhs .gt. 0)
     &   call daxpy
     &      ( nrhs, (-one), rhs, 1, y, 1 )
      call dload
     &     ( nrhsLP, zero, rhsLP, 1 )

      do i = nnCon+1, m
         j = n + i
         if (bl(j) .lt. bu(j)) then
            bl(j)    = - infBnd
            bu(j)    = + infBnd
            rhsLP(i) =   y(i)
         end if
      end do

      end ! subroutine s5LG

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5LPit
     &   ( iExit, feasbl, incres, needpi, Elastc,
     &     m1, m, nb, nDegen, LUreq,
     &     kp, jBq, jSq, jBr, jSr, jq,
     &     featol, pivot, step, tolinc,
     &     hElast, hEstat, hfeas, hs, kBS,
     &     bl, bu, blBS, buBS,
     &     x, xBS, pBS, y1,
     &     iw, leniw, rw, lenrw )

      implicit
     &     none
      logical
     &     feasbl, incres, needpi, Elastc
      integer
     &     iExit, m1, m, nb, nDegen, LUreq, kp, jBq, jSq, jBr, jSr, jq,
     &     leniw, lenrw, hfeas(m1), kBS(m1), hElast(nb), hEstat(nb),
     &     hs(nb), iw(leniw)
      double precision
     &     featol, pivot, step, tolinc,
     &     bl(nb), bu(nb), x(nb), blBS(m1), buBS(m1), xBS(m1),
     &     pBS(m1), y1(m1), rw(lenrw)

*     ==================================================================
*     s5LPit computes one step of the primal simplex method.
*     jq is the variable entering the basis and djq is its reduced cost.
*
*      iExit       Status
*      -----       ------
*        0         Normal exit
*        1         LP is unbounded
*
*     10 Sep 1991: First version based on Minos routine m5lpit.
*     02 Aug 1996: First min sum version added by PEG.
*     12 Jul 1997: Thread-safe version.
*     01 Aug 2003: snPRNT adopted.
*     27 Dec 2003: Current version of s5LPit.
*     =================================================================
      character
     &     str*50
      logical
     &     hitlow, move, onbnd, unbndd
      integer
     &     inform, infpiv, jEs, jr, jrstat, js, LUmod
      double precision
     &     bound, exact, bigdx, infBnd, stepP, stepmx, tolpiv
*     ------------------------------------------------------------------
      integer            xBStox
      parameter         (xBStox = 1)
      double precision   zero,            one
      parameter         (zero   = 0.0d+0, one  = 1.0d+0)
      parameter         (LUmod  = 216) ! number of LU mods
*     ------------------------------------------------------------------
      tolpiv    = rw( 60) ! excludes small elements of pBS.
      infBnd    = rw( 70) ! definition of an infinite bound
      bigdx     = rw( 72) ! unbounded step.

*     s5step assumes that the first m1 components of xBS can move.

      jSq       =     jq             ! Entering superbasic
      jSr       =     jq             ! leaving  superbasic
      hfeas(m1) =      0
      kBS  (m1) =     jq
      xBS  (m1) =  x(jq)
      blBS (m1) =  bl(jq)
      buBS (m1) =  bu(jq)

*     ==================================================================
*     Set hEstat(jq) and the elastic parts of blBS and buBS.
*     ==================================================================
      if ( Elastc ) then

*        If the new superbasic is an elastic variable
*        and it wants to move infeasible, set its elastic state.

         if (hElast(jq) .gt. 0) then
            js  = hs(jq)
            if ( incres ) then
               if (js .eq. 1  .or.  js .eq. 4) then
                  hEstat(jq) =   2
                  buBS(m1)   =   infBnd
                  if ( feasbl ) then
                     blBS(m1) = bu(jq)
                  end if
               end if
            else
               if (js .eq. 0  .or.  js .eq. 4) then
                  hEstat(jq) =   1
                  blBS(m1)   = - infBnd
                  if ( feasbl ) then
                     buBS(m1) = bl(jq)
                  end if
               end if
            end if
         end if
      end if ! Elastic mode

*     ==================================================================
*     Select a variable to be dropped from B.
*     s5step  uses the (m+1)th element of  blBS, buBS, xBS and pBS.
*     ==================================================================
      if (incres) then
         call dscal
     &      ( m, (- one), pBS, 1 )
         pBS(m1) =   one
      else
         pBS(m1) = - one
      end if

      stepmx = bigdx
      call s5step
     &   ( m1, nDegen,
     &     featol, infBnd, stepmx, tolinc, tolpiv,
     &     hfeas, blBS, buBS, xBS, pBS,
     &     hitlow, move, onbnd, unbndd,
     &     infpiv, kp, bound, exact, step, stepP )

      if (.not. unbndd) then
*        ---------------------------------------------------------------
*        Update the basic variables xBS and copy them into x.
*        ---------------------------------------------------------------
         jr      = kBS(kp)

         call daxpy
     &      ( m1, step, pBS, 1, xBS, 1 )
         call s5BSx
     &      ( xBStox, m1, nb, kBS, x, xBS )

c$$$!        10 Mar 2004: Care is needed to prevent the
c$$$!        new nonbasic variable jr from ending up slightly inside
c$$$!        its bound.  EXPAND normally ensures that x(jr) will be
c$$$!        ON or slightly OUTSIDE its bound, but now we realise that
c$$$!        rounding error might make it slightly INSIDE.
c$$$
c$$$         if (onbnd) then
c$$$            x(jr) = bound
c$$$         else if (hitlow) then
c$$$            x(jr) = min( x(jr), bl(jr) )
c$$$         else
c$$$            x(jr) = max( x(jr), bu(jr) )
c$$$         end if

         if (onbnd) x(jr) = bound

         if (kp .eq. m1) then
*           ------------------------------------------------------------
*           Variable jq reaches its opposite bound.
*           ------------------------------------------------------------
            if (incres) then
               hs(jq)  = 1
            else
               hs(jq)  = 0
            end if
            hfeas(kp)  = 0
            pivot      = zero
            if (.not. feasbl  .or.  Elastc) needpi = .true.

         else
*           ------------------------------------------------------------
*           Variable jq replaces the kp-th variable of  B.
*           It could be a fixed variable, whose new state must be 4.
*           ------------------------------------------------------------
            needpi =   .true.
            jBq    =   jq
            jBr    =   jr
            hs(jq) =   3
            pivot  = - pBS(kp)

            jEs    = hEstat(jr)
            hEstat(jr) = 0

            if (jEs .eq. 0) then
               if (blBS(kp) .eq. buBS(kp)) then
                  jrstat = 4
               else if (hitlow) then
                  jrstat = 0
               else
                  jrstat = 1
               end if

            else if (jEs .eq. 1) then
               if (bl(jr) .eq. bu(jr)) then
                  jrstat =   4
               else if (onbnd) then
                  jrstat =   0
               else if (x(jr) .lt. bu(jr)) then
                  jrstat = - 1
               else
                  jrstat =   1
               end if

            else !   jEs .eq. 2
               if (bl(jr) .eq. bu(jr)) then
                  jrstat =   4
               else if (onbnd) then
                  jrstat =   1
               else if (x(jr) .gt. bl(jr)) then
                  jrstat = - 1
               else
                  jrstat =   0
               end if
            end if

            hs(jr)   = jrstat
            kBS(kp)  = jq
            xBS(kp)  =  xBS(m1)
            blBS(kp) = blBS(m1)
            buBS(kp) = buBS(m1)

*           Update the LU factors.

            iw(LUmod)  = iw(LUmod) + 1
            call s2Bmod
     &         ( inform, kp, m, y1, iw, leniw, rw, lenrw )
            if (inform .eq. -1) LUreq = 5 ! Singular after LU mod
            if (inform .eq.  2) LUreq = 6 ! Unstable LU mod
            if (inform .eq.  7) LUreq = 7 ! Insufficient free memory.
         end if ! kp ne m1

      else

*        The solution is apparently unbounded.

         if (incres) then
            write(str, 1000) jq
         else
            write(str, 1100) jq
         end if
         call snPRNT( 21, str, iw, leniw )
         iExit = 1              ! Unbounded direction
      end if

      return

 1000 format(' Variable', i6, '  can increase indefinitely')
 1100 format(' Variable', i6, '  can decrease indefinitely')

      end ! subroutine s5LPit

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5pric
     &   ( Elastc, feasbl, incres, gotg, subopt,
     &     itn, m, n, nb, leng, neg, nnL,
     &     nS, nFreez, nonOpt, weight, sgnObj, piNorm,
     &     jq, djq, kPrc, toldj,
     &     ne, nlocA, locA, indA, Acol,
     &     hElast, hs, bl, bu, g, pi, rc, x, xFreez,
     &     iw, leniw, rw, lenrw )

      implicit
     &     none
      logical
     &     Elastc, feasbl, gotg, incres
      integer
     &     m, n, nb, ne, nlocA, nFreez, leng, neg, nnL, nS,
     &     nonOpt, itn, jq, kPrc, leniw, lenrw, subopt,
     &     hElast(nb), hs(nb), locA(nlocA), indA(ne), iw(leniw)
      double precision
     &     weight, piNorm, sgnObj, djq, toldj(3), Acol(ne), bl(nb),
     &     bu(nb), g(leng), pi(m), rc(nb), x(nb), xFreez(nb), rw(lenrw)

*     ==================================================================
*     s5pric  selects a nonbasic variable to enter the basis,
*     using the reduced gradients  dj = g(j) - pi'*a(j).
*
*     This version does partial pricing on both structurals and slacks.
*     Dynamic tolerances are used if partial price is in effect.
*
*     Partial pricing here means sectional pricing, because the three
*     blocks of  (A  -I)  are sliced up into nParPr sections
*     of equal size.  (The last section of each may be a little bigger,
*     since nParPr is unlikely to divide evenly into  n  or  m.)
*
*     input    g      = gradient for nonlinear variables.
*              pi     = pricing vector.
*              kPrc   = the no. of the section where  s5pric  last found
*                       a useful dj.
*                       (kPrc = 0 at the start of each major iteration.)
*              toldj(1-2) hold the current told if partial pricing, for
*                       phase 1 and 2 respectively.  told is used to
*                       determine if a dj is significant.
*              toldj(3) holds the specified optimality tolerance.
*              biggst   keeps track of the biggest dj found during the
*                       last scan of all sections of ( A  -I ).
*
*     output   kPrc   = the last section scanned.
*              nonOpt = the no. of useful djs found in that section.
*              jq     = best column found.
*              djq    = best dj.
*              toldj(1-2) save the current told if partial pricing.
*              incres   says if variable jq should increase or decrease.
*
*     In the code below,
*     the next section of  A  contains nPr1 structurals (j1+1 thru k1),
*     the next section of -I  contains nPr2 slacks      (j2+1 thru k2).
*     If  nParPr  is rather large, either nPr1 or nPr2 could be zero,
*     but not both.
*
*     If subopt > 0,  variables that haven't moved are not
*     priced, thereby limiting the number of superbasic variables.
*
*     ------------------------------------------------------------------
*     09 Aug 1992: First version of s5pric based on Minos 5.4 m5pric.
*     29 Jul 1996: Multiple pricing removed.
*     05 Aug 1996: First version with Elastic mode.
*     12 Jul 1997: Thread-safe version.
*     22 Dec 1999: subopt implemented.
*     31 Mar 2001: Free variables with hs(j)=-1 eligible for the basis.
*     01 Apr 2001: Current version of s5pric.
*     ==================================================================
      character
     &     str*100
      integer
     &     j, je, jj, js, jSlk, j1, j2, jfree, k1, k2,
     &     kPsav, lprDbg, lvldj, np, nParPr, nPrc, nParP, nPr1, nPr2
      double precision
     &     b1, b2, d, d1, d2, dj, dj1, dj2, djmax, infBnd, told, tolmin
*     ------------------------------------------------------------------
*-->  parameter         (zero = 0.0d+0, reduce = 0.25d+0)
      double precision   zero,          reduce
      parameter         (zero = 0.0d+0, reduce = 0.2d+0 )
*     ------------------------------------------------------------------
      infBnd    = rw( 70) ! definition of an infinite bound
      lprDbg    = iw( 85) ! > 0    => private debug print
      nParPr    = iw( 94) ! # of partial pricing sections

      djmax     = - infBnd
      djq       =   zero

      jq        = 0
      jfree     = 0
      nFreez    = 0
      nonOpt    = 0
      nPrc      = 0
      nParP     = nParPr
      nPr1      = n  / nParP
      nPr2      = m  / nParP
      if (max( nPr1, nPr2 ) .le. 0) nParP = 1

*     Set the tolerance for a significant dj.

      tolmin    = toldj(3) * piNorm
      if ( feasbl ) then
         lvldj  = 2
      else
         lvldj  = 1
      end if
      told      = toldj(lvldj)
      if (nParP .eq. 1) told = tolmin

*     Set pointers to the next section of  A  and  -I.
*     nPrc counts how many sections have been scanned in this call.
*     kPrc keeps track of which one to start with.

  100 nPrc = nPrc + 1
      kPrc = kPrc + 1
      if (kPrc .gt. nParP) kPrc = 1

      nPr1 = n  / nParP
      j1   =      (kPrc - 1)*nPr1
      k1   = j1 + nPr1
      if (kPrc .eq. nParP) k1 = n
      nPr1 = max( 0, k1-j1 )

      nPr2 = m  / nParP
      j2   = n  + (kPrc - 1)*nPr2
      k2   = j2 + nPr2
      if (kPrc .eq. nParP) k2 = nb
      nPr2 = max( 0, k2-j2 )

*     ------------------------------------------------------------------
*     Main loops for partial pricing (or full pricing).
*     Compute reduced costs rc(*)
*     for the kPrc-th section of structurals
*     and the kPrc-th section of slacks.
*     ------------------------------------------------------------------
      call s5rc
     &   ( j1+1, k1, gotg, m, n, leng, neg, sgnObj,
     &     ne, nlocA, locA, indA, Acol,
     &     hs, g, pi, rc )

      do j = j2+1, k2
         rc(j)  = pi(j-n)
      end do

*     ------------------------------------------------------------------
*     Main loop for pricing structural and slack reduced costs.
*     dj is rc(j), the reduced cost.
*     d  is -dj or +dj, depending on which way x(j) can move.
*     We are looking for the largest d (which will be positive).
*     ------------------------------------------------------------------
      np   = nPr1 + nPr2
      j    = j1
      jSlk = nPr1 + 1

      do jj = 1, np
         if (jj .eq. jSlk) j = j2
         j       = j + 1
         js      = hs(j)

         if (js .le. 1) then
            dj     = rc(j)

            if      (js .eq. 0) then
*              xj  is allowed to increase.
               d      = - dj
            else if (js .eq. 1) then
*              xj  is allowed to decrease.
               d      =   dj
            else
*              js is -1.
*              xj  is free to move either way.
*              Remember him as jfree in case he is the only one.
               d      = abs( dj )
               jfree  = j
            end if

            if (subopt .gt. 0) then
               if (x(j) .eq. xFreez(j)) then
                  if (d  .gt. told) then
                     nFreez = nFreez + 1
                     d      = zero
                  end if
               end if
            end if

*           See if this dj is significant.
*           Also see if it is the biggest dj so far.

            if (d  .gt. told) nonOpt = nonOpt + 1
            if (djmax .lt. d) then
               djmax  = d
               djq    = dj
               jq     = j
               kPsav  = kPrc
            end if
         end if
      end do

      if ( Elastc ) then
*        ---------------------------------------------------------------
*        Scan this section again, looking for nonbasic elastics.
*        ---------------------------------------------------------------
*        Compute reduced costs rc(j) for fixed nonbasic columns.

         call s5rcE
     &      ( j1+1, k1, gotg, m, n, leng, neg, sgnObj,
     &        ne, nlocA, locA, indA, Acol,
     &        hElast, hs, g, pi, rc )

         j    = j1
         do jj = 1, np
            if (jj .eq. jSlk) j = j2
            j      = j + 1
            je     = hElast(j)

            if (je .gt. 0) then
               js  = hs(j)
               dj  = rc(j)

               if      (js .eq. 0) then
*                 ------------------------------------------------------
*                 Nonbasic at its lower bound.
*                 An elastic xj can decrease through the bound.
*                 ------------------------------------------------------
                  if (je .eq. 1  .or.  je .eq. 3) then
                     dj  =   dj - weight
                     d   =   dj
                  end if

               else if (js .eq. 1) then
*                 ------------------------------------------------------
*                 Nonbasic at its upper bound.
*                 The default is to allow xj to decrease.
*                 However, an elastic xj can increase through the bound.
*                 ------------------------------------------------------
                  if (je .eq. 2  .or.  je .eq. 3) then
                     dj  =   dj + weight
                     d   = - dj
                  end if

               else if (js .eq. 4) then
*                 ------------------------------------------------------
*                 Fixed elastic variable.
*                 xj is free to move either way.
*                 ------------------------------------------------------
                  if (je .eq. 2) then
                     dj1 =   zero
                     d1  =   zero
                  else
                     dj1 =   dj - weight
                     d1  =   dj1
                  end if

                  if (je .eq. 1) then
                     dj2 =   zero
                     d2  =   zero
                  else
                     dj2 =   dj + weight
                     d2  = - dj2
                  end if

                  if (d1 .ge. d2) then
*                    xj  is allowed to decrease.
                     dj =   dj1
                     d  =   d1
                  else
*                    xj  is allowed to increase.
                     dj =   dj2
                     d  =   d2
                  end if
               else
                  d  = zero
                  dj = zero
               end if

               if (subopt .gt. 0) then
                  if (x(j) .eq. xFreez(j)) then
                     if (d  .gt. told) then
                        nFreez = nFreez + 1
                        d      = zero
                     end if
                  end if
               end if

*              See if this dj is significant.
*              Also see if it is the biggest dj so far.

               if (d  .gt. told) nonOpt = nonOpt + 1
               if (djmax .lt. d) then
                  djmax  = d
                  djq    = dj
                  jq     = j
                  kPsav  = kPrc
               end if
            end if
         end do
      end if

*     ------------------------------------------------------------------
*     End of loop looking for biggest dj in the kPrc-th section.
*     ------------------------------------------------------------------
      if (nonOpt .eq. 0) then
         if (nParP .gt. 1) then
*           ============================================================
*           No significant dj has been found.  (All are less than told.)
*           Price the next section, if any remain.
*           ============================================================
            if (nPrc .lt. nParP) go to 100

*           ============================================================
*           All sections have been scanned.  Reduce told
*           and grab the best dj if it is bigger than tolmin.
*           ============================================================
            if (djmax .gt. tolmin) then
               nonOpt = 1
               kPrc   = kPsav
               told   = max( reduce * djmax, tolmin  )
               toldj(lvldj) = told
               if (lprDbg .ge. 1) then
                 write(str, 1000) itn, told, piNorm, weight
                 call snPRNT( 23, str, iw, leniw )
               end if
            end if
         end if
      end if

*     -----------------------------------------------------------------
*     Finish if we found someone nonoptimal (nonOpt .gt. 0)
*     or if there's a nonbasic floating free
*     between its bounds                    (jfree  .eq. 0)
*     or if the problem is nonlinear        (nnL    .gt. 0)
*     or there are some superbasics         (nS     .gt. 0).
*     -----------------------------------------------------------------
      incres = djq .lt. zero
      if (nonOpt .gt. 0) go to 900
      if (jfree  .eq. 0) go to 900
      if (nS     .gt. 0) go to 900
      if ( feasbl ) then
         if (nnL .gt. 0) go to 900
      end if

*     -----------------------------------------------------------------
*     jfree > 0 and nS = 0.
*     We prefer not to end an LP problem (or an infeasible problem)
*     with some nonbasic variables floating free between their bounds
*     (hs(j) = -1).  Their dj's will be essentially zero
*     but we might as well jam them into the basis.
*
*     First, we try leaving true free variables alone -- they will
*     probably be zero.
*     -----------------------------------------------------------------
      do j = 1, nb
         if (hs(j) .eq. -1) then
            b1     = bl(j)
            b2     = bu(j)
            if (b1 .gt. -infBnd  .or.  b2 .lt. infBnd) then

*              We just found a floating variable with finite bounds.
*              Ask for a move towards the bound nearest zero.

               incres = abs(b1) .ge. abs(b2)
               nonOpt = 1
               jq     = j
               djq    = rc(j)
               go to 900
            end if
         end if
      end do

*     -----------------------------------------------------------------
*     If we are about to declare the problem infeasible we move nonzero
*     free variables towards zero.
*     -----------------------------------------------------------------
      if (.not. feasbl  .and.  nonOpt .eq. 0) then
         do j = 1, nb
            if (hs(j) .eq. -1  .and.  x(j) .ne. zero) then

*              We just found a true free variable.
*              Ask for a move towards zero.

               incres = x(j) .le. zero
               nonOpt = 1
               jq     = j
               djq    = rc(j)
               go to 900
            end if
         end do
      end if

*     Exit.

  900 return

 1000 format(' Itn', i7, ': toldj =', 1p, e8.1,
     &       '    Norm pi =', e8.1, '    weight = ', e8.1)

      end ! subroutine s5pric

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5rc
     &   ( j1, j2, gotg, m, n, leng, neg, sgnObj,
     &     ne, nlocA, locA, indA, Acol,
     &     hs, g, pi, rc )

      implicit
     &     none
      logical
     &     gotg
      integer
     &     j1, j2, m, n, nlocA, leng, neg, ne,
     &     locA(nlocA), indA(ne), hs(n)
      double precision
     &     sgnObj, Acol(ne), g(leng), pi(m), rc(n)

*     ==================================================================
*     s5rc   computes reduced costs rc(j) for nonbasic columns of A
*     in the range j = j1 to j2.  It is called by s5pric.
*
*     The loop for computing dj for each column could conceivably be
*     optimized on some machines.  However, there are seldom more than
*     5 or 10 entries in a column.
*
*     Note that we could skip fixed variables by passing in the bounds
*     and testing if bl(j) .eq. bu(j), but these are relatively rare.
*     But see comment for 08 Apr 1992 in m5pric.
*
*     31 Jan 1992: First version of s5rc.
*     08 Apr 1992: Internal values of hs(j) are now used, so fixed
*                  variables (hs(j) = 4) are skipped as we would like.
*     03 Apr 1999: Linear objective stored as row 0 of A.
*     11 Apr 1999: Current version.
*     ==================================================================
      integer
     &     i, j, l
      double precision
     &     dj
*     ------------------------------------------------------------------
      double precision   zero
      parameter         (zero = 0.0d+0)
*     ------------------------------------------------------------------
      do j = j1, j2
         if (hs(j) .le. 1) then
            dj    = zero
            do l  = locA(j), locA(j+1)-1
               i  = indA(l)
               dj = dj  +  pi(i) * Acol(l)
            end do
            rc(j) = - dj
         end if
      end do

*     Include the nonlinear objective gradient if relevant.

      if ( gotg ) then
         if (j1 .le. neg) then
            do j = j1, min( j2, neg )
               if (hs(j) .le. 1) rc(j) = rc(j) + sgnObj*g(j)
            end do
         end if
      end if

      end ! subroutine s5rc

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5rcE
     &   ( j1, j2, gotg, m, n, leng, neg, sgnObj,
     &     ne, nlocA, locA, indA, Acol,
     &     hElast, hs, g, pi, rc )

      implicit
     &     none
      logical
     &     gotg
      integer
     &     j1, j2, m, n, nlocA, leng, neg, ne,
     &     hElast(n), hs(n), locA(nlocA), indA(ne)
      double precision
     &     sgnObj, Acol(ne), g(leng), pi(m), rc(n)

*     ==================================================================
*     s5rcE  computes reduced costs rc(j) in the range j = j1 to j2
*     for fixed nonbasic columns that have one or more elastic bounds.
*     It is called by s5pric.
*
*     07 Feb 1998: First version based on s5rc.
*     03 Apr 1999: Objective stored as row 0 of A.
*     26 Jul 1999: Current version.
*     ==================================================================
      integer
     &     i, j, l
      double precision
     &     dj
*     ------------------------------------------------------------------
      double precision   zero
      parameter         (zero = 0.0d+0)
*     ------------------------------------------------------------------
      do j = j1, j2
         if (hs(j) .eq. 4  .and.  hElast(j) .gt. 0) then
            dj    = zero
            do l  = locA(j), locA(j+1)-1
               i  = indA(l)
               dj = dj  +  pi(i) * Acol(l)
            end do
            rc(j) = - dj
         end if
      end do

*     Include the nonlinear gradient term if present.

      if ( gotg ) then
         if (j1 .le. neg) then
            do j = j1, min( j2, neg )
               if (hs(j) .eq. 4  .and.  hElast(j) .gt. 0) then
                  rc(j) = rc(j) + sgnObj*g(j)
               end if
            end do
         end if
      end if

      end ! subroutine s5rcE

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5setE
     &   ( nBS, nb, nElast, featol, infBnd,
     &     hElast, hEstat, kBS,
     &     bl, bu, blBS, buBS, xBS )

      implicit
     &     none
      integer
     &     nb, nBS, nElast, hElast(nb), hEstat(nb), kBS(nBS)
      double precision
     &     featol, infBnd, bl(nb), bu(nb), blBS(nBS), buBS(nBS),
     &     xBS(nBS)

*     ==================================================================
*     s5setE  resets the upper and lower bounds on the elastic variables
*     for Elastic Phase 1 and 2.
*
*     s5setE is called in elastic mode after a factorize.
*     The bounds blBS and buBS are redefined for the basic elastics.
*
*     23 Aug 1996: First version of s5setE.
*     12 Jun 2000: Elastic mode cleaned up.
*     10 Jun 2001: Current version.
*     ==================================================================
      integer
     &     j, jE, jEs, k
      double precision
     &     res
*     ------------------------------------------------------------------
      nElast = 0

      do k   = 1, nBS
         j   = kBS(k)
         jE  = hElast(j)
         jEs = hEstat(j)

         if (jEs .gt. 0) then
*           ------------------------------------------------------------
*           x(j) is predicted to violate an elastic bound.
*           ------------------------------------------------------------
            if (jEs .eq. 1) then

*              Elastic lower bound.

               res = bl(j) - xBS(k)
               if (res .gt. -featol) then

*                 x(j) is indeed beyond its elastic lower bound.

                  nElast  =   nElast + 1
                  blBS(k) = - infBnd
                  buBS(k) =   bl(j)

               else

*                 x(j) satisfies its elastic lower bound.
*                 Check if x(j) violates an opposite elastic bound.

                  jEs = 0
                  if (jE .eq. 2  .or.  jE .eq. 3) then
                     res = xBS(k) - bu(j)
                     if (res .gt. featol) then
                        nElast  =   nElast + 1
                        jEs     =   2
                        blBS(k) =   bu(j)
                        buBS(k) =   infBnd
                     end if
                  end if
               end if

            else if (jEs .eq. 2) then

*              Elastic upper bound.

               res = xBS(k) - bu(j)
               if (res .gt. -featol) then

*                 x(j) is indeed beyond its elastic upper bound.

                  nElast  =   nElast + 1
                  buBS(k) =   infBnd
                  blBS(k) =   bu(j)
               else

*                 x(j) satisfies its elastic upper bound.
*                 Check if x(j) violates an opposite elastic bound.

                  jEs = 0
                  if (jE .eq. 1  .or.  jE .eq. 3) then
                     res = bl(j) - xBS(k)
                     if (res .gt. featol) then
                        nElast  =   nElast + 1
                        jEs     =   1
                        blBS(k) = - infBnd
                        buBS(k) =   bl(j)
                     end if
                  end if
               end if
            end if

         else if (jE .gt. 0) then
*           ------------------------------------------------------------
*           x(j) is predicted to satisfy its bounds.
*           Check that this is the case for any elastic bounds.
*           ------------------------------------------------------------
            if (jE .eq. 1  .or.  jE .eq. 3) then
               res = bl(j) - xBS(k)
               if (res .gt. featol) then
                  nElast  =   nElast + 1
                  jEs     =   1
                  blBS(k) = - infBnd
                  buBS(k) =   bl(j)
               end if
            end if

            if (jE .ge. 2) then ! jE .eq. 2  .or.  jE .eq. 3
               res = xBS(k) - bu(j)
               if (res .gt. featol) then
                  nElast  =   nElast + 1
                  jEs     =   2
                  blBS(k) =   bu(j)
                  buBS(k) =   infBnd
               end if
            end if
         end if

         hEstat(j) = -jEs

      end do

*     ------------------------------------------------------------------
*     Check all elements of  hEstat(j)  in case the factorize
*     changed kBS.
*     ------------------------------------------------------------------
      do j = 1, nb
         if (hEstat(j) .eq. 0) then
*           Relax
         else if (hEstat(j) .lt. 0) then
            hEstat(j) = - hEstat(j)
         else
            hEstat(j) =   0
         end if
      end do

      end ! subroutine s5setE

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5setp
     &   ( iExit, m, chkpi, piNorm, rhs, pi, iw, leniw, rw, lenrw )

      implicit
     &     none
      logical
     &     chkpi
      integer
     &     iExit, leniw, lenrw, m, iw(leniw)
      double precision
     &     piNorm, pi(m), rhs(m), rw(lenrw)

*     ==================================================================
*     s5setp solves  B' pi = rhs.  Beware -- rhs is altered by s2Bsol.
*     If a new x has just been computed, the norm is computed by dnormj.
*
*     On Exit:
*      iExit  = -1 if pi contains a NaN or Inf.
*      iExit  =  0 if pi was computed successfully.
*      iExit  >  0 if there was an unexpected error in the solver.
*
*     08 Aug 1996: First version of s5setE.
*     16 Nov 2001: Infinity norm used for piNorm (no longer dnrm1s).
*     28 Dec 2003: Current version.
*     ==================================================================
      external
     &     dnormi, dnormj
      double precision
     &     dnormi, dnormj, flMax
*     ------------------------------------------------------------------
      integer            WithBt
      parameter         (WithBt = 2)
      double precision   one
      parameter         (one    = 1.0d+0)
*     ------------------------------------------------------------------
      flMax = rw(  8) ! est. of the largest pos. real

      iExit = 0
      call s2Bsol
     &   ( iExit, WithBt, m, rhs, pi,  iw, leniw, rw, lenrw )
      if (iExit .ne. 0) return

      if ( chkpi ) then
         piNorm = dnormj( m, pi, 1 )
         if (piNorm .lt. flMax) then ! false if pi = inf, nan
            piNorm = max( piNorm, one )
         else
            iExit  = -1
         end if
      else
         piNorm = max( dnormi( m, pi, 1 ), one )
      end if

      end ! subroutine s5setp

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5setx
     &   ( iExit, Task, itn,
     &     m, n, nb, nBS, rowerr,
     &     ne, nlocA, locA, indA, Acol,
     &     kBS, xBS, nrhs0, nrhs, rhs, x, y, y1,
     &     iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     Task, iExit, itn, nrhs0, m, n, nb, nBS, ne, nlocA, nrhs,
     &     leniw, lenrw, locA(nlocA), indA(ne), kBS(nBS), iw(leniw)
      double precision
     &     rowerr, Acol(ne), rhs(nrhs0), xBS(nBS), x(nb), y(m),
     &     y1(m), rw(lenrw)

*     ==================================================================
*     s5setx performs the following functions:
*      Task            Function
*      ====           ========
*       0 (Resetx)    the basic components of x are computed to satisfy
*                     Ax - s = b; that is  (A -I)*x = b. Then a row
*                     check is performed to see how well  (A -I)*x = b
*                     is satisfied.  y is set to be the row residuals,
*                     y = b - Ax + s,  and the row error is norm(y).
*
*       1 (GetRes)    just get the row error.
*
*     The row error is a measure of how well x satisfies (A -I)*x = b.
*
*     18 Nov 1991: First version of s5setx based on Minos routine m5setx.
*     12 Jul 1997: Thread-safe version.
*     25 Jul 2003: Realized dx can sometimes be a NAN (or INF) but
*                  norm(NAN) > 1/eps is always false.
*     26 Jul 2003: Current version of s5setx.
*     ==================================================================
      character
     &     str*110
      external
     &     jdamax, dnormi, dnormj
      logical
     &     bigres, goodx
      integer
     &     imax, lprDbg, jdamax
      double precision
     &     dnormi, dnormj, eps0, rmax, tolrow, xNorm, dxNorm
*     ------------------------------------------------------------------
      integer            Normal,     WithB
      parameter         (Normal = 0, WithB  = 1)
      integer            Resetx
      parameter         (Resetx = 0)
      integer            xtoxBS,     xBStox
      parameter         (xtoxBS = 0, xBStox = 1)
      double precision   zero,            one
      parameter         (zero   = 0.0d+0, one = 1.0d+0)
*     ------------------------------------------------------------------
      eps0      = rw(  2) ! eps**(4/5)       IEEE DP  3.00e-13
      tolrow    = rw( 61) ! tolerance for the row error.
      lprDbg    = iw( 85) ! > 0    => private debug print

      iExit    = 0

      call s5BSx
     &   ( xtoxBS, nBS, nb, kBS, x, xBS )
      xNorm     = dnormi( nBS, xBS, 1 )
      dxNorm    = zero
      goodx     = .true.

*     ------------------------------------------------------------------
*     Compute row residuals  y  =  rhs - (A -I)*x.
*     The slack columns are done separately.
*     ------------------------------------------------------------------
      if (nrhs .gt. 0)
     &   call dcopy ( nrhs, rhs, 1, y        , 1 )
      if (nrhs .lt. m)
     &   call dload ( m-nrhs, zero, y(nrhs+1), 1 )

      call s2Aprd
     &   ( Normal, eps0,
     &     ne, nlocA, locA, indA, Acol,
     &     (-one), x, n, one, y, m )
      call daxpy
     &   ( m, one, x(n+1), 1, y, 1 )

*     ------------------------------------------------------------------
*     Do a row check, perhaps after recomputing the basic x.
*     -----------------------------------------------------------------
      if (Task .eq. Resetx) then
*        ================================================
*        Extract xBS, the basics and superbasics, from x.
*        See if iterative refinement is worth doing.
*        ================================================
         rowerr = dnormj( m, y, 1 )
         bigres = rowerr .gt. eps0

         if (bigres) then
*           ------------------------------------------------------------
*           Compute a correction to basic x from  B*y1 = y.
*           Extract the basic and superbasic variables from x.
*           Set basic x = x + y1.
*           Store the new basic variables in x.
*           ------------------------------------------------------------
            call s2Bsol
     &         ( iExit, WithB, m, y, y1, iw, leniw, rw, lenrw )
            if (iExit .ne. 0) return
            dxNorm = dnormj( m, y1, 1 )
            goodx  = dxNorm*eps0 .le. one ! false if dxNorm = inf, nan

            if ( goodx ) then
               call daxpy
     &            ( m, one, y1, 1, xBS, 1 )
               call s5BSx
     &            ( xBStox, m, nb, kBS, x, xBS )

*              Compute  y  =  rhs  -  (A -I)*x  again for the new x.

               if (nrhs .gt. 0)
     &            call dcopy ( nrhs, rhs, 1, y, 1 )
               if (nrhs .lt. m)
     &            call dload ( m-nrhs, zero, y(nrhs+1), 1 )

               call s2Aprd
     &            ( Normal, eps0,
     &              ne, nlocA, locA, indA, Acol,
     &              (-one), x, n, one, y, m )
               call daxpy
     &            ( m, one, x(n+1), 1, y, 1 )
            else
               iExit = 11      ! big dx (may be nan or inf)
            end if
         end if ! bigres
      end if ! Task .eq. Reset

*     Find the norm of xBS, the basic and superbasic x.
*     Find the maximum row residual.

      imax   = jdamax( m, y, 1 )
      if (imax .gt. 0) then
         rmax   = abs( y(imax) )
         rowerr = rmax / (one + xNorm )
      else
         imax   = -imax
         rmax   =  dnormj( m, y, 1 ) ! = flmax!
         rowerr =  rmax
      end if

      bigres = rowerr .gt. tolRow
      if (bigres) iExit = 10

      if (iExit .gt. 0  .or.  lprDbg .ge. 2) then
         write(str, 1000) itn, rmax, imax, xNorm, dxNorm
         call snPRNT( 21, str, iw, leniw )
         write(str, 1001) itn, rmax, imax
         call snPRNT( 22, str, iw, leniw )
      end if

      return

 1000 format(  ' Itn', i7, ': Row check',
     &         '.  Max residual =', 1p, e8.1, ' on row', i5,
     &         '.  Norm x =', e8.1, '.  Norm dx =', e8.1 )
 1001 format(  ' Itn', i7, ': Row check',
     &         '.  Max residual =', 1p, e8.1, ' on row', i5 )

      end ! subroutine s5setx

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5step
     &   ( nBS, nDegen,
     &     featol, infBnd, stepmx, tolinc, tolpiv,
     &     hfeas, blBS, buBS, xBS, pBS,
     &     hitlow, move, onbnd, unbndd,
     &     infpiv, kp, bound, exact, alpha, alphap )

      implicit
     &     none
      logical
     &     hitlow, move, onbnd, unbndd
      integer
     &     infpiv, kp, nBS, nDegen, hfeas(nBS)
      double precision
     &     featol, infBnd, stepmx, tolinc, tolpiv, bound, exact,
     &     alpha, alphap, blBS(nBS), buBS(nBS), xBS(nBS), pBS(nBS)

*     ==================================================================
*     s5step  finds a step  alpha  such that the point  xBS + alpha*pBS
*     reaches one of the bounds on  xBS.
*
*     In this  version of  s5step, when x  is infeasible, the  number of
*     infeasibilities  will never  increase.   If the  number stays  the
*     same,  the sum of  infeasibilities will  decrease.  If  the number
*     decreases by one or more,  the sum of infeasibilities will usually
*     decrease also,  but occasionally it  will increase after  the step
*     alpha is taken.  (Convergence  is still assured because the number
*     has decreased.)
*
*     Two possible steps are computed as follows:
*
*     alphaf = the maximum step that can be taken without violating
*              one of the bounds that are currently satisfied.
*
*     alphai = the maximum (nonzero) step that has the property of
*              reaching a bound that is currently violated,
*              subject to the pivot being reasonably close to the
*              maximum pivot among infeasible variables.
*              (alphai is not defined if x is feasible.)
*
*     alphai is  needed occasionally  when infeasible, to  prevent going
*     unnecessarily far when alphaf is quite large.  It will always come
*     into  effect when  x  is about  to  become feasible.   The sum  of
*     infeasibilities  will decrease initially  as alpha  increases from
*     zero, but may start increasing for larger steps.  choosing a large
*     alphai allows several elements of x to become feasible at the same
*     time.
*
*     In the end, we take  alpha = alphaf  if x is feasible, or if
*     alphai > alphap (where alphap is the perturbed step from pass 1).
*     Otherwise,  we take  alpha = alphai.
*
*     Input parameters
*     ----------------
*     nBS    is  m + 1  for s5lpit,  m + nBS  for s5QPit.
*     stepmx defines what should be treated as an unbounded step.
*     infBnd provides insurance for detecting unboundedness.
*            if alpha reaches a bound as large as infBnd, it is
*            classed as an unbounded step.
*     tolpiv is a tolerance to exclude negligible elements of pBS.
*     featol is the current feasibility tolerance used by s5QP.
*            Typically in the range 0.5*tolx to 0.99*tolx,
*            where tolx is the featol specified by the user.
*     tolinc is used to determine stepmn (see below), the minimum
*            positive step.
*     hfeas  is set by  s5Inf  as follows:
*            hfeas(j) = -2  if x(j) .lt. bl(j) - featol
*                     =  0  if x(j)  is feasible
*                     = +2  if x(j) .gt. bu(j) + featol
*     blBS   the lower bounds on the basic and superbasic variables.
*     BLBU   the upper bounds on ditto.
*     xBS    the values of       ditto.
*     pBS    the search direction.
*
*
*     Output parameters
*     -----------------
*     hitlow  = true  if a lower bound restricted alpha.
*             = false otherwise.
*     move    = true  if exact ge stepmn (defined at end of code).
*     onbnd   = true  if alpha =  exact.
*                     this means that the step alpha moves xBS(kp)
*                     exactly onto one of its bounds, namely bound.
*             = false if the exact step would be too small
*                     ( exact lt stepmn ).
*               (with these definitions,  move = onbnd).
*     unbndd  = true  if alpha = stepmx.  kp may possibly be zero.
*               the parameters hitlow, move, onbnd, bound and exact
*               should not be used.
*     infpiv  = the number of indices such that |pBS(i)| < tolpiv and
*               xBS(i) + alphap*pBS(i) is infeasible.
*     kp      = the index (if any) such that xBS(kp) reaches a bound.
*     bound   = the bound value blBS(kp) or buBS(kp) corresponding
*               to hitlow.
*     exact   = the step that would take xBS(kp) exactly onto bound.
*     alpha   = an allowable, positive step.
*               if unbndd is true,  alpha = stepmx.
*               otherwise,          alpha = max( stepmn, exact ).
*     alphap  = the perturbed step from pass 1.
*
*     07 Nov 1991: First version based on Minos routine m5chzr.
*     27 Dec 2003: infpiv added to monitor unwanted infeasibilities.
*     27 Dec 2003: Current version of s5step.
*     ==================================================================
      logical
     &     blockf, blocki
      integer
     &     j, jtype, jhiti, jhitf
      double precision
     &     alphai, delta, pivot, pivabs, pivmxi, pivmxf, res, stepmn
*     ------------------------------------------------------------------
      double precision   zero,          gamma
      parameter         (zero = 0.0d+0, gamma = 0.001d+0)
*     ------------------------------------------------------------------
*     First pass.
*     For  feasible variables,  find the  step alphap  that  reaches the
*     nearest perturbed (expanded) bound.   alphap will be slight larger
*     than the step to the nearest true bound.
*     For infeasible variables, find the maximum pivot pivmxi.
*     ------------------------------------------------------------------
      delta  = featol
      alphap = stepmx
      pivmxi = zero
      jhitf  = 0

      do 200, j  = 1, nBS
         pivot   = pBS(j)
         pivabs  = abs( pivot )
         if (pivabs .gt. tolpiv) then
            jtype  = hfeas(j)
            if (pivot  .gt. zero  ) go to 150

*           x  is decreasing.
*           Test for smaller alphap if lower bound is satisfied.

            if (jtype .lt. 0) go to 200
            res    = xBS(j) - blBS(j) + delta
            if (alphap*pivabs .gt. res) then
               alphap = res / pivabs
               jhitf  = j
            end if

*           Test for bigger pivot if upper bound is violated.

            if (jtype .gt. 0) pivmxi = max( pivmxi, pivabs )
            go to 200

*           x  is increasing.
*           Test for smaller alphap if upper bound is satisfied.

  150       if (jtype .gt. 0) go to 200
            res    = buBS(j) - xBS(j) + delta
            if (alphap*pivabs .gt. res) then
               alphap = res / pivabs
               jhitf  = j
            end if

*           Test for bigger pivot if lower bound is violated.

            if (jtype .lt. 0) pivmxi = max( pivmxi, pivabs )
         end if
  200 continue

*     ------------------------------------------------------------------
*     Second pass.
*     For  feasible  variables,  recompute steps  without  perturbation.
*     Choose  the largest  pivot element  subject to  the step  being no
*     greater than alphap.
*     For infeasible variables, find the largest step subject to the
*     pivot element being no smaller than gamma * pivmxi.
*     ------------------------------------------------------------------
      alphai = zero
      pivmxf = zero
      pivmxi = gamma * pivmxi
      jhiti  = 0
      infpiv = 0

      do 400, j = 1, nBS
         pivot  = pBS(j)
         pivabs = abs( pivot )
         jtype  = hfeas(j)
         if (pivabs .gt. tolpiv) then
            if (pivot  .gt. zero  ) go to 350

*           x  is decreasing.
*           Test for bigger pivot if lower bound is satisfied.

            if (jtype    .lt.     0   ) go to 400
            if (pivabs   .le.   pivmxf) go to 340
            res    = xBS(j) - blBS(j)
            if (alphap*pivabs .lt. res) go to 340
            pivmxf = pivabs
            jhitf  = j

*           Test for bigger alphai if upper bound is violated.

  340       if (jtype    .eq.     0   ) go to 400
            if (pivabs   .lt.   pivmxi) go to 400
            res    = xBS(j) - buBS(j)
            if (alphai*pivabs .ge. res) go to 400
            alphai = res / pivabs
            jhiti  = j
            go to 400

*           x  is increasing.
*           Test for bigger pivot if upper bound is satisfied.

  350       if (jtype    .gt.     0   ) go to 400
            if (pivabs   .le.   pivmxf) go to 360
            res    = buBS(j) - xBS(j)
            if (alphap*pivabs .lt. res) go to 360
            pivmxf = pivabs
            jhitf  = j

*           Test for bigger alphai if lower bound is violated.

  360       if (jtype    .eq.     0   ) go to 400
            if (pivabs   .lt.   pivmxi) go to 400
            res    = blBS(j) - xBS(j)
            if (alphai*pivabs .ge. res) go to 400
            alphai = res / pivabs
            jhiti  = j

         else if (jtype .eq. 0  .and.  pivabs .gt. zero) then

*           Small pivot.
*           Check that the perturbed step is feasible.

            if (pivot  .lt. zero) then ! x  is decreasing.
               res    = xBS(j) - blBS(j)
               if (alphap*pivabs .gt. res) infpiv = infpiv + 1
            else                       ! x  is increasing.
               res    = buBS(j) - xBS(j)
               if (alphap*pivabs .gt. res) infpiv = infpiv + 1
            end if
         end if
  400 continue

*     ------------------------------------------------------------------
*     See if a feasible and/or infeasible variable blocks.
*     ------------------------------------------------------------------
      blockf = jhitf .gt. 0
      blocki = jhiti .gt. 0
      unbndd = .not. ( blockf  .or.  blocki )
      if ( unbndd ) go to 900

      if ( blockf ) then
*        ---------------------------------------------------------------
*        A variable hits a bound for which it is currently feasible.
*        The step alphaf is not used, so no need to get it,
*        but we know that alphaf .le. alphap, the step from pass 1.
*        ---------------------------------------------------------------
         kp     = jhitf
         pivot  = pBS(kp)
         hitlow = pivot .lt. zero
      end if

*     If there is a choice between alphaf and alphai, it is probably
*     best to take alphai (so that the infeasible variable jhiti can be
*     kicked out of the basis).
*     However, we can't if alphai is bigger than alphap.

      if (blocki  .and.  alphai .le. alphap) then
         kp     = jhiti
         pivot  = pBS(kp)
         hitlow = pivot .gt. zero
      end if

*     ------------------------------------------------------------------
*     Try to  step exactly onto bound,  but make sure the  exact step is
*     sufficiently positive.   (exact will be alphaf  or alphai.)  Since
*     featol increases by tolinc each  iteration, we know that a step as
*     large as stepmn  (below) will not cause any  feasible variables to
*     become infeasible  (where feasibility  is measured by  the current
*     featol).
*     ------------------------------------------------------------------
      if ( hitlow ) then
         bound = blBS(kp)
      else
         bound = buBS(kp)
      end if
      unbndd = abs( bound ) .ge. infBnd
      if ( unbndd ) go to 900

      stepmn = tolinc / abs( pivot )
      exact  = (bound - xBS(kp)) / pivot
      alpha  = max( stepmn, exact )
      onbnd  = alpha .eq. exact
      move   = exact .ge. stepmn
      if (.not. move) nDegen = nDegen + 1
      return

*     ------------------------------------------------------------------
*     Unbounded.
*     ------------------------------------------------------------------
  900 alpha  = stepmx
      move   = .true.
      onbnd  = .false.

      end ! subroutine s5step

