*
*===program crint======================================================*
*
C      OPTIONS/ EXTEND_SOURCE
C     SUBROUTINE CRINT
*     KEEP_PHI, KEEP_KSTAR are switch to store phi and K*0 and its decay daughter chain information in output        
      SUBROUTINE DT_PRODUCEEVENT(ENERGY_SL, NPARTICLES, KEEP_PHI,
     & KEEP_KSTAR)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL ENERGY_SL
      INTEGER INIT,KEEP_PHI, KEEP_KSTAR
      REAL ne,etest,prob,slump
      SAVE

* Call the init sub routine in the first event
      DATA INIT /0/

      INCLUDE 'inc/dtflka'

      INCLUDE 'inc/dtevt1'

*     event flag
      INCLUDE 'inc/dtevno'

      LPRi = 20
      LOUt = 6
      IF(INIT.EQ.0) THEN
         OPEN (UNIT = 50, file = "my.input")    
	 LINP = 50
         CALL DT_DTUINI(NEVTS,EPN,NPMASS,NPCHAR,NTMASS,NTCHAR,IDP,IEMU)
*        Init called, make sure it's not called again
         INIT = 1
      ENDIF

*-----------------------------------------------------------------------
*     generation of one event
      NEVENT = 1
      KKMAT = -1

*   If an energy-range has been defined with the ENERGY input-card the
*   laboratory energy ELAB can be set to any value within that range,..
C        ELAB = DT_RNDM(EPN)*(EPN-0.5D7)+0.5D7

*   ..otherwise it has to coincide with EPN.
C        ELAB = EPN

      ELAB = ENERGY_SL

*   sampling of one event

*     TEST

      CALL DT_KKINC(NPMASS,NPCHAR,NTMASS,NTCHAR,IDP,ELAB,KKMAT,IREJ)

      IF (IREJ.NE.0) RETURN

c     Return the number of particles produced
*     KEEP_PHI, KEEP_KSTAR are switch to store phi and K*0 and its decay daughter chain information in output        
c     Fill the particle info 
      CALL DT_GETPARTICLES(NPARTICLES, KEEP_PHI, KEEP_KSTAR)

      END


      SUBROUTINE DT_GETPARTICLES(NPARTICLES, KEEP_PHI, KEEP_KSTAR)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER pid,qch,q_sum,Ntpc,Nfinal,NACCEPT,IPART,RES
      DOUBLE PRECISION yrap,pt,mass,mt,etot
      DOUBLE PRECISION pt_cut_tpc
      PARAMETER(pt_cut_tpc=0.050)

      SAVE
*
* COMMON /DTEVT1/ :
*                   NHKK         number of entries in common block
*                   NEVHKK       number of the event
*                   ISTHKK(i)    status code for entry i
*                   IDHKK(i)     identifier for the entry
*                                (for particles: identifier according
*                                 to the PDG numbering scheme)
*                   JMOHKK(1,i)  pointer to the entry of the first mother
*                                of entry i
*                   JMOHKK(2,i)  pointer to the entry of the second mother
*                                of entry i
*                   JDAHKK(1,i)  pointer to the entry of the first daughter
*                                of entry i
*                   JDAHKK(2,i)  pointer to the entry of the second daughter
*                                of entry i
*                   PHKK(1..3,i) 3-momentum
*                   PHKK(4,i)    energy
*                   PHKK(5,i)    mass
*
* event history

      INCLUDE 'inc/dtevt1'

* extended event history
      INCLUDE 'inc/dtevt2'

      DOUBLE PRECISION SLPX, SLPY, SLPZ, SLE, SLM
      INTEGER SLPID, SLCHARGE
      COMMON /DPMJETPARTICLE/ SLPX(NMXHKK), SLPY(NMXHKK), SLPZ(NMXHKK),
     &       SLE(NMXHKK), SLM(NMXHKK), SLPID(NMXHKK), SLCHARGE(NMXHKK)

* Declare and dimension the new arrays for phi and K*0 before the COMMON
      INTEGER SLMOTH1, SLMOTH2, SLSTATUS
      DIMENSION SLMOTH1(NMXHKK), SLMOTH2(NMXHKK), SLSTATUS(NMXHKK)
      COMMON /DPMJETMOTHERS/ SLMOTH1, SLMOTH2, SLSTATUS

* Switches for each resonance (set in main program or defaults)
C     >> Set Counter to Zero
      LOGICAL KEEPPARTICLE
      Nfinal=0
      
      DO 42 I=1, NHKK
c      I = IPART

* --- Decide if particle is kept ---
     
      KEEPPARTICLE = .FALSE.

* --- Always keep stable particles ---
      IF (ISTHKK(I).EQ.1 .OR. ISTHKK(I).EQ.-1 .OR. ISTHKK(I).EQ.1001)
     & KEEPPARTICLE = .TRUE.
     
* --- Keep selected Phi and K*0 resonances according to switches ---
      IF (KEEP_PHI.EQ.1 .AND. ISTHKK(I).EQ.2 .AND. IDHKK(I).EQ.333)
     & KEEPPARTICLE = .TRUE.
      IF (KEEP_KSTAR.EQ.1 .AND. ISTHKK(I).EQ.2 .AND. (IDHKK(I).EQ.313
     & .OR. IDHKK(I).EQ.-313)) KEEPPARTICLE = .TRUE.
    
      IF (.NOT.KEEPPARTICLE) GOTO 42
         
C	>> Find Particle Charge, qch
        IF((ABS(ISTHKK(I)).eq.1).and.(IDHKK(I).ne.80000))THEN
C         >> final state ptcles except nuclei

          qch=IPHO_CHR3(IDHKK(I),1)/3
        ELSEIF(IDHKK(I).eq.80000)THEN
C         >> final state nuclei
          qch=IDXRES(I)
        ELSE
C         >> not a final state particle, qch not interesting
          qch=-999
        ENDIF

	Nfinal = Nfinal + 1
	SLPX(Nfinal) = PHKK(1,I)
        SLPY(Nfinal) = PHKK(2,I)
        SLPZ(Nfinal) = PHKK(3,I)
        SLE(Nfinal) = PHKK(4,I)
        SLM(Nfinal) = PHKK(5,I)
        SLPID(Nfinal) = IDHKK(I)
        SLCHARGE(Nfinal) = qch
* Now assign mother(phi and K*0) info to the same index
        SLSTATUS(Nfinal) = ISTHKK(I)
        SLMOTH1(Nfinal) = JMOHKK(1,I)
        SLMOTH2(Nfinal) = JMOHKK(2,I)

 42     CONTINUE
        NPARTICLES = Nfinal
  
      END

      SUBROUTINE DT_USRHIS(MODE)
c Dummy to make the linker happy
      END

