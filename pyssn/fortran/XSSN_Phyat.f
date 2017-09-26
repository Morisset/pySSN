c  gfortran XSSN_Phyat.f
c  gfortran -fcheck=all XSSN_Phyat.f
c  ./a.out
c
c-------------------------------------------------------
      Program XSSN_Phyat
c
c  Compute recombination lines for nebulae in physical conditions to be specified.
c  Output is in format corresponding to code pySSN (Spectral Synthesis for Nebulae).
c  Daniel Péquignot, Oct. 2015 - Sept. 2017.
c
      IMPLICIT NONE
c  ***  EXPLICIT: IONS ***
c (IONS is generic name for the list of ions to be changed in liste_phyat..)
      CHARACTER IONS*(*)
      PARAMETER (IONS='ions_rec.dat') !name imposed for now
c
c  *** REQUIRED FILES: ***
c  ***********************
c   HI: e1bx.d, HeII: e2bx.d, H-like: all r1....d to r8....d (Storey & Hummer 1995).
c   BA5.out (for H-like branchings, Hoang Binh, Dy, 1990),
c   HeI: HeI_Porter13.dat, all Te_1_0_50_Results.txt 
c                             (Te = 5000-25000 (1000), Bauman et al. 21007)
c   CII: CII_diel.dat
c   OII_Sto: line_...._dr (Storey 2013),
c   OUTPUTCOND (eg, outputcond.dat), read in ions_rec.dat
c   PHYAT (eg, liste_phyat_NF.dat), phyat.old (previous saved),
c   IONS (ATT: name is imposed for now : ions_rec.dat),
c   all .lab (input) and .res (output) files for ions listed in IONS 
c   (note that some ions require several .lab and as many .res).
c  ***********************
c
      DOUBLE PRECISION hc,hc8,hcsk,Ryd_inf
      DOUBLE PRECISION csaha
      DOUBLE PRECISION hplank,clight,mel,kboltz,pi
      DOUBLE PRECISION int_ref,lambda_ref
      DOUBLE PRECISION Te,Ne
      DOUBLE PRECISION vzero
      CHARACTER lambda_refc*6,nomion*7,ligne*30
      CHARACTER OUTPUTCOND*30, PHYAT*30
      INTEGER lmax_del1,lmax_del2,lmax_del3
      INTEGER istor,istow,ibr,iprint
      INTEGER ifre1_ref  !VERIF USAGE ?
      INTEGER in_niv,out_niv,iwr_phyat,in_ions
      COMMON/MOC1/in_niv,out_niv,iwr_phyat
      COMMON/MOC2/iprint
      COMMON/MOCC1/hc,hc8,hcsk,Ryd_inf
      COMMON/MOCC2/csaha,kboltz
      COMMON/MOCH3/ibr
      COMMON/MOCS3/istor,istow
      COMMON/MOCV1/vzero
c Memo print for imult = 0 --> imult > 0 of one given ion :
c Print/Write correct if imult=0 is listed 1st in Main for the given ion
      COMMON/MOCX1/lambda_refc
      COMMON/MOCX2/int_ref,lambda_ref,ifre1_ref
      INTEGER yn,ntot,k,kmax,indice
      COMMON/MOCR1/wlmin,wlmax,wlim,Irelmin_H,Ionab,Ionab_min,iwlim
      DOUBLE PRECISION int_ref_mod_H
      COMMON/MOCR2/int_ref_mod_H
      DOUBLE PRECISION Ionab,Ionab_min
      DOUBLE PRECISION wlim(20),Irelmin_H(21)
      DOUBLE PRECISION wlmin,wlmax
      INTEGER iwlim
      CHARACTER bidon*1
      INTEGER iw,il,i,im
      INTEGER in_lim
c
c Phys Const: CODATA 2010
        hplank=6.62606957d-27 ! erg.s (29)
        clight=2.99792458d10  ! cm.s-1 (0)
        mel=9.10938291d-28    ! g (40) 
c        mp=1.672621777d-24   ! g (74)
c        mp/mu=1.007276466812 ! (90) 
c        mp/mel=1836.15267245 ! (75)    
c        mn/mp=1.00137841917  ! (45)
        kboltz=1.3806488d-16  ! erg.K-1 (13)
        pi=3.1415926535898d0
        Ryd_inf=109737.31568539d0 !cm-1 (55); (109737.3153 CS99)
        hc=hplank*clight !erg.cm
        hc8=1.d8*hc      ! (hc8=1.98644568d-8)
        hcsk=hc/kboltz   ! 1.438777 cm.K
        csaha=((hplank/mel)*(hplank/kboltz)/2.d0/pi)**1.5 ! 4.1413306d-16
c        m_alph=6.64465675d-24    ! g (29)
c        m_alph/mel=7294.2995361  ! (29) 
c        m_alph/mp=3.97259968933  ! (36) 
c        md=3.34358348d-24        ! g (15)
c
c Control addresses fort files XSSN_Hlike, XSSN_di, etc.:
      in_niv=10  !adr. data .lab (niv_input): atomic levels, etc. (10 arb)
      out_niv=11 !adr. results .res (lambda_output): outputs, incl. X-SSN (11 arb)
c Control address fort ReadSto: r read, w write
      istor=13   ! adr. selected Storey Hummer files (13 arb)
c Control address for ion astrophysical data: 
      in_ions=17 ! adr. ions_rec.dat file (17 arb)
      in_lim=32  ! adr. outputcond.dat file (32 arb)
c**** ESSAI 12/11/15
      istow=0  ! istow=0: no write, stdd ! istow=20 (\rm alfsto.res) (arb)
c Control address fort HBD : write in BA5.out
c****ESSAI 12/11/15
      ibr=0    ! ibr=0: no write, stdd ! ibr=15 (touch BA5.out): write (arb)
c Control write XSSN:   (ATT: update iprint.ne.0 to be checked ?)
      iprint=0 ! iprint=0: input + XSSN output only, stdd ! .ne.0: print all 1/14: OBSO ?
c
      OPEN(unit=in_ions,file=IONS,status='old')
c
      ntot=0
      Do
        read(in_ions,*,end=1000) ligne
        ntot=ntot+1
      Enddo
 1000 continue
c  Total nbre of lines: ntot
c  Total nbre of ions listed (not all will be processed):
      kmax=ntot-5
c
      REWIND(in_ions)
c
      read(in_ions,100) iwr_phyat
 100  format(20i4)
      read(in_ions,101) vzero !nebular source property (*OBSO?*)
         print *,'Intrinsic HW1/e (km/s) vzero =',vzero
 101  format(10d6.0)
c***********************!!! OBSO?
c  vzero is half-width at 1/e intensity of the reference intrinsic gaussian profile of emission lines from the source (vzero is input for prog profil_emis.pro):
c*** Ref NGC7027 21.5 km/s; NGC6302 13.0 km/s; N49 50 km/s; NGC3918 13.0 km/s
c***********************!!!
      read(in_ions,102) PHYAT
        print *,' PHYAT = ', PHYAT
      read(in_ions,102) OUTPUTCOND 
        print *,' OUTPUTCOND = ', OUTPUTCOND
 102  format(a30) 
c***********************!!!
cc      iwr_phyat=0 ! .res is NOT transcribed in PHYAT
ccc      iwr_phyat=1 ! .res IS transcribed (hence, previous data erased)
      IF(iwr_phyat.eq.0) THEN
         print *,' *** New .res will NOT be transcribed'
      ELSE
         print *,' *** New .res will be transcribed in ',PHYAT
      ENDIF
c***********************!!!
c
c  Reading emissivity limits for output of HI recombination lines.
c  For other ions, these limits are scaled according to reference line intensity.
c
      OPEN(unit=in_lim,file=OUTPUTCOND,status='old')
c
c Title
      read(in_lim,1) bidon
 1    format(A1)
c lambda range for X-SSN output: wlmin, wlmax AA; 
      read(in_lim,61) wlmin,wlmax
 61   format(6x,d7.2,6x,d7.2)
      read(in_lim,1) bidon
c minimal output intensity is a function of lambda : iwlim+1 intervals
      read(in_lim,22) iwlim
 22   format(20i4)
       read(in_lim,26) (wlim(iw),iw=1,iwlim)
       read(in_lim,27) (Irelmin_H(iw),iw=1,iwlim+1)
 26   format(5x,20d8.2)
 27   format(21d8.2)
c minimal ionic abundance for first estimate of the Irelmin, see ions_rec.dat:
       read(in_lim,27) Ionab_min
c
      CLOSE(in_lim)
c**************************
c  5th line of file IONS:
      read(in_ions,102) bidon  ! (comment)
c*
      DO k=1,kmax
c*
       read(in_ions,110) nomion,yn,Te,Ne,Ionab
 110   format(a7,i2,3d9.2)
c
       If(k.eq.1) then
c  HI necessary to store Hbeta intensity:
      if(nomion.ne.'HI     ') then
       print *,' *** HI MUST BE 1ST ION IN FILE ions_rec.dat: STOP ***'
       STOP
      endif
      if(yn.eq.0) then
c  Preliminary call to store Hbeta intensity:
        yn=1
        print *,' *** HI in file ions_rec.dat: yn changed from 0 to 1'
      endif
       Endif
c
        IF(yn.ne.0) THEN
         indice=0
c
      if(nomion.eq.'HI     ') then
      Call XSSN_Hlike(Te,Ne,'HI_2.lab','HI_2.res',PHYAT)
        indice=1
      endif
c
      if(nomion.eq.'DI     ') then
      Call XSSN_Hlike(Te,Ne,'DI_2.lab','DI_2.res',PHYAT)
        indice=1
      endif
c
      if(nomion.eq.'HeI    ') then
      Call HeI_BP(Te,Ne,'HeI.lab','HeI.res',PHYAT)
        indice=1
      endif
c
c      if(nomion.eq.'3HeI   ') then
c      Call HeI_BP(Te,Ne,'3HeI.lab','3HeI.res',PHYAT)
c        indice=1
c      endif
c
      if(nomion.eq.'HeII   ') then
      Call XSSN_Hlike(Te,Ne,'HeII_2.lab','HeII_2.res',PHYAT)
        indice=1
      endif
c
      if(nomion.eq.'3HeII  ') then
      Call XSSN_Hlike(Te,Ne,'3HeII_2.lab','3HeII_2.res',PHYAT)
        indice=1
      endif
c
      if(nomion.eq.'LiII   ') then
      Call XSSN_Hlike(Te,Ne,'LiII_3.lab','LiII_3.res',PHYAT) ! en 1er
      Call XSSN_Hlike(Te,Ne,'LiII_1.lab','LiII_1.res',PHYAT)
        indice=1
      endif
c
      if(nomion.eq.'LiIII  ') then
      Call XSSN_Hlike(Te,Ne,'LiIII_2.lab','LiIII_2.res',PHYAT)
        indice=1
      endif
c
c *** CI TBD
c
      if(nomion.eq.'CII    ') then
      Call CII_DSK(Te,Ne,'CII_DSK.lab','CII_DSK.res',PHYAT) !Ref CII -> 1st
      Call XSSN_Hlike(Te,Ne,'CII_2.lab','CII_2.res',PHYAT) !D5/2 3/2 togeth
      Call CII_diel(Te,Ne,'CII_diel.lab','CII_diel.res',PHYAT)
        indice=1
      endif
c SUPPR 1/12: Call XSSN_Hlike(Te,Ne,'CII_2D5.lab','CII_2D5.res',PHYAT) !prev
c SUPPR 1/12: Call XSSN_Hlike(Te,Ne,'CII_2D3.lab','CII_2D3.res',PHYAT) !obso
c
      if(nomion.eq.'CIII   ') then
      Call XSSN_Hlike(Te,Ne,'CIII_3F4.lab','CIII_3F4.res',PHYAT) ! en 1er
      Call XSSN_Hlike(Te,Ne,'CIII_3F3.lab','CIII_3F3.res',PHYAT)
      Call XSSN_Hlike(Te,Ne,'CIII_3F2.lab','CIII_3F2.res',PHYAT)
      Call XSSN_Hlike(Te,Ne,'CIII_1F3.lab','CIII_1F3.res',PHYAT)
      Call XSSN_di(Te,Ne,'CIII_di.lab','CIII_di.res',PHYAT)
      Call CIII_5ga(Te,Ne,'CIII_5ga.lab','CIII_5ga.res',PHYAT)
        indice=1
      endif
c
      if(nomion.eq.'CIV    ') then
      Call XSSN_Hlike(Te,Ne,'CIV_2.lab','CIV_2.res',PHYAT) !D5/2 3/2 regr
        indice=1
      endif
ccc SUPPR 11/11: Call XSSN_Hlike(Te,Ne,'CIV_2D5.lab','CIV_2D5.res',PHYAT)
ccc SUPPR 11/11: Call XSSN_Hlike(Te,Ne,'CIV_2D3.lab','CIV_2D3.res',PHYAT)
c
c *** NI TBD
c
      if(nomion.eq.'NII    ') then
      Call NII_Sto(Te,Ne,'NII_13.lab','NII_13.res',PHYAT)
        indice=1
      endif
c
      if(nomion.eq.'NIII   ') then
      Call XSSN_Hlike(Te,Ne,'NIII_2F7.lab','NIII_2F7.res',PHYAT) ! en 1er
      Call XSSN_Hlike(Te,Ne,'NIII_2F5.lab','NIII_2F5.res',PHYAT)
        indice=1
      endif
c
      if(nomion.eq.'NIV    ') then
      Call XSSN_Hlike(Te,Ne,'NIV_3.lab','NIV_3.res',PHYAT) ! en 1er
      Call XSSN_Hlike(Te,Ne,'NIV_1.lab','NIV_1.res',PHYAT)
        indice=1
      endif
c
      if(nomion.eq.'NV     ') then
      Call XSSN_Hlike(Te,Ne,'NV_2.lab','NV_2.res',PHYAT)
        indice=1
      endif
c
c *** OI TBD
c
      if(nomion.eq.'OII    ') then
      Call OII_Sto(Te,Ne,'OII_24.lab','OII_24.res',PHYAT)
        indice=1
      endif
c
      if(nomion.eq.'OIII   ') then
c *** full OIII TBD
      Call OIII_5g(Te,Ne,'OIII_5g.lab','OIII_5g.res',PHYAT)
        indice=1
      endif
c
      if(nomion.eq.'OIV    ') then
      Call XSSN_Hlike(Te,Ne,'OIV_2.lab','OIV_2.res',PHYAT) ! en 1er
      Call XSSN_di(Te,Ne,'OIV_di.lab','OIV_di.res',PHYAT)
        indice=1
      endif
c
      if(nomion.eq.'OV     ') then
      Call XSSN_Hlike(Te,Ne,'OV_3G5.lab','OV_3G5.res',PHYAT) ! en 1er
      Call XSSN_Hlike(Te,Ne,'OV_3G4.lab','OV_3G4.res',PHYAT)
      Call XSSN_Hlike(Te,Ne,'OV_3G3.lab','OV_3G3.res',PHYAT)
      Call XSSN_Hlike(Te,Ne,'OV_1G4.lab','OV_1G4.res',PHYAT)
        indice=1
      endif
c
      if(nomion.eq.'OVI    ') then
      Call XSSN_Hlike(Te,Ne,'OVI_2.lab','OVI_2.res',PHYAT)
        indice=1
      endif
c
c *** NeI, II TBD
c
ccSUPPR      Call XSSN_Hlike(Te,Ne,'NeIII_3D3.lab','NeIII_3D3.res',PHYAT)
ccSUPPR      Call XSSN_Hlike(Te,Ne,'NeIII_3D2.lab','NeIII_3D2.res',PHYAT)
ccSUPPR      Call XSSN_Hlike(Te,Ne,'NeIII_3D1.lab','NeIII_3D1.res',PHYAT)
ccSUPPR      Call XSSN_Hlike(Te,Ne,'NeIII_5D.lab','NeIII_5D.res',PHYAT)
c  sans 4p et sans split 4d3D :
      if(nomion.eq.'NeIII  ') then
      Call XSSN_Hlike(Te,Ne,'NeIII_5.lab','NeIII_5.res',PHYAT) ! en 1er
      Call XSSN_Hlike(Te,Ne,'NeIII_3.lab','NeIII_3.res',PHYAT)
      Call XSSN_di(Te,Ne,'NeIII_di.lab','NeIII_di.res',PHYAT)
        indice=1
      endif
c
      if(nomion.eq.'NeIV   ') then
      Call XSSN_Hlike(Te,Ne,'NeIV_2.lab','NeIV_2.res',PHYAT)
        indice=1
      endif
c
      if(nomion.eq.'NeV    ') then
      Call XSSN_Hlike(Te,Ne,'NeV_3.lab','NeV_3.res',PHYAT) ! en 1er
      Call XSSN_Hlike(Te,Ne,'NeV_1.lab','NeV_1.res',PHYAT)
        indice=1
      endif
c
      if(nomion.eq.'NeVI   ') then
      Call XSSN_Hlike(Te,Ne,'NeVI_2.lab','NeVI_2.res',PHYAT)
        indice=1
      endif
c
      if(nomion.eq.'NeVII  ') then
      Call XSSN_Hlike(Te,Ne,'NeVII_3H6.lab','NeVII_3H6.res',PHYAT) ! en 1er
      Call XSSN_Hlike(Te,Ne,'NeVII_3H54.lab','NeVII_3H54.res',
     $    PHYAT)
      Call XSSN_Hlike(Te,Ne,'NeVII_1H5.lab','NeVII_1H5.res',PHYAT)
        indice=1
      endif
c
      if(nomion.eq.'NeVIII ') then
      Call XSSN_Hlike(Te,Ne,'NeVIII_2.lab','NeVIII_2.res',PHYAT)
        indice=1
      endif
c
c *** MgI TBD
c
      if(nomion.eq.'MgII   ') then
      Call XSSN_Hlike(Te,Ne,'MgII_2.lab','MgII_2.res',PHYAT)
        indice=1
      endif
c
c *** MgIII, IV, V, VI, VII TBD
c
c *** SiI TBD
c
      if(nomion.eq.'SiII   ') then
      Call XSSN_Hlike(Te,Ne,'SiII_2D5.lab','SiII_2D5.res',PHYAT) ! en 1er
      Call XSSN_Hlike(Te,Ne,'SiII_2D3.lab','SiII_2D3.res',PHYAT)
        indice=1
      endif
c
      if(nomion.eq.'SiIII  ') then
      Call XSSN_Hlike(Te,Ne,'SiIII_3F4.lab','SiIII_3F4.res',PHYAT) ! en 1er
      Call XSSN_Hlike(Te,Ne,'SiIII_3F3.lab','SiIII_3F3.res',PHYAT)
      Call XSSN_Hlike(Te,Ne,'SiIII_3F2.lab','SiIII_3F2.res',PHYAT)
      Call XSSN_Hlike(Te,Ne,'SiIII_1F3.lab','SiIII_1F3.res',PHYAT)
        indice=1
      endif
c
      if(nomion.eq.'SiIV   ') then
      Call XSSN_Hlike(Te,Ne,'SiIV_2D5.lab','SiIV_2D5.res',PHYAT) ! en 1er
      Call XSSN_Hlike(Te,Ne,'SiIV_2D3.lab','SiIV_2D3.res',PHYAT)
        indice=1
      endif
c
c *** SiV, VI, VII TBD
c
c *** SI, II, III, IV TBD
c
      if(nomion.eq.'SV     ') then
      Call XSSN_Hlike(Te,Ne,'SV_3G5.lab','SV_3G5.res',PHYAT) ! en 1er
      Call XSSN_Hlike(Te,Ne,'SV_3G4.lab','SV_3G4.res',PHYAT)
      Call XSSN_Hlike(Te,Ne,'SV_3G3.lab','SV_3G3.res',PHYAT)
      Call XSSN_Hlike(Te,Ne,'SV_1G4.lab','SV_1G4.res',PHYAT)
        indice=1
      endif
c
      if(nomion.eq.'SVI    ') then
      Call XSSN_Hlike(Te,Ne,'SVI_2.lab','SVI_2.res',PHYAT)
        indice=1
      endif
c
c *** SVII TBD
c
       If(indice.eq.0) then
         print *,' *** NOT FOUND:',nomion, ' ***'
       Endif
c
        ENDIF  !end yn.ne.0
c
      ENDDO  !k
c
      CLOSE(in_ions)
c
      STOP
c
      END
c----------------------end XSSN_main---------------------
c--------------------------------------------------------
c
      Subroutine XSSN_Hlike(Te,Ne,niv_input,lambda_output,PHYAT)
c
      IMPLICIT NONE
c
      CHARACTER niv_input*(*),lambda_output*(*),PHYAT*(*)
c
c   niv_input : data: level energies, parameters, etc. Ex: 'NV_2.lab'
c   lambbda_output: outputs, including for X-SSN. Ex: 'NV_2.res'
c *NB: minimun relative intensities of lines according to wlgths, etc. in:
c                outputcond.dat, read by:
c
c I) Compute wavelengths for lines with 'high' (n,l) (more or less 'yrast')  
c   using table cor(n,l) = corrections to hydrogenic energies given in : 
c                 niv_input.
c  Only one energy is assumed for each "l" (degenerescence). Otherwise, the same
c   ion is computred several times, with auxiliary tabulation restricted 
c   to lines corresponding to the changing energy level. Cf logics nmult, imult. 
c  Tables start from l=1 (p), followed by l=2 (d), etc., but l=0 (s) 
c   is considered implicitely in the tables (extended fortran 77). 
c  The 'l' considerds are generaly > 1, except for 'exact' H-likes
c   HI and HeII (with 'l' grouped together or not for each 'n'). 
c  Level energy of a level E(n,l) = chi_ion - chi_Hli/n**2 + cor(n,l) in cm-1.
c  Known energies EN(n,l) ('N' for 'NIST') are read from nmil1 to nmal1
c   providing corN(n,l) for known (considered) levels (otherwise corN(n,l)=0).
c  A preliminary run provides corN(n,l) in lambda_output for these known 
c   levels. These corN(n,l) are adopted as cor(n,l) (transcribed "by hand") 
c   in niv_input.
c  Other cor(n,l), hopefully relatively small when they are truly useful, are 
c   obtained by extrapolation of corN(n,l) and empirically from the observation 
c   of lines by means of programme X-SSN applied to actual nebulae.
c  In a second step, EN(n,l) may be empirically refined, taking cor(n,l) 
c   different from 'known' corN(n,l), again with the help of observation
c   and syhthesis. 
c
c II) Extract transitions listed in X-SSN. 
c
c III) Calls :
c                 Function Int_Hlike
c  Compute emissivity of X-SSN lines in quasi-hydrogenic approximation,
c   using effective recombination coefficients listed by  
c   Storey et Hummer (1995) for n < 51, read by:
c                 Subroutine Readsto, 
c  and branching ratios computed in:
c                 Subroutine HBD,
c  adapted from prog BA5.2 of Hoang Binh, Dy: 1990A&A...238..449H.
c   Possibility to extrapolate H-like strict emissivities n = 51 --> 99.
c *For H and He+, cf tabul PSto H_like strict --> ns = 500 (but without l's).
c     *Cf also:
c                 Subr XSSN_di
c   computing di-electronic lines 'alla Nussbaumer & Storey 1984' with NIST wl's.
c *SEE 'Peculiar cases'.
c
c IV) Outputs .res, including relative intensities in X-SSN format.
c  Possibility to group together several orbitals l (narrow blends), 
c   cf notations 'regr'.
c
c V) The .res are transfered in PHYAT (= liste_phyat... for X-SSN) by:
c                 Subr Write_Phyat
c
c VI) Treatments in separate Subr: HeI, ... OIII_5g, etc. 
c
      INTEGER ni_ref,li_ref,ns_ref,k,km !ATT: normally k=1 on calling Int_Hlike
      INTEGER numin,numax !Br table
      DOUBLE PRECISION Mat
      REAL az,am
      CHARACTER zc*1,bidon*1
      INTEGER iz,izz,iprint,i,im
      DOUBLE PRECISION alfeff,b_SH
      DOUBLE PRECISION Br,A_tp
      DOUBLE PRECISION Ne,Te
      DOUBLE PRECISION Int_Hlike,int_ref_Hli
      DOUBLE PRECISION int_ref,lambda_ref,leakns_ref
cc      INTEGER lc1,lc2    ! check string tables
      INTEGER lmax_del1,lmax_del2,lmax_del3
      INTEGER n,nmi,nma,nma_orig,nma2,nh2,l,lmi,lmi_orig,lma
      INTEGER nmil1,nmal1,lmal1,nmil2,nmal2,lmal2
      INTEGER lmin,lmax,lmin_int,lmax_int,lmin_wr
      INTEGER nmitab,nmatab,ni,ns,n1,n2,ini(99)
      INTEGER inmitab,inmatab,iraies
      INTEGER mult,imult,nmult,nimult(20),limult(20),j
      INTEGER niregr,nsregr,indregr,lregr,hydrostrict
      INTEGER inep1
      INTEGER iextrapol
      INTEGER iwlim,iw,iecr,iecrb
      DOUBLE PRECISION ryd,charge,chi_Hli,chi_ion,x,A
      DOUBLE PRECISION cmult,comult(99,0:99)
      DOUBLE PRECISION E(99,0:99),EN(99,0:99),cor(99,0:99),corN(99,0:99)
      DOUBLE PRECISION wl(99,0:99,99),Emisto(0:99) !trans 'normales' ns l+1 --> ni l
      DOUBLE PRECISION wk(50,0:50,50)           !trans eventuelles ns l-1 --> ni l
      DOUBLE PRECISION wlhydrostr(22,500) !wl pr 10 series H et He+ Sto
      DOUBLE PRECISION emhydr(22,500),int_hydr_ref !emis rel H et He+ Sto
      COMMON/MOCR1/wlmin,wlmax,wlim,Irelmin_H,Ionab,Ionab_min,iwlim
      DOUBLE PRECISION int_ref_mod_H
      COMMON/MOCR2/int_ref_mod_H
      DOUBLE PRECISION Ionab,Ionab_min
      DOUBLE PRECISION wlim(20),Irelmin_H(21)
      DOUBLE PRECISION wlmin,wlmax
      DOUBLE PRECISION Irelmin(21)
      COMMON/MOCR3/Irelmin
      DOUBLE PRECISION va,vb,pente1,ve,vf,vg
      DOUBLE PRECISION vh   !pr test H
      DOUBLE PRECISION vzero
      CHARACTER*17 titre
      CHARACTER numion*4,numion4*4,nomion*7 !format: ' 705', '0705', 'N_V    '
      CHARACTER*3 inmitabc,inmatabc,ni_refc,ns_refc
      CHARACTER Irelminc*5,lambda_refc*6,iraiesc*4,wlminc*8,wlmaxc*8
      CHARACTER*1 nom(21),nom2
      CHARACTER tabsto*9
      CHARACTER zero*1,nc1*1,nc2*2
      CHARACTER ns_alph*2,ni_alph*2,lp1_alph*2,ns_alph3*3
      INTEGER ifre1,ifre2,ifre2b,ifre1_ref !caracts supplem No X-SSN de la raie
c  add leakage NIV 4f-5g1 'by hand' :
      DOUBLE PRECISION BrNIV5g1,wlNIV5g1fu
c  add split NeIII 4d-5f 3D-3F 'by hand' :
      INTEGER lp1
      DOUBLE PRECISION wlNeIII5f(3),coNeIII5f(3)
c  add NeIII 4p-4d (mult 5 ou 3) 'by hand' :
      DOUBLE PRECISION wlNeIII4d5(3),wlNeIII4d3(5)
      DOUBLE PRECISION coNeIII4d5(3),coNeIII4d3(5)
      DOUBLE PRECISION BrNeIII4d(2)
      CHARACTER*10 comNeIII4d(2)
c  add NeIII 4s-4p (mult 5 ou 3) 'by hand' :
      DOUBLE PRECISION wlNeIII4p5(3),wlNeIII4p3(3)
      DOUBLE PRECISION coNeIII4p5(3),coNeIII4p3(3)
      DOUBLE PRECISION BrNeIII4p(2)
      CHARACTER*10 comNeIII4p(2)
c  add NeIII : auxil memory
      DOUBLE PRECISION wlNeIII(5),coNeIII(5),BrNeIII
      CHARACTER*10 comNeIII
c  add CIV NV OVI splits np et nd et qlq k=2 'by hand' :
      DOUBLE PRECISION del_Li(11,2,3),wlp,vc,vd ! del_Li(n,l,i_Li)
      CHARACTER*8 com_Li(11,2,3) !comment A/B com_Li(n,l,i_Li)
c  add CIV NV OVI 2s-2p, 3s-3p, 3p-3d 'by hand' :
      INTEGER i_Li ! 1,2,3 CIV,NV,OVI
      DOUBLE PRECISION wl_Li(2),co_Li(2),IBr_Li(3,3),BsA  ! IBr_Li(ns-1+l,i_Li)= Inverse Branching 
c
      INTEGER in_niv,out_niv,iwr_phyat
      COMMON/MOC1/in_niv,out_niv,iwr_phyat
c
      COMMON/MOC2/iprint
      COMMON/MOCS2/inep1,iz,zc,tabsto
      COMMON/MOCH2/az,am,numin,numax
      COMMON/MOCX1/lambda_refc
      COMMON/MOCX2/int_ref,lambda_ref,ifre1_ref
      COMMON/MOCV1/vzero
c For direct computation of a few non Hli trans (NeIII):
      DOUBLE PRECISION hc,hc8,hcsk,Ryd_inf
      COMMON/MOCC1/hc,hc8,hcsk,Ryd_inf
      COMMON/MOCH1/Br(100,100,100,2),A_tp(100,100,100,2) !Br(nu,lup,nl,k) lup=lu+1
      COMMON/MOCS1/alfeff(50,50),b_SH(50,0:50) !alfeff(nu,lup) lup=lu+1
c
      data zero/'0'/
c
c orbital 's' in pos 21:
      data nom/'p','d','f','g','h','i','k','l','m','n','o','q','r'
     $   ,'t','u','v','w','x','y','z','s'/
c
c split NeIII 4d-5f 3D-3F 'by hand'
      data wlNeIII5f/4229.305d0,4226.011d0,4224.029d0/
      data coNeIII5f/0.4667d0,0.3333d0,0.2000d0/  !LS
c NeIII 4p-4d 5 ou 3 'by hand'
      data wlNeIII4d5/6130.55d0,6134.44d0,6141.32d0/
      data coNeIII4d5/0.4667d0,0.3333d0,0.2000D0/  !LS
      data wlNeIII4d3/7178.85d0,7197.32d0,7211.08d0,7188.09d0,7203.09d0/
      data coNeIII4d3/0.4667d0,0.2500d0,0.1111d0,0.0833d0,0.0833d0/  !LS
      data BrNeIII4d/0.55d0,0.33d0/ !fm Aij_PvH: 5 puis 3, 3=CasB,A=B/28.3
      data comNeIII4d/'          ','A=B/28.3  '/
c NeIII 4s-4p 5 ou 3 'by hand'
      data wlNeIII4p5/7584.67d0,7594.45d0,7600.69d0/
      data coNeIII4p5/0.4667d0,0.3333d0,0.2000d0/  !LS
      data wlNeIII4p3/6895.41d0,6886.94d0,6879.65d0/
      data coNeIII4p3/0.5556d0,0.3333d0,0.1111d0/  !LS
      data BrNeIII4p/0.122d0,0.179d0/ !fm Aij_PvH: 5 puis 3
      data comNeIII4p/'          ','          '/
c IBr_Li(ns-1+l,i_Li)= Inverse Branching 
      Data IBr_Li/1.0d0,146.0d0,3.02d4,    ! CIV 2s-2p 3s-3p 3p-3d NIST
     $            1.0d0,292.0d0,4.45d4,    ! NV 
     $     1.0d0,513.0d0,6.32d4/ ! OVI
cc 5/9/17 ***ATT: see E_levels CIV in Tunklev, M. et al. 1997 PhyS...55..707
cc       wlCIV(1)=1548.202 .667?       wlCIV(2)=1550.774  .333
cc       wlCIV(1)=5801.31       wlCIV(2)=5811.97
cc       wlCIV(1)=20705. .4      wlCIV(2)=20796.  .6
c del_Li(11,2,3) del_Li(n,l,i_Li)
      Data del_Li/0.d0,-107.7d0,-31.6d0,-13.1d0, -7.1d0,-3.9d0, !CIV l=1
     $            -2.4d0,-1.6d0,-1.1d0,-0.7d0,-0.4d0,
     $            0.d0,    0.d0,-10.5d0, -4.0d0,-1.7d0,-1.2d0, !CIV l=2 ext n>6
     $            -0.8d0,-0.5d0,-0.3d0,-0.2d0,-0.1d0,
     $            0.d0,-258.7d0,-76.3d0,-32.8d0,-16.1d0,-9.4d0, !NV l=1 ext n>10
     $            -5.9d0,-3.9d0,-2.8d0,-2.0d0,-1.4d0,
     $            0.d0,    0.d0,-22.0d0, -9.2d0,-4.9d0,-3.0d0, !NV l=2 cor n=6, n>8
     $            -2.0d0,-1.4d0,-1.0d0,-0.7d0,-0.5d0,
     $            0.d0,-532.5d0,-156.6d0,-63.9d0,-30.d0,-16.d0, !OVI l=1 ext n>4
     $            -11.0d0,-7.0d0,-5.0d0,-3.5d0,-2.5d0,
     $            0.d0,    0.d0,-51.1d0,-21.4d0,-10.5d0,-6.5d0, !OVI l=2 cor n>4
     $            -4.5d0,-3.0d0,-2.0d0,-1.5d0,-1.0d0/
       Data com_Li/'','','','','','A=B/2.8','A=B/2.6', !CIV l=1
     $             'A=B/2.5','A=B/2.4','','',
     $             '','','','','','','','','','','', !CIV l=2
     $             '','','','','','A=B/2.8','A=B/2.6', !NV l=1
     $             'A=B/2.5','A=B/2.4','','',
     $             '','','','','','','','','','','', !NV l=2
     $             '','','','','','A=B/2.8','A=B/2.6', !OVI l=1
     $             'A=B/2.5','A=B/2.4','','',
     $             '','','','','','','','','','',''/ !OVI l=2
c
c ryd C12,N14,O16: 109732.30d0,109733.016d0,109733.554d0
c21/3/11 rydberg ion H-like AND normal isotopic composition (NIST) : 
c    H  109677.6228d0/.6155Clegg, He 109722.2755d0/.2772Clegg
c    C  109732.303d0, N  109733.017d0, O 109733.552d0, Ne 109734.332d0
c    Mg 109734.838d0, Si 109735.172d0, infinite 109737.3157d0
c Q: relativistic Dirac corrections () TBD?
c
c Energies NIST EN() :
c Rem : EN=0.d0 for unusefully small 'l' or unknown energies
c
c Correc atm refrac x=lambda_vac in um, dry air 15C p=760mmHg, Allen_3, p.124
      A(x)=1.+(64.328+29498.1/(146.-1./x**2)+255.4/(41.-1./x**2))*1.d-6
c Check:
cc      lc1=LEN(niv_input)
cc      lc2=LEN(lambda_output)
cc      print *,niv_input,'  ',lambda_output,'   (sizes:',lc1,lc2,')'
c
      print *,niv_input,'  ',lambda_output
c
      OPEN(unit=in_niv,file=niv_input,status='old')
c
      read(in_niv,15) titre
 15   format(A17)
c rydberg, charge, ionization potential, atomic mass (uma)
      read(in_niv,20) ryd,charge,chi_ion,Mat
 20   format(4x,d12.3,7x,d4.0,8x,d12.3,4x,d5.0)
c  mult = multiplicity, nmult = nber of lower lev with mult, cmult = coeff mult
c  imult= 1,2,.. ==> 2nd,3rd,.. multiplicity; only lower levels with mult shown
c Rem: imult sert aussi bien pour des sous-niv differents dans une meme multiplicite
c Rem: if nmult>0 for an ion, only imult=0 (TO BE OUT 1st) gives the ref line. 
c  ATT: imult also used in line labeling ==> 0 < imult < 9 
c       (not necessarily with nmult > 0).
c  ifre1,2 (0 to 9, std = 0) : for further differenciation in line labeling
c     where needed.
      read(in_niv,21) mult,nmult,imult,cmult,ifre1,ifre2
 21   format(6x,i4,6x,i4,6x,i4,6x,d7.0,6x,i4,6x,i4)
       if(nmult.gt.9) then  !ATT: nmult < 10
          print *,' ********** nmult =',nmult, ' > 9 STOP *******'
          STOP
       endif
       if(nmult.gt.0) then
c  ATT: provide the nimult's by increasing order ***
        read(in_niv,23) (nimult(j),limult(j),j=1,nmult)
 23   format(10(i4,i3))
c  (check nimult)
          if(nmult.gt.1) then
           do j=1,nmult-1
            if(nimult(j+1).lt.nimult(j)) then
         print *,' ***STOP: in',niv_input,' NIMULT TO BE RE-ORDERED **'
             STOP
            endif
           enddo
          endif
c  (end check)
       else   ! nmult=0
        read(in_niv,1) bidon
 1    format(A1)
        nimult(1)=0
        limult(1)=0
cc 1/14 suppr check suiv car imult.ne.0 possible sans 'multiplicite'
cc      (eg, convention imult=2 pour di-el)
cc         if(imult.ne.0) then
cc          print *,'****** imult =',imult,' put to 0 (cf nmult =',nmult,' )'
cc          imult=0
cc         endif
       endif
c potential H-like cm-1 :
      chi_Hli=ryd*charge*charge
c charge & mass (atomic number = izz) :
      iz=IDINT(charge)
c Conversion to character (note: iz > 8 treated below)
      write(zc,'(I1)') iz
c Check:
      If(iz.lt.1) then
       print *,'**** iz=',iz,' < 1 STOP***'
       STOP
      Endif
      If(iz.gt.8) then
       print *,'*** iz=',iz,' : zc put to 8 for Storey tables.'
       print *,'Scaling factor applied to recombination coefficients'
       zc='8'
c  (Note: Br has little dependence on az and am)
      Endif
	az=REAL(iz)
	am=Mat
c Table branching ratios Br :
      read(in_niv,66) numin,numax
 66   format(6x,i4,6x,i4)
c Check:
      if(numax.le.numin) then
        print *,'numax=',numax, 'le. numin=',numin,'***STOP***'
        STOP
      endif
c
c Label + name ion (form X-SSN): numion=' 705' numion4='0705' nomion='N_V    '
      read(in_niv,60) numion,numion4,nomion
 60   format(9x,a4,11x,a4,10x,a7)
c
         BACKSPACE(in_niv)  ! recule d'une ligne
        read(in_niv,160) izz ! atomic nber
 160   format(9x,i2)
c
c Transition ref emissivity X-SSN :
      read(in_niv,22) ni_ref,li_ref,ns_ref,hydrostrict
 22   format(20i4)
c Ref. int: sum over l fm li_ref to ni_ref-1. Check :
       if(li_ref.lt.0.or.li_ref.gt.ni_ref-1) then
        print *,'**** li_ref=',li_ref,' out of range. STOP**'
        STOP
       endif
c Correction for possible 'leakage' of ns_ref lev
      read(in_niv,25) leakns_ref ! +/- useful now?
 25   format(d9.3)

c Lim. n & l calcul
c  IF mult > 0, choose lma according to effective limult (hence lma < 100)
      read(in_niv,22) nmi,nma,lmi,lma
c nma2 keep initial value read for nma (nma = maximum n considered):
        nma2=nma !tabulation PSto nl-nu H et He+ up to nl=10 and nu=500
c Lim. for energy reading: lab (NIST, etc.) then effectivs (X-SSN)
      read(in_niv,22) nmil1,nmal1,lmal1,nmil2,nmal2,lmal2
c Lim. for tables X-SSN and limites for l grouping:
      read(in_niv,22) nmitab,nmatab,niregr,nsregr
c
c Control interval l of computed transitions:
c  lmin=max0(lmi,lmax-lmax_del1) general tabul
c  lmin=max0(lmi,lmax-lmax_del2) X-SSN tabul, *ATT: lmax_del2.le.lmax_del1
c  Only the largest l's for each ni: 3 (lmax_del=2), 4 (lmax_del=3)...
      read(in_niv,22) lmax_del1,lmax_del2 !general, then X-SSN
c
c 'average' wl adopted for X-SSN lines grouped together in one only:
c  lmax=min0(lma,ni-1-lmax_del3)
c  lmin=lmax
c (Rem: further weighting below)
      lmax_del3=min0(lmax_del2/3,3) !lmax_del3=1: penultimate orbital, etc.
c
c SECU (dimensions, writings)
       if(nma.gt.99.and.
     $    numion.ne.' 101'.and.numion.ne.' 202') then ! 3 nbers EXPLICIT
         print *,'*** ',nomion,': nma=',nma,' LEVELED OFF AT 99 ***'
         nma=99
       endif
        nma_orig=min0(nma,99) !for H-like extrapol ns>50 (except H, He+)
       if(hydrostrict.eq.2.and.
     $    (numion.eq.' 101'.or.numion.eq.' 202')) then
c tabulation PSto nl-nu H and He+ up to nl=10 and nu=500 : nma2=nma higher
        nma=nsregr !(no splitting of multiplets for ns>nsregr)
c (**Rem 1/15 : below, hydrostrict.eq.2 only for HI and HeII. Correct logics?)
       endif
c  Controls possible extrapolation of emissivities:
        iextrapol=0
c For lack of alfeff, emissivities are extrapolees for ns = 50 - 99 in H-like 
c  (except for H et He+, treated by tabulation PSto nl,nu) :
        if(nma.gt.50) then
         nma=50
          if(nsregr.le.40) then !extrapol only if basis sufficient
         iextrapol=1
          else
      print *,'nma=',nma_orig,' nsregr=',nsregr,' >40, nma=50 no extrap'
          endif
        endif
cc No extrapol for HI HeII. HeI treated separately. 
       if(numion.eq.' 101'.or.numion.eq.' 202') iextrapol=0 ! 2 EXPLICIT: H, He
c INIT   (!REM: lmi and lmin = 0 possible; extended gfortran)
      do n=1,nma
        do l=0,n-1
         E(n,l)=0.d0
         EN(n,l)=0.d0
         cor(n,l)=0.d0
         corN(n,l)=0.d0
         comult(n,l)=1.d0
        enddo
        ini(n)=0
      enddo
c
      do n=nmil1,nmal1
      lmax=min0(lmal1,n-1)
       read(in_niv,24) (EN(n,l),l=lmi,lmax)
 24    format(6d12.3)
        do l=lmi,lmax
         E(n,l)=EN(n,l)
         corN(n,l)=E(n,l)-chi_ion+chi_Hli/dble(n*n)
        enddo
      enddo
c
c  Computation done for a defibnit ion and a definite multiplicity. 
c When 2 multiplicities exist, eg, singlet+triplet, the computation is 
c performed several times, applying a multiplicative coefficient cmult  < 1
c to the only lower levels (ni,li) whose energies change according to 
c multiplicity. The ouput params of X-SSN aere chosen so as to prevent 
c line duplications. 
c
      if(nmult.gt.0) then
       do j=1,nmult
        n=nimult(j)
        l=limult(j)
        comult(n,l)=cmult
       enddo
      endif
c
      do n=nmil2,nmal2
      lmax=min0(lmal2,n-1)
       read(in_niv,26) (cor(n,l),l=lmi,lmax)
      enddo
 26    format(10d8.2)
c
      CLOSE(in_niv)
c
      do n=nmi,nma
        do l=lmi,n-1
         E(n,l)=chi_ion-chi_Hli/dble(n*n)+cor(n,l)
        enddo
      enddo
c
c Tabulation of wavelengths wl in AA with
c  lower lev ni fm nmitab to nmatab and upper lev ns fm ni+1 to nma 
c Normaly, nma < 100, except for HI and HeII possibly ->500  by PSto without l's.
c If needed, FIT nmatab to range of useful wl. **Cf also multiplicities. 
c
      DO ni=nmitab,nmatab !tabul wl, wk
c
      lmax=ni-1
c If only the greatest l for each ni, eg: 3 (lmax_del1=2), 4 (lmax_del1=3)...
c HeI: all of the 'l' : take large lmax_del1 (eg, 100) and lmi = 0
      lmin=max0(lmi,lmax-lmax_del1) !REM: lmi and lmin = 0 possible
c IMPLICITLY HERE** : 
c   'l' corresp. to an orbital of lower level ni
c   wl: spontaneous transition is l+1 --> l
       Do l=lmin,lmax  ! l of lower level (ni)
         do ns=ni+1,nma
          va=E(ns,l+1)-E(ni,l)
          if(va.gt.0.d0) then
           wl(ni,l,ns)=1.d8/va
c   air refraction wl > 2000A AND wl < 20000A:
            IF(wl(ni,l,ns).gt.2.d3.and.wl(ni,l,ns).lt.2.d4)THEN
             wl(ni,l,ns)=wl(ni,l,ns)/A(1.d4/va)
            ENDIF
          else !va </= 0 pathol
           wl(ni,l,ns)=1.d0
       print *,nomion,' wl : ni,l,ns',ni,l,ns,'** E(ns,l+1)-E(ni,l)=',va
       print *,'****STOP*****************************'
       STOP
          endif
         enddo
       Enddo
c 13/3/14 for HeI, CIV, etc.
c IMPLICITLY HERE*: 
c   'l' corresp. to an orbital of lower level ni
c   wk: spontaneous transition is l-1 --> l
        If(lmin.lt.lmax) Then
       Do l=lmin+1,lmax  ! l of lower level (ni)
         do ns=ni+1,nma
          va=E(ns,l-1)-E(ni,l)
          if(va.gt.0.d0) then
           wk(ni,l,ns)=1.d8/va
c   air refraction wl > 2000A AND wl < 20000A:
            IF(wk(ni,l,ns).gt.2.d3.and.wk(ni,l,ns).lt.2.d4)THEN
             wk(ni,l,ns)=wk(ni,l,ns)/A(1.d4/va)
            ENDIF
          else !va </= 0 pathol
           wk(ni,l,ns)=1.d0
       print *,nomion,' wk : ni,l,ns',ni,l,ns,'** E(ns,l-1)-E(ni,l)=',va
       print *,'****STOP*****************************'
       STOP
          endif
         enddo
       Enddo
        Endif  !lmin.lt.lmax
c
      ENDDO ! DO ni (tabulation wl, wk)
c
c Tabulation wlhydrostr (without l) for H or He+; 
c  useful for hydrostrict=2, but nalso for int corr
c  (P Storey intensities also for large ns). In parallel with wl(). 
c**** SECU hydrostrict :
      IF((numion.ne.' 101'.and.numion.ne.' 202').and.
     $    (hydrostrict.eq.2)) Then      ! 2 nbers EXPLICIT
       print *,'**** numion=',numion,' : hydrostrict 2 --> 0 ****'
        hydrostrict=0  !ATT 1/15: REALLY only H and He+ ? and why not =1 ? VOIR
      ENDIF
c
cccc    If(hydrostrict.eq.2) then !If removed: systematic computation for int cor
       do ni=nmitab,nmatab
         do ns=ni+1,nma2   !nma2 < 500 defini +haut (nma lu)
          va=chi_Hli/dble(ni**2)-chi_Hli/dble(ns**2)
          wlhydrostr(ni,ns)=1.d8/va
c   air refraction wl > 2000A AND wl < 20000A:
          IF(wlhydrostr(ni,ns).gt.2.d3.and.wlhydrostr(ni,ns).lt.2.d4)
     $    THEN
           wlhydrostr(ni,ns)=wlhydrostr(ni,ns)/A(1.d4/va)
          ENDIF
         enddo
       enddo
cccc    Endif
c
      OPEN(unit=out_niv,file=lambda_output,status='old')
c
      write(out_niv,15) titre
      write(out_niv,40) ryd,charge,chi_ion,chi_Hli,Mat,iz
 40   format(1p,' ryd =',e17.10,' charge =',e8.1,' chi_ion =',e17.10,
     $  /,' chi_Hli =',e17.10,' Mat =',e8.1,'  (iz =',0p,i3,')')
      write(out_niv,41) mult,nmult,imult,cmult,ifre1,ifre2
 41   format('  mult =',i4,' nmult =',i4,' imult =',i4,
     $   ' cmult =',f7.4,' ifre1 =',i4,' ifre2 =',i4,' ni,li-mult :')
      if(nmult.gt.0) then
        write(out_niv,23) (nimult(j),limult(j),j=1,nmult)
      else
        write(out_niv,23) nimult(1),limult(1)  !pour eviter ligne blanche
      endif
      write(out_niv,67) Ne,Te
 67   format(1p,' Ne =',e9.2,'  Te =',e9.2, ' alfeff')
      write(out_niv,68) numin,numax
 68   format(2i4,'  numin numax Br (HBD)')
      write(out_niv,62) numion,numion4,nomion
 62   format(' numion :',a4,' numion4 :',a4,' nomion :',a7,' a4,a4,a7')
      write(out_niv,63) wlmin,wlmax
 63   format(1p,2e10.3,' wlmin, wlmax X-SSN')
c This modified nma corresponds to explicit treatmt  of 'multiplets' (l)
c By def, the larges ns are only possible with grouped l's. 
c nma2 still useful for HI HeII (PSto treatmt of large n's)
      write(out_niv,42) nmi,nma,lmi,lma,nma2
 42   format(5i4,'  nmi,nma,lmi,lma,nma2 = nma lu')
      write(out_niv,72) nmil1,nmal1,lmal1,nmil2,nmal2,lmal2
 72   format(6i4,'  nmil1,nmal1,lmal1,nmil2,nmal2,lmal2 lectures')
      write(out_niv,70) ni_ref,li_ref,ns_ref,hydrostrict
 70   format(4i4,'  ni_ref,l=li_ref->ni_ref-1,ns_ref,hydrostrict=0-2')
      write(out_niv,71) leakns_ref
 71   format(1p,d10.3,'  fact. leakns_ref')
      write(out_niv,64) nmitab,nmatab,niregr,nsregr
 64   format(4i4,'  nmitab,nmatab,niregr,nsregr X-SSN')
      write(out_niv,46) lmax_del1,lmax_del2,lmax_del3
 46   format(3i4,'  lmax_del1,lmax_del2,lmax_del3')
      write(out_niv,43) 
 43   format(4x,'EN(n,l),corN(n,l)  Lab')
c
      do n=nmil1,nmal1
      lmax=min0(lmal1,n-1)
       write(out_niv,44) n,(EN(n,l),l=lmi,lmax)
 44   format(i3,(10f12.4))
       write(out_niv,45) n,(corN(n,l),l=lmi,lmax)
 45   format(i3,(10f12.4)) !(i3,(10(4x,f8.4))) 11/13
      enddo
c (corN repeated to make easier comparison/transcription to cor)
c *Q: change format for convenience?
      write(out_niv,49) 
 49   format(4x,'corN(n,l)  Lab')
      do n=nmil1,nmal1
      lmax=min0(lmal1,n-1)
      write(out_niv,51) n,(corN(n,l),l=lmi,lmax)
 51   format(i3,6x,(10f12.4))  !(i3,6x,(10f8.4))
      enddo
c cor = correction to H-like effectively adopted
      do n=nmil2,nmal2
      lmax=min0(lmal2,n-1)
       write(out_niv,48) n,lmi,lmax,(cor(n,l),l=lmi,lmax)
      enddo
 48   format(3i3,(10f12.4))  !(3i3,(10f8.4))
c**
        IF(iprint.ne.0) THEN ! ** iprint NOT UPDATED (obso ?)
      write(out_niv,83) 
 83   format(4x,'E(n,l),l=lmi,lmax')
      do n=nmil2,nmal2
      lmax=min0(lmal2,n-1)
       write(out_niv,47) n,lmi,lmax,(E(n,l),l=lmi,lmax)
      enddo
 47   format(3i3,(10f12.4))
c
      write(out_niv,73) 
 73   format('  ni,ns,lmin,lmax,   wl(ni,l,ns),l=lmin,lmax')
      do ni=nmitab,nmatab
      lmax=ni-1
      lmin=max0(lmi,lmax-lmax_del1)
       do ns=ni+1,nma
       write(out_niv,50) ni,ns,lmin,lmax,(wl(ni,l,ns),l=lmin,lmax)
       enddo
      enddo
 50   format(4i3,(11f10.3))
        ENDIF !iprint end
c**
c************TRIAL lmi
c      lmi=0
c******
c  **********************
c  Tabulation for X-SSN : 
c  **********************
c
c*****  Tabulation of alfeff:  (n,l, n->50)
      Call ReadSto(Te,Ne)
c**     Computation/tabulation of branching ratio Br:
      Call HBD
c**
c****   Emissivities (ni, ns) of PSto for H and He+. (Possibly 
c   modified LATER for splitting of 1st lines, ie if n<nregr). 
      IF(numion.eq.' 101'.and.hydrostrict.eq.2) Then ! EXPLICIT
c  H0: reading e1bx.d  (n, n->500)
       Call hdatx
         Call hlinex(ns_ref,ni_ref,Te,Ne,int_hydr_ref,1)
       Do ni=nmitab,nmatab
        Do ns=ni+1,nma2
         Call hlinex(ns,ni,Te,Ne,va,1) !iopt=1 for linees
         emhydr(ni,ns)=va/int_hydr_ref !usually 'vb'
        Enddo
       Enddo
      ENDIF
c
      IF(numion.eq.' 202'.and.hydrostrict.eq.2) Then ! nbre EXPLICIT
c  He+: reading e2bx.d  (n, n->500)
       Call hepdatx
         Call heplinex(ns_ref,ni_ref,Te,Ne,int_hydr_ref,1)
       Do ni=nmitab,nmatab
        Do ns=ni+1,nma2
         Call heplinex(ns,ni,Te,Ne,va,1) !iopt=1 pour raies
         emhydr(ni,ns)=va/int_hydr_ref !vb habituel
        Enddo
       Enddo
      ENDIF !end H He+, hydrostrict.eq.2
c
c Labeling treatment of possible 2nd multiplicity of the ion.
c  same nmitab,nmatab,nmult,nimult,limult assumed / 1st mult.
c imult > 0 ==> only sequences (nimult,limult) in output X-SSN
c
c General ref line X-SSN : 
c     IF several times the same ion, lambda_ref is from imult = 0 
c  lambda_ref, int_ref, ifre1_ref, in COMMON, are known for imult > 0
c    (ATT.: imult = 0 computed 1st)
c
          IF(imult.eq.0) THEN
             ifre1_ref=ifre1  !anticipating write imult > 0
c
c Reference H-like transition:
      if(li_ref.eq.0) then
c   k=1 et 2 (H like complet), if all l's considered
       km=2
      else
c   k=1 only, in case only large l's considered
       km=1
      endif
      vb=0.
      DO k=1,km
       DO l=li_ref,ni_ref-1
        vb=vb+Int_Hlike(iz,Mat,ni_ref,l,ns_ref,k,Te)
       ENDDO
      ENDDO
c Emissivity and wavelength of H-like ref line:
      int_ref_Hli=vb
cc9/17 out:      va=chi_Hli*(1./dble(ni_ref**2)-1./dble(ns_ref**2))
c Effective wl of ref line: 
      lambda_ref=wl(ni_ref,ni_ref-1,ns_ref)
c Corr for non-H-like wl and possible 'leakage' du niv sup: (fuite: slmt CII)
      int_ref=int_ref_Hli*leakns_ref*
     $     wlhydrostr(ni_ref,ns_ref)/lambda_ref !int_ref en common
c
c Storage of Hbeta emissivity for normalisations and choice of Irelmin's:
      if(numion.eq.' 101'.and.ifre1.eq.0) then
      int_ref_mod_H=int_ref
      endif
c
c Wavelength ranges:
c
      do iw=1,iwlim+1
         Irelmin(iw)=Irelmin_H(iw)*
     $   int_ref_mod_H/(int_ref*dmax1(Ionab,Ionab_min))
      enddo
c     
          ENDIF  !end imult.eq.0
c
        if(li_ref.eq.0) li_ref=21  !for nom(0) undefined
        nom2=' '
        if(li_ref.ne.ni_ref-1) nom2=nom(ni_ref-1)
c
      write(out_niv,85) iwlim
 85   format(' iwlim ='i4,' then wlim(1->iwlim), Irelmin(1->iwlim+1)')
      write(out_niv,87) (wlim(iw),iw=1,iwlim)
 87   format(5x,1p,13e10.3)
      write(out_niv,86) (Irelmin(iw),iw=1,iwlim+1)
 86   format(1p,13e10.3)
c
      Write(out_niv,105) ni_ref,nom(li_ref),nom2,ns_ref,
     $  wlhydrostr(ni_ref,ns_ref),int_ref_Hli
 105  Format(' X-SSN: Ref =',i2,2a1,'-',i2,' lambda H-like =',
     $  1p,e12.6,' emissivity H-like =',e10.4,' erg.cm3.s-1')
      Write(out_niv,106) lambda_ref,int_ref
 106  Format(' X-SSN: effective lambda reference =',
     $  1p,e12.6,' effective emissivity =',e10.4,' erg.cm3.s-1')
c P Storey for H, He+ with interpolation (Te, Ne):
      If(hydrostrict.eq.2) Then
       Write(out_niv,107) int_hydr_ref
 107   Format(' **Emissivity (interpol) P Storey =',1p,e10.4,
     $   ' (NOT adopted)')
      Endif
c
c  START WRITING FOR SYNTHESIS: 
c   START ref lin or sub-ref
       IF(imult.eq.0) THEN  !ref
c
      Write(out_niv,210)  ! imult = 0
 210  Format(' Ref. for liste_modele.dat :')
      Write(out_niv,110) numion,ifre1_ref,nomion,vzero,int(lambda_ref),
     $  ni_ref,nom(li_ref),nom2,ns_ref  !!!,int_ref
 110  Format(1x,A4,I1,'00000000',1x,A7,'        1.0   0.     1.000e+0',
     $  '  1.                 0   1  ',f5.1,1x,I5,'++',
     $  I2,2A1,'-',I2)   !!!,1p,e10.3,' erg.cm3.s-1')
      Write(out_niv,211)
 211  Format(' Ref. for liste_phyat.dat :')
        Write(out_niv,111) numion4,ifre1_ref,nomion,int_ref,
     $  int(lambda_ref),ni_ref,nom(li_ref),nom2,ns_ref,tabsto,inep1
 111  Format('9',A4,I1,'00000000',1x,A7,'        1.0   0.   ',
     $  1p,e10.3,'  1.               999   1   1.00 ',I5,'++',
     $  I2,2A1,'-',I2,1x,'(',a9,i3,')')
c
       ELSE   ! imult > 0  sub-ref
c
      Write(out_niv,212)
 212  Format(' Sub-Ref. for liste_phyat.dat :')
      Write(out_niv,112) numion,ifre1,imult,ifre2,nomion,
     $  numion,ifre1_ref,int(lambda_ref),
     $  ni_ref,nom(li_ref),nom2,ns_ref
 112  Format(1x,A4,3I1,'000000',1x,A7,'        1.0   0.     1.000e+0',
     $  '  1.     ',A4,I1,'00000000   1   1.00 (',I5,
     $  I2,2A1,'-',I2,')')
c
       ENDIF  ! imult 0 ou > 0
c   END ref line or sub-ref
c
        if(li_ref.eq.21) li_ref=0 !reset after writing
c
c line counter:
               iraies=0
c
                  n1=nmitab
                  n2=nmatab
c  imult > 0 --> 2nd, 3rd 'multiplicity', etc, PROVIDED THAT nmult > 0 AS WELL
               If(imult.gt.0.and.nmult.gt.0) Then
c  ATT: nimult read by increasing ordrer
c  ATT: CHOOSE lma in accordance with effective limult's (thus lma < 100)
c  ATT: this logics produces a 'square' that may strictly contain 
c   the (ni,li) to be kept; +discriminant test later. 
                  n1=nimult(1)
                  n2=nimult(nmult)
c  looking for a 'useful' lmi, given the limult's:
c  **possible CHANGE of lmi
                  lmi=lma
                do j=1,nmult
                  lmi=min0(lmi,limult(j))
                enddo
               Endif
                  print *,'  Practical limits : n1,n2,lmi',n1,n2,lmi
cc#####
c********  ********** ************
	   Do ni=n1,n2 !ni lower lev******* loop ni + external
c********  ********** ************
c#####
c For Write:
        if(ni.lt.10) then
       write(nc1,'(I1)') ni
        ni_alph=zero//nc1
        else   !(assume ni < 100)
       write(ni_alph,'(I2)') ni
        endif
c
        if(iextrapol.eq.1) then  !H-like (H and He+ excluded)
         do ns=1,nma_orig !nma_orig < 100 (potentially reduced, cf above)
          Emisto(ns)=0.d0
         enddo
        endif
c
c lma < n-1 to be used only without grouping, e.g., 
c *********            ****
c   for 2nd multiplicities, such that mult=1 completing mult=3;
c   otherwise, choose lma large (100) to inhibit this possibility.
      lmax=min0(ni-1,lma)
      lmin=max0(lmi,lmax-lmax_del2) !imposed: lmax-lmin < 9 ; 2/15: still this way?
c  useful only when lines are grouped: domaine l regroupes
        lmax_int=lmax
        lmin_int=lmin
c  for write only :
        lmin_wr=lmin_int  !lmin_wr for nom() of smaller l of the sum over l
           if(lmin_wr.le.0) then
c  nom(0) does not exist: 's' in nom(21)
            lmin_wr=21
           endif
c
c************************* (here nma < 50, cf extrapol)
       Do ns=ni+1,nma !ns = niv sup  ****** loop ns intermediary 1
c**************** Rem hydrostrict=2 : nma=nsregr, prevent interference with indregr
c For Write:
        if(ns.lt.10) then
       write(nc1,'(I1)') ns
        ns_alph=zero//nc1
        else   !(assume ns < 100)
       write(ns_alph,'(I2)') ns
        endif
c
c  Re-initialise line grouping index:
       indregr=0
        If(ni.ge.niregr.or.ns.ge.nsregr) Then 
c  One only line per (ni,ns) if ni > niregr-1 or ns > nsregr-1
          indregr=1  !VOIR : ? interference with hydrostrict=1 or 2 ?
c 'Average' wl adopted for grouped lines:
cc **Here, last l (yrast) imposed for the 'average':
c    (perhaps inappropriate for large charges, see alternative)
c ***1/15 wl grouping inappropriate for HI and HeII? 
c       Average over lower sub-levels required, NOT lmin=lmax=ni-1
          lmax=ni-1
         if(charge.gt.5.d0) Then  !*** alternative to be watched at *** :
          lmax=max0(0,ni-1-lmax_del3,(4*(ni-1))/5)
         endif
          lmin=lmax !here, notation lmin lmax for write only
c  and emissivity = sum over l fm lmin_int to lmax_int
        Endif
c
c  Auxil variable vf useful only for 'J' treatment of HI and HeII components:
           vf=0.
c (for HI and HeII, required: lmax= ni-1, lmin= 0 and hydrostrict > 0)
c 1/15 DEF hydrostrict = 1 : int of each line wk added to int of corresponding wl
c
c*************** l associated to ni ***
        DO l=lmax,lmin,-1  ! ******** loop l internal (lmin=0 acceptable in fortran)
c***************
c For Write:
        if(l+1.lt.10) then
       write(nc1,'(I1)') l+1
       lp1_alph=zero//nc1
        else   !(assume l+1 < 100)
       write(lp1_alph,'(I2)') l+1
        endif
c
c********************
c suppressing line duplication:
         IF(imult.eq.0.or.(imult.gt.0.and.nmult.eq.0).or.
     $ (imult.gt.0.and.comult(ni,l).lt.1.d0)) THEN !imult,nmult, **ATT: comult VOIR?
c********************
c**********
c Lambda range:
c**********
          IF(wl(ni,l,ns).gt.wlmin.and.wl(ni,l,ns).lt.wlmax) THEN !wl range
c **Call Function Int_Hlike :
c  (k=1: spont trans lu=ll+1 --> ll +intenses, only considered here)
         k=1  ! Rem : k=2 considered in a few special cases, eg Li_like, cf wk 
c***###***
           If(indregr.eq.0) Then !indregr=0, no line grouping
c***###***
              ve=Int_Hlike(iz,Mat,ni,l,ns,k,Te)
           vb=ve*comult(ni,l)/int_ref
ccc Case hydrostrict.eq.1 or 2: all transitions; k=2 put in same ll (convention)
            if(hydrostrict.gt.0) Then  !(.gt.0 required for HI and HeII; 1 optionnel otherwise)
          if(l.gt.0) then
            k=2  !(int. k=2 added to usual int. k=1)
              ve=Int_Hlike(iz,Mat,ni,l,ns,k,Te)
             vb=vb+ve*comult(ni,l)/int_ref
            k=1
          endif
            endif
c intensity correction for wl departure:
           vb=vb*wlhydrostr(ni,ns)/wl(ni,l,ns)
c  potential leak ns_ref (**obsolete? CII ref now according to DSK):
            if(ni.eq.ni_ref.and.l.ge.li_ref.and.ns.eq.ns_ref) then
             vb=vb*leakns_ref
            endif
c
c  INTENSITIES IN SPECIAL CASES :
c
c HI, DI, HeII and 3HeII: wl, AND THEREFORE INTENSITES redistributed 
C              according to the J's (NOT the l's) of lower level ni
               If(numion.eq.' 101'.or.numion.eq.' 202') Then !HI et HeII
cc            print *,'numion,ni,l,ns,vb,vg,vf',numion,ni,l,ns,vb,vg,vf
c On utilise chaque fois le vb ant et le vb actu. Pourtant une fois sur 2 le 
c    vb ant n'est-il pas associé à une raie de ns diff ??********
                 vg=vb !***? 1/15 logicS: 1st vf=0. 11/15 check: sum multiplets is conserved
                 vb=(l+1.)*(vf/(2.*l+3.)+vg/(2.*l+1.)) !OK. Explain formula?
                 vf=vg
               Endif    !end special case HI and HeII
c
c Eliminate / correct CII lines already computed according to DSK:
c  3d-nf and 4d-nf, remove n=4,9 and correct n>9. 
c  Prefer ref int 4267 by DSK.
                If(numion.eq.' 602') Then ! 7 nbres EXPLICIT
                 If((ni.eq.3.or.ni.eq.4).and.l.eq.2) Then  ! l of lower level (ni)
                  If(ns.ge.4.and.ns.le.9) Then
                   vb=-1.d0
                  Elseif(ni.eq.3) Then  !ns>9  (No correction for ni=4)
                   vb=vb/1.1
                  Endif
                 Endif
                 If((ni.eq.3.or.ni.eq.4).and.l.eq.1) Then  ! l of lower level (ni)
                  If(ns.ge.4.and.ns.le.7) Then
                   vb=-1.d0
                  Endif
                 Endif
                Endif   !end special case CII
c  END INTENSITY TREATMENT OF SPECIAL CASES
c
c  START WRITINGS:
c
c Lower lim of relative intensity depends on wl range:
          iecr=0
c
         iw=1
         if(wl(ni,l,ns).le.wlim(iw).and.vb.ge.Irelmin(iw)) then
          iecr=1
         endif
        do iw=2,iwlim
         if(wl(ni,l,ns).gt.wlim(iw-1).and.wl(ni,l,ns).le.wlim(iw)
     $     .and.vb.ge.Irelmin(iw)) then
          iecr=1
         endif
        enddo
        if(wl(ni,l,ns).gt.wlim(iwlim).and.vb.ge.Irelmin(iwlim+1)) then
          iecr=1
       endif
c
c*****##
             IF(iecr.eq.1) THEN !1st iecr=1
c*****##
c counter of lines:
               iraies=iraies+1
c counter of effectively written ni's:
               ini(ni)=1
                lmin_wr=l
               if(lmin_wr.le.0) then
c  nom(0) does not exist
                lmin_wr=21
               endif
c
c    BEGIN SPECIAL CASES OF NON-H-LIKE WRITE
c
c   Supplementary or split lines
c  Replacing the std wl() by several nearby wl's
c   and write without re-testing modified wl or va.
c
c A) NeIII                 split NeIII 4d-5f 3D-3F :
                If(numion.eq.'1003'.and.mult.eq.3.and.
     $              ni.eq.4.and.l.eq.2.and.ns.eq.5) Then ! split NeIII
                iraies=iraies+2  ! (1 raie --> 3)
              ifre2b=ifre2
                 do i=1,3
              ifre2=i+1
        Write(out_niv,187) numion,ifre1,imult,ifre2,ni_alph,ns_alph,
     $    lp1_alph,nomion,wlNeIII5f(i),vb*coNeIII5f(i),numion,
     $    ifre1,imult,ifre2b,ni,nom(lmin_wr),mult,ns,nom(l+1)
                 enddo
              ifre2=ifre2b
c
c B) NIV(5g1) :
                Elseif  ! end split NeIII  ! begin leak NIV(5g1)
     $             (numion.eq.' 704'.and.mult.eq.1.and.
     $              ni.eq.4.and.l.eq.3.and.ns.eq.5) Then
                iraies=iraies+1  ! (1 line --> 2)
              ifre2b=ifre2
              ifre2=ifre2+1
          BrNIV5g1=0.24
          wlNIV5g1fu=2080.3
        Write(out_niv,187) numion,ifre1,imult,ifre2,ni_alph,ns_alph,
     $    lp1_alph,nomion,wl(ni,l,ns),vb*BrNIV5g1,numion,
     $    ifre1,imult,ifre2b,ni,nom(lmin_wr),mult,ns,nom(l+1)
              ifre2=ifre2+1
        Write(out_niv,187) numion,ifre1,imult,ifre2,ni_alph,ns_alph,
     $    lp1_alph,nomion,wlNIV5g1fu,
     $    vb*(1.-BrNIV5g1)*wl(ni,l,ns)/wlNIV5g1fu,numion,
     $    ifre1,imult,ifre2b,ni,nom(lmin_wr),mult,ns,nom(l+1)
              ifre2=ifre2b
c
c C) Li_like 1 :
                Elseif  ! end NIV(5g1)  ! split Li_like np trans s-p
     $ ((numion.eq.' 604'.or.numion.eq.' 705'.or.numion.eq.' 806').and.
     $ l.eq.0.and.ns.le.9) Then
c
            i_Li=izz-5  ! 1 CIV, 2 NV, 3 OVI
c split trans s-p :
c del_Li(11,2,3) del_Li(n,l,i_Li)
                iraies=iraies+1  ! (1 raie --> 2)
              ifre2b=ifre2
              ifre2=ifre2+1
        Write(out_niv,187) numion,ifre1,imult,ifre2,ni_alph,ns_alph,
     $    lp1_alph,nomion,wl(ni,l,ns),vb*2./3.,numion,
     $    ifre1,imult,ifre2b,ni,nom(lmin_wr),mult,ns,nom(l+1),
     $    com_Li(ns,l+1,i_Li)
              ifre2=ifre2+1
        wlp=1./(1./wl(ni,l,ns)+1.d-8*del_Li(ns,l+1,i_Li))
c        write(out_niv,1187) del_Li(ns,l+1,i_Li),l+1
c 1187   Format(' del_Li(ns,l+1,i_Li)=',f7.3,' l+1=',i2)
        Write(out_niv,187) numion,ifre1,imult,ifre2,ni_alph,ns_alph,
     $    lp1_alph,nomion,wlp,vb/3.,numion,
     $    ifre1,imult,ifre2b,ni,nom(lmin_wr),mult,ns,nom(l+1),
     $    com_Li(ns,l+1,i_Li)
              ifre2=ifre2b
c
c D) Li_like 2 :
                Elseif  ! split Li_like np trans p-d + trans p-s (k=2)
     $ ((numion.eq.' 604'.or.numion.eq.' 705'.or.numion.eq.' 806').and.
     $ l.eq.1.and.ns.le.9) Then
c
            i_Li=izz-5  ! 1 CIV, 2 NV, 3 OVI
c split trans p-d : 
                iraies=iraies+1  ! (1 raie --> 2)
              ifre2b=ifre2
              ifre2=ifre2+1
        Write(out_niv,187) numion,ifre1,imult,ifre2,ni_alph,ns_alph,
     $    lp1_alph,nomion,wl(ni,l,ns),vb*0.6,numion,
     $    ifre1,imult,ifre2b,ni,nom(lmin_wr),mult,ns,nom(l+1)
              ifre2=ifre2+1
         wlp=1/(1/wl(ni,l,ns)+
     $                1.d-8*(del_Li(ns,l+1,i_Li)-del_Li(ni,l,i_Li)))
c        write(out_niv,1188) del_Li(ns,l+1,i_Li),del_Li(ni,l,i_Li),l
c 1188 Format(' * del_Li(ns,l+1,i_Li),del_Li(ni,l,i_Li)=',2f7.3,' l=',i2)
        Write(out_niv,187) numion,ifre1,imult,ifre2,ni_alph,ns_alph,
     $    lp1_alph,nomion,wlp,vb*0.4,numion,
     $    ifre1,imult,ifre2b,ni,nom(lmin_wr),mult,ns,nom(l+1)
c trans p-s and split : accessible only if the p-d int is already suff for write
           k=2
           vc=Int_Hlike(iz,Mat,ni,l,ns,k,Te)*comult(ni,l)/int_ref
           k=1
c corr departure wl H-like :
           vc=vc*wlhydrostr(ni,ns)/wk(ni,l,ns) !(with corr air refr)
c lim int p-s :
              iecrb=iecr
c *Rem: provisional modif of iecr within the 'IF' of the 1st iecr=1:
c
          iecr=0
         iw=1
         if(wk(ni,l,ns).le.wlim(iw).and.vc.ge.Irelmin(iw)) then
          iecr=1
         endif
        do iw=2,iwlim
         if(wk(ni,l,ns).gt.wlim(iw-1).and.wk(ni,l,ns).le.wlim(iw)
     $     .and.vc.ge.Irelmin(iw)) then
          iecr=1
         endif
        enddo
        if(wk(ni,l,ns).gt.wlim(iwlim).and.vc.ge.Irelmin(iwlim+1)) then
          iecr=1
        endif
c
           if(iecr.eq.1) then !iecr=1
               iraies=iraies+2  ! (2 nouvelles raies)
              ifre2=ifre2+1
        Write(out_niv,187) numion,ifre1,imult,ifre2,ni_alph,ns_alph,
     $    lp1_alph,nomion,wk(ni,l,ns),vc*2./3.,numion,
     $    ifre1,imult,ifre2b,ni,nom(lmin_wr),mult,ns,nom(21) !nom(21) pour 's'
              ifre2=ifre2+1
             wlp=1./(1./wk(ni,l,ns)-1.d-8*del_Li(ni,l,i_Li))
c        write(out_niv,1189) del_Li(ni,l,i_Li),l
c 1189   Format(' ** del_Li(ni,l,i_Li)=',f7.3,' l=',i2)
        Write(out_niv,187) numion,ifre1,imult,ifre2,ni_alph,ns_alph,
     $    lp1_alph,nomion,wlp,vc/3.,numion,
     $    ifre1,imult,ifre2b,ni,nom(lmin_wr),mult,ns,nom(21)
           endif
              ifre2=ifre2b
              iecr=iecrb  !end provisional modif
c
c E) Li_like 3 :
                Elseif  ! split Li_like nd trans d-f & d-p (k=2)
     $ ((numion.eq.' 604'.or.numion.eq.' 705'.or.numion.eq.' 806').and.
     $ (l.eq.2.and.ns.le.9)) Then
c
            i_Li=izz-5  ! 1 CIV, 2 NV, 3 OVI
c split trans d-f : 
                iraies=iraies+1  ! (1 raie --> 2)
              ifre2b=ifre2
              ifre2=ifre2+1
        Write(out_niv,187) numion,ifre1,imult,ifre2,ni_alph,ns_alph,
     $    lp1_alph,nomion,wl(ni,l,ns),vb*0.6,numion,
     $    ifre1,imult,ifre2b,ni,nom(lmin_wr),mult,ns,nom(l+1)
              ifre2=ifre2+1
             wlp=1./(1./wl(ni,l,ns)-1.d-8*del_Li(ni,l,i_Li))
c        write(out_niv,1190) del_Li(ni,l,i_Li),l
c 1190   Format(' *** del_Li(ni,l,i_Li)=',f7.3,' l=',i2)
        Write(out_niv,187) numion,ifre1,imult,ifre2,ni_alph,ns_alph,
     $    lp1_alph,nomion,wlp,vb*0.4,numion,
     $    ifre1,imult,ifre2b,ni,nom(lmin_wr),mult,ns,nom(l+1)
c trans d-p and split: accessible only if d-f int already suff for write
           k=2
           vc=Int_Hlike(iz,Mat,ni,l,ns,k,Te)*comult(ni,l)/int_ref
           k=1
c corr ecart wl H-like :
           vc=vc*wlhydrostr(ni,ns)/wk(ni,l,ns) !(with corr air refr)
c lim int d-p :
              iecrb=iecr
          iecr=0
c
         iw=1
         if(wk(ni,l,ns).le.wlim(iw).and.vc.ge.Irelmin(iw)) then
          iecr=1
         endif
        do iw=2,iwlim
         if(wk(ni,l,ns).gt.wlim(iw-1).and.wk(ni,l,ns).le.wlim(iw)
     $     .and.vc.ge.Irelmin(iw)) then
          iecr=1
         endif
        enddo
        if(wk(ni,l,ns).gt.wlim(iwlim).and.vc.ge.Irelmin(iwlim+1)) then
          iecr=1
        endif
c
           if(iecr.eq.1) then !iecr=1
               iraies=iraies+2  ! (2 new lines)
              ifre2=ifre2+1
        Write(out_niv,187) numion,ifre1,imult,ifre2,ni_alph,ns_alph,
     $    lp1_alph,nomion,wk(ni,l,ns),vc*2./3.,numion,
     $    ifre1,imult,ifre2b,ni,nom(lmin_wr),mult,ns,nom(l-1),
     $    com_Li(ns,l-1,i_Li)
              ifre2=ifre2+1
        wlp=1/(1/wk(ni,l,ns)+
     $               1.d-8*(del_Li(ns,l-1,i_Li)-del_Li(ni,l,i_Li)))
c        write(out_niv,1191) del_Li(ns,l-1,i_Li),del_Li(ni,l,i_Li),l
c 1191    Format(' **** del_Li(ns,l-1,i_Li),del_Li(ni,l,i_Li)=',2f7.3,
c     $          ' l=',i2)
        Write(out_niv,187) numion,ifre1,imult,ifre2,ni_alph,ns_alph,
     $    lp1_alph,nomion,wlp,vc/3.,numion,
     $    ifre1,imult,ifre2b,ni,nom(lmin_wr),mult,ns,nom(l-1),
     $    com_Li(ns,l-1,i_Li)
           endif
              ifre2=ifre2b
c
              iecr=iecrb
c
c F) Notation J H_like:
                Elseif  ! end split Li_like
     $     (numion.eq.' 101'.or.numion.eq.' 202') then
c  HI, DI, HeII et 3HeII special impressions: notation J
c
        Write(out_niv,487) numion,ifre1,imult,ifre2,ni_alph,ns_alph,
     $    lp1_alph,nomion,wl(ni,l,ns),vb,numion,ifre1,imult,ifre2,
     $    ni,2*l+1,ns
 487    Format(1x,A4,3I1,3A2,1x,A7,f13.3,' 0.   ',1p,e10.3,'  1.     ',
     $    A4,3I1,'000000   1   1.00 ',I2,'(J=',I1,'/2)-',I2)
c
c G)  General Case:
                Else  ! END SPECIAL CASES: back to general case write
c
        Write(out_niv,187) numion,ifre1,imult,ifre2,ni_alph,ns_alph,
     $    lp1_alph,nomion,wl(ni,l,ns),vb,numion,ifre1,imult,ifre2,
     $    ni,nom(lmin_wr),mult,ns,nom(l+1)
 187    Format(1x,A4,3I1,3A2,1x,A7,f13.3,' 0.   ',1p,e10.3,'  1.     ',
     $    A4,3I1,'000000   1   1.00 ',I2,A1,I1,'-',I2,A1,1x,A8)
c
                Endif  ! END ALTERNATIVE SPECIAL CASE / GENERAL CASE
c
c*****##
             ENDIF              ! end 1st iecr=1
c*****##
c***###***
           Else !indregr.ne.0; end indregr=0
c               --> Grouping: in that case, lmin=lmax in DO loop above
c***###***
              IF(hydrostrict.ne.2) THEN !hydrostrict eq 0 or 1
            vb=0.
            do lregr=lmax_int,lmin_int,-1
            vb=vb+Int_Hlike(iz,Mat,ni,lregr,ns,k,Te)*comult(ni,lregr)
            enddo
ccc  hydrostrict.eq.1 : Case H-like strict: all transitions, including the k=2
c (Rem : hydrostrict.eq.2 and indregr.ne.0 : lecture et interp int tot dirctmt dans PSto)
cc  The sum l=ni-1 to 0 is not imposed by the logics, but is desirable in general:
            If(hydrostrict.eq.1) Then
            k=2
             do lregr=lmax_int,lmin_int,-1
              vb=vb+Int_Hlike(iz,Mat,ni,lregr,ns,k,Te)*comult(ni,lregr)
             enddo
            k=1
            Endif
ccc 
           vb=vb/int_ref
ccc
              ELSE    !hydrostrict.eq.2
            vb=-1.d0   !inhibit indregr, cf Bracket hydrostrict=2
              ENDIF   !hydrostrict end
c
c  (***1/15 ATT origin ni ns l to be clarified here)
          iecr=0
c
         iw=1
         if(wl(ni,l,ns).le.wlim(iw).and.vb.ge.Irelmin(iw)) then
          iecr=1
         endif
        do iw=2,iwlim
         if(wl(ni,l,ns).gt.wlim(iw-1).and.wl(ni,l,ns).le.wlim(iw)
     $     .and.vb.ge.Irelmin(iw))then
          iecr=1
         endif
        enddo
        if(wl(ni,l,ns).gt.wlim(iwlim).and.vb.ge.Irelmin(iwlim+1)) then
          iecr=1
        endif
c
             IF(iecr.eq.1) THEN ! 2nd iecr=1 --> line kept
               iraies=iraies+1
               ini(ni)=1
c
           vb=vb*wlhydrostr(ni,ns)/wl(ni,l,ns) !ATT: with air refraction
c
          Write(out_niv,178) numion,ifre1,imult,ifre2,ni_alph,ns_alph,
     $    lp1_alph,nomion,wl(ni,l,ns),vb,numion,ifre1,imult,ifre2,
     $     ni,nom(lmax_int),nom(lmin_wr),mult,ns
 178      Format(1x,A4,3I1,3A2,1x,A7,f13.3,' 0.   ',1p,e10.3,'  1.     ',
     $      A4,3I1,'000000   1   1.00 ',I2,A1,'_',A1,I1,'-',I2)
             ENDIF             ! end 2e iecr=1
c***###***
           Endif   !end else indregr.ne.0 (grouping)
c***###***
c**********
          ENDIF   !wl range
c**********
c********************
         ENDIF    !imult, nmult, comult
c********************
c***************
        ENDDO     !l=lmax,lmin,-1   loop internal
c***************
c Storage for extrapolation in case nma_orig > 50 :
            if(iextrapol.eq.1) Emisto(ns)=vb
c*************************
       Enddo      !ns = upper level, here ns.le.50  loop intermed 1
c*************************
c
c******  (Bracket 100 > nma_orig > 50 beginning)
c Extrapolation ns > 50 :
         l=lmax  ! (cf 'average' wl)
c***######
       If(iextrapol.eq.1.and.         !begin iextrapol=1
     $     wl(ni,l,50).gt.wlmin.and.wl(ni,l,50).lt.wlmax) then
        pente1=dlog10(Emisto(50)/Emisto(40))/dlog10(50.d0/40.d0)
        print *,ni,' extrapolation ns > 50 : pente1 =',pente1
         do ns=51,nma_orig !storage extrap emissiv
          Emisto(ns)=Emisto(50)*(dble(ns)/50.)**pente1
         enddo
c**********
         Do ns=51,nma_orig  ! 100 > ns > 50   loop intermed 2
c**********
          IF(wl(ni,l,ns).gt.wlmin.and.wl(ni,l,ns).lt.wlmax) THEN !wl range
          vb=Emisto(ns)
c
          iecr=0
c
         iw=1
         if(wl(ni,l,ns).le.wlim(iw).and.vb.ge.Irelmin(iw)) then
          iecr=1
         endif
        do iw=2,iwlim
         if(wl(ni,l,ns).gt.wlim(iw-1).and.wl(ni,l,ns).le.wlim(iw)
     $     .and.vb.ge.Irelmin(iw))then
          iecr=1
         endif
        enddo
        if(wl(ni,l,ns).gt.wlim(iwlim).and.vb.ge.Irelmin(iwlim+1)) then
          iecr=1
        endif
c
             IF(iecr.eq.1) THEN !iecr=1
               iraies=iraies+1
               ini(ni)=1
        Write(out_niv,187) numion,ifre1,imult,ifre2,ni_alph,ns_alph,
     $    lp1_alph,nomion,wl(ni,l,ns),vb,numion,ifre1,imult,ifre2,
     $    ni,nom(lmin_wr),mult,ns,nom(l+1)
             ENDIF              !end iecr=1
c**********
          ENDIF    !wl range
c**********
         Enddo     !ns > 50    loop intermed 2
c**********
       Endif       !end iextrapol=1
c***######
c******  (Bracket 100 > nma_orig > 50 end)
c
c******  (Bracket hydrostrict=2 beginning)
c  H, He+ : PSto (with interpol) if the l's are grouped
c  (transitions inhibited by vb=-1 above, hence no duplication)
       If(hydrostrict.eq.2) Then
         if(ni.lt.niregr) then
           nh2=nsregr
         else
           nh2=ni+1
         endif
c**********
         Do ns=nh2,nma2 !nma2 < 500    loop intermed 3
c********** !!!!!***** ICIC
          va=wlhydrostr(ni,ns) !**1/15: PB: USE EXACT AVERAGE WL, NOt HYDROSTRICT WL
          IF(va.gt.wlmin.and.va.lt.wlmax) THEN !wl range
          vb=emhydr(ni,ns)
c
          iecr=0
c
         iw=1
         if(va.le.wlim(iw).and.vb.ge.Irelmin(iw)) then
          iecr=1
         endif
        do iw=2,iwlim
         if(va.gt.wlim(iw-1).and.va.le.wlim(iw)
     $     .and.vb.ge.Irelmin(iw))then
          iecr=1
         endif
        enddo
        if(va.gt.wlim(iwlim).and.vb.ge.Irelmin(iwlim+1)) then
          iecr=1
        endif
c
             IF(iecr.eq.1) THEN !iecr=1
               iraies=iraies+1
               ini(ni)=1
c   cas part H et He+ ns > 99 possible (PSto) :
        if(ns.lt.10) then
       write(nc1,'(I1)') ns
        ns_alph3=zero//zero//nc1
         elseif(ns.lt.100) then
        write(nc2,'(I2)') ns
       ns_alph3=zero//nc2
          else   !(assume ns < 1000)
        write(ns_alph3,'(I3)') ns
        endif
c  
        if(l+1.lt.10) then
       write(nc1,'(I1)') l+1
       lp1_alph=zero//nc1
        else   !(assume l+1 < 100)
       write(lp1_alph,'(I2)') l+1
        endif
c
                  if(ns.lt.100) then ! ns<100,  ifre2 suppressed here
          Write(out_niv,179) numion,ifre1,imult,ni_alph,ns_alph3,
     $    lp1_alph,nomion,va,vb,numion,ifre1,imult,ifre2,
     $     ni,nom(lmax_int),nom(lmin_wr),mult,ns
 179      Format(1x,A4,2I1,A2,A3,A2,1x,A7,f13.3,' 0.   ',1p,e10.3,
     $      '  1.     ',
     $      A4,3I1,'000000   1   1.00 ',I2,A1,'_',A1,I1,'-',I2)
                  else  ! 99<ns<1000,           ifre2 suppressed here
          Write(out_niv,180) numion,ifre1,imult,ni_alph,ns_alph3,
     $    lp1_alph,nomion,va,vb,numion,ifre1,imult,ifre2,
     $     ni,nom(lmax_int),nom(lmin_wr),mult,ns
 180      Format(1x,A4,2I1,A2,A3,A2,1x,A7,f13.3,' 0.   ',1p,e10.3,
     $      '  1.     ',
     $      A4,3I1,'000000   1   1.00 ',I2,A1,'_',A1,I1,'-',I3)
                  endif
             ENDIF !iecr
c**********
          ENDIF    !wl range
c**********
         Enddo     !ns    loop intermed 3
c**********
       Endif       !hydrostrict=2
c******  (Bracket hydrostrict=2 end)
c
c#####
c********  ********** ************
           Enddo   !ni lower level, ****** loop + external ni
c********  ********** ************
c#####
c
c  BEGINNING SPECIAL CASES : non H-like optical rad trans (eg, ns=ni)
c
c A) ADD NeIII 4s-4p, 4p-4d quint or tripl, with alim 4p, 4d admited H-like
c
         IF(numion.eq.'1003') THEN ! 3 EXPLICIT
           ni=4
           ns=4
           ni_alph='04'
           ns_alph='04'
c
          ifre2b=ifre2
           DO lp1=1,2
              write(nc1,'(I1)') lp1
              lp1_alph=zero//nc1
            l=lp1-1
             lmin_wr=l
             if(lmin_wr.le.0) then
c  nom(0) does not exist
             lmin_wr=21
             endif
c imult=0, 1 for quintet, triplet, but the tabulation must start from 1
           im=3 ! EXPLICIT pratical component number
            if(imult.eq.1.and.l.eq.1) im=5 !only 4p3-4d3
           iraies=iraies+im
c
          if(imult.eq.0) then !5
           if(lp1.eq.1) then
             BrNeIII=BrNeIII4p(imult+1)
             comNeIII=comNeIII4p(imult+1)
            do i=1,im
             wlNeIII(i)=wlNeIII4p5(i)
             coNeIII(i)=coNeIII4p5(i)
            enddo             
           else  !(lp1=2)
             BrNeIII=BrNeIII4d(imult+1)
             comNeIII=comNeIII4d(imult+1)
            do i=1,im
             wlNeIII(i)=wlNeIII4d5(i)
             coNeIII(i)=coNeIII4d5(i)
            enddo             
           endif
          elseif(imult.eq.1) then !3
           if(lp1.eq.1) then
             BrNeIII=BrNeIII4p(imult+1)
             comNeIII=comNeIII4p(imult+1)
            do i=1,im
             wlNeIII(i)=wlNeIII4p3(i)
             coNeIII(i)=coNeIII4p3(i)
            enddo             
           else  !(lp1=2)
             BrNeIII=BrNeIII4d(imult+1)
             comNeIII=comNeIII4d(imult+1)
            do i=1,im
             wlNeIII(i)=wlNeIII4d3(i)
             coNeIII(i)=coNeIII4d3(i)
            enddo             
           endif
          else
            print *,'***STOP NeIII anomaly : imult =', imult
            STOP
          endif
c
             do i=1,im
c alfeff(nu,lup), lup=lu+1 idem lp1+1
          vb=alfeff(ns,lp1+1)*BrNeIII*hc8/wlNeIII(i)/int_ref
c **Rem: no intensity test
          ifre2=ifre2+1
c (comult not necessarily defined here)
          Write(out_niv,1781) numion,ifre1,imult,ifre2,ni_alph,ns_alph,
     $    lp1_alph,nomion,wlNeIII(i),vb*coNeIII(i)*cmult,numion,
     $    ifre1,imult,ifre2b,ni,nom(lmin_wr),mult,ns,nom(l+1),comNeIII
 1781   Format(1x,A4,3I1,3A2,1x,A7,f13.3,' 0.   ',1p,e10.3,'  1.     ',
     $        A4,3I1,'000000   1   1.00 ',I2,A1,I1,'-',I2,A1,1x,A10)
             enddo
          ifre2=ifre2b
           ENDDO  !lp1=1,2
         ENDIF  ! end add NeIII
c
c B) begin ADD Li_like 2s-2p, 3s-3p, 3p-3d
      IF(numion.eq.' 604'.or.numion.eq.' 705'.or.numion.eq.' 806') THEN !EXPLICIT
c
          i_Li=izz-5  ! 1 CIV, 2 NV, 3 OVI
        ifre2b=ifre2
c
      co_Li(1)=2./3. !(3/2 -> 3/2 (I/15) mise dans 5/2 -> 3/2 pour 3p-3d)
      co_Li(2)=1./3.
           DO ni=2,3
          DO k=2,ni
       l=k-2
       ns=ni
        write(nc1,'(I1)') ni
        ni_alph=zero//nc1
        ns_alph=ni_alph
        lmin_wr=l
        if(lmin_wr.le.0) lmin_wr=21
        write(nc1,'(I1)') l+1
        lp1_alph=zero//nc1
            va=E(ns,l+1)-E(ni,l)
          wl_Li(1)=1.d8/va
            va=va+del_Li(ns,l+1,i_Li)
            if(ns-1+l.eq.3) va=va-del_Li(ni,l,i_Li)  ! 3p-3d
          wl_Li(2)=1.d8/va
           IF(wl_Li(1).gt.wlmin.and.wl_Li(1).lt.wlmax) THEN !wl range
         If(ns-1+l.eq.1) then      ! 2s-2p (UV)
c case B
       vb=alfeff(2,1)+alfeff(2,2)
         Elseif(ns-1+l.eq.2) then  ! 3s-3p
c case B with cascade correction from 5p,6p,7p 
c    Br(nu,lup,nl,k), lup=lu+1
       vb=alfeff(3,2)+alfeff(5,2)*Br(5,2,2,1)*
     $  (Br(5,2,4,1)*Br(4,1,3,2)+Br(5,2,4,2)*Br(4,3,3,1))+
     $  alfeff(6,2)*Br(6,2,2,1)*
     $  (Br(6,2,4,1)*Br(4,1,3,2)+Br(6,2,4,2)*Br(4,3,3,1)+
     $  Br(6,2,5,1)*Br(5,1,3,2)+Br(6,2,5,2)*Br(5,3,3,1))+
     $  alfeff(7,2)*Br(7,2,2,1)*
     $  (Br(7,2,4,1)*Br(4,1,3,2)+Br(7,2,4,2)*Br(4,3,3,1)+
     $  Br(7,2,5,1)*Br(5,1,3,2)+Br(7,2,5,2)*Br(5,3,3,1)+
     $  Br(7,2,6,1)*Br(6,1,3,2)+Br(7,2,6,2)*Br(6,3,3,1)+
     $  (Br(7,2,6,1)*Br(6,1,5,2)+Br(7,2,6,2)*Br(6,3,5,1))*
     $  (Br(5,2,4,1)*Br(4,1,3,2)+Br(5,2,4,2)*Br(4,3,3,1)))
          wl_Li(1)=wl_Li(1)/A(1.d4/wl_Li(1))
          wl_Li(2)=wl_Li(2)/A(1.d4/wl_Li(2))
         Else               ! 3p-3d  (3/2 -> 3/2 (I/15) mise dans 5/2 -> 3/2)
c case B with cascade correc from 4p,5p,6p, 
c   but with 2p-3d trans remaining optically thin (branching maintained)
c (case C if 2p-3d optically thick: improbable)
c    Br(nu,lup,nl,k), lup=lu+1
       vb=(alfeff(3,3)+alfeff(4,2)*Br(4,2,2,1)*Br(4,2,3,2)+
     $  alfeff(5,2)*Br(5,2,2,1)*Br(5,2,3,2)+
     $  alfeff(6,2)*Br(6,2,2,1)*(Br(6,2,3,2)+(Br(6,2,5,1)*Br(5,1,4,2)+
     $  Br(6,2,5,2)*Br(5,3,4,1))*Br(4,2,2,1)))/IBr_Li(ns-1+l,i_Li)
c   air refraction wl > 2000A and wl < 20000A:
           IF(wl_Li(1).gt.2.d3.and.wl_Li(1).lt.2.d4)THEN
          wl_Li(1)=wl_Li(1)/A(1.d4/wl_Li(1))
          wl_Li(2)=wl_Li(2)/A(1.d4/wl_Li(2))
           ENDIF
         Endif
       BsA=vb/(alfeff(ns,l+1)/IBr_Li(ns-1+l,i_Li))     ! case B / case A
       vb=vb*hc8/(wl_Li(1)*co_Li(1)+wl_Li(2)*co_Li(2))/int_ref ! case B
c
            iecr=0
c
          do iw=2,iwlim
           if(wl_Li(1).gt.wlim(iw-1).and.wl_Li(1).le.wlim(iw)
     $       .and.vb.ge.Irelmin(iw))then
            iecr=1
           endif
          enddo
           if(wl_Li(1).gt.wlim(iwlim).and.vb.ge.Irelmin(iwlim+1))then
            iecr=1
           endif
c
             IF(iecr.eq.1) THEN !iecr=1
       iraies=iraies+2
        ifre2=0
        do i=1,2
        ifre2=ifre2+1
        Write(out_niv,877) numion,ifre1,imult,ifre2,ni_alph,ns_alph,
     $   lp1_alph,nomion,wl_Li(i),vb*co_Li(i),numion,ifre1,imult,ifre2b,
     $   ni,nom(lmin_wr),mult,ns,nom(l+1),BsA
        enddo
 877    Format(1x,A4,3I1,3A2,1x,A7,f13.3,' 0.   ',1p,e10.3,'  1.     ',
     $  A4,3I1,'000000   1   1.00 ',I2,A1,I1,'-',I2,A1,' A=B/',e8.2)
             ENDIF !iecr
           ENDIF !wl range
          ENDDO     ! ni
           ENDDO    ! k (= l+2)
      ENDIF         ! numion.eq  END Li-like
c
c  END SPECIAL CASES
c
         Write(out_niv,75) iraies
 75      format('Number of X-SSN lines: iraies =',i5)
c
      CLOSE(out_niv)
c
        inmitab=0
       do ni=nmitab,nmatab
         if(ini(ni).ne.0.and.inmitab.eq.0) then
        inmitab=ni
         endif
       enddo
        inmatab=0
       do ni=nmatab,nmitab,-1
         if(ini(ni).ne.0.and.inmatab.eq.0) then
        inmatab=ni
         endif
       enddo
c Conversion to character (note: iz > 8 treated below)
      write(inmitabc,'(I2)') inmitab
      write(inmatabc,'(I2)') inmatab
      write(iraiesc,'(I4)') iraies
      write(Irelminc,'(I5)') int(Irelmin(4)*1.d5)
      write(wlminc,'(I8)') int(wlmin)
      write(wlmaxc,'(I8)') int(wlmax)
          IF(imult.eq.0) THEN
c Print correct if imult=0 is 1st in main programme for ions nmult>0
      write(lambda_refc,'(I5)') int(lambda_ref)
          ENDIF
c
      write(ni_refc,'(I2)') ,ni_ref
      write(ns_refc,'(I2)') ,ns_ref
      print *,'effective ni range:',inmitabc,inmatabc,
     $ '  wl range:',wlminc,wlmaxc
        if(li_ref.eq.0)li_ref=21
      print *,'number of emission lines = ',iraiesc,' / ',
     $ '   Irelmin =',Irelminc,'E-5',lambda_refc,' (',
     $ ni_refc,nom(li_ref),nom(ni_ref-1),'  ',ns_refc,')'
        if(li_ref.eq.21)li_ref=0
c
c  If requested, insert the result .res in PHYAT for SYNTHESIS:
      if(iwr_phyat.ne.0) then
         Call Write_Phyat(lambda_output,PHYAT)
      endif
c
      RETURN
c
      END
c
c-----------------end XSSN_Hlike-------------------------
c------------------------------------------------
c
      Double Precision Function Int_Hlike(iz,Mat,nl,ll,nu,k,Te)
c
c      Units: erg.cm3.s-1
c
c Emissivity of line (nu,lu,nl,ll) [Notation: k=1: lu=ll+1, k=2: lu=ll-1] 
c  in the quasi hydrogenic approximation for charge iz, atomic mass Mat, 
c  based on H-like effective recombination coefficients alfeff from 
c  tables of Storey & Hummer (1995, SH95), read elsewhere by:
c                 Subroutine Readsto, 
c  for given conditions (Ne, Te), and on H-like branching ratios Br_nl,n'l'
c   computed by:
c                 Subroutine HBD,
c  adapted from program BA5.2 by Hoang Binh, Dy: 1990, A&A, 238, 449.
c
c     CONCERNING the use of PARAMETER Te, SEE BELOW
c
c IMPORTANT: the alphaeff of SH95 include all cascades and n-changing collisions, 
c     but NOT the l-changing collisions for given n. This is adequate for 'quasi' 
c     (but not exactly) H-like ions, in which small energy shifts will decrease 
c     the efficiency of nl --> nl' collisions relative to the pure hydrogenic case.
c
c   For PURE H-like (normally HI, HeII) line emissivities, use the departure 
c   coefficients from LTE, b_SH(Ne,Te) tabulated by SH95, Saha and the A_tp_nl,n'l' 
c   (unlike the alfaeff, the b_SH account for the nl --> nl' collisions). 
c
c     The two cases are differenciated locally by means of Te: 
c                 Te </= 0. for any quasi H-like.
c                 Te = Te for pure H-like (Te used in Saha equation)
c
      IMPLICIT NONE
c
      SAVE
c
      INTEGER iz
      INTEGER nu,nl,ll,lup,k
      DOUBLE PRECISION Te
      DOUBLE PRECISION alfeff,b_SH
      DOUBLE PRECISION csaha,cons,kboltz
      DOUBLE PRECISION Mat
      DOUBLE PRECISION Br,A_tp
      DOUBLE PRECISION chi_Hli,Eline
      DOUBLE PRECISION hc,hc8,hcsk,Ryd_inf
c
      COMMON/MOCC1/hc,hc8,hcsk,Ryd_inf
      COMMON/MOCC2/csaha,kboltz
      COMMON/MOCH1/Br(100,100,100,2),A_tp(100,100,100,2) !Br(nu,lup,nl,k) lup=lu+1
      COMMON/MOCS1/alfeff(50,50),b_SH(50,0:50) !alfeff(nu,lup) lup=lu+1
c
      chi_Hli=dble(iz*iz)*hc*Ryd_inf/(1.d0+1.d0/(Mat*1823.d0)) !erg
c
c alfeff(nu,lup), Br(nu,lup,nl,k), b_SH(nu,lup), A_tp(nu,lup,nl,k), lup=lu+1
c
      if(k.eq.1) then
       lup=ll+2 !trans lu=ll+1 --> ll
      elseif(k.eq.2) then
       lup=ll   !trans lu=ll-1 --> ll
      else
       print *,'*** in Int_Hlike, k=',k,' STOP***'
       STOP
      endif
c
      If(lup.eq.0) then
c impossible transition:
         Int_Hlike=0.d0
      Else
      Eline=chi_Hli*(1.d0/dble(nl*nl)-1.d0/dble(nu*nu))
c
        if(Mat.ge.5.d0.and.(2.*dble(iz)+1.).lt.Mat) then
c Quasi H-like (exclude pure H-like CVI, etc.): 
      Int_Hlike=alfeff(nu,lup)*Br(nu,lup,nl,k)*Eline !erg.cm3.s-1
        else
c Pure H-like: (Note: if Mat = 3 or 4: Int_Hlike not called for He I)
      Int_Hlike=csaha*2*(2*dble(lup-1)+1.)/(2.*1)/Te**1.5*b_SH(nu,lup)*
     $ A_tp(nu,lup,nl,k)*Eline*dexp(chi_Hli/nu/nu/kboltz/Te) !erg.cm3.s-1 (*1 for proton)
        endif
      Endif
c
      Return
c
      End
c----------------end Int_Hlike---------------------------
c--------------------------------------------------------
c Ex:  Call Write_Phyat('OVI_2.res','liste_phyat_NF.dat')
c Condition: compatibility of 6 1st digits labeling the "ion" to be replaced
c  CONSEQUENCES: - ifre1 and imult contribute to caracterize the ion.
c                - only ifre2 (7th digit) can be used to label the lines 
c                   (in the same way as the 6 last digits).
c
      Subroutine Write_Phyat(lambda_output,PHYAT)
c
      IMPLICIT NONE
c
      CHARACTER lambda_output*(*),PHYAT*(*)
c
      character*130 ligne,ligne_ref_modele
      integer n,n_ini,n_suiv,ind_new
      integer io
      integer imult,inil,ntot,iraies
      integer label,label_ion
      INTEGER out_niv
      integer in_phy,out_phy,old_phy,ntot_phy,ntot_phynew
c
      out_niv=11 ! .res results (lambda_output), including for X-SSN (11 arb)
      in_phy=14 ! initial PHYAT of X-SSN before replacement (14 arb)
      old_phy=15 ! backup PHYAT before replacement (15 arb)
      out_phy=16 ! intermediate PHYAT after replacement (16 arb), erase initial PHYAT
c
 101  format(A130)
c
      OPEN(out_niv,file=lambda_output,status='old')
cc ex :     OPEN(out_niv,file='OIV_di.res',status='old')
c
      ntot=0
      Do n=1,3
        read(out_niv,*) ligne
        ntot=ntot+1
      Enddo
       read(out_niv,193) imult
 193   format(32X,I4)
        ntot=ntot+1
      Do
        read(out_niv,*,end=1000) ligne
        ntot=ntot+1
      Enddo
 1000 continue
c  Total number of lines in .res : ntot
c
      BACKSPACE(out_niv)  !    back twice
      BACKSPACE(out_niv)  ! ... to get last line !?
      read(out_niv,195) iraies !get emission line number
 195   format(31X,I5)
c
      REWIND(out_niv)
c
      inil=ntot-iraies-1 !inil 1st emission line for PHYAT (liste_phyat..)
       if(imult.eq.0) then
        Do n=1,inil-3
         read(out_niv,*) ligne
        Enddo
        read(out_niv,101) ligne_ref_modele !for liste_modele
        read(out_niv,*) ligne  !title
       else
        Do n=1,inil-1
         read(out_niv,*) ligne
        Enddo
       endif
c 1st line: read ion label (i6)
c  made of numion, ifre1 and imult for liste_phyat :
       read(out_niv,190) label_ion
       print *,' label_ion =',label_ion
 190   format(1X,I6)
c
      OPEN(unit=in_phy,file=PHYAT,status='old')
c
c Spotting ion in PHYAT:
      n_ini=0
      n_suiv=0
      ind_new=0
      n=0
      DO
         n=n+1
       Read(in_phy,190,IOSTAT=io) label
      IF (io.gt.0) THEN
       print *, '***Pb reading PHYAT, line =',n,' STOP***'
        STOP
      ELSEIF (io.eq.0) then
       If(ind_new.eq.0) THEN
        if(label.eq.label_ion) then
         ind_new=1
c Beginning ion :
         n_ini=n
        endif
       Elseif(ind_new.eq.1) THEN
        if(label.ne.label_ion) then
         ind_new=2
c Beginning next ion:
         n_suiv=n
        endif
       Endif
      ELSE !end file (io.lt.0)
c  Total nber of lines in PHYAT before modifying : ntot_phy
       ntot_phy=n-1
       if(ind_new.eq.0) then !anomaly:
        print *, '***end phyat : label_ion =',label_ion,' not found ***'
       elseif(ind_new.eq.1) then !end last ion
        n_suiv=n
       endif
        GOTO 99
      ENDIF
      ENDDO
 99   CONTINUE
        print *,' ntot_phy, n_ini, n_suiv =',ntot_phy,n_ini,n_suiv
c
      REWIND(in_phy)
c
c  Safeguard previous PHYAT: .dat --> .old
c   suppr previous .old if any:
      OPEN(unit=old_phy,file='phyat.old',
     $     iostat=io,status='unknown')
      if(io.eq.0) then !really opened, thus can be deleted:
       CLOSE(old_phy,status='delete')
      endif
c
      OPEN(unit=old_phy,file='phyat.old',status='new')
c
       Do n=1,ntot_phy
        read(in_phy,101) ligne
        write(old_phy,101) ligne
       Enddo
c
      CLOSE(old_phy)
c
c  Positioning on .res :
      REWIND(out_niv)
       Do n=1,inil-1
        read(out_niv,*) ligne
       Enddo
c
c Replacement in .new:
c  suppr possible previous .new:
      OPEN(unit=out_phy,file='phyat.new',
     $     iostat=io,status='unknown')
      if(io.eq.0) then !really opened, thus can be deleted:
       CLOSE(out_phy,status='delete')
      endif
c
      OPEN(out_phy,file='phyat.new',status='new')
c
      REWIND(in_phy)
c
      ntot_phynew=0
c Before ion :
       Do n=1,n_ini-1
        read(in_phy,101) ligne
        write(out_phy,101) ligne
      ntot_phynew=ntot_phynew+1
       Enddo
c Modified ion:
       Do n=inil,ntot-1
        read(out_niv,101) ligne
        write(out_phy,101) ligne
      ntot_phynew=ntot_phynew+1
       Enddo
c Jump (remove) ion before modif:
       Do n=n_ini,n_suiv-1
        read(in_phy,*) ligne
       Enddo
c After ion:
       Do n=n_suiv,ntot_phy
        read(in_phy,101) ligne
        write(out_phy,101) ligne
      ntot_phynew=ntot_phynew+1
       Enddo
        print *,' new total ntot_phynew =',ntot_phynew
c
      CLOSE(out_niv) !closing .res
c
      REWIND(out_phy)
c
      CLOSE(in_phy,status='delete')
      OPEN(in_phy,file=PHYAT,status='new')
c
c  Replacement PHYAT : .new --> .dat :
c 
       Do n=1,ntot_phynew
        read(out_phy,101) ligne
        write(in_phy,101) ligne
       Enddo
c
      CLOSE(in_phy) 
c
      CLOSE(out_phy,status='delete') !suppr .new
c
      RETURN
c
      END
c---------------end Write_Phyat -------------------------
c--------------------------------------------------------
c Ex: Call XSSN_di(...,'CIII_di.lab','CIII_di.res',PHYAT)
c
      Subroutine XSSN_di(Te,Ne,niv_input,lambda_output,PHYAT)
c
c     Di-electronic transitions.
c
c  THE REFERENCE LINE MUST ALREADY BE KNOWN BY A CALL OF THE SAME ION 
c  (XSSN_Hlike OR ANY OTHER SUBR)
c
c Using  alphaeff fits by H. Nussbaumer & P.J. Storey, A&A 1984 and NIST levels
c   data read in, eg, CIII_di.lab
c
      IMPLICIT NONE
c
      CHARACTER niv_input*(*),lambda_output*(*),PHYAT*(*)
c
c   niv_input: atomic data, etc. Ex: 'CIII_di.lab'
c   lambbda_output: output, including for X-SSN. Ex: 'CIII_di.res'
c
      DOUBLE PRECISION hc,hc8,hcsk,Ryd_inf
      DOUBLE PRECISION int_ref,lambda_ref
      DOUBLE PRECISION vzero
      DOUBLE PRECISION Te,Ne,t4
      DOUBLE PRECISION wl(99,99),relint(99,99),relref(99),va
      DOUBLE PRECISION cdi(99,5) !coeff alla NS84
      COMMON/MOCR1/wlmin,wlmax,wlim,Irelmin_H,Ionab,Ionab_min,iwlim
      DOUBLE PRECISION int_ref_mod_H
      COMMON/MOCR2/int_ref_mod_H
      DOUBLE PRECISION Ionab,Ionab_min
      DOUBLE PRECISION wlim(20),Irelmin_H(21)
      DOUBLE PRECISION wlmin,wlmax
      DOUBLE PRECISION Irelmin(21)
      COMMON/MOCR3/Irelmin
      DOUBLE PRECISION cmult
      INTEGER imult,ifre1,ifre1_ref,ifre2,mult,nmult
      INTEGER iwlim,iw,iecr,icomment
      INTEGER ncomp(99),kmax,k,nc,l,iraies
      CHARACTER*17 titre
      CHARACTER numion*4,numion4*4,nomion*7 !format: ' 705', '0705', 'N_V    '
      CHARACTER Irelminc*5,lambda_refc*6,iraiesc*4,wlminc*8,wlmaxc*8
      CHARACTER bidon*1
      CHARACTER commentraie(99,30)*28,comment(99)*28,comw*28,comwb*28
c
      INTEGER in_niv,out_niv,iwr_phyat
      COMMON/MOC1/in_niv,out_niv,iwr_phyat
c
      COMMON/MOCC1/hc,hc8,hcsk,Ryd_inf
      COMMON/MOCV1/vzero
      COMMON/MOCX1/lambda_refc
      COMMON/MOCX2/int_ref,lambda_ref,ifre1_ref
c
      print *,niv_input,'  ',lambda_output
c
      t4=Te/1.d4
      iraies=0 !counter of X-SSN emission lines
c
        OPEN(unit=in_niv,file=niv_input,status='old')
c
      comwb='                            ' ! 28 blanks
c
      read(in_niv,15) titre
 15   format(A17)
 1    format(a1)
c Stdd di_electronics : imult = 2
      read(in_niv,121) mult,nmult,imult,cmult,ifre1,ifre2
 121  format(6x,i4,6x,i4,6x,i4,6x,d7.0,6x,i4,6x,i4)
      read(in_niv,1) bidon
      read(in_niv,60) numion,numion4,nomion
 60   format(9x,a4,11x,a4,10x,a7)
      read(in_niv,1) bidon   !notations
c Usage comments: 1=mult only (comment), 2=each emission line (commentraie)
        Read(in_niv,22) icomment
        Read(in_niv,22) kmax  !nber of diel multiplets
 22   format(20i4)
        Do k=1,kmax
         Read(in_niv,32) ncomp(k),wl(k,1),(cdi(k,l),l=1,5),comment(k)
 32     Format(i2,f10.3,5f8.4,a28)
c (int_ref from previous computation of radiative lines)
         relref(k)=1.d-12/t4**1.5*dexp(-cdi(k,5)/t4)*hc8/wl(k,1)*
     $      (cdi(k,1)/t4+cdi(k,2)+cdi(k,3)*t4+cdi(k,4)*t4*t4)/int_ref
         relint(k,1)=1.d0
         if(ncomp(k).gt.1) then
          do nc=1,ncomp(k) !(wl(k,1) erased if several components are read)
           Read(in_niv,33) wl(k,nc),relint(k,nc),commentraie(k,nc)
 33     Format(f10.3,f6.3,a28)
          enddo
         endif
        Enddo
c
        CLOSE(in_niv)
c
c  Wavelength ranges:
c
      do iw=1,iwlim+1
         Irelmin(iw)=Irelmin_H(iw)*
     $   int_ref_mod_H/(int_ref*dmax1(Ionab,Ionab_min))
      enddo
c     
	   OPEN(unit=out_niv,file=lambda_output,status='old')
c
      write(out_niv,15) titre
      write(out_niv,116)
 116   format(' bbb',/,' bbb')
      write(out_niv,141) mult,nmult,imult,cmult,ifre1,ifre2
 141   format('  mult =',i4,' nmult =',i4,' imult =',i4,
     $   ' cmult =',f7.4,' ifre1 =',i4,' ifre2 =',i4)
      write(out_niv,67) Ne,Te
 67   format(1p,' Ne =',e9.2,'  Te =',e9.2)
      write(out_niv,62) numion,numion4,nomion
 62   format(' numion :',a4,' numion4 :',a4,' nomion :',a7,' a4,a4,a7')
      write(out_niv,63) wlmin,wlmax
 63   format(1p,2e10.3,' wlmin, wlmax X-SSN')
      write(out_niv,85) iwlim
 85   format(' iwlim ='i4,' puis wlim(1->iwlim), Irelmin(1->iwlim+1)')
      write(out_niv,87) (wlim(iw),iw=1,iwlim)
 87   format(5x,1p,13e10.3)
      write(out_niv,86) (Irelmin(iw),iw=1,iwlim+1)
 86   format(1p,13e10.3)
      write(out_niv,88) kmax
 88   format(' kmax =',i4,' nbre de multiplets')
c
c  **********************
c  Tabulation for X-SSN :
c  **********************
c
      Write(out_niv,106) lambda_ref,int_ref
 106  Format(' X-SSN: lambda reference effective =',
     $  1p,e12.6,' effective emissivity =',e9.3,' erg.cm3.s-1')
c
 212  Format(' Sub-Ref. for liste_phyat.dat :')
      Write(out_niv,112) numion,ifre1,imult,ifre2,nomion,
     $  numion,ifre1_ref,int(lambda_ref)
 112  Format(1x,A4,3I1,'000000',1x,A7,'        1.0   0.     1.000e+0', ! 3I1=290
     $  '  1.     ',A4,I1,'00000000   1   1.00 ',I5,'++')
c
        Do k=1,kmax
         Do nc=1,ncomp(k)
c wavelength range:
          IF(wl(k,nc).gt.wlmin.and.wl(k,nc).lt.wlmax) THEN !wl range
c
          va=relref(k)*relint(k,nc)
c
          iecr=0
c
        iw=1
         if(wl(k,nc).le.wlim(iw).and.va.ge.Irelmin(iw)) then
          iecr=1
         endif
        do iw=2,iwlim
           if(wl(k,nc).gt.wlim(iw-1).and.wl(k,1).le.wlim(iw)) then
            if(va.ge.Irelmin(iw)) then
             iecr=1
            endif
           endif
        enddo
        if(wl(k,nc).gt.wlim(iwlim).and.va.ge.Irelmin(iwlim+1)) then
          iecr=1
        endif
c
            IF(iecr.eq.1) THEN !iecr=1
c Counter X-SSN emission lines:
               iraies=iraies+1
c Comments written according to icomment:
         comw=comwb ! reinit 28 blanks
         if(icomment.eq.1.and.nc.eq.1) comw=comment(k)
         if(icomment.eq.2) comw=commentraie(k,nc)
         if(icomment.eq.2.and.ncomp(k).eq.1) comw=comment(k)
c
             If(k.lt.10.and.nc.lt.10) then ! suppose k < 100, nc < 100
       Write(out_niv,40) numion,ifre1,imult,ifre2,k,nc,nomion,wl(k,nc),
     $     va,numion,ifre1,imult,ifre2,comw
 40    Format(1x,A4,3I1,'00',I1,'00',I1,1x,A7,f13.3,' 0.   ',
     $      1p,e10.3,'  1.     ',A4,3I1,'000000   1   1.00 ',a28)
             Elseif(k.ge.10.and.nc.lt.10) then
       Write(out_niv,41) numion,ifre1,imult,ifre2,k,nc,nomion,wl(k,nc),
     $     va,numion,ifre1,imult,ifre2,comw
 41    Format(1x,A4,3I1,'0',I2,'00',I1,1x,A7,f13.3,' 0.   ',
     $      1p,e10.3,'  1.     ',A4,3I1,'000000   1   1.00 ',a28)
             Elseif(k.lt.10.and.nc.ge.10) then
       Write(out_niv,42) numion,ifre1,imult,ifre2,k,nc,nomion,wl(k,nc),
     $     va,numion,ifre1,imult,ifre2,comw
 42      Format(1x,A4,3I1,'00',I1,'0',I2,1x,A7,f13.3,' 0.   ',
     $      1p,e10.3,'  1.     ',A4,3I1,'000000   1   1.00 ',a28)
             Elseif(k.ge.10.and.nc.ge.10) then
       Write(out_niv,43) numion,ifre1,imult,ifre2,k,nc,nomion,wl(k,nc),
     $     va,numion,ifre1,imult,ifre2,comw
 43      Format(1x,A4,3I1,'0',I2,'0',I2,1x,A7,f13.3,' 0.   ',
     $      1p,e10.3,'  1.     ',A4,3I1,'000000   1   1.00 ',a28)
             Endif
            ENDIF !iecr
c
          ENDIF !wl range
         Enddo  !nc
        Enddo  !k
c
         Write(out_niv,75) iraies
 75      format('Number of X-SSN lines: iraies =',i5)
c
      write(iraiesc,'(I4)') iraies
      write(Irelminc,'(I5)') int(Irelmin(4)*1.d5)
      write(wlminc,'(I8)') int(wlmin)
      write(wlmaxc,'(I8)') int(wlmax)
      print *,'  wl range:',wlminc,wlmaxc
      print *,'number of emission lines = ',iraiesc,' / ',
     $ '   Irelmin =',Irelminc,'E-5',lambda_refc
c
	   CLOSE(out_niv)
c
c  If requested, insert the result .res in PHYAT for SYNTHESIS:
      if(iwr_phyat.ne.0) then
         Call Write_Phyat(lambda_output,PHYAT)
      endif
c
        Return
c
        End
c
c-----------------end XSSN_di--------------------
c--------------------------------------------------------
c
c Usage: Call HeI_BP(Te,Ne,HeI.lab,HeI.res,PHYAT)

      Subroutine HeI_BP(Te,Ne,niv_input,lambda_output,PHYAT)
c
c ****** Opening files: ******
c   1/ niv_input = HeI.lab (unit in_niv = 10)
c   2/ lambda_output = HeI.res (unit out_niv = 11)
c   3/ HeI_Porter13.dat (unit 15)
c   4/ TeK_1_0_50_Results.txt, where TeK = 5000 (1000) 25000 (unit 16)
c   5/ outputcond.dat (unit in_lim = 32)
c
c ***Bauman, Porter et al 2005, 2007 --> J-resolved, n levels, see update?!
c *** assumed: n=50, recombination topoff = 1, intercombination = 0
c      --> HeI_Bauman07.dat (the TeK_1_0_50_Results.txt), BUT only Ne=1cm-3 for now (?)
c ***Porter et al 2013 update + erratum 
c      --> HeI_Porter13.dat = 44 mult = f(Te,Ne), used to scale high level lines
c          Te=5000-25000(1000) lgNe= 1-14(1)

c  1/ Options for Wlgh's: from Bauman07 (WL_B07=1) OR recalculated by X-SSN (WL_B07=0 STANDARD)
c  2/ Only triplets involving 2^3P are split into 2 lines (2^3P12,2^3P0)
c  3/ Use 40 multiplets n < 6 Porter_2013 at given Te, Ne.
c  4/ Assume that mults from high levels behave like mults from n=5 
c    --> 8 ref mult: 2p-5s 2s-5p 2p-5d 3d-5f for singlet and triplet
c        behaviour g,h,i,k assumed the same as 5f
c    --> corresponding scaling coeff I_P13/I_B07 applied to HeI_Bauman07.dat
c
      IMPLICIT NONE
c
      CHARACTER niv_input*(*),lambda_output*(*),PHYAT*(*)
c
      INTEGER ncounter  ! proviso pr check
c      DOUBLE PRECISION emi(10000),wli(10000),ES(10000),EI(10000)
      DOUBLE PRECISION hc,hc8,hcsk,Ryd_inf
      DOUBLE PRECISION ryd,charge,chi_Hli,chi_ion,x,A
      DOUBLE PRECISION vzero
      COMMON/MOCR1/wlmin,wlmax,wlim,Irelmin_H,Ionab,Ionab_min,iwlim
      DOUBLE PRECISION int_ref_mod_H
      COMMON/MOCR2/int_ref_mod_H
      DOUBLE PRECISION Ionab,Ionab_min
      DOUBLE PRECISION wlim(20),Irelmin_H(21)
      DOUBLE PRECISION wlmin,wlmax
      DOUBLE PRECISION Irelmin(21)
      COMMON/MOCR3/Irelmin
      CHARACTER bidon*1
      CHARACTER zero*1,nc1*1,nc2*2,nc3*3
      CHARACTER ns_alph*2,ni_alph*2,lp1_alph*2,ns_alph3*3
      CHARACTER titre*27
      CHARACTER*3 inmitabc,inmatabc,ni_refc,ns_refc
      CHARACTER Irelminc*5,lambda_refc*6,iraiesc*4,wlminc*8,wlmaxc*8
      CHARACTER numion*4,numion4*4,nomion*7
      CHARACTER*1 nom(21),nom2
      CHARACTER tabsto*9
      DOUBLE PRECISION Ne,Te
      DOUBLE PRECISION cmult,comult(99,99)
      DOUBLE PRECISION E(99,99),EN(99,99),cor(99,99),corN(99,99)
      DOUBLE PRECISION wl(99,99,99)      !trans 'normales' ns l+1 --> ni l
      DOUBLE PRECISION wk(99,99,99)      !trans eventuelles ns l-1 --> ni l
      DOUBLE PRECISION wl0(99),wk0(99)   !trans avec 2^3P0 pour HeI
      INTEGER ls !added for HeI
      INTEGER it,in
c 14300 <--> 4.2um (include dominant trans 5-4) for testrun_50
      DOUBLE PRECISION wl_m(44),I_m(44),wvlng_b(14300),j_int_b(14300)
      CHARACTER*1 Ls_m(44),Li_m(44),Ls_b(14300),Li_b(14300)
      DOUBLE PRECISION int_ref,lambda_ref,leakns_ref
      DOUBLE PRECISION Int1_m(7,4,8,8),Int3_m(7,4,8,8) !ns_m,ni_m,lp1s_m,lp1i_m
      DOUBLE PRECISION Int1(50,14,8,8),corr_I(3,8)  !corr_I(mult,lp1s)
      DOUBLE PRECISION Int3(50,14,8,8)              !ns,ni,lp1s,lp1i
      DOUBLE PRECISION Int3pi(50,14,8,8),Int3ps(50,14,8,8)
      DOUBLE PRECISION wl1(50,14,8,8)
      DOUBLE PRECISION wl3(50,14,8,8)
      DOUBLE PRECISION wl3ps(50,14,8,8),wl3pi(50,14,8,8)
      DOUBLE PRECISION wlmaxscale
      DOUBLE PRECISION wlmax_Bau,wlmax_Por
      DOUBLE PRECISION rap
      DOUBLE PRECISION va,vb,vc
      INTEGER WL_B07,i_rwd
      INTEGER i,j,k,l
      INTEGER ns_m(44),mult_m(44),ni_m(44),ntot_m,n_ref_m
      INTEGER lp1s_m(44),lp1i_m(44),lp3s_m(44),lp3i_m(44)
      INTEGER ns_b(14300),mults_b(14300),Js_b(14300)
      INTEGER ni_b(14300),multi_b(14300),Ji_b(14300)
      INTEGER lp1s_b(14300),lp1i_b(14300)
      INTEGER levs_b(14300),levi_b(14300),ntot,Te_m,dlgNe_m
      INTEGER ns_max,ni_max,lp1s_max,lp1i_max
c
      INTEGER lmax_del1,lmax_del2,lmax_del3
      INTEGER n,nmi,nma,nma_orig,nma2,nh2,lmi,lmi_orig,lma ! l ailleurs
      INTEGER lp1,lp1i,lp1s
      CHARACTER*1 lp1i_alph,lp1s_alph
      INTEGER nmil1,nmal1,lmal1,nmil2,nmal2,lmal2
      INTEGER lmin,lmax,lmin_int,lmax_int,li_wr,ls_wr
      INTEGER nmitab,nmatab,ni,ns,n1,n2,ini(99)
      INTEGER inmitab,inmatab,iraies
c
      INTEGER mult,imult,imult_ref,nmult,nimult(20),limult(20) ! imult_ref pour HeI
      INTEGER niregr,nsregr,indregr,lregr,hydrostrict
      INTEGER ifre1,ifre2,ifre2_3P0,ifre1_ref,ifre2_ref
      INTEGER ni_ref,ns_ref,li_ref   !HeI 4471 : for synth impr
      INTEGER numin,numax !Br table
      DOUBLE PRECISION Mat
      REAL az,am
      INTEGER iwlim,iw,iecr,inep1
c
      INTEGER in_niv,out_niv,iwr_phyat
      COMMON/MOC1/in_niv,out_niv,iwr_phyat
      COMMON/MOCH2/az,am,numin,numax
      COMMON/MOCX1/lambda_refc
      COMMON/MOCX2/int_ref,lambda_ref,ifre1_ref
      COMMON/MOCV1/vzero
      COMMON/MOCC1/hc,hc8,hcsk,Ryd_inf

      CHARACTER typemis*19,Tec4*4,Tec5*5
      INTEGER Tei

c Bauman et al 2007 /w topoff, no intercomb, n=50, Ne=1/cm3
c   Temperature 5000 (1000) 25000 added in front of typemis to define table to be used
      DATA typemis/'_1_0_50_Results.txt'/
      DATA ntot_m/44/     !number of mult in Porter13
      DATA n_ref_m/10/    !mult 4471 in Porter13
      DATA ni_ref/2/ns_ref/4/li_ref/1/ !4471 : for synth impr
c lambda_ref NIST:
c 4471.4703652  6.8275e+05   1s.2p 3P* 2    1s.4d 3D  1 }
c 4471.4740693  6.1440e+06   1s.2p 3P* 2    1s.4d 3D  2 }
c 4471.4743096  2.4579e+07   1s.2p 3P* 2    1s.4d 3D  3 } --> 3P12
c 4471.4856503  1.0241e+07   1s.2p 3P* 1    1s.4d 3D  1 }
c 4471.4893544  1.8432e+07   1s.2p 3P* 1    1s.4d 3D  2 }
c 4471.6832433  1.3655e+07   1s.2p 3P* 0    1s.4d 3D  1   --> 3P0
c EXPLICIT:
      DATA wlmax_Bau/4.160d+4/ ! Bauman07 table: 4.160d+4 for dimension = 14300 ***
      DATA wlmax_Por/2.170d+4/ ! Porter13 : 2.170d+4 to reach f scaling n=7 (1.280d+4 n=5)
c
c orbital 's' in 21st pos :
      DATA nom/'p','d','f','g','h','i','k','l','m','n','o','q','r'
     $   ,'t','u','v','w','x','y','z','s'/
c
      DATA zero/'0'/
c
c Correc atm refrac x=lambda_vac in um, dry air 15C p=760mmHg, Allen_3, p.124
      A(x)=1.+(64.328+29498.1/(146.-1./x**2)+255.4/(41.-1./x**2))*1.d-6
c
      wlmaxscale=dmax1(wlmax_Bau,wlmax_Por)
c
c  Wavelength ranges:
c
      do iw=1,iwlim+1
         Irelmin(iw)=Irelmin_H(iw)*
     $   int_ref_mod_H/(int_ref*dmax1(Ionab,Ionab_min))
      enddo
c
        OPEN(unit=in_niv,file=niv_input,status='old')
c
c beginning HeI.lab read and write only once:
c
c### CHOICE origin wavelengths : WL_B07=0 from HeI.dat (stdd), WL_B07=1 from Bauman07
      read(in_niv,15) titre,WL_B07
 15   format(a27,i2)
c###
 1    format(a1)
c rydberg, charge, ionization potential, atomic mass (uma)
      read(in_niv,20) ryd,charge,chi_ion,Mat
 20   format(4x,d12.3,7x,d4.0,8x,d12.3,4x,d5.0)
c potential H-like cm-1 :
      chi_Hli=ryd*charge*charge
c
      read(in_niv,21) mult,nmult,imult,cmult,ifre1,ifre2
 21   format(6x,i4,6x,i4,6x,i4,6x,d7.0,6x,i4,6x,i4)
        read(in_niv,1) bidon
c Table branching ratios Br : (UNUSEFUL here)
      read(in_niv,66) numin,numax
 66   format(6x,i4,6x,i4)
c Labeling + name ion (X-SSN format):
      read(in_niv,60) numion,numion4,nomion
 60   format(9x,a4,11x,a4,10x,a7)
c
c Transition ref emissivity X-SSN :
      read(in_niv,22) ni_ref,li_ref,ns_ref,hydrostrict !hydrostrict=0
 22   format(20i4)
c Int de ref : somme des l de li_ref a ni_ref-1. Check :
       if(li_ref.lt.0.or.li_ref.gt.ni_ref-1) then
        print *,'**** li_ref=',li_ref,' out of range. STOP'
        STOP
       endif
c  correction for possible 'leak' fm level ns_ref
      read(in_niv,25) leakns_ref ! +/- useful ? usually 1.0
 25   format(d9.3)

      read(in_niv,22) nmi,nma,lmi,lma
c nma2 conserve la valeur initialmt lue de nma = n maxi considere :
        nma2=nma !tabulation PSto nl-nu H et He+ jusqu'a nl=10 et nu=500
c Lim. lecture energies : lab (NIST etc) puis effectives (X-SSN)
      read(in_niv,22) nmil1,nmal1,lmal1,nmil2,nmal2,lmal2
c Lim. tables X-SSN et limites pour regroupements des l
      read(in_niv,22) nmitab,nmatab,niregr,nsregr
      read(in_niv,22) lmax_del1,lmax_del2 !general puis X-SSN
c * Suite lecture in_niv +bas ***
c
      OPEN(unit=out_niv,file=lambda_output,status='old')
c
c beginning WRITE out_niv = HeI.res
       if(WL_B07.eq.0) then
        write(out_niv,115) titre,WL_B07
 115  format(A27,' wl computed from HeI.dat (WL_B07 =',i2,')')
       else
        write(out_niv,116) titre,WL_B07
 116  format(A27,' wl from table Bauman07 (WL_B07 =',i2,')')
       endif
      write(out_niv,401) ryd,charge,chi_ion,chi_Hli,Mat
 401  format(1p,' ryd =',e17.10,' charge =',e8.1,' chi_ion =',e17.10,
     $  /,' chi_Hli =',e17.10,' Mat =',e8.1)
      write(out_niv,41) mult,nmult,imult,cmult,ifre1,ifre2
 41   format('  mult =',i4,' nmult =',i4,' imult =',i4,
     $   ' cmult =',f7.4,' ifre1 =',i4,' ifre2 =',i4,' ni,li-mult :')
        write(out_niv,1) bidon  !to prevent blanck line
      write(out_niv,167) Ne,Te
 167  format(1p,' Ne =',e9.2,'  Te =',e9.2)
      write(out_niv,62) numion,numion4,nomion
 62   format(' numion :',a4,' numion4 :',a4,' nomion :',a7,' a4,a4,a7')
      write(out_niv,63) wlmin,wlmax
 63   format(1p,2e10.3,' wlmin, wlmax X-SSN')
      write(out_niv,85) iwlim    !iwlim for Irelmin
 85   format(' iwlim ='i4,' puis wlim(1->iwlim), Irelmin(1->iwlim+1)')
      write(out_niv,87) (wlim(iw),iw=1,iwlim)
 87   format(5x,1p,13e10.3)
      write(out_niv,86) (Irelmin(iw),iw=1,iwlim+1)
 86   format(1p,13e10.3)
      write(out_niv,42) nmi,nma,lmi,lma,nma2
 42   format(5i4,'  nmi,nma,lmi,lma,nma2 = nma lu')
      write(out_niv,72) nmil1,nmal1,lmal1,nmil2,nmal2,lmal2 !only nmal1 depends on mult
 72   format(6i4,'  nmil1,nmal1,lmal1,nmil2,nmal2,lmal2 lectures')
      write(out_niv,70) ni_ref,li_ref,ns_ref,hydrostrict
 70   format(4i4,'  ni_ref,l=li_ref->ni_ref-1,ns_ref,hydrostrict=0-2')
      write(out_niv,71) leakns_ref
 71   format(1p,d10.3,'  fact. leakns_ref')
      write(out_niv,64) nmitab,nmatab,niregr,nsregr
 64   format(4i4,'  nmitab,nmatab,niregr,nsregr X-SSN')
      write(out_niv,46) lmax_del1,lmax_del2,lmax_del3
 46   format(3i4,'  lmax_del1,lmax_del2,lmax_del3')
c
      if(WL_B07.eq.0) then
      print *,niv_input,' ',lambda_output,' (wavelengths by X-SSN)'
      else
      print *,niv_input,' ',lambda_output,' (wavelengths from Bauman07)'
      endif
c
c##
c  Pre-reading just for Write in case data in HeI.lab will be adopted: 
      IF(WL_B07.eq.0) THEN  !wl's: 0 <--> using data in HeI.lab (stdd); else: Bauman07
c##
            imult=-1
             DO mult=3,1,-2
            imult=imult+1
      do n=1,nma
        do lp1=1,n
         E(n,lp1)=0.d0
         EN(n,lp1)=0.d0
         cor(n,lp1)=0.d0
         corN(n,lp1)=0.d0
        enddo
      enddo
c READ   reading triplet first
c  EXPLICIT nmal1:
        if(mult.eq.3) nmal1=10
        if(mult.eq.1) nmal1=15
      do n=nmil1,nmal1
      lmax=min0(lmal1,n-1)
       read(in_niv,24) (EN(n,lp1),lp1=lmi+1,lmax+1)
 24    format(6d12.3)
        do lp1=lmi+1,lmax+1
         E(n,lp1)=EN(n,lp1)
         corN(n,lp1)=E(n,lp1)-chi_ion+chi_Hli/dble(n*n)
        enddo
      enddo
      do n=nmil2,nmal2
      lmax=min0(lmal2,n-1)
       read(in_niv,26) (cor(n,lp1),lp1=lmi+1,lmax+1)
      enddo
 26    format(10d8.2)
c WRITE
      do n=nmil1,nmal1
      lmax=min0(lmal1,n-1)
       write(out_niv,44) n,(EN(n,lp1),lp1=lmi+1,lmax+1)
 44   format(i3,(10f12.4))
       write(out_niv,45) n,(cor(n,lp1),lp1=lmi+1,lmax+1)
 45   format(i3,(10f12.4)) !(i3,(10(4x,f8.4))) 11/13
      enddo
c (corN repeated for easier comparison/transcription to cor)
      write(out_niv,49) 
 49   format(4x,'corN(n,l)  Lab')
      do n=nmil1,nmal1
      lmax=min0(lmal1,n-1)
      write(out_niv,51) n,(corN(n,lp1),lp1=lmi+1,lmax+1)
 51   format(i3,6x,(10f12.4))  !(i3,6x,(10f8.4))
      enddo
      do n=nmil2,nmal2
      lmax=min0(lmal2,n-1)
       write(out_niv,48) n,lmi,lmax,(cor(n,lp1),lp1=lmi+1,lmax+1)
      enddo
 48   format(3i3,(10f12.4))  !(3i3,(10f8.4))
c
             ENDDO
      ENDIF
c
      REWIND(in_niv)
      do i=1,12
        read(in_niv,1) bidon
      enddo
c
c End pre-reading. See effective reading of E's, etc.  below
c
c******************debut
c Special HeI   Bauman - Porter
c******************
c  Look for nearest log(Ne) without interpolation:
c  The first log(Ne) (in=1) is 1 ; steps 1 ; last is 14 (in=14)
      in=IDNINT(dlog10(Ne))
c
      If(in.lt.1) then
       in=1
       elseif(in.gt.14) then
       in=14
      Endif
c
c  Look for nearest Te without interpolation:
c  The first Te (it=1) is 5000; steps 1000; last 25000 (it=21)
      it=IDNINT(Te/1000.)-4
      If(it.lt.1) then
       it=1
       elseif(it.gt.21) then
       it=21
      Endif
c
        OPEN(unit=15,status='old',file='HeI_Porter13.dat')
c
          do i=1,9 !header = 9 lines, incl 2 blanks (made explicit)
           read(15,1) bidon
          enddo
c  Low multiplets: approximate wl(A) and designation
c Notation: _m for multiplets from HeI_Porter13.dat
          do i=1,ntot_m       !(no intercombination transition)
      Read(15,1010) wl_m(i),ns_m(i),mult_m(i),Ls_m(i),ni_m(i),Li_m(i)
c      print *, wl_m(i),ns_m(i),mult_m(i),Ls_m(i),ni_m(i),Li_m(i)
          enddo
 1010 format(x,f5.0,2x,i2,3x,i1,x,a1,2x,i2,5x,a1)
          do i=1,4  ! 4 lines of which 1 long
           read(15,1) bidon
          enddo
c
        If(in.eq.1) then
         if(it.gt.1) then
          do i=1,it-1  ! long lines
           read(15,1) bidon
          enddo
         endif
        Else
         do j=1,in-1
          do i=1,23  ! 21 long + 2 blank
           read(15,1) bidon
          enddo
         enddo
         if(it.gt.1) then
          do i=1,it-1  ! long
           read(15,1) bidon
          enddo
         endif
        Endif
       Read(15,*) Te_m,dlgNe_m,(I_m(i),i=1,ntot_m) !erg cm^3 s^-1 !ntot_m=44
        print *,Te_m,dlgNe_m
c        print *,(I_m(i),i=1,ntot_m)
c
        CLOSE(15)  ! Porter13
c
c Conversion to J+1:
      Do i=1,ntot_m
         lp1s_m(i)=0
         lp1i_m(i)=0
         lp3s_m(i)=0
         lp3i_m(i)=0
         If(mult_m(i).eq.1) Then  !mult=1
        if(Ls_m(i).eq.'S') lp1s_m(i)=1
        if(Ls_m(i).eq.'P') lp1s_m(i)=2
        if(Ls_m(i).eq.'D') lp1s_m(i)=3
        if(Ls_m(i).eq.'F') lp1s_m(i)=4
        if(Li_m(i).eq.'S') lp1i_m(i)=1
        if(Li_m(i).eq.'P') lp1i_m(i)=2
        if(Li_m(i).eq.'D') lp1i_m(i)=3
         Else                     !mult=3
        if(Ls_m(i).eq.'S') lp3s_m(i)=1
        if(Ls_m(i).eq.'P') lp3s_m(i)=2
        if(Ls_m(i).eq.'D') lp3s_m(i)=3
        if(Ls_m(i).eq.'F') lp3s_m(i)=4
        if(Li_m(i).eq.'S') lp3i_m(i)=1
        if(Li_m(i).eq.'P') lp3i_m(i)=2
        if(Li_m(i).eq.'D') lp3i_m(i)=3
         Endif
      Enddo
c Emissivity of reference multiplet: # 4471-- 4d^{3}D - 2p^{3}P
       int_ref=10.**I_m(n_ref_m) !alias Int3_m(4,2,3,2) !n_ref_m=10
c
c Normalized emissivity of multiplets in Porter13: Int1_m(7,4,4,3),..
      Do l=1,7  !ns_m
       Do k=1,4  !ni_m
        Do j=1,4  !lp1s_m
         Do i=1,3  !lp1i_m=lp1s_m +1 or -1
          Int1_m(l,k,j,i)=0.
          Int3_m(l,k,j,i)=0.
         Enddo
        Enddo
       Enddo
      Enddo
      Do i=1,ntot_m
       if(mult_m(i).eq.1) then
         Int1_m(ns_m(i),ni_m(i),lp1s_m(i),lp1i_m(i))=10.**I_m(i)/int_ref
       else
         Int3_m(ns_m(i),ni_m(i),lp3s_m(i),lp3i_m(i))=10.**I_m(i)/int_ref
       endif
      Enddo
c
c  Nearest tabulated Te (no interpolation) converted to character Tec. 
c   The first Tec is 5000 (it=1) ; steps 1000 ; last is 25000.
c  Opening file from Bauman07 for this Tec:
c
              IF(it.le.5) then
      write(Tec4,'(I4)') (it+4)*1000 ! conversion to character
c Bauman07
        OPEN(unit=16,status='old',file=Tec4//typemis)
c
              ELSE
      write(Tec5,'(I5)') (it+4)*1000 ! conversion to character
c
        OPEN(unit=16,status='old',file=Tec5//typemis)
c
              ENDIF
c
          do i=1,5 !header = 5 lines
           read(16,1) bidon
          enddo
c
      i=0
      Do
      i=i+1
      Read(16,1050,end=1000) ns_b(i),mults_b(i),Ls_b(i),Js_b(i),
     $             ni_b(i),multi_b(i),Li_b(i),Ji_b(i),
     $       wvlng_b(i),j_int_b(i),levs_b(i),levi_b(i)
c**
        if(wvlng_b(i).gt.wlmaxscale) go to 1000
c**
      Enddo
 1050 format(i2,x,i1,a1,i1,7x,i2,x,i1,a1,i1,
     $    d18.6,d14.6,56x,i4,x,i4)
 1000 continue
      ntot=i-1
cc      print *,'ntot=',ntot,' last line :',wvlng_b(ntot)
c
        CLOSE(16)  ! Bauman07
c
c  Init intensities: Int1(50,14,8,8)
      Do l=1,50  !ns_b(i)
       Do k=1,14  !ni_b(i)  !Note: 1st ni=7: wl=41680A N°=14350 > 14300
        Do j=1,8  !lp1s_b(i) !table Bauman07 up to lp1s=8
         Do i=1,8  !lp1i_b(i)=lp1s_b(i) +1 or -1 !table Bauman07 up to lp1i=8
          Int1(l,k,j,i)=0. !singlet
          Int3(l,k,j,i)=0. !triplet incl lp1s or lp1i=2 with 3P12 involved
          Int3ps(l,k,j,i)=0. !triplet lp1s=2 with upper 3P0 involved j=2 & i=1 or 3
          Int3pi(l,k,j,i)=0. !triplet lp1i=2 with lower 3P0 involved i=2 & j=1 or 3
         Enddo
        Enddo
       Enddo
      Enddo
c Conversion to J+1:
      Do i=1,ntot
         lp1s_b(i)=0
         lp1i_b(i)=0
        if(Ls_b(i).eq.'S') lp1s_b(i)=1
        if(Ls_b(i).eq.'P') lp1s_b(i)=2
        if(Ls_b(i).eq.'D') lp1s_b(i)=3
        if(Ls_b(i).eq.'F') lp1s_b(i)=4
        if(Ls_b(i).eq.'G') lp1s_b(i)=5
        if(Ls_b(i).eq.'H') lp1s_b(i)=6
        if(Ls_b(i).eq.'I') lp1s_b(i)=7
        if(Ls_b(i).eq.'K') lp1s_b(i)=8
        if(Li_b(i).eq.'S') lp1i_b(i)=1
        if(Li_b(i).eq.'P') lp1i_b(i)=2
        if(Li_b(i).eq.'D') lp1i_b(i)=3
        if(Li_b(i).eq.'F') lp1i_b(i)=4
        if(Li_b(i).eq.'G') lp1i_b(i)=5
        if(Li_b(i).eq.'H') lp1i_b(i)=6
        if(Li_b(i).eq.'I') lp1i_b(i)=7
        if(Li_b(i).eq.'K') lp1i_b(i)=8
      Enddo
       ns_max=0
       ni_max=0
       lp1s_max=1
       lp1i_max=1
c  Emissivities and wl's Bauman (wl's substituted later if WL_B07=0):
      Do i=1,ntot
c In the data set used, all intercombination lines have flux = 0.
cc  multiplets involving ^3P are split in 2 components, P0 and P12
cc  Here, wl() is the last wvlng_b(i) of multiplet or sub-multiplet
cc  (rigorously, if over one component in the sum, a weighted average 
cc  for wl should be computed; wl's will be recomputed if WL_B07=0 stdd)
       If(mults_b(i).eq.1.and.multi_b(i).eq.1) Then
        Int1(ns_b(i),ni_b(i),lp1s_b(i),lp1i_b(i))=
     $       Int1(ns_b(i),ni_b(i),lp1s_b(i),lp1i_b(i))+j_int_b(i)
        wl1(ns_b(i),ni_b(i),lp1s_b(i),lp1i_b(i))=wvlng_b(i)
       Elseif(mults_b(i).eq.3.and.multi_b(i).eq.3) Then
        if(lp1s_b(i).eq.2) then      ! P upper level involved
           if(Js_b(i).eq.0) then ! Js_b(i) = 0
        Int3ps(ns_b(i),ni_b(i),lp1s_b(i),lp1i_b(i))=
     $       Int3ps(ns_b(i),ni_b(i),lp1s_b(i),lp1i_b(i))+j_int_b(i)
        wl3ps(ns_b(i),ni_b(i),lp1s_b(i),lp1i_b(i))=wvlng_b(i)
           else                ! Js_b(i) = 1 or 2 grouped together
        Int3(ns_b(i),ni_b(i),lp1s_b(i),lp1i_b(i))=
     $       Int3(ns_b(i),ni_b(i),lp1s_b(i),lp1i_b(i))+j_int_b(i)
        wl3(ns_b(i),ni_b(i),lp1s_b(i),lp1i_b(i))=wvlng_b(i)
           endif
        elseif(lp1i_b(i).eq.2) then  ! P lower level involved
           if(Ji_b(i).eq.0) then ! Ji_b(i) = 0
        Int3pi(ns_b(i),ni_b(i),lp1s_b(i),lp1i_b(i))=
     $       Int3pi(ns_b(i),ni_b(i),lp1s_b(i),lp1i_b(i))+j_int_b(i)
        wl3pi(ns_b(i),ni_b(i),lp1s_b(i),lp1i_b(i))=wvlng_b(i)
           else                ! Ji_b(i) = 1 or 2 grouped together
        Int3(ns_b(i),ni_b(i),lp1s_b(i),lp1i_b(i))=
     $       Int3(ns_b(i),ni_b(i),lp1s_b(i),lp1i_b(i))+j_int_b(i)
        wl3(ns_b(i),ni_b(i),lp1s_b(i),lp1i_b(i))=wvlng_b(i)
           endif
        else         ! P not involved, no splitting
        Int3(ns_b(i),ni_b(i),lp1s_b(i),lp1i_b(i))=
     $       Int3(ns_b(i),ni_b(i),lp1s_b(i),lp1i_b(i))+j_int_b(i)
        wl3(ns_b(i),ni_b(i),lp1s_b(i),lp1i_b(i))=wvlng_b(i)
        endif
       Endif
        ns_max=MAX0(ns_max,ns_b(i))
        ni_max=MAX0(ni_max,ni_b(i))
        lp1s_max=MAX0(lp1s_max,lp1s_b(i))
        lp1i_max=MAX0(lp1i_max,lp1i_b(i))
      Enddo
      print *,'max ns,ni,lp1s,lp1i:',ns_max,ni_max,lp1s_max,lp1i_max
c
c  Ref scaling lines:
c  corr_I(mult,lp1s) for ns > 5. I_Port(Te,Ne)/I_Baum(Te) for ns=5
      corr_I(1,1)=Int1_m(5,2,1,2)*int_ref/Int1(5,2,1,2) !s
      corr_I(1,2)=Int1_m(5,2,2,1)*int_ref/Int1(5,2,2,1) !p
      corr_I(1,3)=Int1_m(5,2,3,2)*int_ref/Int1(5,2,3,2) !d
      corr_I(1,4)=Int1_m(5,3,4,3)*int_ref/Int1(5,3,4,3) !f
      corr_I(1,5)=corr_I(1,4)                !g
      corr_I(1,6)=corr_I(1,4)                !h
      corr_I(1,7)=corr_I(1,4)                !i
      corr_I(1,8)=corr_I(1,4)                !k
      corr_I(3,1)=Int3_m(5,2,1,2)*int_ref/
     $     (Int3(5,2,1,2)+Int3pi(5,2,1,2))            !s
      corr_I(3,2)=Int3_m(5,2,2,1)*int_ref/
     $     (Int3(5,2,2,1)+Int3ps(5,2,2,1))            !p
      corr_I(3,3)=Int3_m(5,2,3,2)*int_ref/
     $     (Int3(5,2,3,2)+Int3pi(5,2,3,2))            !d
      corr_I(3,4)=Int3_m(5,3,4,3)*int_ref/Int3(5,3,4,3) !f
      corr_I(3,5)=corr_I(3,4)                           !g
      corr_I(3,6)=corr_I(3,4)                           !h
      corr_I(3,7)=corr_I(3,4)                           !i
      corr_I(3,8)=corr_I(3,4)                           !k
c
c Normalized emissivity of multiplets ns > 5 from scaled Bauman07:
         if(ns_max.gt.5) then
      Do l=6,ns_max  !ns > 5
       Do k=2,ni_max  !ni
        Do j=1,lp1s_max  !lp1s
         Do i=1,lp1i_max  !lp1i=lp1s +1 or -1
          Int1(l,k,j,i)=Int1(l,k,j,i)*corr_I(1,j)/int_ref !singl
          Int3(l,k,j,i)=Int3(l,k,j,i)*corr_I(3,j)/int_ref !tripl incl lp1s=2 or lp1i=2 with 3P12 involved
          Int3ps(l,k,j,i)=Int3ps(l,k,j,i)*corr_I(3,j)/int_ref !tripl lp1s=2 with upper 3P0 involved, lp1i=1 or 3
          Int3pi(l,k,j,i)=Int3pi(l,k,j,i)*corr_I(3,j)/int_ref !tripl lp1i=2 with lower 3P0 involved, lp1s=1 or 3
         Enddo
        Enddo
       Enddo
      Enddo
         endif
c
c Normalized emissivity of multiplets ns < 6 from Porter13:
      ncounter=0
      Do l=2,5  !ns < 6 ! ns=2 ok for 2s-2p, cf max0/min0
       Do k=2,MAX0(2,MIN0(l-1,3))  !ni ; no 5-4 trans in Porter13
         If(l.gt.2) then !test is useful for 2s-2p    2 2 2 1
          j=1  !lp1s
          i=2  !lp1i
          Int1(l,k,j,i)=Int1_m(l,k,j,i) !singl
           rap=Int3pi(l,k,j,i)/(Int3pi(l,k,j,i)+Int3(l,k,j,i))
          Int3(l,k,j,i)=Int3_m(l,k,j,i)*(1.-rap) !tripl lp1i=2 /w 3P12 involved
          Int3pi(l,k,j,i)=Int3_m(l,k,j,i)*rap !tripl lp1i=2 /w lower 3P0 involved, lp1s=1 or 3
      ncounter=ncounter+1
         Endif
        Do j=2,l  !lp1s
        Do i=j-1,j+1,2  !lp1i=lp1s -1 or +1
      ncounter=ncounter+1
             If(i.le.k) then
         Int1(l,k,j,i)=Int1_m(l,k,j,i) !singl
        if(j.eq.2) then ! lp1s=2
           rap=Int3ps(l,k,j,i)/(Int3ps(l,k,j,i)+Int3(l,k,j,i))
          Int3(l,k,j,i)=Int3_m(l,k,j,i)*(1.-rap) !tripl lp1s=2 /w 3P12 involved
          Int3ps(l,k,j,i)=Int3_m(l,k,j,i)*rap !tripl lp1s=2 /w upper 3P0 involved, lp1i=1 or 3
        elseif(i.eq.2) then ! lp1i=2
           rap=Int3pi(l,k,j,i)/(Int3pi(l,k,j,i)+Int3(l,k,j,i))
          Int3(l,k,j,i)=Int3_m(l,k,j,i)*(1.-rap) !tripl lp1i=2 /w 3P12 involved
          Int3pi(l,k,j,i)=Int3_m(l,k,j,i)*rap !tripl lp1i=2 /w lower 3P0 involved, lp1s=1 or 3
        else
          Int3(l,k,j,i)=Int3_m(l,k,j,i) !tripl excluding lp1s=2 and lp1i=2
         endif
             Endif !(i.le.k)
         Enddo
        Enddo
       Enddo
      Enddo
c
c******************fin
c End Special HeI : wl & Int  Bauman - Porter
c******************
c
c  BEGINING SYNTH
c
      do n=1,nma
        ini(n)=0
      enddo
c###########
c          ifre1 stable
c          ifre2 stable
           ifre2_3P0=9 !2nd value ifre2 !EXPLICIT (by convention)
            imult=-1    !2 imult values
c###########
             DO mult=3,1,-2  ! triplet first
c###########
                imult=imult+1   ! imult=0 for triplet et =1 for singlet
c *Note : for correct labeling of HeI, imult=1 will appear in the place of ifre2
c    (ifre2_3P0=9 only used for transitions from 2^3P, will not interfere)
c###########
c  Choice wl, using EITHER HeI.dat (stdd below) OR Bauman07 (above)
c#################
          IF(WL_B07.eq.0) THEN  !wl's: 0 <--> using HeI.dat (stdd); else: Bauman07
c#################
c INIT wl Traitement triplet puis singlet (lectures successives)
      do n=1,nma
        do lp1=1,n  ! 9/15 l changed to lp1=l+1
         E(n,lp1)=0.d0
         EN(n,lp1)=0.d0
         cor(n,lp1)=0.d0
         corN(n,lp1)=0.d0
        enddo
      enddo
c READ   reading triplet first
c  EXPLICIT nmal1:
        if(mult.eq.3) nmal1=10
        if(mult.eq.1) nmal1=15 ! for np singlet (only nmal1 diff fm trip to sing)
      do n=nmil1,nmal1
      lmax=min0(lmal1,n-1)
       read(in_niv,24) (EN(n,lp1),lp1=lmi+1,lmax+1)
        do lp1=lmi+1,lmax+1
         E(n,lp1)=EN(n,lp1)
         corN(n,lp1)=E(n,lp1)-chi_ion+chi_Hli/dble(n*n)
        enddo
      enddo
      do n=nmil2,nmal2
      lmax=min0(lmal2,n-1)
       read(in_niv,26) (cor(n,lp1),lp1=lmi+1,lmax+1)
      enddo
c
      do n=nmi,nma
        do lp1=lmi+1,n
         E(n,lp1)=chi_ion-chi_Hli/dble(n*n)+cor(n,lp1)
        enddo
      enddo
c
c Tabulation of wavelengths wl in air in \AA with 
c  lower niv ni fm nmitab to nmatab and upper niv ns fm ni+1 to nma 
c
c    Special case 2s-2p:
       ni=2
       ns=2
       lp1=1
       va=E(ns,lp1+1)-E(ni,lp1)
       wl(ni,lp1,ns)=(1.d8/va)/A(1.d4/va)
           if(mult.eq.3) then
c *EXPLICIT: 1.03cm-1 between 2p3P0 et 2p3P12
       wl0(ns)=(1.d8/(va+1.03))/A(1.d4/(va+1.03)) !2^3P0, here 2^3P = niv sup
           endif
c
      DO ni=nmitab,nmatab !tabul wl, wk
c
      lmax=ni-1
c (If only the greatest l for each ni, eg: 3 (lmax_del1=2), 4 (lmax_del1=3)...)
c HeI: all 'l' : take large lmax_del1 (eg, 100) and lmi = 0
      lmin=max0(lmi,lmax-lmax_del1)
c IMPLICITLY HERE** : 
c   'l' corresp. to an orbital of lower level ni
c   spontaneous transition is l+1 --> l
       Do lp1=lmin+1,lmax+1  ! l of lower level (ni) lp1=l+1
         Do ns=ni+1,nma
          va=E(ns,lp1+1)-E(ni,lp1)
          if(va.gt.0.d0) then
           wl(ni,lp1,ns)=1.d8/va
          If(ni.eq.2.and.lp1.eq.2.and.mult.eq.3) wl0(ns)=1.d8/(va-1.03) ! 2^3P0-n^3D1
c   air refraction wl > 2000A AND wl < 20000A:
            IF(wl(ni,lp1,ns).gt.2.d3.and.wl(ni,lp1,ns).lt.2.d4)THEN
             wl(ni,lp1,ns)=wl(ni,lp1,ns)/A(1.d4/va)
          If(ni.eq.2.and.lp1.eq.2.and.mult.eq.3) 
     $                                wl0(ns)=wl0(ns)/A(1.d4/(va-1.03)) ! 2^3P0-n^3D1
            ENDIF
          else !va </= 0 pathol
           wl(ni,lp1,ns)=1.d0
       print *,nomion,' wl : ni,lp1,ns',ni,lp1,ns,
     $                '** E(ns,lp1+1)-E(ni,lp1)=',va
       print *,cor(ns,lp1+1),cor(ni,lp1)
       print *,E(ns,lp1+1),E(ni,lp1)
       print *,'**STOP******************************'
       STOP
          endif
         Enddo
       Enddo
c IMPLICITLY HERE** : 
c   'l' corresp. to an orbital of lower level ni
c   spontaneous transition is l-1 --> l
           If(lmin.lt.lmax) Then
       Do lp1=lmin+2,lmax+1  ! l of lower level (ni)
         Do ns=ni+1,nma
          va=E(ns,lp1-1)-E(ni,lp1)
          if(va.gt.0.d0) then
           wk(ni,lp1,ns)=1.d8/va
          If(ni.eq.2.and.lp1.eq.2.and.mult.eq.3) wk0(ns)=1.d8/(va-1.03) ! 2^3P0-n^3D1
c   air refraction wl > 2000A AND wl < 20000A:
            IF(wk(ni,lp1,ns).gt.2.d3.and.wk(ni,lp1,ns).lt.2.d4)THEN
             wk(ni,lp1,ns)=wk(ni,lp1,ns)/A(1.d4/va)
          If(ni.eq.2.and.lp1.eq.2.and.mult.eq.3) 
     $                                wk0(ns)=wk0(ns)/A(1.d4/(va-1.03)) ! 2^3P0-n^3D1
            ENDIF
          else !va </= 0 pathol
           wk(ni,lp1,ns)=1.d0
       print *,nomion,' wk : ni,lp1,ns',ni,lp1,ns,
     $                '** E(ns,lp1-1)-E(ni,lp1)=',va
       print *,cor(ns,lp1+1),cor(ni,lp1)
       print *,E(ns,lp1+1),E(ni,lp1)
       print *,'**STOP******************************'
       STOP
          endif
         Enddo
       Enddo
           Endif
c
      ENDDO ! DO ni (tabulation wl, wk X-SSN)
c
c  Storage:
c
c    Special case 2s^3S-2p^3P :
        ni=2
        ns=2
        lp1=1
           if(mult.eq.1) then
          wl1(ns,ni,lp1+1,lp1)=wl(ni,lp1,ns)
           else   ! mult.eq.3
          wl3(ns,ni,lp1+1,lp1)=wl(ni,lp1,ns)
          wl3ps(ns,ni,lp1+1,lp1)=wl0(ns)
           endif
c
      DO ni=nmitab,nmatab
c
      lmax=ni-1
      lmin=max0(lmi,lmax-lmax_del1)
c IMPLICITLY HERE** : 
c   'l' corresp. to orbital of lower level ni
c   spontaneous transition is l+1 --> l
       Do lp1=lmin+1,lmax+1  ! l of lower level (ni)
         do ns=ni+1,nma
           if(mult.eq.1) then
          wl1(ns,ni,lp1+1,lp1)=wl(ni,lp1,ns)
           else
          wl3(ns,ni,lp1+1,lp1)=wl(ni,lp1,ns)
          If(ni.eq.2.and.lp1.eq.2) then ! 2^3P0-n^3D1
            wl3pi(ns,ni,lp1+1,lp1)=wl0(ns)
          Endif
           endif
         enddo
c   spontaneous transition is l-1 --> l
         if(lp1-1.gt.0) then
         do ns=ni+1,nma
           if(mult.eq.1) then
          wl1(ns,ni,lp1-1,lp1)=wk(ni,lp1,ns)
           else
          wl3(ns,ni,lp1-1,lp1)=wk(ni,lp1,ns)
          If(ni.eq.2.and.lp1.eq.2) then ! 2^3P0-n^3S1
            wl3pi(ns,ni,lp1-1,lp1)=wk0(ns)
          Endif
           endif
         enddo
         endif
       Enddo
      ENDDO  ! ni
c
c#################
           ENDIF ! end WL_B07.eq.0 (substitution of wl)
c#################
c
c ***** Reference line only once:
          if(imult.eq.0) THEN
             imult_ref=imult  !anticipating possible write imult > 0
       lambda_ref=(wl3(4,2,3,2)*Int3(4,2,3,2)+wl3pi(4,2,3,2)*
     $  Int3pi(4,2,3,2))/(Int3(4,2,3,2)+Int3pi(4,2,3,2))
      Write(out_niv,906) wl3(4,2,3,2),wl3pi(4,2,3,2),lambda_ref,int_ref
 906  Format(/,' X-SSN: lambda reference =',1p,e12.6,' + ',e12.6,
     $ ' = ',e12.6,' emissivity =',e10.4,' erg.cm3.s-1')
c
c  START WRITING FOR SYNTHESIS: 
c   START ref line or sub-ref
c
      Write(out_niv,210)  ! imult = 0 ****
 210  Format(' Ref. for liste_modele.dat :')
      Write(out_niv,1101) numion,ifre1,nomion,vzero,int(lambda_ref),
     $  ni_ref,nom(li_ref),ns_ref,Te_m,dlgNe_m   !suppr nom2
 1101 Format(1x,A4,I1,'00000000',1x,A7,'        1.0   0.     1.000e+0',
     $  '  1.                 0   1  ',f5.1,1x,I5,'++',
     $  I2,A1,'^3P -',I2,i7,i2)

      Write(out_niv,211)
 211  Format(' Ref. for liste_phyat.dat :')
      Write(out_niv,1111) numion4,ifre1,nomion,int_ref,
     $  int(lambda_ref),ni_ref,nom(li_ref),ns_ref,Te_m,dlgNe_m  !suppr nom2
 1111 Format('9',A4,I1,'00000000',1x,A7,'        1.0   0.   ',
     $  1p,e10.3,'  1.               999   1   1.00 ',I5,'++',
     $  I2,A1,'^3P -',I2,i7,i2)
c line counter: for triplet and singlet together
               iraies=0
          endif  !end imult=0
c *****
c  Cas part 2s-2p :
        ni=2
        ns=2
        lp1=1
        lp1s=2
        ni_alph='02'
        ns_alph='02'
        li_wr=21
        ls_wr=1
        nom(li_wr)='s'
        nom(ls_wr)='p'
           if(mult.eq.1) then
          va=wl1(ns,ni,lp1+1,lp1)
          vb=Int1(ns,ni,lp1+1,lp1)
           else   ! mult.eq.3 (imult.eq.0)
          va=wl3(ns,ni,lp1+1,lp1)
          vb=Int3(ns,ni,lp1+1,lp1)
           endif
         if(va.gt.wlmin.and.va.lt.wlmax) then
            iraies=iraies+1
c ifre1,ifre1,imult: see below
       Write(out_niv,187) numion,ifre1,ifre1,imult,ni_alph,ns_alph,
     $    lp1i,lp1s,nomion,va,vb,numion,ifre1,ifre1,imult_ref,
     $    ni,nom(li_wr),mult,ns,nom(ls_wr)
            If(mult.eq.3) then ! 2^3P0
               iraies=iraies+1
          va=wl3ps(ns,ni,lp1+1,lp1)
          vb=Int3ps(ns,ni,lp1+1,lp1)
       Write(out_niv,1871) numion,ifre1,ifre1,ifre2_3P0,ni_alph,ns_alph,
     $    lp1i,lp1s,nomion,va,vb,numion,ifre1,ifre1,imult_ref,
     $    ni,nom(li_wr),mult,ns,nom(ls_wr)
            Endif ! 2^3P0
         endif ! wl range
c  end Cas part 2s-2p
c
c********  ********** ************
	 DO ni=nmitab,nmatab !ni lower lev ******* loop ni + external
c********  ********** ************
c For Write:
        if(ni.lt.10) then
       write(nc1,'(I1)') ni
        ni_alph=zero//nc1
        else   !(assume ni < 100)
       write(ni_alph,'(I2)') ni
        endif
        lmax=ni-1
        lmin=lmi  !normally lmi=0 for HeI
c
c*************** l+1 associated to ni *** ATT : order of loops modified for HeI
          DO lp1i=lmax+1,lmin+1,-1  ! ******** loop l intermed external
c***************
c*************** ls associated to ns *** ATT : loop added for HeI
           DO lp1s=lp1i+1,lp1i-1,-2  ! ******** loop ls intermed internal
c***************
            IF(lp1s.ge.1) THEN
c************************* (here nma < 51)
            DO ns=ni+1,nma !ns = niv sup  ****** loop ns internal
c***************
c For Write:
        if(ns.lt.10) then
       write(nc1,'(I1)') ns
        ns_alph=zero//nc1
        else   !(assume ns < 100)
       write(ns_alph,'(I2)') ns
        endif
c
c**********
c wl range:
c**********
        if(mult.eq.3) then
           va=wl3(ns,ni,lp1s,lp1i) ! wl
           vb=Int3(ns,ni,lp1s,lp1i) ! relative emissivity
        elseif(mult.eq.1) then
           va=wl1(ns,ni,lp1s,lp1i) ! wl
           vb=Int1(ns,ni,lp1s,lp1i) ! relative emissivity
        else
           print *,'**** HeI: mult=', mult,' .ne.1 or 3, STOP ****'
           STOP
        endif
             IF(va.gt.wlmin.and.va.lt.wlmax) THEN ! wl range
c
c  START WRITINGS LINES
c
          iecr=0
c
         iw=1
         if(va.le.wlim(iw).and.vb.ge.Irelmin(iw)) then
          iecr=1
         endif
        do iw=2,iwlim
         if(va.gt.wlim(iw-1).and.va.le.wlim(iw)
     $     .and.vb.ge.Irelmin(iw)) then
          iecr=1
         endif
        enddo
        if(va.gt.wlim(iwlim).and.vb.ge.Irelmin(iwlim+1)) then
          iecr=1
        endif
c
c*****##
              IF(iecr.eq.1) THEN
c*****##
c line counter:
               iraies=iraies+1
c counter of effectively written ni:
               ini(ni)=1
                li_wr=lp1i-1
               if(li_wr.le.0) then
c  nom(0) does not exist
                li_wr=21
               endif
                ls_wr=lp1s-1
               if(ls_wr.le.0) then
                ls_wr=21
               endif
c
c To correctly spot HeI (triplet (imult=0) and singlet (imult=1) are treated together) :
c  1/ ifre1=0 always (position 5)
c  2/ position 6 (normaly for imult) is also always =0. ifre1 repeated.
c  3/ imult=1 (singlet) is put in the place of ifre2 (position 7)
c  4/ ifre2 takes the value 9 for the 2nd transitions from 2^3P (position 7)
c
       Write(out_niv,187) numion,ifre1,ifre1,imult,ni_alph,ns_alph,
     $    lp1i,lp1s,nomion,va,vb,numion,ifre1_ref,ifre1_ref,imult_ref, !imult replacing ifre2
     $    ni,nom(li_wr),mult,ns,nom(ls_wr)
 187    Format(1x,A4,3I1,2A2,2I1,1x,A7,f13.3,' 0.   ',
     $    1p,e10.3,'  1.     ',   
     $    A4,3I1,'000000   1   1.00 ',I2,A1,I1,'-',I2,A1,1x,A8)
               If(mult.eq.3.and.ni.eq.2.and.lp1i.eq.2) then ! duplication mult. 2^3P-...
               iraies=iraies+1
           va=wl3pi(ns,ni,lp1s,lp1i)
           vb=Int3pi(ns,ni,lp1s,lp1i)
       Write(out_niv,1871) numion,ifre1,ifre1,ifre2_3P0,ni_alph,ns_alph, !ifre2=9
     $    lp1i,lp1s,nomion,va,vb,numion,ifre1_ref,ifre1_ref,imult_ref,
     $    ni,nom(li_wr),mult,ns,nom(ls_wr)
 1871   Format(1x,A4,3I1,2A2,2I1,1x,A7,f13.3,' 0.   ',
     $    1p,e10.3,'  1.     ',   
     $    A4,3I1,'000000   1   1.00 ',I2,A1,I1,'-',I2,A1,1x,A8)
               Endif
c
c*****##
              ENDIF ! iecr=1
c*****##
             ENDIF ! wl range
c***************
            ENDDO ! ns = niv sup  ****** loop ns, + internal
c***************
            ENDIF ! (ls+1.ge.1)
c***************
           ENDDO ! ******** loop ls+1 associated to ns, intermed internal, special HeI
c***************
          ENDDO ! ******** loop l+1 associated to ni, intermed external
c********  ********** ************
	 ENDDO ! ni = niv inf ******* loop ni, + external
c###########
             ENDDO ! loop mult 3 & 1
c###########
c
         Write(out_niv,75) iraies
 75      format('Number of X-SSN lines: iraies =',i5)
c
        CLOSE(in_niv) !HeI.dat
c
        CLOSE(out_niv) !HeI.res
c
        inmitab=0
       do ni=nmitab,nmatab
         if(ini(ni).ne.0.and.inmitab.eq.0) then
        inmitab=ni
         endif
       enddo
        inmatab=0
       do ni=nmatab,nmitab,-1
         if(ini(ni).ne.0.and.inmatab.eq.0) then
        inmatab=ni
         endif
       enddo
c Conversion to character
      write(inmitabc,'(I2)') inmitab
      write(inmatabc,'(I2)') inmatab
      write(wlminc,'(I8)') int(wlmin)
      write(wlmaxc,'(I8)') int(wlmax)
      write(iraiesc,'(I4)') iraies
      write(Irelminc,'(I5)') int(Irelmin(4)*1.d5)
      write(lambda_refc,'(I5)') int(lambda_ref)
      write(ni_refc,'(I2)') ,ni_ref
      write(ns_refc,'(I2)') ,ns_ref
      print *,'effective ni range:',inmitabc,inmatabc,
     $ '  wl range:',wlminc,wlmaxc
      print *,'number of emission lines = ',iraiesc,' / ',
     $ '   Irelmin =',Irelminc,'E-5',lambda_refc,' (',
     $ ni_refc,'3',nom(li_ref),nom(ni_ref-1),'  ',ns_refc,')'
c
c  If requested, insert the result .res in PHYAT for SYNTHESIS:
      if(iwr_phyat.ne.0) then
         Call Write_Phyat(lambda_output,PHYAT)
      endif
c
      Return
c
      End
c
c-----------------end HeI_BP--------------------
c--------------------------------------------------------
      Subroutine CII_DSK(Te,Ne,niv_input,lambda_output,PHYAT)
c
c CII extracted and adapted from Davey Storey Kisielius 2000 AA Sup 142, 85
c       s, p, d, f levels, fits t4 = 0.35 - 2.0.
c
c See also complements H-like CII_2
c
      IMPLICIT NONE
c
      CHARACTER niv_input*(*),lambda_output*(*),PHYAT*(*)
c
      CHARACTER*17 titre
      CHARACTER*9 AsB(100)
      CHARACTER numion*4,numion4*4,nomion*7
      CHARACTER*1 orbs(100),orbi(100),SPDF(4),OE(2),nom2,bidon
      CHARACTER Irelminc*5,lambda_refc*6,iraiesc*4,wlminc*8,wlmaxc*8
      INTEGER doub(3,11,4,2),nc,doubs,doubi,doudou
      INTEGER i,nmax,mp,n,Pp,Lp,JJ,ntrans,i_ref,iraies
      INTEGER mps(100),ns(100),Ls(100),Ps(100)
      INTEGER mpi(100),ni(100),Li(100),Pi(100)
      INTEGER iwlim,iw,iecr
      INTEGER imult,ifre1,ifre1_ref,ifre2,mult,nmult
      DOUBLE PRECISION cmult
      DOUBLE PRECISION int_ref,lambda_ref
      DOUBLE PRECISION vzero
      DOUBLE PRECISION hc,hc8,hcsk,Ryd_inf
      COMMON/MOCR1/wlmin,wlmax,wlim,Irelmin_H,Ionab,Ionab_min,iwlim
      DOUBLE PRECISION int_ref_mod_H
      COMMON/MOCR2/int_ref_mod_H
      DOUBLE PRECISION Ionab,Ionab_min
      DOUBLE PRECISION wlim(20),Irelmin_H(21)
      DOUBLE PRECISION wlmin,wlmax
      DOUBLE PRECISION Irelmin(21)
      COMMON/MOCR3/Irelmin
      DOUBLE PRECISION Ne,Te
      DOUBLE PRECISION t,emi_ref,emi,x,R,vb
      DOUBLE PRECISION E(3,11,4,2,2),Ei(2),Es(2),DE(4),wl(4),CO(4)
      DOUBLE PRECISION a(100),b(100),c(100),d(100),f(100),wlDSK(100)
c
      INTEGER in_niv,out_niv,iwr_phyat
      COMMON/MOC1/in_niv,out_niv,iwr_phyat
c
      COMMON/MOCC1/hc,hc8,hcsk,Ryd_inf
      COMMON/MOCX1/lambda_refc
      COMMON/MOCX2/int_ref,lambda_ref,ifre1_ref
      COMMON/MOCV1/vzero
c
      data SPDF/'S','P','D','F'/
      data OE/'e','o'/
c
      data mult/2/nmult/0/imult/0/ifre1/1/ifre2/0/ !imult=0 -> raie Ref CII
      data cmult/1.d0/  !pour write
c
c Correc atm refrac x=lambda_vac in um, dry air 15C p=760mmHg, Allen_3, p.124
      R(x)=1.+(64.328+29498.1/(146.-1./x**2)+255.4/(41.-1./x**2))*1.d-6
c
c EXPLICIT:  Ref. multiplet 3d-4f in DSK list:
      i_ref=36
c
      t=Te/1.d4
c
      print *,niv_input,'  ',lambda_output
c
      iraies=0 !line counter
c
           OPEN(unit=in_niv,file=niv_input,status='old')
c
      read(in_niv,15) titre
 15   format(A17)
 1    format(a1)
      read(in_niv,121) mult,nmult,imult,cmult,ifre1,ifre2
 121  format(6x,i4,6x,i4,6x,i4,6x,d7.0,6x,i4,6x,i4)
c Label + name ion (form X-SSN):
      read(in_niv,60) numion,numion4,nomion
 60   format(9x,a4,11x,a4,10x,a7)
      nom2=' ' ! (CII: li_ref=ni_ref-1)
c
      read(in_niv,1) bidon  !titles
      read(in_niv,1) bidon  !titles
c
      read(in_niv,20) nmax
 20   format(15i5)
c
c  Init
      Do mp=1,3  !parent mult 1 or 3
      Do n=1,nmax  !principal
      Do Pp=1,2 !Pp=P+1 w/ P=0 even, P=1 odd
       Do Lp=1,4 !Lp=L+1 w/ L= 0, 1, 2, 3 for 2S, 2P, 2D, 2F (no quartet)
        Do JJ=1,2 !JJ=1,2 for J=Jmin,Jmax=Jmin+1 (doublet)
         E(mp,n,Lp,Pp,JJ)=0.d0
        Enddo
        doub(mp,n,Lp,Pp)=1 ! doub=2,3 for doublet 2P,2D explicit
       Enddo
      Enddo
      Enddo
      Enddo
c 1-el
      mp=1      !parent multiplicity 1 for one-el levels
      Pp=1      !Pp indifferent for mp=1, 1 arbitrary
      Do n=2,5  !principal
       Do JJ=1,2 !JJ=1,2 for J=Jmin,Jmax=Jmin+1 (doublet only)
        read(in_niv,21)(E(mp,n,Lp,Pp,JJ),Lp=1,4) !Lp=L+1 Lp --> Pp
       Enddo
       Do Lp=1,4
        If(E(mp,n,Lp,Pp,2).gt.0.d0) Then
         IF(Lp.eq.2) THEN
          doub(mp,n,Lp,Pp)=2     !2P explicit
         ELSE
          doub(mp,n,Lp,Pp)=3     !2D explicit (2F excluded)
         ENDIF
        Endif
       Enddo
      Enddo
 21   format(6d12.2)
      JJ=1
      Do n=6,nmax !principal; JJ degeneracy: JJ=2 left empty
       read(in_niv,21)(E(mp,n,Lp,Pp,JJ),Lp=1,4) !Lp=L+1
      Enddo
c 2-el
      read(in_niv,1) bidon  !titles
c
      mp=3      !parent multiplicity 3 for di-el levels
      Do n=2,3  !principal
      Do Pp=1,2 !Pp=P+1 w/ P=0 even, P=1 odd
       Do Lp=1,3 !Lp=L+1
        Do JJ=1,2 !JJ=1,2 for J=Jmin,Jmax=Jmin+1 (doublet only)
         read(in_niv,222) E(mp,n,Lp,Pp,JJ)
        Enddo
       Enddo
       Do Lp=1,3
        If(E(mp,n,Lp,Pp,2).gt.0.d0) Then
         IF(Lp.eq.2) THEN
          doub(mp,n,Lp,Pp)=2     !2P explicit
         ELSE
          doub(mp,n,Lp,Pp)=3     !2D explicit (no 2F)
         ENDIF
        Endif
       Enddo
      Enddo
      Enddo
 222   format(24x,d12.2)
c
      Do i=1,3
      read(in_niv,1) bidon  !titles
      Enddo
      read(in_niv,20) ntrans ! number of transitions selected from DSK
c
      Do i=1,ntrans
      read(in_niv,25) mps(i),ns(i),orbs(i),Ls(i),Ps(i),
     $                mpi(i),ni(i),orbi(i),Li(i),Pi(i),
     $        wlDSK(i),a(i),b(i),c(i),d(i),f(i),AsB(i)
      Enddo
 25   format(1x,i1,3x,i2,a1,1x,i1,i1,
     $       2x,i1,3x,i2,a1,1x,i1,i1,f7.1,5f8.4,a9)
c
      CLOSE(in_niv)
c
c  Wavelength ranges:
c
      do iw=1,iwlim+1
         Irelmin(iw)=Irelmin_H(iw)*
     $   int_ref_mod_H/(int_ref*dmax1(Ionab,Ionab_min))
      enddo
c
          OPEN(unit=out_niv,file=lambda_output,status='old')
c
      write(out_niv,15) titre
      write(out_niv,116)
 116   format(' bbb',/,' bbb')
      write(out_niv,141) mult,nmult,imult,cmult,ifre1,ifre2
 141   format('  mult =',i4,' nmult =',i4,' imult =',i4,
     $   ' cmult =',f7.4,' ifre1 =',i4,' ifre2 =',i4)
      write(out_niv,67) Ne,Te
 67   format(1p,' Ne =',e9.2,'  Te =',e9.2, ' alfeff')
      write(out_niv,62) numion,numion4,nomion
 62   format(' numion :',a4,' numion4 :',a4,' nomion :',a7,' a4,a4,a7')
      write(out_niv,63) wlmin,wlmax
 63   format(1p,2e10.3,' wlmin, wlmax X-SSN')
      write(out_niv,85) iwlim
 85   format(' iwlim ='i4,' puis wlim(1->iwlim), Irelmin(1->iwlim+1)')
      write(out_niv,87) (wlim(iw),iw=1,iwlim)
 87   format(5x,1p,13e10.3)
      write(out_niv,86) (Irelmin(iw),iw=1,iwlim+1)
 86   format(1p,13e10.3)
c
c Reference multiplet for X-SSN output:
      i=i_ref
c Emissivity (absolute) erg.cm3.s-1
      emi_ref=1.d-14*a(i)*t**f(i)*
     $  (1.+b(i)*(1.-t)+c(i)*(1.-t)**2+d(i)*(1.-t)**3)*
     $  hc*1.d7/wlDSK(i)  !wlDSK(i) in nm
c for common:
      int_ref=emi_ref
      lambda_ref=wlDSK(i_ref)*1.d1
c Print correct si imult=0 mis en 1er dans main pour ions nmult>0
      write(lambda_refc,'(I5)') int(lambda_ref)
c
      Write(out_niv,210)
 210  Format(' Ref. for liste_modele.dat :')
      Write(out_niv,110) numion,ifre1_ref,nomion,
     $  vzero,int(wlDSK(i_ref)*10.),
     $  ni(i_ref),orbi(i_ref),nom2,ns(i_ref),orbs(i_ref)
 110  Format(1x,A4,I1,'00000000',1x,A7,'        1.0   0.     1.000e+0',
     $  '  1.                 0   1  ',f5.1,1x,I5,'++',
     $  I2,2A1,'-',I2,A1)
      Write(out_niv,211)
 211  Format(' Ref. for liste_phyat.dat :')
c  for liste_phyat.dat:
      Write(out_niv,111) numion4,ifre1_ref,nomion,emi_ref,
     $  int(wlDSK(i_ref)*10.),
     $  ni(i_ref),orbi(i_ref),nom2,ns(i_ref),orbs(i_ref),Te
 111  Format('9',A4,I1,'00000000',1x,A7,'        1.0   0.   ',
     $  1p,e10.3,'  1.               999   1   1.00 ',I5,'++',
     $  I2,2A1,'-',I2,A1,1x,'(Te=',f7.0,')')
c
      DO i=1,ntrans
c
c  Fit alfDSK t=.35-2.
c Emissivity relative to reference:
      emi=1.d-14*a(i)*t**f(i)*
     $  (1.+b(i)*(1.-t)+c(i)*(1.-t)**2+d(i)*(1.-t)**3)*
     $  hc*1.d7/wlDSK(i)/emi_ref ! ATT: some slighly wrong wlDSK (nm) changed
c
      If(mps(i).eq.1) Then
      doubs=doub(mps(i),ns(i),Ls(i)+1,1)
      Else
      doubs=doub(mps(i),ns(i),Ls(i)+1,Ps(i)+1)
      Endif
c
      If(mpi(i).eq.1) Then
      doubi=doub(mpi(i),ni(i),Li(i)+1,1)
      Else
      doubi=doub(mpi(i),ni(i),Li(i)+1,Pi(i)+1)
      Endif
c
      doudou=doubs*doubi
c
      If(mps(i).eq.1) Then
      Es(1)=E(mps(i),ns(i),Ls(i)+1,1,1)
      Else
      Es(1)=E(mps(i),ns(i),Ls(i)+1,Ps(i)+1,1)
      Endif
c
      If(mps(i).eq.1) Then
      Es(2)=E(mps(i),ns(i),Ls(i)+1,1,2)
      Else
      Es(2)=E(mps(i),ns(i),Ls(i)+1,Ps(i)+1,2)
      Endif
c
      If(mpi(i).eq.1) Then
      Ei(1)=E(mpi(i),ni(i),Li(i)+1,1,1)
      Else
      Ei(1)=E(mpi(i),ni(i),Li(i)+1,Pi(i)+1,1)
      Endif
c
      If(mpi(i).eq.1) Then
      Ei(2)=E(mpi(i),ni(i),Li(i)+1,1,2)
      Else
      Ei(2)=E(mpi(i),ni(i),Li(i)+1,Pi(i)+1,2)
      Endif
c
c  Nature of transition and number of distinct components within multiplet. 
c  Explicitly doublet levels (doub() > 1): 
c   mp=1: only levels 3p,4p,5p,3d,4d(?)
c   mp=3: all levels, except obviously 2p2 2S and 2p.3p 2S
c  Line ratios close to LS coupling (NIST): 
c     S-P and 'P only' 2/3, 1/3; 'D only' 3/5, 2/5;
c     P-D 9/15,5/15,1/15; P-P 5/9,2/9,1/9,1/9; D-D 14/25,9/25,1/25,1/25
c
      If(doudou.eq.1) Then
c 1-1            ! 1 for 2S or degenerate 2P,2D,2F
        nc=1
        CO(1)=1.d0
        DE(1)=Es(1)-Ei(1)
c
      Elseif(doudou.eq.2) Then
        nc=2
        CO(1)=2.d0/3.d0
        CO(2)=1.d0/3.d0
       if(doubs.eq.2) then
c 1-P
        DE(1)=Es(2)-Ei(1) !dominant comp first
        DE(2)=Es(1)-Ei(1)
       else
c P-1
        DE(1)=Es(1)-Ei(2) !dominant comp first
        DE(2)=Es(1)-Ei(1)
       endif
c
      Elseif(doudou.eq.3) Then
        nc=2
        CO(1)=3.d0/5.d0
        CO(2)=2.d0/5.d0
       if(doubs.eq.3) then
c 1-D
        DE(1)=Es(2)-Ei(1) !dominant comp first
        DE(2)=Es(1)-Ei(1)
       else
c D-1
        DE(1)=Es(1)-Ei(2) !dominant comp first
        DE(2)=Es(1)-Ei(1)
       endif
c
      Elseif(doudou.eq.6) Then
        nc=3
        CO(1)=9.d0/15.d0
        CO(2)=5.d0/15.d0
        CO(3)=1.d0/15.d0
       if(doubs.eq.3) then
c P-D
        DE(1)=Es(2)-Ei(2) !dominant comp first
        DE(2)=Es(1)-Ei(1)
        DE(3)=Es(1)-Ei(2)
        
       else
c D-P
        DE(1)=Es(2)-Ei(2) !dominant comp first
        DE(2)=Es(1)-Ei(1)
        DE(3)=Es(2)-Ei(1)
       endif
c
      Elseif(doudou.eq.4) Then
        nc=4
        CO(1)=5.d0/9.d0
        CO(2)=2.d0/9.d0
        CO(3)=1.d0/9.d0
        CO(4)=1.d0/9.d0
c P-P
        DE(1)=Es(2)-Ei(2) !dominant comp first
        DE(2)=Es(1)-Ei(1)
        DE(3)=Es(2)-Ei(1)
        DE(4)=Es(1)-Ei(2)
c
      Elseif(doudou.eq.9) Then
        nc=4
        CO(1)=14.d0/25.d0
        CO(2)= 9.d0/25.d0
        CO(3)= 1.d0/25.d0
        CO(4)= 1.d0/25.d0
c D-D
        DE(1)=Es(2)-Ei(2) !dominant comp first
        DE(2)=Es(1)-Ei(1)
        DE(3)=Es(2)-Ei(1)
        DE(4)=Es(1)-Ei(2)
c
      Else
c Pb!
        print *,'***CII: i,doudou=',i,doudou,' STOP***'
        STOP
      Endif
c
      DO n=1,nc
       if(DE(n).gt.0.d0) then
        wl(n)=1.d8/DE(n)
c   air refraction wl > 2000A AND wl < 20000A:
        if(wl(n).gt.2.d3.and.wl(n).lt.2.d4)THEN
         wl(n)=wl(n)/R(1.d4/DE(n))
        endif
       else
         wl(n)=1.d0
         print *,'**CII: i=',i,' nc,n =',nc,n,' ** DE(n)=',DE(n)
       endif
c lambda range:
          IF(wl(n).gt.wlmin.and.wl(n).lt.wlmax) THEN !wl range
            vb=emi*CO(n) ! int relative
c
c Limite Intensite relative admise fonction du domaine lambda
          iecr=0
         iw=1
         if(wl(n).le.wlim(iw).and.vb.ge.Irelmin(iw)) then
          iecr=1
         endif
        do iw=2,iwlim
         if(wl(n).gt.wlim(iw-1).and.wl(n).le.wlim(iw)
     $     .and.vb.ge.Irelmin(iw))then
          iecr=1
         endif
        enddo
        if(wl(n).gt.wlim(iwlim).and.vb.ge.Irelmin(iwlim+1))then
          iecr=1
        endif
             IF(iecr.eq.1) THEN
c counter des raies ecrites:
               iraies=iraies+1
c
        Write(out_niv,119) numion,ifre1,imult,ifre2,i+10,n,nomion,
     $    wl(n),vb,numion,ifre1_ref,
     $    mpi(i),ni(i),orbi(i),SPDF(Li(i)+1),OE(Pi(i)+1),
     $    mps(i),ns(i),orbs(i),SPDF(Ls(i)+1),OE(Ps(i)+1),AsB(i)
 119    Format(1x,A4,3I1,'00',I2,'0',I1,1x,A7,f13.3,   ! 3I1=100
     $    ' 0.   ',1p,e10.3,'  1.     ',A4,I1,'00000000   1   1.00 ',
     $    '(',I1,')',I1,A1,'2',A1,A1,'-','(',I1,')',I2,A1,'2',A1,A1,A9)
c
             ENDIF  !iecr
          ENDIF !wl range
c
      ENDDO   !n
c
      ENDDO   !i
c
         Write(out_niv,75) iraies
 75      format('Number of X-SSN lines: iraies =',i5)
c
           CLOSE(out_niv)
c
      write(iraiesc,'(I4)') iraies
      write(Irelminc,'(I5)') int(Irelmin(4)*1.d5)
      write(wlminc,'(I8)') int(wlmin)
      write(wlmaxc,'(I8)') int(wlmax)
      print *,'effective ni range:   2  5',
     $ '  wl range:',wlminc,wlmaxc
      print *,'number of emission lines = ',iraiesc,' / ',
     $ '   Irelmin =',Irelminc,'E-5',lambda_refc
c
c  If requested, insert the result .res in PHYAT for SYNTHESIS:
      if(iwr_phyat.ne.0) then
         Call Write_Phyat(lambda_output,PHYAT)
      endif
c
      Return
c
      End
c
c-------------------end CII_DSK--------------------------
c--------------------------------------------------------
c     Call CII_diel(Te,Ne,'CII_diel.lab','CII_diel.res',PHYAT)
c
      Subroutine CII_diel(Te,Ne,niv_input,lambda_output,PHYAT)
c
c  ***Sochi, Taha ; Storey, Peter 2013
c  ***Use the file CII_diel.dat
c
      IMPLICIT NONE
c
      CHARACTER niv_input*(*),lambda_output*(*),PHYAT*(*)
c
      DOUBLE PRECISION hc,hc8,hcsk,Ryd_inf
      DOUBLE PRECISION int_ref,lambda_ref
      DOUBLE PRECISION vzero
      COMMON/MOCR1/wlmin,wlmax,wlim,Irelmin_H,Ionab,Ionab_min,iwlim
      DOUBLE PRECISION int_ref_mod_H
      COMMON/MOCR2/int_ref_mod_H
      DOUBLE PRECISION Ionab,Ionab_min
      DOUBLE PRECISION wlim(20),Irelmin_H(21)
      DOUBLE PRECISION wlmin,wlmax
      DOUBLE PRECISION Irelmin(21)
      COMMON/MOCR3/Irelmin
      CHARACTER*17 titre
      CHARACTER Irelminc*5,lambda_refc*6,iraiesc*4,wlminc*8,wlmaxc*8
      CHARACTER numion*4,numion4*4,nomion*7
      DOUBLE PRECISION Ne,Te
      DOUBLE PRECISION cmult
      DOUBLE PRECISION vb
      INTEGER i,j,imaxl,it
      CHARACTER i_alph*4
      CHARACTER zero*1,nc1*1,nc2*2,nc3*3
      INTEGER imult,ifre1,ifre2,mult,nmult,ifre1_ref
      INTEGER iwlim,iw,iecr,iraies,inep1
c
      CHARACTER*1 text
      CHARACTER Stat(6187)*2,ExTh(6187)*2
      CHARACTER Conf1a(6187)*4,Conf1b(6187)*7,Term1(6187)*3
      CHARACTER Conf2a(6187)*4,Conf2b(6187)*7,Term2(6187)*3
      INTEGER J1(6187),J2(6187)
      DOUBLE PRECISION WlVac(6187),WlAir(6187),RC(27,6187),wl
c
      INTEGER n,ntot,inil,ir_DSK,nsuprs
      CHARACTER*130 ligne
      DOUBLE PRECISION wl_DSK(300)
c
      INTEGER in_niv,out_niv,iwr_phyat
      COMMON/MOC1/in_niv,out_niv,iwr_phyat
c
      COMMON/MOCC1/hc,hc8,hcsk,Ryd_inf
      COMMON/MOCV1/vzero
      COMMON/MOCX1/lambda_refc
      COMMON/MOCX2/int_ref,lambda_ref,ifre1_ref
c
       data imaxl/6187/
       data out_niv/11/
       data zero/'0'/
c
      print *,niv_input,'  ',lambda_output
c
        OPEN(unit=in_niv,file=niv_input,status='old')
c
      read(in_niv,15) titre
 15   format(A17)
 1    format(a1)
c Stdd di_electronics : imult = 2                 ****VOIR pr CII?***
      read(in_niv,121) mult,nmult,imult,cmult,ifre1,ifre2
 121  format(6x,i4,6x,i4,6x,i4,6x,d7.0,6x,i4,6x,i4)
      read(in_niv,60) numion,numion4,nomion
 60   format(9x,a4,11x,a4,10x,a7)
c
        CLOSE(in_niv)
c
c  Look for nearest log(Te) column without interpolation.
c   (Note: Sochi's table is for N_e=1.e10 cm-3 ***SEE EFFECT?)
c  The first Te (it=1) is 10.**2.0; steps .1dex; closest in log: 
      it=IDNINT(10.*dlog10(Te))-19
      If(it.lt.1) then
       it=1
       elseif(it.gt.27) then
       it=27
      Endif
c
        OPEN(unit=15,status='old',file='CII_diel.dat')
c
          do i=1,16 !initial text = 16 lines, incl 3 blanks (made explicit)
           read(15,*) text
          enddo
c
          do i=1,imaxl
      Read(15,1010)Stat(i),ExTh(i),Conf1a(i),Conf1b(i),Term1(i),J1(i),
     $  Conf2a(i),Conf2b(i),Term2(i),J2(i),WlVac(i),WlAir(i),
     $  (RC(j,i),j=1,27)
          enddo
 1010 format(7x,a2,4x,a2,4x,a4,a7,4x,a3,2x,i2,15x,a4,a7,4x,a3,2x,i2,3x,
     $  f11.2,4x,f11.2,21x,27e12.3)
c
        CLOSE(15)
c
c  Wavelength ranges:
c
      do iw=1,iwlim+1
         Irelmin(iw)=Irelmin_H(iw)*
     $   int_ref_mod_H/(int_ref*dmax1(Ionab,Ionab_min))
      enddo
c
         OPEN(unit=out_niv,file=lambda_output,status='old')
c
      write(out_niv,15) titre
      write(out_niv,116)
 116   format(' bbb',/,' bbb')
      write(out_niv,141) mult,nmult,imult,cmult,ifre1,ifre2
 141   format('  mult =',i4,' nmult =',i4,' imult =',i4,
     $   ' cmult =',f7.4,' ifre1 =',i4,' ifre2 =',i4)
      write(out_niv,67) Ne,Te
 67   format(1p,' Ne =',d9.2,'  Te =',d9.2)
      write(out_niv,62) numion,numion4,nomion
 62   format(' numion :',a4,' numion4 :',a4,' nomion :',a7,' a4,a4,a7')
      write(out_niv,63) wlmin,wlmax
 63   format(1p,2d10.3,' wlmin, wlmax X-SSN')
      write(out_niv,85) iwlim
 85   format(' iwlim ='i4,' puis wlim(1->iwlim), Irelmin(1->iwlim+1)')
      write(out_niv,87) (wlim(iw),iw=1,iwlim)
 87   format(5x,1p,13e10.3)
      write(out_niv,86) (Irelmin(iw),iw=1,iwlim+1)
 86   format(1p,13d10.3)
c
c  *******************************
c  Tabulation for X-SSN : CII diel
c  *******************************
c Normalisation by radiative ref line:
      Write(out_niv,106) lambda_ref,int_ref
 106  Format(' X-SSN: lambda reference effective =',
     $  1p,e12.6,' effective emissivity =',e9.3,' erg.cm3.s-1')
c
c   pour liste_phyat.dat :
 212  Format(' Sub-Ref. for liste_phyat.dat :')
      Write(out_niv,112) numion,ifre1,imult,ifre2,nomion,
     $  numion,ifre1_ref,int(lambda_ref)
 112  Format(1x,A4,3I1,'000000',1x,A7,'        1.0   0.     1.000e+0',
     $  '  1.     ',A4,I1,'00000000   1   1.00 (',I5,
     $  ' 3d - 4) pure diel')
c
c For suppression lines also in CII_DSK.res:
            OPEN(unit=21,file='CII_DSK.res',status='old')
c   comptage des lignes:
      ntot=0
      Do
        read(21,*,end=1000) ligne
        ntot=ntot+1
      Enddo
 1000 continue
c  Nbre total de lignes sur CII_DSK.res : ntot
c
      BACKSPACE(21)  ! recule d'une ligne 2 fois
      BACKSPACE(21)  ! ... pour derniere ligne !?
      read(21,195) ir_DSK !recup nbre de raies
 195   format(31X,I5)
c
      REWIND(21)
c
      inil=ntot-ir_DSK !ici : inil 1ere raie apres raie reference
        Do n=1,inil-1
         read(21,*) ligne
        Enddo
        Do n=1,ntot-inil
       read(21,390) wl_DSK(n)
        Enddo
 390    format(23X,f12.3)
c
               iraies=0
               nsuprs=0 !counter raies suppr car deja dans CII_DSK.res
        DO i=1,imaxl
          vb=RC(it,i)*1.e+6*hc8/WlVac(i)/int_ref
              iecr=0
             IF(WlVac(i).gt.wlmin.and.WlVac(i).lt.wlmax) THEN
c Adopted minimum rel emiss depends on lambda range:
         iw=1
         if(WlVac(i).le.wlim(iw).and.vb.ge.Irelmin(iw)) then
              iecr=1
         endif
        do iw=2,iwlim
         if(WlVac(i).gt.wlim(iw-1).and.WlVac(i).le.wlim(iw)
     $     .and.vb.ge.Irelmin(iw))then
              iecr=1
         endif
        enddo
        if(WlVac(i).gt.wlim(iwlim).and.vb.ge.Irelmin(iwlim+1))then
          iecr=1
        endif
             ENDIF
             If(iecr.eq.1) Then !check for lines not in CII_DSK.res
          Do n=1,ntot-inil
           if(WlVac(i).le.2.d3) then
         if(dabs(wl_DSK(n)-WlVac(i)).lt.0.03d0) iecr=0 !tolerance 0.03A
           else
         if(dabs(wl_DSK(n)-WlAir(i)).lt.0.03d0) iecr=0
           endif
          Enddo
             if(iecr.eq.0) nsuprs=nsuprs+1
             Endif
c
             If(iecr.eq.1) Then
c
c emission line counter: 
               iraies=iraies+1
c
        if(i.lt.10) then
        write(nc1,'(I1)') i
        i_alph=zero//zero//zero//nc1
         elseif(i.lt.100) then
        write(nc2,'(I2)') i
       i_alph=zero//zero//nc2
         elseif(i.lt.1000) then
        write(nc3,'(I3)') i
       i_alph=zero//nc3
          else   !(assume i < 10000)
        write(i_alph,'(I4)') i
        endif
c
       if(WlVac(i).le.2.d3)then !**ATT: REVOIR pour lambda > 2.d4 ? **
         wl=WlVac(i)
       else
         wl=WlAir(i)
       endif
        Write(out_niv,190) numion,ifre1,imult,ifre2,i_alph,       !ns_alph,ni_alph,
     $      nomion,wl,vb,numion,ifre1,imult,ifre2,
     $      Stat(i),ExTh(i),Conf1a(i),Conf1b(i),Term1(i),J1(i),
     $      Conf2a(i),Conf2b(i),Term2(i),J2(i)
 190    Format(1x,A4,3I1,'0',A4,'0',1x,A7,f13.3,' 0.   ',1p,E10.3,
     $   '  1.     ',A4,3I1,'000000   1   1.00 ',
     $   a2,a2,x,a4,a7,a3,i2,'/2',x,a4,a7,a3,i2,'/2')
             Endif
c
        ENDDO                     !i
c
            CLOSE(21) !CII_DSK.res
c
         Write(out_niv,75) iraies
 75      format('Number of X-SSN lines: iraies =',i5)
c  (Control)
cc         write(out_niv,90) inot
cc         write(out_niv,90) (i,wlnot(i),i=1,inot)
cc 90      format(6(i5,f11.2))
c
      write(lambda_refc,'(I5)') int(lambda_ref)
      write(iraiesc,'(I4)') iraies
      write(Irelminc,'(I5)') int(Irelmin(4)*1.d5)
      write(wlminc,'(I8)') int(wlmin)
      write(wlmaxc,'(I8)') int(wlmax)
      print *,'effective ni range: 2    ',
     $ '  wl range:',wlminc,wlmaxc
      print *,'number of emission lines = ',iraiesc,' / ',
     $ '   Irelmin =',Irelminc,'E-5',lambda_refc
      print *,' ***',nsuprs,' suppressed lines, already in CII_DSK.res'
c
           CLOSE(out_niv)
c
c  If requested, insert the result .res in PHYAT for SYNTHESIS:
      if(iwr_phyat.ne.0) then
         Call Write_Phyat(lambda_output,PHYAT)
      endif
c
      return
c
      end
c
c---------------------------end CII_diel------------------
c--------------------------------------------------------
c  Call CIII_5ga('CIII_5ga.lab','CIII_5ga.res',PHYAT)
c
      Subroutine CIII_5ga(Te,Ne,niv_input,lambda_output,PHYAT)
c
c Based on data provided by P.J. Storey, 16 Sept 1998, private communication.
c data read in CIII_5ga.lab
c
c   UPPER LEVELS 5g ARE AUTOIONIZING.
c   Eff rec coefs alf(k,l) are computed by means of Saha equation, 
c   then emissivities relint(k,l) of lines 4f-5g, are obtained relative 
c   to emissivities absref(k) (in erg.cm3.s-1) of sub-reference lines. 
c   
c   Levels 5g (1--12 equivalent to 13--24 in PJS) : 
c
      IMPLICIT NONE
c
      CHARACTER niv_input*(*),lambda_output*(*),PHYAT*(*)
c
      DOUBLE PRECISION int_ref,lambda_ref,absref_new
      DOUBLE PRECISION vzero
	Double Precision Te,Ne
        DOUBLE PRECISION widthmax,ECIII,pstatplus,cons
	Double Precision alf_ref
c  (dim 12 correspond to kmax)
	Double Precision pstat(12)
	Double Precision relint(12,9),absref(12),relref(12)
	Double Precision vred(12,9)
	Double Precision wvlgth(12,9),width(12,9),beca(12,9),dv(12,9)
	Double Precision ener(12,9)
      COMMON/MOCR1/wlmin,wlmax,wlim,Irelmin_H,Ionab,Ionab_min,iwlim
      DOUBLE PRECISION int_ref_mod_H
      COMMON/MOCR2/int_ref_mod_H
      DOUBLE PRECISION Ionab,Ionab_min
      DOUBLE PRECISION wlim(20),Irelmin_H(21)
      DOUBLE PRECISION wlmin,wlmax
      DOUBLE PRECISION Irelmin(21)
      COMMON/MOCR3/Irelmin
      DOUBLE PRECISION hc,hc8,hcsk,Ryd_inf
      DOUBLE PRECISION csaha,kboltz
      INTEGER imult,ifre1,ifre1_ref,ifre2,mult,nmult
      DOUBLE PRECISION cmult
      INTEGER iwlim,iw,iecr
      INTEGER iv
        Integer nlmax(12),lref(12),i_tl(12,9)
        Integer kmax,kref,k,l,iwr,iraies
        Integer i_touteslarg
      CHARACTER*17 titre
      CHARACTER numion*4,numion4*4,nomion*7
	Character bidon*1
	Character*7 wltheo(12,9)
      CHARACTER Irelminc*5,lambda_refc*6,iraiesc*4,wlminc*8,wlmaxc*8
c
      INTEGER in_niv,out_niv,iwr_phyat
      COMMON/MOC1/in_niv,out_niv,iwr_phyat
c
      COMMON/MOCC1/hc,hc8,hcsk,Ryd_inf
      COMMON/MOCC2/csaha,kboltz
      COMMON/MOCV1/vzero
      COMMON/MOCX1/lambda_refc
      COMMON/MOCX2/int_ref,lambda_ref,ifre1_ref
c
c  ECIII(cm-1) from Moore or NIST : 
        Data ECIII/386241.d0/   ! +/-2 Moore
c  Nber of upper levels:
        Data kmax/12/
c  No of the level of the 'general' ref line:
        Data kref/9/
c  statistical weight pstat=2J+1 of upper levels (5g) : 
	Data pstat/9.d0,7.d0,9.d0,11.d0,11.d0,9.d0,11.d0,
     $            13.d0,9.d0,7.d0,5.d0,7.d0/
c  pstatplus= statistical weight GS of next ion: 
	Data pstatplus/1.d0/
c  Nber of lines from each level:
	Data nlmax/5,7,5,3,4,5,3,1,3,7,7,7/
c  No of the sub-reference line of each level:
	Data lref/1,2,2,1,2,3,3,1,3,6,6,6/
c
        print *,niv_input,'  ',lambda_output
        print *, 'kmax =', kmax
c
      iraies=0 !Line counter
c
        OPEN(unit=in_niv,file=niv_input,status='old')
c
      read(in_niv,15) titre
 15   format(A17)
 1      Format(a1)
      read(in_niv,121) mult,nmult,imult,cmult,ifre1,ifre2
 121  format(6x,i4,6x,i4,6x,i4,6x,d7.0,6x,i4,6x,i4)
      read(in_niv,60) numion,numion4,nomion
 60   format(9x,a4,11x,a4,10x,a7)
c Maximum natural width:
      read(in_niv,69) widthmax
 69   format(d8.0)
c
        Do k=1,7
         Read(in_niv,1) bidon  !titles
        Enddo
c Ex:
c  13  1  1414  1303  4361.97    31.67  0  0.999988  2.7380E+08  411080.0
        Do k=1,kmax
         Do l=1,nlmax(k)
          Read(in_niv,32) wvlgth(k,l),width(k,l),i_tl(k,l),
     $    beca(k,l),dv(k,l),ener(k,l)
         Enddo
        Enddo
 32     Format(19x,2F9.2,I3,2x,F8.6,2x,D10.4,F10.1)
c
        CLOSE(in_niv)
c
c  Init
        do k=1,kmax
          do l=1,9
            relint(k,l)=0.d0
          enddo
        enddo
c
        Do k=1,kmax
         Do l=1,nlmax(k)
           wltheo(k,l)='       '
          If(i_tl(k,l).eq.0) Then
           wltheo(k,l)='wl theo'
          Endif
         Enddo
        Enddo
c 
        cons=csaha/Te**1.5/(2.*pstatplus)
c (* wvlgth not corrected back for air refraction)
        do k=1,kmax
          do l=1,nlmax(k)
            relint(k,l)=cons*pstat(k)*dv(k,l)*beca(k,l)*
     $      dexp(-hcsk*(ener(k,l)-ECIII)/Te)/wvlgth(k,l)
          enddo
        enddo
c
        do k=1,kmax
          absref(k)=relint(k,lref(k))
          do l=1,nlmax(kmax)
            relint(k,l)=relint(k,l)/absref(k)
          enddo
        enddo
c
        do k=1,kmax
          do l=1,nlmax(kmax)
c  Approximate reduced widths of lines
c  ('width' is for 'gamma' = FWHM lorentz):
            vred(k,l)=width(k,l)/wvlgth(k,l)*3.d5/2./vzero
            vred(k,l)=dsqrt(vred(k,l)**2+1.)
          enddo
        enddo
c
c  Absolute emissivities of sub-ref lines in erg.cm3.s-1:
        do k=1,kmax
           absref(k)=absref(k)*hc8
        enddo
c  Sub-ref relative to ref:
        do k=1,kmax
           relref(k)=absref(k)/absref(kref)
        enddo
c  Recomb coeff of ref line (cm3/s) :
        k=kref
        l=lref(kref)
        alf_ref=cons*pstat(k)*dv(k,l)*beca(k,l)*
     $  dexp(-hcsk*(ener(k,l)-ECIII)/Te)
c
c  Wavelength ranges:
c
      do iw=1,iwlim+1
         Irelmin(iw)=Irelmin_H(iw)*
     $   int_ref_mod_H/(int_ref*dmax1(Ionab,Ionab_min))
      enddo
c
          OPEN(unit=out_niv,file=lambda_output,status='old')
c
      write(out_niv,15) titre
      write(out_niv,116)
 116   format(' bbb',/,' bbb')
      write(out_niv,141) mult,nmult,imult,cmult,ifre1,ifre2
 141   format('  mult =',i4,' nmult =',i4,' imult =',i4,
     $   ' cmult =',f7.4,' ifre1 =',i4,' ifre2 =',i4)
      write(out_niv,67) Ne,Te
 67   format(1p,' Ne =',e9.2,' Te =',e9.2)
      write(out_niv,62) numion,numion4,nomion
 62   format(' numion :',a4,' numion4 :',a4,' nomion :',a7,' a4,a4,a7')
      write(out_niv,169) widthmax
 169  format(' MAXIMUM NATURAL WIDTH widthmax =',e9.2)
      write(out_niv,63) wlmin,wlmax
 63   format(1p,2e10.3,' wlmin, wlmax X-SSN')
      write(out_niv,85) iwlim
 85   format(' iwlim ='i4,' puis wlim(1->iwlim), Irelmin(1->iwlim+1)')
      write(out_niv,87) (wlim(iw),iw=1,iwlim)
 87   format(5x,1p,13e10.3)
      write(out_niv,86) (Irelmin(iw),iw=1,iwlim+1)
 86   format(1p,13e10.3)
	 Write(out_niv,3)
 3     Format(6x,'Di-electronic recombination lines of C III 4f''-5g''')
	 Write(out_niv,7) Te
 7	Format(1p,'   Te=',e9.3)
	 Write(out_niv,5) 
 5      Format('    Walelengths, then relative emissivities:')
	   Do k=1,kmax
	 Write(out_niv,8) k,wvlgth(k,lref(k)),absref(k)
 8      Format(1x,'k=',i2,'   ref line :',f10.3, 
     $     ' emissivity =',1p,d10.3,' erg.cm3.s-1')
	 Write(out_niv,6) (wvlgth(k,l),l=1,nlmax(k))
	 Write(out_niv,2) (relint(k,l),l=1,nlmax(k))
	   Enddo
 2      Format(1p,6e12.3)
 6	Format(6f12.3)
         Write(out_niv,4)
 4      Format('   Relative intensites of sub-ref lines:')
 	 Write(out_niv,6) (wvlgth(k,lref(k)),k=1,kmax)
 	 Write(out_niv,2) (relref(k),k=1,kmax)
         Write(out_niv,9) wvlgth(kref,lref(kref)),alf_ref,absref(kref)
 9      Format('   Rec coef of ',f8.3,' = ',1p,e10.4,' cm3/s',
     $    '  emissivity =',e9.3,' erg.cm3.s-1')
c
c  **********************
c  Tabulation for X-SSN :
c  **********************
c  Use general ref line of CIII: 
      Write(out_niv,106) lambda_ref,int_ref
 106  Format(' X-SSN: effective lambda reference =',
     $  1p,e12.6,' effective emissivity =',e9.3,' erg.cm3.s-1')
c
c  'general' ref of CIII 4f-5g autoio normalized by CIII reference:
      absref_new=absref(kref)/int_ref
c
c   for liste_phyat.dat :
      Write(out_niv,212)
 212  Format(' Sub-Ref. for liste_phyat.dat :')
         Write(out_niv,1112) numion,ifre1,imult,ifre2,nomion,
     $      absref_new,numion,ifre1_ref,int(wvlgth(kref,lref(kref)))
 1112    Format(1x,A4,3I1,'000000',1x,A7,'        1.0   0.   ', ! 3I1=280
     $      1p,e10.3,'  1.     ',A4,I1,'00000000   1   1.00 ',I5,'++')
c
         iw=1
         if(4.4d3.le.wlim(iw)) iv=iw
        do iw=2,iwlim
         if(4.4d3.gt.wlim(iw-1).and.4.4d3.le.wlim(iw)) iv=iw
        enddo
        if(4.4d3.gt.wlim(iwlim)) iv=iwlim
c
       Do k=1,kmax
         iwr=0
         do l=1,nlmax(k)
          if(relint(k,l).gt.Irelmin(iv).and.width(k,l).lt.widthmax)then
           iwr=iwr+1
          endif
         enddo
             IF(iwr.gt.0) THEN
         If(k.lt.10) then
c Ref lev :
               iraies=iraies+1
           Write(out_niv,110) numion,ifre1,imult,ifre2,k,nomion,
     $      relref(k),numion,ifre1,imult,ifre2,int(wvlgth(k,lref(k)))
 110    Format(1x,A4,3I1,'000',I1,'00',1x,A7,'        1.0   0.   ',   !k<10
     $      1p,e10.3,'  1.     ',A4,3I1,'000000   1   1.00 ',I5,'+')
           Do l=1,nlmax(k)
          if(relint(k,l).gt.Irelmin(iv).and.width(k,l).lt.widthmax)then
cc   (selected) Reduced width(k,l) for liste_phyat
c line counter:
             iraies=iraies+1
             Write(out_niv,40) numion,ifre1,imult,ifre2,k,l,nomion,
     $         wvlgth(k,l),relint(k,l),numion,ifre1,imult,ifre2,k,
     $         vred(k,l),wltheo(k,l)
 40     Format(1x,A4,3I1,'000',I1,'0',I1,1x,A7,f13.3,' 0.   ',
     $    1p,e10.3,'  1.     ',A4,3I1,'000',I1,'00   1 ',0p,f6.2,' ',a7)
            endif
           Enddo
c
         Else   !k>10
               iraies=iraies+1
           Write(out_niv,111) numion,ifre1,imult,ifre2,k,nomion,
     $      relref(k),numion,ifre1,imult,ifre2,int(wvlgth(k,lref(k)))
 111    Format(1x,A4,3I1,'00',I2,'00',1x,A7,'        1.0   0.   ',   !k>10
     $      1p,e10.3,'  1.     ',A4,3I1,'000000   1   1.00 ',I5,'+')
           Do l=1,nlmax(k)
          if(relint(k,l).gt.Irelmin(iv).and.width(k,l).lt.widthmax)then
               iraies=iraies+1
             Write(out_niv,41) numion,ifre1,imult,ifre2,k,l,nomion,
     $         wvlgth(k,l),relint(k,l),numion,ifre1,imult,ifre2,k,
     $         vred(k,l),wltheo(k,l)
 41     Format(1x,A4,3I1,'00',I2,'0',I1,1x,A7,f13.3,' 0.   ',
     $    1p,e10.3,'  1.     ',A4,3I1,'00',I2,'00   1 ',0p,f6.2,' ',a7)
            endif
           Enddo
c
         Endif   !
             ENDIF     !iwr>0
c
        Enddo    !k=1,kmax
c
         Write(out_niv,75) iraies
 75      format('Number of X-SSN lines: iraies =',i5)
c
c   EXPLICIT, STDD: i_touteslarg = 0
         i_touteslarg=0
c         i_touteslarg=1 --> ATT: Pb insertion in liste_phyat
        IF(i_touteslarg.ne.0) THEN
           Write(out_niv,225)
 225       Format(//,3x,'All widths:')
         Do k=1,kmax
          If(k.lt.10) then
           Do l=1,nlmax(k)
             Write(out_niv,40) numion,ifre1,imult,ifre2,k,l,nomion,
     $         wvlgth(k,l),relint(k,l),numion,k,vred(k,l),wltheo(k,l)
           Enddo
          Else
           Do l=1,nlmax(k)
             Write(out_niv,41) numion,ifre1,imult,ifre2,k,l,nomion,
     $         wvlgth(k,l),relint(k,l),numion,k,vred(k,l),wltheo(k,l)
           Enddo
          Endif
         Enddo    !k=1,kmax
        ENDIF
c
      write(iraiesc,'(I4)') iraies
      write(Irelminc,'(I5)') int(Irelmin(4)*1.d5)
      write(wlminc,'(I8)') int(wlmin)
      write(wlmaxc,'(I8)') int(wlmax)
      print *,'effective ni range:   4f  4f',
     $ '  wl range:',wlminc,wlmaxc
      print *,'number of emission lines = ',iraiesc,' / ',
     $ '   Irelmin =',Irelminc,'E-5',lambda_refc,' *',absref_new
c
	   CLOSE(out_niv)
c
c  If requested, insert the result .res in PHYAT for SYNTHESIS:
      if(iwr_phyat.ne.0) then
         Call Write_Phyat(lambda_output,PHYAT)
      endif
c
        Return
c
        End
c---------------------------end CIII_5ga-----------
c-------------------------------------------------------
c
      Subroutine NII_Sto(Te,Ne,niv_input,lambda_output,PHYAT)
c
c  Call NII_Sto('NII_13.lab','NII_13.res',PHYAT)
c       Adapted from Fang, Storey, Liu 2011, A&A 530, A18
c     *** PB in paper *** Erratum 2013, A&A , A
      IMPLICIT NONE
c
      CHARACTER niv_input*(*),lambda_output*(*),PHYAT*(*)
c
      DOUBLE PRECISION a0,a1,a2,a3,a4,a5,a6,a7,a8,alflow,alfit
cc      DOUBLE PRECISION b0,b1,b2,b3,b4,b5,b6,alfhig
      DOUBLE PRECISION vzero
      DOUBLE PRECISION norm,va,vb,vc,vd
      DOUBLE PRECISION lambda(218),ERC(4,218,7),T_lin
      DOUBLE PRECISION A(4,218),lp(4,218),parl(4,218),parh(4,218),pl,ph
      DOUBLE PRECISION ac(4,55,6),a0_cor(4,55),del1(4,55),lambdf(55)
      DOUBLE PRECISION bc(4,55,7),del2(4,55)
      DOUBLE PRECISION cc(4,55,8),del3(4,55),del4(4,55)
      DOUBLE PRECISION dd(4,55,8),del5(4,55),del6(4,55)
      DOUBLE PRECISION rap(7),ERCf(7),TeTab(7),t4,dlt4
      DOUBLE PRECISION moyrap(7),moyrap2(7),moyfit(7),moyfit2(7)
      DOUBLE PRECISION flog
      DOUBLE PRECISION alf(218),ERC1,ERC0,lambda_ref,int_ref
      INTEGER imult,ifre1,ifre2,mult_wr,nmult
      DOUBLE PRECISION cmult
      INTEGER ine,irai,ite,item5,j,I_detail,iraies,ine_synth,I_2011
      INTEGER irai1(55),if_rai(218),indrai,irai_ref
      INTEGER nu(218),lu(218),nl,ll,k
      INTEGER iwlim,iw,iecr
      COMMON/MOCR1/wlmin,wlmax,wlim,Irelmin_H,Ionab,Ionab_min,iwlim
      DOUBLE PRECISION int_ref_mod_H
      COMMON/MOCR2/int_ref_mod_H
      DOUBLE PRECISION Ionab,Ionab_min
      DOUBLE PRECISION wlim(20),Irelmin_H(21)
      DOUBLE PRECISION wlmin,wlmax
      DOUBLE PRECISION Irelmin(21)
      COMMON/MOCR3/Irelmin
      CHARACTER*17 titreb
      CHARACTER Irelminc*5,lambda_refc*6,iraiesc*4,wlminc*8,wlmaxc*8
      CHARACTER numion*4,numion4*4,nomion*7,bidon*1
      CHARACTER Tr*23,Trans*32,n_lambda*1,luC*1  !2011: Tr*22,Trans*33
      CHARACTER Tranf*26,mult*4
      CHARACTER titre*80,titr1*80,titr2*80,titr3*80,titr4*80,titra0*80
      DIMENSION Tr(218),Trans(218),n_lambda(218),luC(218)
      DIMENSION Tranf(55),mult(55)
      DIMENSION titre(4),titr1(4),titr2(4),titr3(4),titr4(4)
c ***For readsto, Hli, HBD (tabul compa Hli) :
      INTEGER numin,numax !Br table
      REAL az,am  !HBD
      DOUBLE PRECISION hc,hc8,hcsk,Ryd_inf
      DOUBLE PRECISION Br
      DOUBLE PRECISION Mat
      DOUBLE PRECISION Int_Hlike
      INTEGER inep1,iz
      DOUBLE PRECISION Ne,Te,lgNe,dlgNe,t,lgt
      CHARACTER zc*1,tabsto*9
c
      INTEGER in_niv,out_niv,iwr_phyat
      COMMON/MOC1/in_niv,out_niv,iwr_phyat
c
      COMMON/MOCC1/hc,hc8,hcsk,Ryd_inf
      COMMON/MOCV1/vzero
      COMMON/MOCH2/az,am,numin,numax
      COMMON/MOCS2/inep1,iz,zc,tabsto
c
c Tables Fang et al:
      Data TeTab/125.d0,500.d0,1000.d0,5000.d0,10000.d0,
     $  15000.d0,20000.d0/
c
c 218 raies Tables 3-6 alf eff bruts : Ne=2 3 4 5,(Te=125 500 1000 5000 10000 15000 20000)
c   *** beg 2011: obsolete
c 55 lines Fit Fang alpha eff : Te < 1e4, tables 7 9 11 13 polynomial lgNe=2 3 4 5
c   log({alpha}_eff_)+15 = a0+a1*t+a2*t^2^+a3*t^3^+a4t^4^+a5*t^5,
c   where t is the log_10_ of electronic temperature, log(Te) (K).
c  t --> dlog_10(Te) (K).
ccc      alflow(t,a0,a1,a2,a3,a4,a5)=1.d-15*10.**
ccc     $  (a0+t*(a1+t*(a2+t*(a3+t*(a4+t*a5)))))
c
c 55 raies Fit Fang alpha eff : Te > 1e4, tables 8 10 12 14 'expon" lgNe=2 3 4 5
c   log({alpha}_eff_)+15 = (b0+b1*t+b2*t^2^+b3*t^3^+b4*t^4^)*t^b5^*exp(b6*t),
c   where t is the reduced electronic temperature Te/10^4 (K).
c  t --> t4=Te/1.d4
ccc      alfhig(t,b0,b1,b2,b3,b4,b5,b6)=1.d-15*10.**
ccc     $  ((b0+t*(b1+t*(b2+t*(b3+t*b4))))*t**b5*dexp(b6*t))
c   *** end 2011 obsolete
c
c ***2013: 55 raies t --> t4=Te/1.d4 (Case A or B=stdd)
      alfit(t,a1,a2,a3,a4,a5,a6,a7,a8)=1.d-15*10**(a1+t*(a2+t*a3)+
     $  (a4+t*(a5+t*a6))*dlog10(t)+a7*(dlog10(t))**2+a8/t)
c
      print *,niv_input,'  ',lambda_output
c
	   OPEN(unit=in_niv,file=niv_input,status='old')
c
      read(in_niv,15) titreb
 15   format(A17)
      read(in_niv,121) mult_wr,nmult,imult,cmult,ifre1,ifre2
 121  format(6x,i4,6x,i4,6x,i4,6x,d7.0,6x,i4,6x,i4)
c
c Id + name ion (form X-SSN):
      read(in_niv,60) numion,numion4,nomion
 60   format(9x,a4,11x,a4,10x,a7)
c I_detail=0 : stdd = only X-SSN output, .ne.0 --> more tables
c I_2011=0 : stdd = use new tables 2013 (see erratum Fang et al),
c            otherwise, ***ATT: change data set before activating .ne.0
      Read(in_niv,23) I_detail,I_2011
 23   Format(20I4)
c
c  Coeffs for alternative power law fits, D. Pequignot, 4 param:
      Read(in_niv,144) pl,ph
 144  Format(2F5.1)
c
c  alf eff raw Cas B Fang Tables 3-6 : 2011 ou 2013 (avc chgmt formats 11 et 311)
c
      Do ine=1,4
        read(in_niv,10) titre(ine)
c (Rem : n_lambda = blank or asterisk, read but not used)
         IF(ine.eq.1) then
       Do irai=1,218
        Read(in_niv,311) nu(irai),luC(irai),Trans(irai),
     $      lambda(irai),n_lambda(irai),(ERC(ine,irai,ite),ite=1,7)
       Enddo
         ELSE
       Do irai=1,218
        Read(in_niv,11) Tr(irai),Trans(irai),
     $      lambda(irai),n_lambda(irai),(ERC(ine,irai,ite),ite=1,7)
       Enddo
         ENDIF
      Enddo
 10   Format(A80)
c2011 11   Format(A22,A33,D9.2,1X,A1,1X,7D8.3)
c2011 311  Format(I1,A1,20X,A33,D9.2,1X,A1,1X,7D8.3)
c 6/13 new format for 2013 data:
 11   Format(A23,A32,D9.2,1X,A1,7D8.2)   !stdd with I_2011=0
 311  Format(I1,A1,21X,A32,D9.2,1X,A1,7D8.2)   !stdd with I_2011=0
c
c (Rem: lu(irai) defined but not used for now)
       Do irai=1,218
        if(luC(irai).eq.'s') lu(irai)=0
        if(luC(irai).eq.'p') lu(irai)=1
        if(luC(irai).eq.'d') lu(irai)=2
        if(luC(irai).eq.'f') lu(irai)=3
        if(luC(irai).eq.'g') lu(irai)=4
        if(luC(irai).eq.'h') lu(irai)=5
        if(luC(irai).eq.'i') lu(irai)=6
        if(luC(irai).eq.'k') lu(irai)=7
        if(luC(irai).eq.'l') lu(irai)=8
       Enddo
c
      If(I_2011.ne.0) Then !***ATT: change data set before activating
c
c****************reading OLD TABLES (2011 fit) begin **********
c
c  Coeff fits Fang 55 lines Tables 7 9 11 13 : low
c
      Do ine=1,4
        read(in_niv,10) titr1(ine)
       Do k=1,55
        Read(in_niv,12) Tranf(k),lambdf(k),
     $      (ac(ine,k,j),j=1,6),del1(ine,k)
       Enddo
      Enddo
 12   Format(A26,D8.2,3D9.4,3D8.4,D6.3)
c
c  Coeff fits Fang 55 lines Tables 8 10 12 14 : high
c
      Do ine=1,4
        read(in_niv,10) titr2(ine)
       Do k=1,55
        Read(in_niv,13) Tranf(k),lambdf(k),
     $      (bc(ine,k,j),j=1,7),del2(ine,k)
       Enddo
      Enddo
 13   Format(A26,D8.2,4D9.4,3D8.4,D6.3)
c
c  Corr Tables 7-13 proposed by D. Pequignot
c   to ensure continuity of Fang fit at 1.d4K
c   (small departure for low Te's)
        read(in_niv,10) titra0
      Do ine=1,4
       Read(in_niv,72) (a0_cor(ine,k),k=1,55)
      Enddo
 72   Format(8f10.5)
      Do irai=1,218
       Do ine=1,4
        Read(in_niv,143) A(ine,irai),lp(ine,irai),
     $    parl(ine,irai),parh(ine,irai)
       Enddo
      Enddo
 143   Format(14X,D9.2,3F9.5)
c
c****************reading OLD TABLES (2011 fit) end **********
c
      Else ! I_2011=0 (stdd: 2013)

c****************reading NEW TABLES (2013 fit) begin **********
c
c     'Accurate new fit coeffs for 55 optical lines : CASE B
c   Tables 7 8 9 10 for log(ne)= 2 3 4 5 (mult=No multiplet)
c   (del3, del4 = av. and max. deviation in %)
c
      Do ine=1,4
        read(in_niv,10) titr3(ine)
       Do k=1,55
        Read(in_niv,122) Tranf(k),mult(k),lambdf(k),
     $      (cc(ine,k,j),j=1,8),del3(ine,k),del4(ine,k)
       Enddo
      Enddo
 122  Format(A26,3X,A4,D8.2,7D9.4,D11.6,2D7.3)
c
c     'Accurate new fit coeffs for 55 optical lines : CASE A
c   Tables 11 12 13 14 for log(ne)= 2 3 4 5 (mult=No multiplet)
c   (del5, del6 = av. and max. deviation in %)
c
      Do ine=1,4
        read(in_niv,10) titr4(ine)
       Do k=1,55
        Read(in_niv,122) Tranf(k),mult(k),lambdf(k),
     $      (dd(ine,k,j),j=1,8),del5(ine,k),del6(ine,k)
       Enddo
      Enddo
c
      Endif                     ! I_2011
c
c****************reading NEW TABLES (2013 fit) end **********
c
	   CLOSE(in_niv)
c
	   OPEN(unit=out_niv,file=lambda_output,status='old')
c
      write(out_niv,15) titreb
      write(out_niv,116)
 116   format(' bbb',/,' bbb')
      write(out_niv,141) mult_wr,nmult,imult,cmult,ifre1,ifre2
 141   format('  mult =',i4,' nmult =',i4,' imult =',i4,
     $   ' cmult =',f7.4,' ifre1 =',i4,' ifre2 =',i4)
      write(out_niv,67) Ne,Te
 67   format(1p,' Ne =',e9.2,' Te =',e9.2)
      write(out_niv,62) numion,numion4,nomion
 62   format(' numion :',a4,' numion4 :',a4,' nomion :',a7,' a4,a4,a7')
      write(out_niv,88) I_detail,I_2011
 88   format(2i4,'  I_detail, I_2011')
c
c **2011 Intro correction a0 for continuite at 1.e4K : obsolete
      Do ine=1,4
       do k=1,55
        ac(ine,k,1)=ac(ine,k,1)+a0_cor(ine,k)
       enddo
      Enddo
c
c  I_detail=0: only X-SSN output
c
c  Correspondance between num Table fit (55 raies) and ERC (218 lines)
c   to possibly choose the fit type:
        IF(I_detail.ne.0) THEN
       Write(out_niv,34)
 34   Format(' k, irai1, lambdf, lambda')
        ENDIF
      Do k=1,55
       indrai=0
       do irai=1,218
        if((lambdf(k).gt.lambda(irai)-0.05).and.
     $     (lambdf(k).lt.lambda(irai)+0.05)) then
         irai1(k)=irai
         indrai=indrai+1
        endif
       enddo
       if(indrai.eq.0) then
        print *,lambdf(k),' no corresp, **TO BE CORRECTED****'
        irai1(k)=300 ! arb > 218 PROVISIONAL, TO BE CORRECTED
       endif
       if(indrai.gt.1) then
        print *,lambdf(k),' ',indrai,' corresp, **TO BE CORRECTED****'
       endif
        IF(I_detail.ne.0) THEN
       Write(out_niv,33) k,irai1(k),lambdf(k),
     $     lambda(irai1(k))
        ENDIF
      Enddo
 33   Format(2I4,2F9.2)
        IF(I_detail.ne.0) THEN
       Write(out_niv,233)
 233   Format(' label in fit table then corresp label in ERC tab')
       Write(out_niv,133) (k,k=1,55)
       Write(out_niv,133) (irai1(k),k=1,55)
 133   Format(30I4)
        ENDIF
c Spotting lines w/ Fang fit
      Do irai=1,218
       if_rai(irai)=0 !0 not in fit table of Fang; fit 'power law w/ corr'
      Enddo
      Do k=1,55
       if_rai(irai1(k))=k ! > 0 Fang fit; k = line label in Fang fit table
      Enddo
c
        IF(I_detail.ne.0) THEN
        Write(out_niv,234)
234   Format('label in ERC tab then corresp label in fit tab')
       Write(out_niv,133) (irai,irai=1,218)
       Write(out_niv,133) (if_rai(irai),irai=1,218)

c  Check Fang fit
       Write(out_niv,176)
 176    Format(' *** Check Fang 2013 fit : ')
      Do ine=1,4
       Write(out_niv,29) ine,(TeTab(ite),ite=1,7)
        Do k=1,55
         irai=irai1(k)
         Do ite=1,7
            t4=TeTab(ite)*1.d-4
            ERCf(ite)=alfit(t4,cc(ine,k,1),cc(ine,k,2),
     $         cc(ine,k,3),cc(ine,k,4),cc(ine,k,5),
     $         cc(ine,k,6),cc(ine,k,7),cc(ine,k,8))
         Enddo
         Do ite=1,7
          rap(ite)=ERC(ine,irai,ite)/(1.d15*ERCf(ite))-1.d0
         Enddo
         Write(out_niv,28) k,lambda(irai),(rap(ite),ite=1,7)
        Enddo
      Enddo
c
c  Fit power law w/ correc low & high Te
      Write(out_niv,77) !cc
 77   Format('*** COMPARING alf Tab to power law without and with corr') !cc
 167   Format(I8,26X,7F10.0)
        ENDIF
c
      flog=dlog10(2.d0)
       vb=dlog10(TeTab(5)/TeTab(3))
c
      Do ine=1,4
c
        IF(I_detail.ne.0) THEN
        Write(out_niv,167) ine,(TeTab(ite),ite=1,7) !cc
        ENDIF
         Do ite=1,7
          moyrap(ite)=0.d0
          moyrap2(ite)=0.d0
          moyfit(ite)=0.d0
          moyfit2(ite)=0.d0
         Enddo
c simple power law:
       Do irai=1,218
         A(ine,irai)=ERC(ine,irai,5)
         lp(ine,irai)=dlog10(A(ine,irai)/ERC(ine,irai,3))/vb
         Do ite=1,7
          rap(ite)=ERC(ine,irai,ite)/
     $     (A(ine,irai)*(TeTab(ite)/TeTab(5))**lp(ine,irai))
          moyrap(ite)=moyrap(ite)+rap(ite)-1.d0
          moyrap2(ite)=moyrap2(ite)+dabs(rap(ite)-1.d0)
         Enddo
        IF(I_detail.ne.0) THEN
         Write(out_niv,168) irai,Tr(irai),lambda(irai),
     $     (rap(ite),ite=1,7),lp(ine,irai)             !cc
        ENDIF
c low Te corr:
         parl(ine,irai)=dlog10(rap(1))/(3.*flog)**pl
         Do ite=1,3
          rap(ite)=rap(ite)/10.**(parl(ine,irai)*
     $     ((-dlog10(TeTab(ite)/1.d4)-1.)**pl))
          moyfit(ite)=moyfit(ite)+rap(ite)-1.d0
          moyfit2(ite)=moyfit2(ite)+dabs(rap(ite)-1.d0)
         Enddo
c ite=4 (5000K) not corrected:
          moyfit(4)=moyfit(4)+rap(4)-1.d0
          moyfit2(4)=moyfit2(4)+dabs(rap(4)-1.d0)
c high Te corr:
         parh(ine,irai)=dlog10(rap(7))/flog**ph
         Do ite=5,7
          rap(ite)=rap(ite)/10.**(parh(ine,irai)*
     $     (dlog10(TeTab(ite)/1.d4))**ph)
          moyfit(ite)=moyfit(ite)+rap(ite)-1.d0
          moyfit2(ite)=moyfit2(ite)+dabs(rap(ite)-1.d0)
         Enddo
        IF(I_detail.ne.0) THEN
         Write(out_niv,172) irai,Tr(irai),lambda(irai),
     $     (rap(ite),ite=1,7),parl(ine,irai),parh(ine,irai) !cc
        ENDIF
       Enddo
c
         Do ite=1,7
          moyrap(ite)=moyrap(ite)/218.d0
          moyrap2(ite)=moyrap2(ite)/218.d0
          moyfit(ite)=moyfit(ite)/218.d0
          moyfit2(ite)=moyfit2(ite)/218.d0
         Enddo
c
        IF(I_detail.ne.0) THEN
       Write(out_niv,162) ine,(moyrap(ite),ite=1,7)
       Write(out_niv,163) ine,(moyrap2(ite),ite=1,7)
       Write(out_niv,160) ine,(moyfit(ite),ite=1,7)
       Write(out_niv,161) ine,(moyfit2(ite),ite=1,7)
        ENDIF
c
      Enddo
c
 162  Format(I10,5X,' moy ecart rap a 1.=',7F10.4)
 163  Format(I10,' moy abs(ecart rap a 1.)=',7F10.4)
 160  Format(I10,5X,' moy ecart fit a 1.=',7F10.4)
 161  Format(I10,' moy abs(ecart fit a 1.)=',7F10.4)
 168  Format(I3,1X,A22,F9.2,1p,7e10.3,0p,' pente=',F10.5)
 172  Format(I3,1X,A22,F9.2,1p,7e10.3,0p,' parl,a=',2F10.5)
c
c  Table of coeffs from D. Pequignot fit for the 218 lines and 4 Ne :
       Write(out_niv,25) pl,ph
 25    Format(' Coeffs fit D. Pequignot w/ pl, ph =',2F5.1,
     $     ' : A,lp,parl,parh')
c
        IF(I_detail.ne.0) THEN
 44   Format('  ine=',i3)
      Do irai=1,218
       Do ine=1,4
        Write(out_niv,43) irai,lambda(irai),A(ine,irai),lp(ine,irai),
     $    parl(ine,irai),parh(ine,irai)
       Enddo
      Enddo
 43   Format(i4,F10.2,1p,e9.2,0p,3F9.5)
c  Check Pequignot fit:
       Write(out_niv,76)
 76    Format(' *** Check Pequignot fit:')
      Do ine=1,4
       Write(out_niv,29) ine,(TeTab(ite),ite=1,7)
 29    Format(I5,8X,7F7.0)
        Do irai=1,218
         Do ite=1,7
           t4=TeTab(ite)*1.d-4
           dlt4=dlog10(t4)
           If(t4.lt.0.1d0) Then
            ERCf(ite)=1.d-15*A(ine,irai)*t4**lp(ine,irai)*
     $       10.**(parl(ine,irai)*(-dlt4-1.)**pl)
           Elseif(t4.lt.1.d0) Then
            ERCf(ite)=1.d-15*A(ine,irai)*t4**lp(ine,irai)
           Else
            ERCf(ite)=1.d-15*A(ine,irai)*t4**lp(ine,irai)*
     $       10.**(parh(ine,irai)*dlt4**ph)
           Endif
         Enddo
         Do ite=1,7
          rap(ite)=ERCf(ite)/(1.d-15*ERC(ine,irai,ite))
         Enddo
         Write(out_niv,28) irai,lambda(irai),(rap(ite),ite=1,7)
        Enddo
      Enddo
 28   Format(I3,F10.2,7F7.3)
        ENDIF
c
c  ****************
c  Table for X-SSN: 
c  ****************
c
      t=Te/1.d4
      lgt=dlog10(t)
      lgNe=dlog10(Ne)
c
        ine_synth=0
        do ine=1,4 ! ine=1 --> lgNe=2., etc.
         if(lgNe.lt.dble(ine+1).and.ine_synth.eq.0) then
          ine_synth=ine
         endif
        enddo
         IF(ine_synth.eq.1.or.ine_synth.eq.0) then
c Without interpol :
          ine=1
          if(ine_synth.eq.0) ine=4
c VOIR cas t exotique ?
        do irai=1,218
         if(if_rai(irai).eq.0) then !218-55 fit Peq
           If(t.lt.0.1d0) Then
            alf(irai)=1.d-15*A(ine,irai)*t**lp(ine,irai)*
     $       10.**(parl(ine,irai)*(-lgt-1.)**pl)
           Elseif(t4.lt.1.d0) Then
            alf(irai)=1.d-15*A(ine,irai)*t**lp(ine,irai)
           Else
            alf(irai)=1.d-15*A(ine,irai)*t**lp(ine,irai)*
     $       10.**(parh(ine,irai)*lgt**ph)
           Endif
         else  !55 Fang fit (Cas B avec cc())
            k=if_rai(irai)
            alf(irai)=alfit(t,cc(ine,k,1),cc(ine,k,2),
     $         cc(ine,k,3),cc(ine,k,4),cc(ine,k,5),
     $         cc(ine,k,6),cc(ine,k,7),cc(ine,k,8))
         endif
        enddo
         ELSE
            ine=ine_synth
            dlgNe=lgNe-dble(ine+1)
c Lin interpol in log(Ne):
        do irai=1,218
         if(if_rai(irai).eq.0) then !218-55 fit Peq
           If(t.lt.0.1d0) Then
            ERC1=1.d-15*A(ine,irai)*t**lp(ine,irai)*
     $       10.**(parl(ine,irai)*(-lgt-1.)**pl)
            ERC0=1.d-15*A(ine-1,irai)*t**lp(ine-1,irai)*
     $       10.**(parl(ine-1,irai)*(-lgt-1.)**pl)
           Elseif(t.lt.1.d0) Then
            ERC1=1.d-15*A(ine,irai)*t**lp(ine,irai)
            ERC0=1.d-15*A(ine-1,irai)*t**lp(ine-1,irai)
           Else
            ERC1=1.d-15*A(ine,irai)*t**lp(ine,irai)*
     $       10.**(parh(ine,irai)*lgt**ph)
            ERC0=1.d-15*A(ine-1,irai)*t**lp(ine-1,irai)*
     $       10.**(parh(ine-1,irai)*lgt**ph)
           Endif
         else  !55 Fang fit
            k=if_rai(irai)
            ERC1=alfit(t,cc(ine,k,1),cc(ine,k,2),
     $         cc(ine,k,3),cc(ine,k,4),cc(ine,k,5),
     $         cc(ine,k,6),cc(ine,k,7),cc(ine,k,8))
            ERC0=alfit(t,cc(ine-1,k,1),cc(ine-1,k,2),
     $         cc(ine-1,k,3),cc(ine-1,k,4),cc(ine-1,k,5),
     $         cc(ine-1,k,6),cc(ine-1,k,7),cc(ine-1,k,8))
         endif
          alf(irai)=ERC1+dlgNe*(ERC1-ERC0)
        enddo
         ENDIF

c Ref line X-SSN :  3  192  5679.56  5679.56   (28 185)_2011 
       irai_ref=192   ! (185_2011)
c Wavelength of ref line : 
       lambda_ref=lambda(irai_ref)
c Absolute emissivity of ref :
       int_ref=alf(irai_ref)*hc8/lambda_ref
c
c  Wavelength ranges:
c
      do iw=1,iwlim+1
         Irelmin(iw)=Irelmin_H(iw)*
     $   int_ref_mod_H/(int_ref*dmax1(Ionab,Ionab_min))
      enddo
      write(out_niv,63) wlmin,wlmax
 63   format(1p,2e10.3,' wlmin, wlmax X-SSN')
      write(out_niv,85) iwlim
 85   format(' iwlim ='i4,' puis wlim(1->iwlim), Irelmin(1->iwlim+1)')
      write(out_niv,87) (wlim(iw),iw=1,iwlim)
 87   format(5x,1p,13e10.3)
      write(out_niv,86) (Irelmin(iw),iw=1,iwlim+1)
 86   format(1p,13e10.3)
c
      Write(out_niv,210)
 210  Format(' Ref. for liste_modele.dat :')
      Write(out_niv,110) numion,ifre1,imult,ifre2,nomion,
     $  vzero,int(lambda_ref),int_ref
 110  Format(1x,A4,3I1,'000000',1x,A7,'        1.0   0.     1.000e+0',
     $  '  1.                 0   1  ',f5.1,1x,I5,'++',
     $  1p,e10.3,' erg.cm3.s-1')
      Write(out_niv,211)
 211  Format('Ref. for liste_phyat.dat :')
c   for liste_phyat.dat:
      Write(out_niv,111) numion4,ifre1,imult,ifre2,nomion,int_ref,
     $  int(lambda_ref),t,lgNe

 111  Format('9',A4,3I1,'000000',1x,A7,'        1.0   0.   ',
     $  1p,e10.3,'  1.               999   1   1.00 ',I5,
     $  '++ (t=',f6.2,' lgNe=',f6.2,')')
c
c counter of written lines:
               iraies=0
c
      Do irai=1,218
        IF(lambda(irai).ge.wlmin.and.lambda(irai).le.wlmax) Then
c
c VOIR air refract?
         vb=alf(irai)*hc8/lambda(irai)/int_ref
          iecr=0
         iw=1
         if(lambda(irai).le.wlim(iw).and.vb.ge.Irelmin(iw)) then
          iecr=1
         endif
        do iw=2,iwlim
         if(lambda(irai).gt.wlim(iw-1).and.lambda(irai).le.wlim(iw)
     $     .and.vb.ge.Irelmin(iw))then
          iecr=1
         endif
        enddo
        if(lambda(irai).gt.wlim(iwlim).and.vb.ge.Irelmin(iwlim+1))then
          iecr=1
        endif
c
             If(iecr.eq.1) Then
c
               iraies=iraies+1
          Write(out_niv,100) numion,ifre1,imult,ifre2,1000+irai,
     $     nomion,lambda(irai),vb,numion,ifre1,imult,ifre2,Trans(irai)
 100      Format(1x,A4,3I1,I4,'00',1x,A7,
     $      f13.3,' 0.   ',1p,e10.3,
     $      '  1.     ',A4,3I1,'000000   1   1.00 ',A32)  !2011: A33
             Endif
c
        ENDIF !end wlmin/max
      Enddo !end irai ecr
c
         Write(out_niv,75) iraies
 75      format('Number of X-SSN lines: iraies =',i5)
c
      write(lambda_refc,'(I5)') int(lambda_ref)
      write(iraiesc,'(I4)') iraies
      write(Irelminc,'(I5)') int(Irelmin(4)*1.d5)
      write(wlminc,'(I8)') int(wlmin)
      write(wlmaxc,'(I8)') int(wlmax)
      print *,'effective ni range: 2    ',
     $ '  wl range:',wlminc,wlmaxc
      print *,'number of emission lines = ',iraiesc,' / ',
     $ '   Irelmin =',Irelminc,'E-5',lambda_refc
c
	   CLOSE(out_niv)
c
c  If requested, insert the result .res in PHYAT for SYNTHESIS:
      if(iwr_phyat.ne.0) then
         Call Write_Phyat(lambda_output,PHYAT)
      endif
c
        Return
c
        End
c--------------------end NII_Sto-----------------------
c------------------------------------------------------
c
      Subroutine OII_Sto(Te,Ne,niv_input,lambda_output,PHYAT)
c
c  ***Based on Pete Storey 2012 private
c  ***Use data sets of the form lines_*.*_*.*_dr
c
      IMPLICIT NONE
c
      CHARACTER niv_input*(*),lambda_output*(*),PHYAT*(*)
c
c***Storey begin 1
      character*1 ul
      character*2 dr,suf
      character*3 ktem(9),kdens(16)
      character*4 endk
      character*5 lines
      character*16 name
      character*80 text
cc      real*8 dens(50),emms(50,50,20),flux(20),f(50),
      double precision  dens(50),emms(50,50,20),flux(20),f(50),
     :  drange,ddens,fv,xv,ss,dd,dummy,wave(20),
     :  temlog(9),logtemp,temp,escale,xx(10),
     :  e(10),fun,fluxerr(20),f1,f2,xmin,fdp,dx,logdens(16),
     :  fluxn(20),emm(20),tls(20,20),dls(20,20),
     :  em(20,20,5000),wl(20,20,10000)
      integer ndens,idens,nline,iline,ifd,ier,j,inorm,
     :  ierr,ntem,item,istop,maxit,icon,iprint,i,n,nf(20),
     :  ifixt,irerr,nrec(20,20),kc,ic
       character*6 up1(10000),lo1(10000)
       character*3 up2(10000),lo2(10000)
c***Storey end 1
c
      double precision emi(10000),wli(10000),ES(10000),EI(10000)
cc      double precision wlnot(10000)
      DOUBLE PRECISION vzero
      DOUBLE PRECISION int_ref,lambda_ref
      COMMON/MOCR1/wlmin,wlmax,wlim,Irelmin_H,Ionab,Ionab_min,iwlim
      DOUBLE PRECISION int_ref_mod_H
      COMMON/MOCR2/int_ref_mod_H
      DOUBLE PRECISION Ionab,Ionab_min
      DOUBLE PRECISION wlim(20),Irelmin_H(21)
      DOUBLE PRECISION wlmin,wlmax
      DOUBLE PRECISION Irelmin(21)
      COMMON/MOCR3/Irelmin
      CHARACTER*17 titre
      CHARACTER Irelminc*5,lambda_refc*6,iraiesc*4,wlminc*8,wlmaxc*8
      CHARACTER numion*4,numion4*4,nomion*7
      DOUBLE PRECISION Ne,Te,lgNe,dlgNe,t,lgTe
      DOUBLE PRECISION cmult
      INTEGER imult,ifre1,ifre2,mult,nmult
      INTEGER iwlim,iw,iecr,iraies,inep1   !,inot
c
       double precision ET(1000),emmin,em_rel
       integer ii,it,id,imaxl,ns,ni,nsmax,nimax
       integer is(1000)
       character*3 ns_alph,ni_alph
       character*1 zero,nc1
       character*2 nc2
c
      INTEGER in_niv,out_niv,iwr_phyat
      COMMON/MOC1/in_niv,out_niv,iwr_phyat
      COMMON/MOCV1/vzero
c
      DOUBLE PRECISION hc,hc8,hcsk,Ryd_inf
      COMMON/MOCC1/hc,hc8,hcsk,Ryd_inf
c***Storey begin 2
      common/calc/emms,wave,fluxn,fluxerr,logtemp,ndens,ntem,
     :          nline,ierr,ifixt
      common/mesh/temlog,logdens
c
      data ktem/'2.6','2.8','3.0','3.2','3.4','3.6','3.8',
     :            '4.0','4.2'/
      data kdens/'2.0','2.2','2.4','2.6','2.8',
     :           '3.0','3.2','3.4','3.6','3.8',
     :           '4.0','4.2','4.4','4.6','4.8','5.0'/
      data temlog/2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0,4.2/
      data logdens/2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0,
     :             4.2,4.4,4.6,4.8,5.0/
      data lines/'lines'/, ul/'_'/, dr/'dr'/
c***Storey end 2
c
       data ntem/9/,ndens/16/,imaxl/8793/,out_niv/11/
       data zero/'0'/
c
      print *,niv_input,'  ',lambda_output
c
      t=Te/1.d4
      lgTe=dlog10(Te)
      lgNe=dlog10(Ne)
cc data mal lu ?
cc      lines='lines'
cc      ul='_'
cc      dr='dr'
ccc data mal lu ?
cc      ktem(8)='4.0'
cc      ktem(9)='4.2'
cc      kdens(6)='3.0'
c
        OPEN(unit=in_niv,file=niv_input,status='old')
c
      read(in_niv,15) titre
 15   format(A17)
 1    format(a1)
c Stdd di_electronic: imult = 2
      read(in_niv,121) mult,nmult,imult,cmult,ifre1,ifre2
 121  format(6x,i4,6x,i4,6x,i4,6x,d7.0,6x,i4,6x,i4)
c Label + name ion (form X-SSN):
      read(in_niv,60) numion,numion4,nomion
 60   format(9x,a4,11x,a4,10x,a7)
c
        CLOSE(in_niv)
c
c  Look for nearest table without interpolation.
c  The first Te (it=1) is 10.**2.6; steps .2dex; closest in log: 
      it=IDNINT(5.*dlog10(Te))-12
      If(it.lt.1) then
       it=1
       elseif(it.gt.9) then
       it=9
      Endif
c  The first Ne (id=1) is 10.**2 ;steps .2dex; closest in log: 
      id=IDNINT(5.*dlog10(Ne))-9
      If(id.lt.1) then
       id=1
       elseif(id.gt.16) then
       id=16
      Endif
      print *,'OII it,id =',it,id
c
* read theoretical emissivities from files "lines_...." and store 
* wavelengths and emissivities

cc      do item=1,ntem
cc        do idens=1,ndens
c (All tables of levels and lambdas identical : read one only)
      do item=it,it
        do idens=id,id
          name=lines//ul//ktem(item)//ul//kdens(idens)//ul//dr
          tls(item,idens)=temlog(item)
          dls(item,idens)=logdens(idens)
         write(6,*) 'open and read file ',name
c
        OPEN(unit=15,status='old',file=name)
c
          do i=1,13 !initial text = 14 lines, but 1 blank line skipped
           read(15,*) text
          enddo
c
          do i=1,imaxl
           IF(i.ne.260.and.i.ne.261) THEN
      read(15,1010)up1(i),up2(i),lo1(i),lo2(i),
     $             emi(i),wli(i),ES(i),EI(i)
      emi(i)= emi(i)*hc8/wli(i)
           ELSE
c To remove '***********' for too long wl:
      read(15,2010)up1(i),up2(i),lo1(i),lo2(i),
     $  emi(i)       ,ES(i),EI(i)
      wli(i)=99999999.99
      emi(i)= 0.
           ENDIF
          enddo
 1010 format(2x,a6,1x,a3,5x,a6,1x,a3,31x,e13.4,f11.2,1x,f11.2,f11.2)
 2010 format(2x,a6,1x,a3,5x,a6,1x,a3,31x,e13.4,  11X,1x,f11.2,f11.2)
c
        close(15)
c
        enddo
      enddo
c
c Table of energy levels: ET
c Levels are numbered by increasing energies
c Only the first 4 metastable levels are always lower levels
      ET(1)=EI(imaxl)
      ET(2)=EI(imaxl-1)
      ET(3)=EI(imaxl-2)
      ET(4)=EI(imaxl-6)
      do ns=1,4
       is(ns)=0 ! 0 for irrelevant
      enddo
      ns=5
      ET(ns)=ES(imaxl)
      is(ns)=imaxl
      do i=imaxl-1,1,-1  !i = line number
        if(ES(i).ne.ET(ns)) then
         ns=ns+1 !next upper energy level
         ET(ns)=ES(i)
         is(ns)=i ! number i of the 1st line from this energy level
        endif
      enddo
      nsmax=ns !total number of (upper) levels
      print *,'nsmax =',nsmax
c  for DO loop i=:
       is(nsmax+1)=0
c
         OPEN(unit=out_niv,file=lambda_output,status='old')
c
c  Sum multiplet V1 '4649' = reference emissivity
c  ***TBD: interpolate for emissivity of ref multiplet?***
      lambda_ref=4649.d0   !EXPLICIT
      int_ref=emi(8713)+emi(8729)+emi(8730)+    !EXPLICIT
     $   emi(8736)+emi(8737)+emi(8738)+emi(8745)+emi(8746)
c
c  Wavelength ranges:
c
      do iw=1,iwlim+1
         Irelmin(iw)=Irelmin_H(iw)*
     $   int_ref_mod_H/(int_ref*dmax1(Ionab,Ionab_min))
      enddo
c
      write(out_niv,15) titre
      write(out_niv,116)
 116   format(' bbb',/,' bbb')
      write(out_niv,141) mult,nmult,imult,cmult,ifre1,ifre2
 141   format('  mult =',i4,' nmult =',i4,' imult =',i4,
     $   ' cmult =',f7.4,' ifre1 =',i4,' ifre2 =',i4)
      write(out_niv,67) Ne,Te
 67   format(1p,' Ne =',d9.2,'  Te =',d9.2)
      write(out_niv,62) numion,numion4,nomion
 62   format(' numion :',a4,' numion4 :',a4,' nomion :',a7,' a4,a4,a7')
      write(out_niv,63) wlmin,wlmax
 63   format(1p,2d10.3,' wlmin, wlmax X-SSN')
      write(out_niv,85) iwlim
 85   format(' iwlim ='i4,' puis wlim(1->iwlim), Irelmin(1->iwlim+1)')
      write(out_niv,87) (wlim(iw),iw=1,iwlim)
 87   format(5x,1p,13e10.3)
      write(out_niv,86) (Irelmin(iw),iw=1,iwlim+1)
 86   format(1p,13d10.3)
c
c  START WRITING OII FOR SYNTHESIS: 
      Write(out_niv,210)
 210  Format(' Ref. for liste_modele.dat :')
      Write(out_niv,110) numion,ifre1,imult,ifre2,nomion,vzero,
     $  int(lambda_ref),int_ref
 110  Format(1x,A4,3I1,'000000',1x,A7,'        1.0   0.     1.000e+0',
     $  '  1.                 0   1  ',f5.1,1x,I5,'++',
     $  1p,e10.3,' erg.cm3.s-1')
      Write(out_niv,211)
 211  Format(' Ref. for liste_phyat.dat :')
      Write(out_niv,111) numion4,ifre1,imult,ifre2,nomion,int_ref,
     $  int(lambda_ref),name
 111  Format('9',A4,3I1,'000000',1x,A7,'        1.0   0.   ',
     $  1p,e10.3,
     $  '  1.               999   1   1.00 ',I5,'++ (',a16,')')

c counter des raies ecrites :
               iraies=0
cc               inot=0
c
      DO ns=5,nsmax              !ns
        DO i=is(ns),is(ns+1)+1,-1  !i
          ni=0
          do ii=1,ns               !ii
           if(EI(i).eq.ET(ii)) then
            ni=ii
           endif
          enddo                    !ii
c secu :
          if(ni.eq.0) then
            print *,'***STOP OII: no ni found for ns, i =',ns,i
            STOP
          endif
         em_rel=emi(i)/int_ref !relative emissivity of line
c
              iecr=0
             IF(wli(i).gt.wlmin.and.wli(i).lt.wlmax) THEN
c
c Adopted minimum rel emiss depends on lambda range:
         iw=1
         if(wli(i).le.wlim(iw).and.em_rel.ge.Irelmin(iw)) then
              iecr=1
         endif
        do iw=2,iwlim
         if(wli(i).gt.wlim(iw-1).and.wli(i).le.wlim(iw)
     $     .and.em_rel.ge.Irelmin(iw))then
              iecr=1
         endif
        enddo
        if(wli(i).gt.wlim(iwlim).and.em_rel.ge.Irelmin(iwlim+1))then
          iecr=1
        endif
             ENDIF
c (Check)
cc          If(iecr.eq.0) then
cc             inot=inot+1
cc             wlnot(inot)=wli(i)
cc          Endif
             If(iecr.eq.1) Then
c counter of written lines: 
               iraies=iraies+1
c
        if(ns.lt.10) then
       write(nc1,'(I1)') ns
        ns_alph=zero//zero//nc1
         elseif(ns.lt.100) then
        write(nc2,'(I2)') ns
       ns_alph=zero//nc2
          else   !(assume ns < 1000)
        write(ns_alph,'(I3)') ns
        endif
c
        if(ni.lt.10) then
       write(nc1,'(I1)') ni
        ni_alph=zero//zero//nc1
         elseif(ni.lt.100) then
        write(nc2,'(I2)') ni
       ni_alph=zero//nc2
          else   !(assume ni < 1000)
        write(ni_alph,'(I3)') ni
        endif
c
        Write(out_niv,190) numion,ifre1,imult,ifre2,ns_alph,ni_alph,
     $      nomion,wli(i),em_rel,numion,ifre1,imult,ifre2,
     $      up1(i),up2(i),lo1(i),lo2(i)
 190    Format(1x,A4,3I1,2A3,1x,A7,f13.3,' 0.   ',1p,e10.3,
     $   '  1.     ',A4,3I1,'000000   1   1.00 ',A6,A3,'-',A6,A3)
             Endif
c
        ENDDO                     !i
      ENDDO                       !ns
c
         Write(out_niv,75) iraies
 75      format('Number of X-SSN lines: iraies =',i5)
c  (Check)
cc         write(out_niv,90) inot
cc         write(out_niv,90) (i,wlnot(i),i=1,inot)
cc 90      format(6(i5,f11.2))
c
      write(lambda_refc,'(I5)') int(lambda_ref)
      write(iraiesc,'(I4)') iraies
      write(Irelminc,'(I5)') int(Irelmin(4)*1.d5)
      write(wlminc,'(I8)') int(wlmin)
      write(wlmaxc,'(I8)') int(wlmax)
      print *,'effective ni range: 2    ',
     $ '  wl range:',wlminc,wlmaxc
      print *,'number of emission lines = ',iraiesc,' / ',
     $ '   Irelmin =',Irelminc,'E-5',lambda_refc
c
           CLOSE(out_niv)
c
c  If requested, insert the result .res in PHYAT for SYNTHESIS:
      if(iwr_phyat.ne.0) then
         Call Write_Phyat(lambda_output,PHYAT)
      endif
c
      return
c
      end
c
c---------------------------end OII_Sto------------------
c--------------------------------------------------------
c
      Subroutine OIII_5g(Te,Ne,niv_input,lambda_output,PHYAT)
c
      Implicit double precision (A-H,O-Z)
c
      CHARACTER niv_input*(*),lambda_output*(*),PHYAT*(*)
c
      DOUBLE PRECISION vzero
      DOUBLE PRECISION Te,Ne
      CHARACTER*17 titre
      CHARACTER Irelminc*5,lambda_refc*6,iraiesc*4,wlminc*8,wlmaxc*8
      CHARACTER numion*4,numion4*4,nomion*7
      INTEGER imult,ifre1,ifre2,mult,nmult
      COMMON/MOCR1/wlmin,wlmax,wlim,Irelmin_H,Ionab,Ionab_min,iwlim
      DOUBLE PRECISION int_ref_mod_H
      COMMON/MOCR2/int_ref_mod_H
      DOUBLE PRECISION Ionab,Ionab_min
      DOUBLE PRECISION wlim(20),Irelmin_H(21)
      DOUBLE PRECISION wlmin,wlmax
      DOUBLE PRECISION Irelmin(21)
      COMMON/MOCR3/Irelmin
      DOUBLE PRECISION cmult
c
      INTEGER in_niv,out_niv,iwr_phyat
      COMMON/MOC1/in_niv,out_niv,iwr_phyat
c
      DOUBLE PRECISION hc,hc8,hcsk,Ryd_inf
      COMMON/MOCC1/hc,hc8,hcsk,Ryd_inf
      COMMON/MOCV1/vzero
c 5g-6h :
      DOUBLE PRECISION alfeff,b_SH,asa
      DOUBLE PRECISION wl6h(6),rel6h(6)
      INTEGER inep1,iz
      CHARACTER zc*1,tabsto*9
      INTEGER m,k1,k2
      COMMON/MOCS1/alfeff(50,50),b_SH(50,0:50) !alfeff(n,lp) lp=l+1 pr 5g-6h
      COMMON/MOCS2/inep1,iz,zc,tabsto !a ver
c
c***  ADAPTED FROM SUBR. O3_5g.f IN NEBU  ***
c
c According to R. Kisielius & P.J. Storey, 1999, A&AS 137, 157
c
c   Calcul. of eff rec coeff alf(k) on to 5g levels k of de OIII
c   then emissivities relint(k,l) of 4f-5g lines, relatives 
c   to absolute reference emissivities absref(k) in erg/s/ion. 
c   
c   Lev 5g (1--12 equivlt to 13--24 in KS) : 
c   (2P1/2)5g  1G4 3G3 3H4 3H5
c   (2P3/2)5g  3G4 3G5 3F4 3F3 3H6 1H5 3F2 1F3
c
c  (dim 12 refers to kmax)
	Double Precision pstat(12),g1mix(12),g2mix(12)
	Double Precision a1(12),b1(12),c1(12),d1(12)
	Double Precision a2(12),b2(12),c2(12),d2(12)
	Double Precision alfdir1(12),alfdir2(12),alf(12)
	Double Precision relint(12,9),absref(12),relref(12)
	Double Precision wl1(6),wl2(8),wl3(6),wl4(4),wl5(7),wl6(4)
	Double Precision wl7(7),wl8(9),wl9(1),wl10(3),wl11(8),wl12(7)
	Double Precision br1(6),br2(8),br3(6),br4(4),br5(7),br6(4)
	Double Precision br7(7),br8(9),br9(1),br10(3),br11(8),br12(7)
	Double Precision branch(12,9),wvlgth(12,9)
        Integer nlmax(12),lref(12)
c
c  EOIV(cm-1) from Feuchtgruber et al (1997, ApJ) : 
        Data EOIV/386.245/
c  A21(s-1) see Mendoza (1983) : 
        Data A21/5.166D-4/
c
c  Nber of upper lev:
        Data kmax/12/
c  Label upper lev of 'general' ref line:
        Data kref/1/
c  pstat=2J+1 of upper lev (5g): 
	Data pstat/9.d0,7.d0,9.d0,11.d0,9.d0,11.d0,9.d0,7.d0,
     $          13.d0,11.d0,5.d0,7.d0/
c  Nber of lines from each upper lev:
	Data nlmax/6,8,6,4,7,4,7,9,1,3,8,7/
c  Label subref line of each lev:
	Data lref/ 1,2,2,1,4,2,6,8,1,3,7,7/
c  Mixing coeffs(g2=1-g1 see below), Table 3 :
	Data g1mix/0.85027d0,0.85015d0,0.96220d0,0.96187d0,
     $          0.03873d0,0.03883d0,0.06558d0,0.06779d0,
     $          0.03393d0,0.03395d0,0.00153d0,0.00185d0/
c  Coeffs for alfdir1, Table 4 :
	Data a1/11.8009d0,9.1260d0,13.1900d0,16.0049d0,
     $      1.0613d0,1.3070d0,1.2819d0,1.0075d0,
     $      1.5039d0,1.2806d0,0.4670d0,0.7113d0/
	Data b1/-0.5620d0,-0.5640d0,-0.5469d0,-0.5518d0,
     $      -1.3123d0,-1.3069d0,-1.2250d0,-1.2138d0,
     $      -1.3207d0,-1.3284d0,-1.4083d0,-1.4063d0/
	Data c1/-0.2070d0,-0.2020d0,-0.2135d0,-0.2059d0,
     $      0.1673d0,0.1578d0,0.1855d0,0.1904d0,
     $      0.1283d0,0.1386d0,-0.0680d0,-0.0680d0/
	Data d1/0.0255d0,0.0241d0,0.0266d0,0.0245d0,
     $      -0.0245d0,-0.0231d0,-0.0311d0,-0.0315d0,
     $      -0.0176d0,-0.0197d0,0.0117d0,0.0121d0/
c  Coeffs for alfdir2, Table 4 :
	Data a2/1.0104d0,0.7891d0,0.2551d0,0.3136d0,
     $      6.5977d0,8.0183d0,6.3031d0,4.9188d0,
     $      9.3645d0,7.9356d0,3.6569d0,5.1570d0/
	Data b2/-0.5473d0,-0.5456d0,-0.5411d0,-0.5413d0,
     $      -0.5445d0,-0.5486d0,-0.5445d0,-0.5439d0,
     $      -0.5432d0,-0.5413d0,-0.5487d0,-0.5444d0/
	Data c2/-0.2022d0,-0.2064d0,-0.2139d0,-0.2127d0,
     $      -0.2134d0,-0.2083d0,-0.2146d0,-0.2148d0,
     $      -0.2167d0,-0.2197d0,-0.2106d0,-0.2164d0/
	Data d2/0.0235d0,0.0248d0,0.0270d0,0.0264d0,
     $      0.0264d0,0.0252d0,0.0267d0,0.0270d0,
     $      0.0271d0,0.0279d0,0.0256d0,0.0271d0/
c  Wavelengths (A, air):
        Data wl1/4366.997d0,4378.300d0,4381.370d0,
     $      4484.716d0,4534.308d0,4537.992d0/
        Data wl2/4366.997d0,4369.573d0,4378.300d0,
     $      4477.926d0,4484.716d0,4534.308d0,4539.002d0,4573.660d0/
        Data wl3/4365.242d0,4376.536d0,4379.604d0,4476.082d0,
     $      4482.866d0,4536.098d0/
        Data wl4/4379.604d0,4482.866d0,4522.002d0,4536.098d0/
        Data wl5/4300.735d0,4311.698d0,4314.675d0,4408.283d0,
     $      4414.863d0,4452.815d0,4466.483d0/
        Data wl6/4314.675d0,4414.863d0,4452.815d0,4466.483d0/
        Data wl7/4292.481d0,4306.367d0,4399.611d0,
     $      4406.165d0,4443.967d0,4454.026d0,4457.580d0/
        Data wl8/4294.970d0,4303.401d0,4306.367d0,4399.611d0,
     $      4406.165d0,4454.026d0,4457.580d0,4458.555d0,4491.991d0/
        Data wl9/4434.604d0/
        Data wl10/4396.960d0,4434.604d0,4448.159d0/
        Data wl11/4277.146d0,4279.618d0,4287.988d0,4383.502d0,
     $      4437.517d0,4442.013d0,4471.037d0,4475.200d0/
        Data wl12/4277.146d0,4279.618d0,4290.933d0,
     $      4390.009d0,4437.517d0,4442.013d0,4475.200d0/
c  Branching ratios: 
        Data br1/8.029d-1,3.142d-2,7.798d-3,
     $      2.879d-2,1.273d-1,1.365d-3/
        Data br2/4.471d-2,7.917d-1,2.718d-3,
     $      2.995d-2,1.209d-3,6.034d-3,1.209d-1,2.477d-3/
        Data br3/2.239d-2,7.260d-1,2.031d-2,2.221d-1,
     $      4.774d-3,4.340d-3/
        Data br4/7.545d-1,2.321d-1,2.393d-3,1.095d-2/
        Data br5/8.400d-3,2.070d-1,5.550d-3,6.818d-1,
     $      3.588d-2,1.720d-3,5.962d-2/
        Data br6/2.376d-1,6.537d-1,7.657d-2,3.214d-2/
        Data br7/1.244d-1,7.362d-2,1.845d-3,
     $      5.775d-2,2.727d-3,7.370d-1,2.357d-3/
        Data br8/1.187d-1,7.532d-2,2.596d-3,6.241d-2,
     $      3.242d-3,3.625d-2,1.913d-3,6.881d-1,1.090d-2/
        Data br9/1.000d-0/
        Data br10/3.619d-2,1.797d-2,9.452d-1/
        Data br11/2.259d-3,4.650d-2,3.286d-3,3.235d-3,
     $      6.380d-3,1.107d-1,7.647d-1,6.286d-2/
        Data br12/4.692d-2,3.029d-3,2.292d-3,
     $      2.890d-3,9.009d-2,2.865d-3,8.509d-1/
c 5g-6h :
      Data wl6h/8238.5d0,8244.1d0,8227.6d0,8250.8d0,8268.8d0,8286.7d0/
c
      print *,niv_input,'  ',lambda_output
c
        dldem2=dlog10(Ne)-2.d0
        tred=Te/1.D4
c
        OPEN(unit=in_niv,file=niv_input,status='old')
c
      read(in_niv,15) titre
 15   format(A17)
      read(in_niv,121) mult,nmult,imult,cmult,ifre1,ifre2
 121  format(6x,i4,6x,i4,6x,i4,6x,d7.0,6x,i4,6x,i4)
c Label + name (X-SSN form):
      read(in_niv,60) numion,numion4,nomion
 60   format(9x,a4,11x,a4,10x,a7)
c
        CLOSE(in_niv)
c
        Do k=1,kmax
           g2mix(k)=1.d0-g1mix(k)
        Enddo
c  Rec coeff by cascades on to OIII(5g), approx H-like (Storey):
c   0.5<t4<2.0. Extrapolate to low t4? VOIR
        alfcasc_5g=0.69808D-12*tred**(-0.94619)*
     $     (1.+tred*(-0.20470+tred*0.02837))*
     $     (1.+dldem2*(-8.6442d-3+dldem2*(-1.3913d-3-5.993d-4*dldem2)))
c  Fit eff coll strength [OIV] 2P-2P Blum & Pradhan (1992) : 
cc        Omega=3.6292-0.5728*dlog10(tred)-1.4231/tred+0.2101/tred**2
cccc
c 10/01 : Extract from nebula/omegt.f : 
c O IV f.s. Blum Pradhan 92 ApJS 80, 425 t4 0.1-->4.
c  + Zhang Graziani Pradhan AA94 ~*0.94 uniform! + large t4.
        t4=tred
	if(t4.lt.0.200d0) then
	Omega=1.424+1.13*t4
	else
	if(t4.lt.1.093d0) then
	Omega=2.27+0.775*(t4-1.0)
	else
	if(t4.lt.4.0d0) then
	Omega=2.39+0.100*(t4-1.5)
	else
	if(t4.lt.16.d0) then
	Omega=2.64
	else
	Omega=2.64*(t4/16.)**(-0.35)
	endif
	endif
	endif
	endif
cccc
        C21=8.63D-6/4./dsqrt(Te)*Omega
        C12=8.63D-6/2./dsqrt(Te)*Omega*dexp(-hcsk*EOIV/Te)
        X=Ne*C12/(Ne*C21+A21)
        popO4_1=1./(1.+X)
        popO4_2=X/(1.+X)
        do k=1,kmax
c  Coeff rec directes et totales :
        alfdir1(k)=1.d-15*a1(k)*tred**b1(k)*(1.+tred*(c1(k)+tred*d1(k)))
        alfdir2(k)=1.d-15*a2(k)*tred**b2(k)*(1.+tred*(c2(k)+tred*d2(k)))
        alf(k)=popO4_1*(alfdir1(k)+g1mix(k)*alfcasc_5g*pstat(k)/36.)+
     $         popO4_2*(alfdir2(k)+g2mix(k)*alfcasc_5g*pstat(k)/72.)
        enddo
        do k=1,kmax
           do l=1,9
              relint(k,l)=0.d0
           enddo
        enddo
           k=1
           do l=1,nlmax(k)
              branch(k,l)=br1(l)
              wvlgth(k,l)=wl1(l)
           enddo
           k=2
           do l=1,nlmax(k)
              branch(k,l)=br2(l)
              wvlgth(k,l)=wl2(l)
           enddo
           k=3
           do l=1,nlmax(k)
              branch(k,l)=br3(l)
              wvlgth(k,l)=wl3(l)
           enddo
           k=4
           do l=1,nlmax(k)
              branch(k,l)=br4(l)
              wvlgth(k,l)=wl4(l)
           enddo
           k=5
           do l=1,nlmax(k)
              branch(k,l)=br5(l)
              wvlgth(k,l)=wl5(l)
           enddo
           k=6
           do l=1,nlmax(k)
              branch(k,l)=br6(l)
              wvlgth(k,l)=wl6(l)
           enddo
           k=7
           do l=1,nlmax(k)
              branch(k,l)=br7(l)
              wvlgth(k,l)=wl7(l)
           enddo
           k=8
           do l=1,nlmax(k)
              branch(k,l)=br8(l)
              wvlgth(k,l)=wl8(l)
           enddo
           k=9
           do l=1,nlmax(k)
              branch(k,l)=br9(l)
              wvlgth(k,l)=wl9(l)
           enddo
           k=10
           do l=1,nlmax(k)
              branch(k,l)=br10(l)
              wvlgth(k,l)=wl10(l)
           enddo
           k=11
           do l=1,nlmax(k)
              branch(k,l)=br11(l)
              wvlgth(k,l)=wl11(l)
           enddo
           k=12
           do l=1,nlmax(k)
              branch(k,l)=br12(l)
              wvlgth(k,l)=wl12(l)
           enddo
c
           do k=1,kmax
           do l=1,nlmax(k)
            Ener=1.d8/wvlgth(k,l)
            refr=1.+6432.8d-8+2949810.d0/(146.d8-Ener**2)+
     $           25540.d0/(41.d8-Ener**2)
            wvlvac=wvlgth(k,l)*refr
            relint(k,l)=alf(k)*branch(k,l)/wvlvac
           enddo
           enddo
c
        do k=1,kmax
           absref(k)=relint(k,lref(k))
           do l=1,nlmax(k)
c  Emissivites relatives aux raies de ref :
              relint(k,l)=relint(k,l)/absref(k)
           enddo
        enddo
c  Emissivites absolues des raies de ref en erg/s/ion :
        do k=1,kmax
           absref(k)=absref(k)*hc8
        enddo
c  Emissivites relatives des raies de ref entre elles :
        do k=1,kmax
           relref(k)=absref(k)/absref(kref)
        enddo
c  Emissivite totale relative a la ref generale (kref) :
        emitot_rel=0.d0
        do k=1,kmax
          va=0.
           do l=1,nlmax(k)
            va=va+relint(k,l)
           enddo
          emitot_rel=emitot_rel+relref(k)*va
        enddo
c
c  Emissivity [O IV]26um /electron /ion :
c       wl in air in A :
c        refr=1.+6432.8d-8+2949810/(146d8-EOIV**2)+
c     $           25540.d0/(41.d8-EOIV**2)
c10/11   wlOIV in vacuum
        wlOIV=1.d8/EOIV
        absOIV=A21*popO4_2*hc8/wlOIV/Ne
c
c  Wavelength ranges: **inactive for OIII_5g
c
      do iw=1,iwlim+1
         Irelmin(iw)=Irelmin_H(iw)*
     $   int_ref_mod_H/(int_ref*dmax1(Ionab,Ionab_min))
      enddo
c
	   OPEN(unit=out_niv,file=lambda_output,status='old')
c
      write(out_niv,15) titre
      write(out_niv,116)
 116   format(' bbb',/,' bbb')
      write(out_niv,141) mult,nmult,imult,cmult,ifre1,ifre2
 141   format('  mult =',i4,' nmult =',i4,' imult =',i4,
     $   ' cmult =',f7.4,' ifre1 =',i4,' ifre2 =',i4)
      write(out_niv,67) Ne,Te
 67   format(1p,' Ne =',e9.2,'  Te =',e9.2)
      write(out_niv,62) numion,numion4,nomion
 62   format(' numion :',a4,' numion4 :',a4,' nomion :',a7,' a4,a4,a7')
      write(out_niv,63) wlmin,wlmax
 63   format(1p,2e10.3,' wlmin, wlmax X-SSN')
      write(out_niv,85) iwlim
 85   format(' iwlim ='i4,' puis wlim(1->iwlim), Irelmin(1->iwlim+1)')
      write(out_niv,87) (wlim(iw),iw=1,iwlim)
 87   format(5x,1p,13e10.3)
      write(out_niv,86) (Irelmin(iw),iw=1,iwlim+1)
 86   format(1p,13e10.3)
	   Write(out_niv,3)
 3         Format(6x,'Recombination lines O III 4f-5g')
	   Write(out_niv,7) Te,Ne,popO4_1,popO4_2
 7	Format(1p,'  Te=',e9.3,' Ne=',e9.3,'  PopOIV = ',e9.3,1x,e9.3)
	   Write(out_niv,9)
 9      Format('   Recombination coefficients for O III 5g')
	      Write(out_niv,1) (alf(k),k=1,kmax)
 1	Format(1p,6e12.3)
	   Write(out_niv,5) 
 5      Format('   Wavelengths, then relative emissivities:')
	   Do k=1,kmax
	      Write(out_niv,8) k,wvlgth(k,lref(k)),absref(k)
 8      Format(1x,'k=',i2,'   raie ref :',f10.3, 
     $     ' emissivity=',1p,e10.3)
	      Write(out_niv,6) (wvlgth(k,l),l=1,nlmax(k))
 6	Format(6f12.3)
	      Write(out_niv,1) (relint(k,l),l=1,nlmax(k))
	   Enddo
              Write(out_niv,4)
 4            Format('   Relative intensites of Ref lines:')
 	      Write(out_niv,6) (wvlgth(k,lref(k)),k=1,kmax)
 	      Write(out_niv,1) (relref(k),k=1,kmax)
c
      Write(out_niv,300) wlOIV/1.d4,wvlgth(kref,lref(kref)),
     $       absOIV/absref(kref),int(wvlgth(kref,lref(kref))),emitot_rel
 300    Format('   I([OIV]',0p,F8.4,'um) / I(OIII',F9.3,'A) =',
     $         1p,e10.3,' I(OIII 4f-5g)/I(',0p,I5,') =',1p,e10.3)
c
c  **********************
c  Tabulation for X-SSN :
c  **********************
c General ref line of 4f-5g (here the 1st one: see relref):
      Write(out_niv,210)
 210  Format(' Ref. for liste_modele.dat :')
c   for liste_modele.dat :
       Write(out_niv,1110) numion,ifre1,imult,ifre2,nomion,vzero,
     $     int(wvlgth(kref,lref(kref))),absref(kref)
 1110  Format(1x,A4,3I1,'000000',1x,A7,
     $      '        1.0   0.     1.000e+0',
     $      '  1.                 0   1  ',f5.1,1x,I5,'++',
     $      1p,e10.3,' erg.cm3.s-1')
       Write(out_niv,211)
 211  Format(' Ref. for liste_phyat.dat :')
c   for liste_phyat.dat :

      Write(out_niv,1100) numion4,ifre1,imult,ifre2,nomion,absref(kref),
     $    int(wvlgth(kref,lref(kref))),Te/1.d4,dlog10(Ne),emitot_rel
 1100   Format('9',A4,3I1,'000000',1x,A7,
     $      '        1.0   0.   ',1p,e10.3,
     $      '  1.               999   1   1.00 ',I5,'++ t=',
     $      0p,f4.2,' lgNe=',f4.1,' 4f-5g/ref=',f5.2)
c Line counter:
          iraies=0
        Do k=1,kmax
          If(k.lt.10) then
c Raie ref niv :
          iraies=iraies+1
           Write(out_niv,110) numion,ifre1,imult,ifre2,k,nomion,
     $      relref(k),numion,ifre1,imult,ifre2,int(wvlgth(k,lref(k)))
 110    Format(1x,A4,3I1,'000',I1,'00',1x,A7,'        1.0   0.   ',
     $      1p,e10.3,'  1.     ',A4,3I1,'000000   1   1.00 ',I5,'+')
            Do l=1,nlmax(k)
          iraies=iraies+1
c Other lines (assumed: < 10 per level)
              Write(out_niv,10) numion,ifre1,imult,ifre2,k,l,nomion,
     $      wvlgth(k,l),relint(k,l),numion,ifre1,imult,ifre2,k
 10     Format(1x,A4,3I1,'000',I1,I1,'0',1x,A7,f13.3,' 0.   ',
     $      1p,e10.3,'  1.     ',A4,3I1,'000',I1,'00   1   1.00 ')
            Enddo
          Else
          iraies=iraies+1
           Write(out_niv,111) numion,ifre1,imult,ifre2,k,nomion,
     $      relref(k),numion,ifre1,imult,ifre2,int(wvlgth(k,lref(k)))
 111    Format(1x,A4,3I1,'00',I2,'00',1x,A7,'        1.0   0.   ',
     $      1p,e10.3,'  1.     ',A4,3I1,'000000   1   1.00 ',I5,'+')
            Do l=1,nlmax(k)
          iraies=iraies+1
              Write(out_niv,11) numion,ifre1,imult,ifre2,k,l,nomion,
     $      wvlgth(k,l),relint(k,l),numion,ifre1,imult,ifre2,k
 11     Format(1x,A4,3I1,'00',I2,I1,'0',1x,A7,f13.3,' 0.   ',
     $      1p,e10.3,'  1.     ',A4,3I1,'00',I2,'00   1   1.00 ')
            Enddo
          Endif
        Enddo
c
c  Adding transitions 5g-6h, approximately transposed from 4f-5g:
c*****  Tabulation of alfeff:
c charge :
      iz=3
      write(zc,'(I1)') iz ! conversion to character 
c
      Call ReadSto(Te,Ne)
c
      ifre2b=6
      nomion='O_III6h'
      asa=alfeff(6,5+1)/alfeff(5,4+1) ! 6h / 5g
c Six main lines. Neglected: ~ 10 lines 5-10% of main lines. 
c  sub-ref :
          iraies=iraies+1
           Write(out_niv,410) numion,ifre1,imult,ifre2b,nomion,
     $      numion,ifre1,imult,ifre2,int(wvlgth(kref,lref(kref)))
 410    Format(1x,A4,3I1,'000000',1x,A7,'        1.0   0.     ',
     $  '1.000e+0  1.     ',A4,3I1,'000000   1   1.00 (',I5,'+) 5g-6h')
       Do m=1,6
        k1=2*m-1
        k2=2*m
        if(k1.eq.7) k1=10 ! EXPLICIT ech de niv  ! 2/16: ?
        if(k2.eq.10) k2=7
        rel6h(m)=relref(k1)*wvlgth(k1,lref(k1))
        rel6h(m)=rel6h(m)+relref(k2)*wvlgth(k2,lref(k2))
        rel6h(m)=rel6h(m)*asa/wl6h(m)
          iraies=iraies+1
              Write(out_niv,310) numion,ifre1,imult,ifre2b,m,nomion,
     $      wl6h(m),rel6h(m),numion,ifre1,imult,ifre2b
 310     Format(1x,A4,3I1,'000',I1,'00',1x,A7,f13.3,' 0.   ',
     $      1p,e10.3,'  1.     ',A4,3I1,'000000   1   1.00 ')
      Enddo
c
         Write(out_niv,75) iraies
 75      format('Number of X-SSN lines: iraies =',i5)
c
      write(lambda_refc,'(I5)') int(wvlgth(kref,lref(kref)))
      write(iraiesc,'(I4)') iraies
      Irelminc='  0.0'
      print *,'number of emission lines = ',iraiesc,' / ',
     $ '   Irelmin =',Irelminc,'E-5',lambda_refc
c
	   CLOSE(out_niv)
c
c  If requested, insert the result .res in PHYAT for SYNTHESIS:
      if(iwr_phyat.ne.0) then
         Call Write_Phyat(lambda_output,PHYAT)
      endif
c
        Return
c
        End
c-----------------------end OIII_5g---------------------
c--------------------------------------------------------
c
      Subroutine ReadSto(Te,Ne)
c
c H-like effective recomb coeff fm tables Storey, Hummer 1995
c   ATT: maxi n = 50
c
      IMPLICIT NONE
c
      INTEGER iz
      INTEGER it,ine,inep1,i,j,n,l,lp,istor,istow
      DOUBLE PRECISION Ne,Te
      DOUBLE PRECISION alfeff,b_SH
      DOUBLE PRECISION t2
      CHARACTER*1 r,b,bidon,zc
      CHARACTER*2 d
      CHARACTER*9 tabsto
      CHARACTER*72 titre1,titre2
      CHARACTER*4 T2T(12)
c
      COMMON/MOCS1/alfeff(50,50),b_SH(50,0:50) !alfeff(n,lp) lp=l+1
      COMMON/MOCS2/inep1,iz,zc,tabsto
      COMMON/MOCS3/istor,istow
c
cc      Data T2T/'0005','0010','0030','0050','0075','0100','0125',
cc     $  '0150','0200','0300','0500','1000'/
c 2/15 Remplace data ci-dessus :
      T2T(1)='0005'
      T2T(2)='0010'
      T2T(3)='0030'
      T2T(4)='0050'
      T2T(5)='0075'
      T2T(6)='0100'
      T2T(7)='0125'
      T2T(8)='0150'
      T2T(9)='0200'
      T2T(10)='0300'
      T2T(11)='0500'
      T2T(12)='1000'
c
      r='r'
      b='b'
      d='.d'
c
c Look for closest table without interpolation.
c  The first Ne (ine=1) is 1.d2 ; closest integer in log: 
      ine=IDNINT(dlog10(Ne))-1
      If(ine.lt.1) then
       ine=1
       elseif(ine.gt.13) then
       ine=13
      Endif
c  inep1 = entier de log(Ne) pour sortie liste_phyat.dat
       inep1=ine+1
c
      if(iz.le.8) then
      t2=Te/100.
      else
c Scaling for charge > 8: Te (Ne not done)
      t2=Te/100.*(8./dfloat(iz))**2
      endif
c
      If(t2.lt.7.1d0) then
       it=1
       elseif(t2.lt.17.3d0) then
       it=2
       elseif(t2.lt.38.7d0) then
       it=3
       elseif(t2.lt.61.2d0) then
       it=4
       elseif(t2.lt.86.6d0) then
       it=5
       elseif(t2.lt.112.d0) then
       it=6
       elseif(t2.lt.137.d0) then
       it=7
       elseif(t2.lt.173.d0) then
       it=8
       elseif(t2.lt.245.d0) then
       it=9
       elseif(t2.lt.387.d0) then
       it=10
       elseif(t2.lt.707.d0) then
       it=11
       else
       it=12
      Endif
      IF(iz.eq.1) it=min0(it,10)
      IF(iz.eq.3) it=max0(it,2)
      IF(iz.ge.4.and.iz.le.6) it=max0(it,3)
      IF(iz.ge.7) it=max0(it,4)
c
      tabsto=r//zc//b//T2T(it)//d
      print *, 'reading ',tabsto
c
c Here file without '' as already alphanumeric
      OPEN(unit=istor,file=tabsto,status='old')
c
      if(ine.gt.1) then
       do i=1,ine-1
        do j=1,1189
         read(istor,15)bidon
        enddo
       enddo
      endif
       do j=1,276
        read(istor,15)bidon
       enddo
 15   format(A1)
      read(istor,16)titre1
 16   format(A72)
      do n=50,2,-1
        read(istor,30)(alfeff(n,lp),lp=1,n)
        read(istor,15)bidon
      enddo
 30   format(4x,d10.3,3x,d10.3,3x,d10.3,3x,d10.3,3x,d10.3,3x,d10.3)
c
       do j=1,263
        read(istor,15)bidon
       enddo
      read(istor,16)titre2
      do n=50,2,-1
        read(istor,30)(b_SH(n,lp),lp=1,n)
        read(istor,15)bidon
      enddo

      if(iz.gt.8) then
c Scaling for charge > 8:
      do n=50,2,-1
       do lp=1,n
        alfeff(n,lp)=alfeff(n,lp)*dfloat(iz)/8.
       enddo
      enddo
c****** 2/16 Scaling  b_SH ?? *******
cc      print *,' *** Intsto: no scaling z > 8 for b_SH, except T/z**2, z = 8 used ***'
      endif
c
      CLOSE(istor)
c
c Verif:
        IF(istow.gt.0) THEN
      OPEN(unit=istow,file='alfsto.res',status='new')
c
      write(istow,16)titre1
      do n=50,2,-1
c l = lp-1
      write(istow,31)(lp-1,alfeff(n,lp),lp=1,n)
      enddo
 31   format(i3,1p,e10.3,i3,e10.3,i3,e10.3,i3,e10.3,i3,e10.3,i3,e10.3)
c
      write(istow,16)titre2
      do n=50,2,-1
c l = lp-1
      write(istow,31)(lp-1,b_SH(n,lp),lp=1,n)
      enddo
c
      CLOSE(istow)
c
        ENDIF
c
      Return
c
      END
c---------------------end ReadSto---------------------------
c-----------------------------------------------------------
c
      Subroutine HBD
c
c H-like Branching ratios Br adapted from program BA5.2 (here HBD)
c   by Hoang-Binh, Dy 1990, A&A 238, 449.
c   Some changes/comments by Daniel Pequignot (4/4/11)
c
c    *	Program for computing exact hydrogenic oscillator strengths, 
c    *    einstein coefficients, and lifetimes
c    *	nu= principal quantum number (n) of upper state
c    *	nl= principal quantum number (n') of lower state
c    *	lu= orbital quantum number (l) of upper state
c    *	ll= orbital quantum number (l') of lower state
c    *	z= Z= nuclear charge (IDEM REAL az)
c    *	am= Mat= nuclear mass in atomic units
c    *	ain= R(n,l;n'l')=radial integral
c    *	ain2= ain**2
c    *	os= f(n',l';n,l)= absorption oscillator strength
c    *	ein= A(n,l;n',l')= Einstein coefficient
c    *	seins= A(n,n')= total Einstein coefficient
c    *	ssein= A(n)= sum of A(n,n') over n' down to n'=1
c    *	tl= 1/A(n)= lifetime of level n in s
c    *	subroutine fb5z(z,am,nl,nu,lu,ain,ain2,os,ein,einl)
c    *	subroutine fa5z(z,am,nl,nu,lu,ain,ain2,os,ein,einl)
c
	DIMENSION AS(100,100),BrS(100,100) !sum AS(nu,lu+1) & BrS(nu,lu+1)
        REAL az,am  !HBD
        DOUBLE PRECISION Br,A_tp
        INTEGER numin,numax,ibr
c
        COMMON/MOCH1/Br(100,100,100,2),A_tp(100,100,100,2) !Br(nu,lup,nl,k) lup=lu+1
        COMMON/MOCH2/az,am,numin,numax
        COMMON/MOCH3/ibr
c Local notation:
        z=az
c
          IF(ibr.gt.0) THEN
cc	iread=2
c	OPEN(file='BA5.out',unit=ibr)
c
	OPEN(unit=ibr,file='BA5.out',status='old')
c
ccc	open(file='BA5.in',unit=iread)
c 4/4/11 Rem: this nl is not used
c	read (iread,*) nu,nl,z,am
c	close(iread)
c
c1	format(2x,' lu     ','ll    R**2',
c     2   9x,'f(nl,ll;nu,lu)',3x,'A_tp(nu,lu;nl,ll)',x//)
c2	format('BA5.out'/)
3	format(' Z=',1p,e9.2,'  Mat=',1p,e9.2)
c4	format(/)
5	format('nu='i4,2x,'nl=',i4,3x, 2(1p,e11.4,3x))
6	format('A_tp(nu,nl)=', 2(1p,e11.4,3x))
c7	format('nu='i4,2x,'A_tp(nu)=',1p,e11.4,2x,'Lifetime=',1p,e11.4)
c55	format (2(i4,3x), 2(1p,e11.4,3x),3x,2(1p,e11.4,3x))
c
c	write(ibr,2)
	write(ibr,3) z,am
c	write(ibr,4)
 51    format(2x,' lu   ll  A_tp(nu,lu;nl,ll)  CASE A')
 52    format(2x,' lu   ll  A_tp(nu,lu;nl,ll)  CASE B')
 53    format(2x,' lu   ll  A_tp(nu,lu;nl,ll)  CASE C')
 54    format(2x,' lu   ll  A_tp(nu,lu;nl,ll) ***numin=',I3,' > CASE C')
        if(numin.eq.2)then
	write(ibr,51)
        elseif(numin.eq.3)then
	write(ibr,52)
        elseif(numin.eq.4)then
	write(ibr,53)
        else
	write(ibr,54) numin
        endif
          ENDIF
c
c * Stock:    Br(nu,lup,nl,k)
c * index lup=lu+1 for orbital lu (prevent index 0 in fortran)
c * k=1: lu -> ll=lu-1, k=2: lu -> ll=lu+1
c
c  Initialization:
	DO nu=1,numax
	 DO lup=1,numax
	  AS(nu,lup)=0.
	  BrS(nu,lup)=0.
	  DO nl=1,numax
	   DO k=1,2
	   A_tp(nu,lup,nl,k)=0.
	   Br(nu,lup,nl,k)=0.
	   ENDDO
	  ENDDO
	 ENDDO
	ENDDO
cccccc
	DO nu=numin,numax
cccc Anu
c	num1=nu-1
c	ssein=0.
c	
	 DO nl=numin-1,nu-1  !CASE A, B, C: numin-1 = 1, 2, 3
c	
          IF(ibr.gt.0) THEN
	write(ibr,5) nu,nl
c	write(ibr,1)
          ENDIF
c	luf=nu-1
c	sein=0.
c	
	 Do iu=1,nl
	lu=iu
	ll=lu-1 !k=1. Rem: ll just for write
	call fb5z(z,am,nl,nu,lu,ain,ain2,os,ein,einl)
c	sein=sein+einl
cc Stock: index lu+1 for lu (prevent index 0 in fortran)
	A_tp(nu,lu+1,nl,1)=ein
	AS(nu,lu+1)=AS(nu,lu+1)+ein
          IF(ibr.gt.0) THEN
c	write(ibr,55) lu,ll,ain2,os,ein
	write(ibr,56) lu,ll,ein
 56	format (i5,i5,1p,e14.4)
          ENDIF
	 Enddo
c
cc	write(ibr,4)
c
c	llf= nl-1
c	
	 Do ll = 1, nl-1
	lu=ll-1 !k=2
	call fa5z(z,am,nl,nu,lu,ain,ain2,os,ein,einl)
c 4/4/11 Rem: einl unuseful (it assumes (n,l) population = (2l+1)/n**2)
c	sein=sein+einl
	A_tp(nu,lu+1,nl,2)=ein
	AS(nu,lu+1)=AS(nu,lu+1)+ein
          IF(ibr.gt.0) THEN
cc	write(ibr,55) lu,ll,ain2,os,ein
	write(ibr,56) lu,ll,ein
          ENDIF
c
ccc 4/4/11 Rem: next 2 statements apparently unuseful
ccc	lu=ll+1
ccc	if(lu .gt. ll) go to 33
	 Enddo
c
cc	write(ibr,6) sein
cc	write(ibr,4)
c 4/4/11 Rem: sums weighted according to statistical weights
cc	ssein=ssein+sein
	 ENDDO  !nl
c
cc	tl=1./ssein
cc	write(ibr,7) nu,ssein,tl
cccccc
	ENDDO  !nu=numin,numax A
cccc Anu
c
c27/11/15 essai recalcul des AS autrement --> AS inchanges (compa BA5.out_old)
c	DO nu=numin,numax
c	DO lup=1,nu
c           AS(nu,lup)=0.
c	 DO nl=numin-1,nu-1  !CASE A, B, C: numin-1 = 1, 2, 3
c           AS(nu,lup)=AS(nu,lup)+A_tp(nu,lup,nl,2)
c           AS(nu,lup)=AS(nu,lup)+A_tp(nu,lup,nl,1)
c         ENDDO
c       ENDDO
c       ENDDO
cc*
          IF(ibr.gt.0) THEN
        write(ibr,156)
 156    format(1x,'Sum transition probabilities AS(nu,lu), lu=0,nu-1')
        DO nu=numin,numax !AS
         write(ibr,157) nu,(AS(nu,lup),lup=1,nu) !lu=0,nu-1
 157     format(i4,1p,(10e11.4))
	ENDDO  !nu=numin,numax AS
          ENDIF
cc*
        DO nu=numin,numax !Branching
cccccc Br
        k=1  ! trans lu -> ll=lu-1
         Do lup=2,nu !lup=lu+1
          DO nl=lup-1,nu-1
	va=A_tp(nu,lup,nl,k)/AS(nu,lup)
	Br(nu,lup,nl,k)=va
        BrS(nu,lup)=BrS(nu,lup)+va
	  ENDDO
	 Enddo
        k=2  ! trans lu -> ll=lu+1
         Do lup=1,nu-2 !lup=lu+1
          DO nl=lup+1,nu-1
	va=A_tp(nu,lup,nl,k)/AS(nu,lup)
	Br(nu,lup,nl,k)=va
        BrS(nu,lup)=BrS(nu,lup)+va
	  ENDDO
	 Enddo
c verif BrS:
         Do lup=1,nu !lup=lu+1
          If(nu.ne.2.or.lup.ne.1) then !2 photon emission for level 2s
           if(BrS(nu,lup).lt.1.0-1.e-6.or.BrS(nu,lup).gt.1.0+1.e-6) then
         print *,'***HBD: nu lu BrS=',nu,lup-1,BrS(nu,lup),' .ne. 1 ***'
cc	 write(ibr,58) nu,lup-1,BrS(nu,lup)
cc 58	format (' *** nu lu BrS=',i5,i5,1p,e14.4)
           endif
          Endif
	 Enddo
cccccc
	ENDDO  !nu=numin,numax Br
cccc
          IF(ibr.gt.0) THEN
	write(ibr,162)
 162    format('  nu  lu   k    Br: k=1:nl=lu,nu-1, k=2:nl=lu+2,nu-1')
        DO nu=numin,numax
         Do lup=1,nu
          if(lup.gt.1) then
         k=1
	write(ibr,62) nu,lup-1,k,(Br(nu,lup,nl,k),nl=lup-1,nu-1)
          endif
          if(lup.lt.nu-1) then
         k=2
	write(ibr,62) nu,lup-1,k,(Br(nu,lup,nl,k),nl=lup+1,nu-1)
          endif
	 Enddo
 62     format(3i4,1p,(10e11.4))
	ENDDO
c
	CLOSE(ibr)
c
          ENDIF
c
	Return
	End
c----------------------------------------
c
	Subroutine fb5z(z,am,nl,nu,lu,ain,ain2,os,ein,einl)
c    *	recurrence on b

	dimension h(1100)

	n=nu
	np=nl
	xu=nu
	xl=nl
	ll=lu-1
	lp=ll

c 4/4/11 update CODATA (but unuseful in simple precision):
c	cl=2.9979250e+10
	cl=2.9979246e+10 ! exact : 2.99792458d10
c	ryi=109737.312
	ryi=109737.316  !109737.31568539
c	em=5.4859e-04
	em=5.4857991e-04   !5.4857990946d-4 u 
	rkay=ryi/(1.+em/am)

	l=lp+1
	a01=-n+l+1.
	a02=a01-2.
	c0=2.*l

	b=a01
	e1=0.
	c = c0
	x=-4.*xu*xl/(xu-xl)**2

	h(1)=1.
	h(2)= 1.-b*x/c
	i1=(np-l)
	if(i1 .eq. 0) go to 40
	if(i1 .eq. 1) go to 41

	do 4 i=2,i1
	j=i+1
	a=-i+1.
	h1=-a*(1.-x)*h(j-2) /(a-c)
	h2=(a*(1.-x)+(a+b*x-c))*h(j-1)/(a-c)
	h(j)= h1+h2
	if(abs(h(j)) .gt. 1.e+25) go to 30
	go to 4

30	continue
	h(j)=h(j)/1.e+25
	h(j-1)=h(j-1)/1.e+25
	e1=e1+25.

4	continue
	p1=h(i1+1)
	go to 50

40	continue
	p1=h(1)
	go to 50

41	continue
	p1=h(2)
	go to 50

50	continue
	b=a02
	e2=0.
	h(1)=1.
	h(2)= 1.-b*x/c
	i1=(np-l)
	if(i1 .eq. 0) go to 42
	if(i1 .eq. 1) go to 43

	do 5 i=2,i1
	j=i+1
	a=-i+1.
	h1=-a*(1.-x)*h(j-2) /(a-c)
	h2=(a*(1.-x)+(a+b*x-c))*h(j-1)/(a-c)
	h(j)= h1+h2
	if(abs(h(j)) .gt. 1.e+25) go to 31
	go to 5

31	continue
	h(j)=h(j)/1.e+25
	h(j-1)=h(j-1)/1.e+25
	e2=e2+25.

5	continue
	p2=h(i1+1)

	go to 51

42	continue
	p2=h(1)
	go to 51

43	continue
	p2=h(2)
	go to 51

51	continue
	cc4=n-np
	cc5=n+np
	ff= p1*(1.-cc4**2/cc5**2*p2/p1*10.**(e2-e1))
	alof=alog10(abs(ff))+e1

c    *	cal of c1, c2, c3, c4, c5

	i2=(2*l-1)
	s1=0.

	do 6 i=1,i2
	ai=i
	s1=s1+ alog10(ai)

6	continue

	c1= - (alog10(4.)+s1)

	s=0.

	i3=(n+l)
	si3=0.

	do 7 i=1,i3
	ai=i
	si3=si3+ alog10(ai)

7	continue
	s=s+si3

	i4=(np+l-1)
	si4=0.

	do 8 i=1,i4
	ai=i
	si4=si4+ alog10(ai)

8	continue
	s=s+si4

	i5=n-l-1
	si5=0.
	if(i5 .eq. 0) go to 2

	do 9 i=1,i5
	ai=i
	si5=si5+ alog10(ai)

9	continue
	s=s-si5

2	continue
	i6=i1
	si6=0.
	if(i6 .eq. 0)go to 3

	do 12 i=1,i6
	ai=i
	si6=si6+ alog10(ai)

12	continue
	 s=s-si6

3	continue
	c2=s/2.

	cc3=4.*n*np
	c3= (l+1.)*alog10(cc3)

	cc4=n-np
	c4= (n+np-2.*l-2.)*alog10(cc4)

	cc5=n+np
	c5=(-n-np)*alog10(cc5)

	c=c1+c2+c3+c4+c5

	ali =alof+c
	ain = 10.**ali/z
	ain2=ain**2

	xu=nu
	xl=nl
	euplo=rkay*z**2/xu**2/xl**2*(xu-xl)*(xu+xl)
	os=1./3.*(xu+xl)*(xu-xl)/(xu*xl)**2*max(lu,ll)/
     2	(2.*ll+1.)*ain2*z**2
c 4/4/11 Q: origin of this numerical value 0.66704?
	ein=0.66704*(2.*ll+1.)/(2.*lu+1.)*(euplo)**2*os
	einl=ein*(2.*lu+1.)/xu**2

	return
	end
c------------------------------------------------------
c
	Subroutine fa5z(z,am,nl,nu,lu,ain,ain2,os,ein,einl)
c    *	recurrence on a

	dimension h(1100)
c 4/4/11 update CODATA (but unuseful in simple precision):
c	cl=2.9979250e+10
	cl=2.9979246e+10 ! exact : 2.99792458d10
c	ryi=109737.312
	ryi=109737.316  !109737.31568539
c	em=5.4859e-04
	em=5.4857991e-04 !5.4857990946d-4 u 
	rkay=ryi/(1.+em/am)

	xu=nu
	xl=nl
	x=-4.*xu*xl/(xu-xl)**2
	n=nl
	ll=lu+1
	l=ll
	b0= 0.-nu+l+0.
	c0=2.*l
	b=b0
	c = c0
	cc4=1.*(n-nu)

	h(1)=1.
	h(2)= 1.-b*x/c
	e1=0.
	e2=e1
	i1=(n-l+1)

	do 4 i=2,i1
	j=i+1
	a=1.-i+0.
	h1=-a*(1.-x)*h(j-2) /(a-c)
	h2=(a*(1.-x)+(a+b*x-c))*h(j-1)/(a-c)
	h(j)= h1+h2
	if(abs(h(j)) .gt. 1.e+25) go to 30
	go to 4

30	continue
	h(j)=h(j)/1.e+25
	h(j-1)=h(j-1)/1.e+25
	h(j-2)=h(j-2)/1.e+25
	e1=e1+25.
	e2=e1

4	continue
	p1=h(i1-1)
	p2 = h(i1+1)

	cc4=1.*(n-nu)
	cc4= abs(cc4)
	cc5=n+nu
	ff= p1*(1.-cc4**2/cc5**2*p2/p1*10.**(e2-e1))
	alof=alog10(abs(ff))+e1

	i2=(2*l-1)
	s1=0.

	do 6 i=1,i2
	ai=i
	s1=s1+ alog10(ai)
6	continue

	c1=-(alog10(4.)+s1)
 	s=0.

	i3=(n+l)
	si3=0.

	do 7 i=1,i3
	ai=i
	si3=si3+ alog10(ai)
7	continue

	s=s+si3

	i4=(nu+l-1)
	si4=0.

	do 8 i=1,i4
	ai=i
	si4=si4+ alog10(ai)
8	continue

	s=s+si4

	i5=n-l-1
	si5=0.
	if(i5 .eq. 0) go to 2

	do 9 i=1,i5
	ai=i
	si5=si5+ alog10(ai)
9	continue

	s=s-si5

2	continue

	i6=nu-l
	si6=0.
	if(i6 .eq. 0) go to 3

	do 12 i=1,i6
	ai=i
	si6=si6+ alog10(ai)
12	continue

	s=s-si6

3	continue
	c2=s/2.

	cc3=4.*n*nu
	c3= (l+1.)*alog10(cc3)

	cc4=cc4
	c4= (n+nu-2.*l-2.)*alog10(cc4)

	cc5=n+nu
	c5=(-n-nu)*alog10(cc5)

	c=c1+c2+c3+c4+c5

	ali =alof+c
	ain = 10.**ali/z
	ain2=ain**2

	xu=nu
	xl=nl
	euplo=rkay*z**2/xu**2/xl**2*(xu-xl)*(xu+xl)
	os=1./3.*(xu+xl)*(xu-xl)/(xu*xl)**2*max(lu,ll)/
     2	(2.*ll+1.)*ain2*z**2
c 4/4/11 Q: origin of this numerical value 0.66704?
	ein=0.66704*(2.*ll+1.)/(2.*lu+1.)*(euplo)**2*os
	einl=ein*(2.*lu+1.)/xu**2

	return
	end	
c----------------------------------------------------
c	END OF BA5
c----------------------------------------------------
c************************************************************
c**** Read + Interpolate H and He+ Storey (4 subr) ****

      subroutine hdatx
c
c  subroutine to read e1bx.d
c
      implicit none
      real*8 densx,tempx,ex,a,a2s,r
      integer ntempx,ndensx,ia,ib,j,i,ne,ntop,ndum,nlu,nll

      common/hdatax/densx,tempx,ex,a,a2s,ntempx,ndensx,ntop,nll,nlu
      
      dimension densx(15),tempx(15),ex(5000,15,15),r(15,15),a(15,15),
     $          a2s(15,15)

      write(6,*) 'reading dataset e1bx.d'
c
      OPEN(unit=15,file='e1bx.d',status='old')
c
c
      read(15,*) ntempx,ndensx
      print *, ntempx,ndensx
      do 101 ia=1,ntempx
           do 100 ib=1,ndensx
c                write(6,*) ia,ib
                read(15,25) densx(ib),tempx(ia),ntop,ndum,nlu,nll
25              format(1x,e10.3,5x,e10.3,5x,4i5)
c                 write(6,25) densx(ib),tempx(ia),ntop,ndum,nlu,nll
                ne=(2*ntop-nlu-nll)*(nlu-nll+1)/2
c                write(6,*) 'ne= ',ne
                read(15,30) (ex(j,ia,ib),j=1,ne)
30              format((8e10.3))
100        continue  
101   continue
c      read(15,*) ((a(i,j),i=1,ndens),j=1,ntemp)
c      read(15,*) ((a2s(i,j),i=1,ndens),j=1,ntemp)
c
      CLOSE(15)
c
      return
      end
c
c**********************************************************************
c
      subroutine hepdatx
c
c  subroutine to read e2bx.d
c
      implicit none
      real*8 densx,tempx,ex,a,a2s,r
      integer ntop,ndum,nlu,nll,ntempx,ndensx,ia,ib,j,i,ne

      common/hepdatax/densx,tempx,ex,a,a2s,ntempx,ndensx,ntop,nll,nlu
      
      dimension densx(15),tempx(15),ex(5000,15,15),r(15,15),a(15,15),
     $          a2s(15,15)

      write(6,*) 'reading dataset e2bx.d'
c
      OPEN(unit=15,file='e2bx.d',status='old')
c
c
      read(15,*) ntempx,ndensx
      do 101 ia=1,ntempx
           do 100 ib=1,ndensx
c                write(6,*) ia,ib
                read(15,25) densx(ib),tempx(ia),ntop,ndum,nlu,nll
25              format(1x,e10.3,5x,e10.3,5x,4i5)
c                 write(6,25) densx(ib),tempx(ia),ntop,ndum,nlu,nll
                ne=(2*ntop-nlu-nll)*(nlu-nll+1)/2
c                write(6,*) 'ne= ',ne
                read(15,30) (ex(j,ia,ib),j=1,ne)
30              format((8e10.3))
100        continue  
101   continue
c      read(15,*) ((a(i,j),i=1,ndens),j=1,ntemp)
c      read(15,*) ((a2s(i,j),i=1,ndens),j=1,ntemp)
c
      CLOSE(15)
c
      return
      end
c
c*********************************************************************
c
      subroutine hlinex(nu,nl,xt,xd,fh,iopt)

c
c  Interpolate in density/temperature for specified line emissivity (iopt=1)
c  or 2s recombination coefficient (iopt=2)
c
      implicit real*8(a-h,o-z)

      common/hdatax/densx,tempx,ex,a,a2s,ntempx,ndensx,ntop,nll,nlu
      
      dimension densx(15),tempx(15),ex(5000,15,15),r(15,15),x(15),y(15),
     $          ni(5),cx(5),cy(5),ri(5),f(2,15),a(15,15),a2s(15,15)

      data max/4/ni/2,3,4,5,6/       ! interpolation parameters
c
c          interpolation variables
c
      do 102 i=1,ndensx
           x(i)=log10(densx(i))
102   continue
      do 103  i=1,ntempx
           y(i)=sqrt(tempx(i))
           f(1,i)=1.0          ! f is emissivity smoothing function in temp 
           f(2,i)=y(i)
103   continue

c      if(iopt.eq.2) then
c           do 106 it=1,ntemp
c                do 105 id=1,ndens
c                     r(it,id)=a2s(id,it)
c105             continue 
c106        continue
c      ns=1
c      goto 109

c      endif 

c
      nus=0
      nls=0
c
c          set keys to locate transitions of interest
c
      if((nus+nls).eq.0) then
           ns=2
           ks=999
      else
           ns=1
           ks=(((ncut-nus)*(ncut+nus-1))/2)+nls
      endif
      k=(2*ntop-nll-nl+1)*(nl-nll)/2+ntop-nu+1
c      write(6,*) 'hlinex, k=',k
c      k=(((ncut-nu)*(ncut+nu-1))/2)+nl
c
c          calculate desired intensity ratio (or emissivity if nus=nls=0)
c
c     write(6,*) 'calc r k,ks,ntemp,ndens,ncut',k,ks,ntemp,ndens,ncut
      do 108 it=1,ntempx
           do 107 id=1,ndensx
                if(ns.eq.1) then
                     r(it,id)=ex(k,it,id)/ex(ks,it,id)
                else
                     r(it,id)=ex(k,it,id)
                endif
107        continue
108   continue
c
c          interpolate in r-table
c
 109  nt=1

      xkeep=xd

      if(xt.lt.tempx(1).or.xt.gt.tempx(ntempx).or.xd.lt.densx(1).
     &     or.xd.gt.densx(ndensx)) then
         write(*,85)
85       format(' requested temp/dens not in table')
         write(6,*) 'xt=',xt,'  xd=',xd
         write(6,*) ntempx,tempx(1),tempx(ntempx)
         write(6,*) ndensx,densx(1),densx(ndensx)
ctemp
         write(6,*) 'xd, densx(1)',xd,densx(1)
         if(xd.lt.densx(1)) then
           xd=densx(1)
           write(6,*) 'carry on with xd =',densx(1)
           go to 555
         endif
ctemp
         stop
      endif

 555  continue

      xp=log10(xd)           ! interpolate in log(dens)
      yp=sqrt(xt)             ! interpolate in temp**0.5
c      write(6,*) 'xp,yp',xp,yp
c
c          find interpolation box 
c
      i=1
86    if(xp.ge.x(i).and.xp.le.x(i+1)) then
          goto 88
      else
          i=i+1
          if(i.eq.ndensx) then
                     stop 'dens overflow'
          endif
          goto 86
      endif
88    i0=i
      j=1
90    if(yp.ge.y(j).and.yp.le.y(j+1)) then
         goto 92
      else
         j=j+1
         if(j.eq.ntempx) then
              stop 'temp overflow'
         endif
         goto 90
      endif
92    j0=j
c
c          interpolate to orders 2,3,4,5 in both directions
c
      int=max
      nint=ni(int)             ! interpolation order
      nint1=nint-1
      nof=nint1/2
c
c          shift i0 to nearest box boundary in each direction if nint is odd
c
      if(nint.eq.3.or.nint.eq.5.or.nint.eq.7) then  ! note ODD order
          if((xp-x(i0)).gt.(x(i0+1)-xp)) then
              is=i0+1-nof
          else
              is=i0-nof
          endif
          if((yp-y(j0)).gt.(y(j0+1)-yp)) then
              js=j0+1-nof
          else
              js=j0-nof
          endif
      else
          is=i0-nof
          js=j0-nof
      endif
c
c          ensure that interpolation box lies in table
c
      if(is.lt.1) then
          is=1
      endif
      if((is+nint1).gt.ndensx) then
          is=ndensx-nint1
      endif
      if(js.lt.1) then
          js=1
      endif
      if((js+nint1).gt.ntempx) then
          js=ntempx-nint1
      endif
c
c          nint**2-point interpolation
c
      do 111 k=1,nint
           i=is+k-1
           cx(k)=1.0
           do 110 kp=1,nint
              if(kp.ne.k) then
                 ip=is+kp-1
                 cx(k)=cx(k)*(xp-x(ip))/(x(i)-x(ip))
              endif
110        continue
111   continue
      do 113 k=1,nint
           j=js+k-1
           cy(k)=1.0
           do 112 kp=1,nint
              if(kp.ne.k) then
                 jp=js+kp-1
                 cy(k)=cy(k)*(yp-y(jp))/(y(j)-y(jp))
              endif
112        continue
113   continue
      rint=0.0
      do 115 kx=1,nint
         do 114 ky=1,nint
            if((js+ky-1).gt.ntempx.or.(is+kx-1).gt.ndensx)
     .                         then
               stop 'final loop error'
            endif
            rrr=r(js+ky-1,is+kx-1)*f(ns,js+ky-1) ! smoothing ftn  
            if(nt.ne.0) then
               rrr=log(rrr)
            endif
            rint=rint+cx(kx)*cy(ky)*rrr
114      continue
115   continue
      ri(int)=rint
      if(nt.ne.0) then
         ri(int)=exp(ri(int))
      endif
      if(ns.eq.2) then
         ri(int)=ri(int)/yp ! remove smoothing function = temp**.5
      endif
      fh=ri(max)

c restore xd if necessary

      if(xd.ne.xkeep) xd=xkeep

      return
      end
c
c***********************************************************************
c
      subroutine heplinex(nu,nl,xt,xd,fh,iopt)
c
c  Interpolate in density/temperature for specified line emissivity (iopt=1)
c  or 2s recombination coefficient (iopt=2)
c
      implicit real*8(a-h,o-z)

      common/hepdatax/densx,tempx,ex,a,a2s,ntempx,ndensx,ntop,nll,nlu
      
      dimension densx(15),tempx(15),ex(5000,15,15),r(15,15),x(15),y(15),
     $          ni(5),cx(5),cy(5),ri(5),f(2,15),a(15,15),a2s(15,15)

      data max/4/ni/2,3,4,5,6/       ! interpolation parameters
c
c          interpolation variables
c
      do 102 i=1,ndensx
           x(i)=log10(densx(i))
102   continue
      do 103  i=1,ntempx
           y(i)=sqrt(tempx(i))
           f(1,i)=1.0          ! f is emissivity smoothing function in temp 
           f(2,i)=y(i)
103   continue

      if(iopt.eq.2) then
           do 106 it=1,ntempx
                do 105 id=1,ndensx
                     r(it,id)=a2s(id,it)
105             continue 
106        continue
      ns=1
      goto 109

      endif 


c
      nus=0
      nls=0

c
c          set keys to locate transitions of interest
c
      if((nus+nls).eq.0) then
           ns=2
           ks=999
      else
           ns=1
           ks=(((ncut-nus)*(ncut+nus-1))/2)+nls
      endif
      k=(2*ntop-nll-nl+1)*(nl-nll)/2+ntop-nu+1
c
c          calculate desired intensity ratio (or emissivity if nus=nls=0)
c
c     write(6,*) 'calc r k,ks,ntemp,ndens,ncut',k,ks,ntemp,ndens,ncut
      do 108 it=1,ntempx
           do 107 id=1,ndensx
                if(ns.eq.1) then
                     r(it,id)=ex(k,it,id)/ex(ks,it,id)
                else
                     r(it,id)=ex(k,it,id)
                endif
107        continue
108   continue
c
c          interpolate in r-table
c
 109  nt=1

      xkeep=xd

      if(xt.lt.tempx(1).or.xt.gt.tempx(ntempx).or.xd.lt.densx(1).
     &     or.xd.gt.densx(ndensx)) then
         write(*,85)
85       format(' requested temp/dens not in table')
         write(6,*) 'xt=',xt,'  xd=',xd
         write(6,*) ntempx,tempx(1),tempx(ntempx)
         write(6,*) ndensx,densx(1),densx(ndensx)
ctemp
         write(6,*) 'xd, densx(1)',xd,densx(1)
         if(xd.lt.densx(1)) then
           xd=densx(1)
           write(6,*) 'carry on with xd =',densx(1)
           go to 555
         endif
ctemp
         stop
      endif

 555  continue

      xp=log10(xd)           ! interpolate in log(dens)
      yp=sqrt(xt)             ! interpolate in temp**0.5
*      write(6,*) 'xp,yp',xp,yp
c
c          find interpolation box 
c
      i=1
86    if(xp.ge.x(i).and.xp.le.x(i+1)) then
          goto 88
      else
          i=i+1
          if(i.eq.ndensx) then
                     stop 'dens overflow'
          endif
          goto 86
      endif
88    i0=i
      j=1
90    if(yp.ge.y(j).and.yp.le.y(j+1)) then
         goto 92
      else
         j=j+1
         if(j.eq.ntempx) then
              stop 'temp overflow'
         endif
         goto 90
      endif
92    j0=j
c
c          interpolate to orders 2,3,4,5 in both directions
c
      int=max
      nint=ni(int)             ! interpolation order
      nint1=nint-1
      nof=nint1/2
c
c          shift i0 to nearest box boundary in each direction if nint is odd
c
      if(nint.eq.3.or.nint.eq.5.or.nint.eq.7) then  ! note ODD order
          if((xp-x(i0)).gt.(x(i0+1)-xp)) then
              is=i0+1-nof
          else
              is=i0-nof
          endif
          if((yp-y(j0)).gt.(y(j0+1)-yp)) then
              js=j0+1-nof
          else
              js=j0-nof
          endif
      else
          is=i0-nof
          js=j0-nof
      endif
c
c          ensure that interpolation box lies in table
c
      if(is.lt.1) then
          is=1
      endif
      if((is+nint1).gt.ndensx) then
          is=ndensx-nint1
      endif
      if(js.lt.1) then
          js=1
      endif
      if((js+nint1).gt.ntempx) then
          js=ntempx-nint1
      endif
c
c          nint**2-point interpolation
c
      do 111 k=1,nint
           i=is+k-1
           cx(k)=1.0
           do 110 kp=1,nint
              if(kp.ne.k) then
                 ip=is+kp-1
                 cx(k)=cx(k)*(xp-x(ip))/(x(i)-x(ip))
              endif
110        continue
111   continue
      do 113 k=1,nint
           j=js+k-1
           cy(k)=1.0
           do 112 kp=1,nint
              if(kp.ne.k) then
                 jp=js+kp-1
                 cy(k)=cy(k)*(yp-y(jp))/(y(j)-y(jp))
              endif
112        continue
113   continue
      rint=0.0
      do 115 kx=1,nint
         do 114 ky=1,nint
            if((js+ky-1).gt.ntempx.or.(is+kx-1).gt.ndensx)
     .                         then
               stop 'final loop error'
            endif
            rrr=r(js+ky-1,is+kx-1)*f(ns,js+ky-1) ! smoothing ftn  
            if(nt.ne.0) then
               rrr=log(rrr)
            endif
            rint=rint+cx(kx)*cy(ky)*rrr
114      continue
115   continue
      ri(int)=rint
      if(nt.ne.0) then
         ri(int)=exp(ri(int))
      endif
      if(ns.eq.2) then
         ri(int)=ri(int)/yp ! remove smoothing function = temp**.5
      endif
      fh=ri(max)

c restore xd if necessary

      if(xd.ne.xkeep) xd=xkeep

      return
      end
c
c*********************************************************************
c
c---------------------------------------------------------------------
