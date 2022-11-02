;; 1. Based on: 
; Model desc: Isoniazid.Model.KDHTB
; Settings for the memory of NONMEM
$SIZES      PD=-1000 LVR=-150 LTH=-200 MAXFCN=10000000 LNP4=-150000
$PROBLEM    |ADVAN13_TRANS1_PZA_DATA|

$INPUT      ID DAT2=DROP TIME OCC WHAT=DROP 
            EVID_RIF AMT_RIF RIFCONC RIF_BLQ MDV_RIF 
			EVID AMT DV INH_BLQ MDV 
            EVID_PZA AMT_PZA PZACONC PZA_BLQ MDV_PZA 
            AMT_EMB
            VPC_TIME
            PATIENT HOSPITALIZED DIED 
			AGE SEXF WT HT FFM FAT BRC 

; DEFINE ETAS FOR OCCASSION
$ABBREVIATED REPLACE ETA(OCC_BIO)=ETA(,7 to 8 by 1)
$ABBREVIATED REPLACE ETA(OCC_KA)=ETA(,9 to 10 by 1)
$ABBREVIATED REPLACE ETA(OCC_MTT)=ETA(,12 to 13 by 1)

; IGNORE=@ will skip any line starting with any non-numerical character
$DATA      alldat_final_simulated.csv IGNORE=@

$SUBROUTINE ADVAN5 TRANS1 ;TOL=9 ATOL=9 ; TOL is the precision to solve differential equations

;Sim_start
$PRIOR      NWPRI NPEXP=1 PLEV=0.9999
;Sim_end

$MODEL      NCOMPS=10 ; NUMBER OF COMPARTMENTS (ABSORPTION COMPATMENT (DEFINED AS FIRST ONE) AND CENTRAL COMPARTMENT DEFIEND AS 2ND COMPARTMENT
            COMP=(TRANSIT1,DEFDOSE) ;1 GUT TRANIST 1 (F1 is associated with first compartment)
            COMP=(TRANSIT2) ;2 GUT TRANIST 2
            COMP=(TRANSIT3) ;3 GUT TRANIST 3
            COMP=(TRANSIT4) ;4 GUT TRANIST 4
            COMP=(TRANSIT5) ;5 GUT TRANIST 5
            COMP=(TRANSIT6) ;6 GUT TRANIST 6
            COMP=(TRANSIT7) ;7 GUT TRANIST 7
            COMP=(ABS) ;8 GUT ABS
            COMP=("CENTRAL",DEFOBS) ;9 CENTRAL CMT
            COMP=(PERI) ; 10 PERIPHERAL CMT
;-----------------------------------------------------------------------------------------------------
;Sim_start
$MIX
NSPOP=2
P(1) = THETA(3)
P(2) = 1 - THETA(3)
;Sim_end

$PK
;-----------------------------------------Allometric scaling----------------------------------------
TVFFM = 43.5
TVWT = 56
TVFAT = 11.3

ALLMCLWT=(WT/TVWT)**0.75
ALLMVWT=WT/TVWT

ALLMCLFFM=(FFM/TVFFM)**0.75
ALLMVFFM=FFM/TVFFM

ALLMCLFAT=(FAT/TVFAT)**0.75
ALLMVFAT=FAT/TVFAT

;-------------------------------For Mixture Modelling-----------------------------------------------
;Sim_start
EST = MIXEST
;EST = POP
;DUMMY_THETA = THETA(3)
IF(MIXNUM.EQ.1) THEN
      TVCL = THETA(4) *ALLMCLFFM
ELSE
      TVCL = THETA(5) *ALLMCLFFM
ENDIF
;IF(POP.EQ.1) THEN
	 ;   TVCL=THETA(4)
	 ;ELSE
	 ;   TVCL=THETA(5)
	 ;ENDIF
;Sim_end
;----------------------------Imputing missing values to median--------------------------------------
; LACTATE_IMP = LACTATE
; IF (LACTATE.EQ.-99) LACTATE_IMP = 1.6

; TIMETOBREAKFAST_IMP = TIMETOBREAKFAST
; IF (TIMETOBREAKFAST.EQ.-999) TIMETOBREAKFAST_IMP = -4.82

; CRP_IMP = CRP
; IF (CRP.EQ.-99) CRP_IMP = 127

; AST_IMP = AST
; IF (AST.EQ.-99) AST_IMP = 37.5

; ALT_IMP = ALT
; IF (ALT.EQ.-99) ALT_IMP = 24

; TOTALPROTEIN_IMP = TOTALPROTEIN
; IF (TOTALPROTEIN.EQ.-999) TOTALPROTEIN_IMP = 80

; ALBUMIN_IMP = ALBUMIN
; IF (ALBUMIN.EQ.-99) ALBUMIN_IMP = 30

; UREA_IMP = UREA
; IF (UREA.EQ.-99) UREA_IMP = 5

; CREATININE_IMP = CREATININE
; IF (CREATININE.EQ.-99) CREATININE_IMP = 68

; CD4_IMP = CD4
; IF (CD4.EQ.-99) CD4_IMP = 71.5

;--------------------------Hospitalized Effect on Variability in CL---------------------------------
; Allometry for liver-------------------------------------------
ALLMCL_WT_HEP = (WT/TVWT)**0.75
ALLMV_WT_HEP = (WT/TVWT)

ALLMCL_FFM_HEP = (WT/TVFFM)**0.75
ALLMV_FFM_HEP = (WT/TVFFM)
;--------------------Typical values (thetas & etas)-------------------------------------------------
TVVRATIO = EXP(THETA(1))
TVQ = EXP(THETA(2))*ALLMCLFFM							 
;population parameters
;TVCL        = THETA(1)
TVV         = THETA(6) *ALLMVFFM
TVBIO       = THETA(7)
TVKA        = THETA(8)
TVMTT       = THETA(9)
TVNN        = THETA(10)
TVV3 = (TVV/TVVRATIO)*ALLMVFFM							  

;BSV
BSVCL       = ETA(1)
BSVV        = ETA(2)
BSVBIO 		= ETA(3)
BSVKA 		= ETA(4)
BSVQ 		=ETA(5)
BSVV3 		=ETA(6)

;BOV
BOVBIO 	= ETA(OCC_BIO) ;7-8
BOVKA       = ETA(OCC_KA) ;9-10
BSVMTT 	= ETA(11)
BOVMTT  = ETA(OCC_MTT) ;12-13

;-----------Individual parameters-------------------------------------------------------------------
VRATIO = TVVRATIO ;						  
CL    = TVCL*EXP(BSVCL)                   ; Clearance
V2     = TVV*EXP(BSVV)                     ; CENTRAL VOL.
BIO   = TVBIO*EXP(BSVBIO + BOVBIO)        ; BIOAVAILABILITY
KA    = TVKA*EXP(BSVKA + BOVKA) ; ABS. RATE
Q	= TVQ*EXP(BSVQ)
V3	= TVV3*EXP(BSVV3)
;LAG   = TVLAG*EXP(BSVLAG + BOVLAG) ; LAG TIME
MTT   = TVMTT*EXP(BSVMTT + BOVMTT) ; MEAN TRANSIT TIME
NN 	= TVNN
;---------------------------------------------------------------------------------------------------
VARCL       = BSVCL ;+ BOVCL
VARBIO      = BSVBIO + BOVBIO
VARAUC      = BSVBIO + BOVBIO - BSVCL ;- BOVCL
VARABS      = BSVKA + BOVKA - BSVMTT - BOVMTT
VARKA       = BSVKA + BOVKA
VARMTT      = BSVMTT + BOVMTT
;===================================================================================================
;--------------------------HEPATIC CL----------------------------------------;
TVQH=THETA(11)*ALLMCL_WT_HEP   ; PLASMA FLOW RATE
TVFU=THETA(12)            ; UNBOUND PLASMA FRACTION OF INH

QH=TVQH
FU=TVFU
CLINT=CL

;--------  Transfer constants for liver model ---------

; define hepatic extraction
EH  = (CLINT*FU)/((CLINT*FU)+QH) ; fraction undergoing first pass extraction
FH  = 1 - EH ; fraction available after 1st pass to go to systemic circulation																		  

;------------------------------------Re-parameterization--------------------------------------------
;===================================================================================================
F1 = BIO
KTR = (NN+1)/MTT
K12	=	KTR  ;Rate between transit CMT
K23	=	KTR   ;Rate between transit CMT
K34	=	KTR   ;Rate between transit CMT
K45	=	KTR   ;Rate between transit CMT
K56	=	KTR   ;Rate between transit CMT
K67	=	KTR   ;Rate between transit CMT
K78	=	KTR   ;Rate between transit CMT
K89	=	KA   ;Rate between transit CMT

K90 = CL/V2 ;(rate constant of elimination)
K9T10 = Q/V2
K10T9 = Q/V3

;----------------------------------------Error model------------------------------------------------
$ERROR
IPRED = A(9)/V2
IRES = DV-IPRED

LLOQ = 0.105

PROP = IPRED*THETA(13)
ADD = 0.2*LLOQ + THETA(14)

; For ADD, in this case we are coding THETA(.) as the additive error on top of 20% of the LLOQ.
; So the lower bound of THETA(.) can be zero. If it goes to zero, we can fix it, and the additive error
; will then be constrained to 20% of LLOQ + the value in THETA(.).
; REMEMBER about this when you report the value of ADD and its uncertainty! NONMEM gives uncertainty on THETA, not ADD
; An alternative approach is to set the lower bound of the THETA for the additive error to 20% of the LLOQ.
; In that case, one does not have to worry about adjusting the precision. On the other hand, this cannot be done if you have different LLOQs within your analysis (e.g. different labs)

; For BLQ==1 (i.e. first BLQ value in a series), we add extra additive error on the concentrations, since the value in DV has been imputed
IF(ICALL.NE.4.AND.INH_BLQ==1) ADD = ADD + (LLOQ/2)

; For BLQ==2 (i.e. the trailing BLQ values in a series), we don't want these to influence the fit,
; we only want them for simulation-based diagnostics such as the VPC.
; So we define a separate error structure for these points. It has no proportional component
; (PROP = 0, as we would not want these points to affect our estimate of proportional error)
; and a FIXED and HUGE additive component (ADD = 1000000000, large with respect to the readings of concentration),
; so that the values do not affect the fit.
; It's also a good idea to repeat the diagnostic plots without the BLQ=2 points
IF(ICALL.NE.4.AND.INH_BLQ==2) THEN
	PROP = 0
	ADD = 100000000
ENDIF

W = SQRT(ADD**2+PROP**2)

IF (W.LE.0.000001) W=0.000001 ; to protect IWRES from overflow

IWRES = IRES/W

Y = IPRED + W*ERR(1)

; To prevent simulation (ICALL==4) of negative values. It set a positive lower bound for Y, so that VPCs in the log-scale can be plotted
IF (ICALL==4.AND.Y<=LLOQ) THEN 
	Y=LLOQ/2
	INH_BLQ = 1
ENDIF

; To calculate time after dose.
IF(AMT.GT.0) THEN
	TIMEDOSE = TIME
	AMOUNTDOSE = AMT
ENDIF

TAD = TIME-TIMEDOSE

;------------------------------------------RETRIEVE AMOUNT IN EACH COMPARTMENT---------------------------------------------------------------------------------------
AA1 = A(8)
AA2 = A(9)
AA3 = A(10)
;--------------------------------------------------------------------------------------------------------------------------------------------------------------------
$THETA  (-10,0.653065,10) ; 1 VRATIO [LOG]
$THETA  (-10,0.35836,10) ; 2 Q [LOG]
;Sim_start
$THETA  (0,0.645125,1) ; 3 PROB_FAST
;$THETA  0 FIX ; 3 PROB_FAST
;Sim_end
$THETA  (0,25.5488,100) ; 4 CL_FAST [L/h]
$THETA  (0,9.76422,40) ; 5 CL_SLOW [L/h]
$THETA  (0,59.0219) ; 6 V [L]
$THETA  1 FIX ; 7 BIO
$THETA  (0,2.42706) ; 8 KA [1/h]
;$THETA  (0,2.20652,50) ; 7 Q
;$THETA  (0,82.0861) ; 8 V3
$THETA  (0,0.441855,3) ; 9 MTT [h]
$THETA  7 FIX ; 10 NN
$THETA  0 FIX ; 11 QH
$THETA  0 FIX ; 12 FU
$THETA  (0,0.138966) ; 13 PROP [%]
$THETA  0 FIX ; 14 ADD [mg/L]
;$THETA  (-5,0.728585,5) ; 15 BXPAR
; PRIORS FOR VRATIO ---------------------------------------------------------
;Sim_start
$THETAP  0.703 FIX ; logVRATIO
;Sim_end
;------------------------------------------------------------------------------------------------------------------------------------------------------
$OMEGA  0.0642373  ;   1 BSV CL
$OMEGA  0  FIX  ;    2 BSV V
$OMEGA  0  FIX  ;  3 BSV BIO
$OMEGA  0  FIX  ;   4 BSV KA
$OMEGA  0  FIX  ;     5 BSVQ
$OMEGA  0  FIX  ;    6 BSVV3
;---------------------------------------------------------------------------------------------------------------------------------------------------------------------
$OMEGA  BLOCK(1)
 0.122052  ; 7,8 BOVBIO
$OMEGA  BLOCK(1) SAME
;---------------------------------------------------------------------------------------------------------------------------------------------------------------------
$OMEGA  BLOCK(1)
 1.48417  ; 9,10 BOVKA
$OMEGA  BLOCK(1) SAME
;---------------------------------------------------------------------------------------------------------------------------------------------------------------------
$OMEGA  0  FIX  ;  11 BSVMTT
$OMEGA  BLOCK(1)
 0.994496  ; 12,13 BOVMTT
$OMEGA  BLOCK(1) SAME
;------------------------------------------------------------------------------------------------------------------------------------------------------
; PRIORS DATA FOR THETAS VRATIO
;Sim_start
$THETAPV  BLOCK(1) FIX
 0.1  ;     VRATIO
;Sim_end
;---------------------------------------------------------------------------------------------------------------------------------------------------------------------
$SIGMA  1  FIX
;-------------------------------------------------------------------------------------------------------------------------------------------------------
;Sim_start
$ESTIMATION MSFO=INH.msf MAXEVAL=9999 PRINT=1 METHOD=1
            INTERACTION NOABORT NSIG=3 NONINFETA=1 ETASTYPE=1 SIGL=9

$COVARIANCE PRINT=E UNCONDITIONAL MATRIX=R PRECOND=1

;Sim_end

