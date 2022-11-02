;; 1. Based on: 
; Model desc: Rifampicin.Model.KDHTB
; Settings for the memory of NONMEM
$SIZES      PD=-1000 LVR=-150 LTH=-200 MAXFCN=10000000 LNP4=-150000
$PROBLEM    |ADVAN2_TRANS1_RIF_DATA|
$INPUT      ID DAT2=DROP TIME OCC WHAT=DROP 
            EVID AMT DV RIF_BLQ MDV 
			EVID_INH AMT_INH INHCONC INH_BLQ MDV_INH 
            EVID_PZA AMT_PZA PZACONC PZA_BLQ MDV_PZA 
            AMT_EMB
            VPC_TIME
            PATIENT HOSPITALIZED DIED 
			AGE SEXF WT HT FFM FAT BRC

; DEFINE ETAS FOR OCCASSION
$ABBREVIATED REPLACE ETA(OCC_BIO)=ETA(,6 to 7 by 1)
$ABBREVIATED REPLACE ETA(OCC_KA)=ETA(,8 to 9 by 1)
$ABBREVIATED REPLACE ETA(OCC_MTT)=ETA(,10 to 11 by 1)

; IGNORE=@ will skip any line starting with any non-numerical character			
$DATA      alldat_final_simulated.csv IGNORE=@

$SUBROUTINE ADVAN5 TRANS1 ;TOL=8 ATOL=8 ; TOL is the precision to solve differential equations

$MODEL      NCOMPS=13 ; NUMBER OF COMPARTMENTS (ABSORPTION COMPATMENT (DEFINED AS FIRST ONE) AND CENTRAL COMPARTMENT DEFIEND AS 2ND COMPARTMENT
            COMP=(TRANSIT1,DEFDOSE) ;1 GUT TRANIST 1 (F1 is associated with first compartment)
            COMP=(TRANSIT2) ;2 GUT TRANIST 2
            COMP=(TRANSIT3) ;3 GUT TRANIST 3
            COMP=(TRANSIT4) ;4 GUT TRANIST 4
            COMP=(TRANSIT5) ;5 GUT TRANIST 5
            COMP=(TRANSIT6) ;6 GUT TRANIST 6
            COMP=(TRANSIT7) ;7 GUT TRANIST 7
            COMP=(TRANSIT8) ;8 GUT TRANIST 8
            COMP=(TRANSIT9) ;9 GUT TRANIST 9
            COMP=(TRANSIT10) ;10 GUT TRANIST 10
            COMP=(TRANSIT11) ;11 GUT TRANIST 11
            COMP=(ABS) ;12 GUT ABS
            COMP=("CENTRAL",DEFOBS) ;13 CENTRAL CMT

;----------------------------------------------------------------------------------------------------------------------------------------------
$PK  
; Allometric scaling 
TVFFM = 43.5
TVWT = 56
TVFAT = 11.3

ALLMCLWT=(WT/TVWT)**0.75
ALLMVWT=WT/TVWT

ALLMCLFFM=(FFM/TVFFM)**0.75
ALLMVFFM=FFM/TVFFM

ALLMCLFAT=(FAT/TVFAT)**0.75
ALLMVFAT=FAT/TVFAT	
				 					 
; LACTATE_IMP = LACTATE
; IF (LACTATE.EQ.-99) LACTATE_IMP = 1.6

; TIMETOBREAKFAST_IMP = TIMETOBREAKFAST
; IF (TIMETOBREAKFAST.EQ.-99) TIMETOBREAKFAST_IMP = -4.82

; CRP_IMP = CRP
; IF (CRP.EQ.-99) CRP_IMP = 127

; AST_IMP = AST
; IF (AST.EQ.-99) AST_IMP = 37.5

; ALT_IMP = ALT
; IF (ALT.EQ.-99) ALT_IMP = 24

; TPROT_IMP = TPROT
; IF (TPROT.EQ.-99) TPROT_IMP = 80

; ALBUMIN_IMP = ALBUMIN
; IF (ALBUMIN.EQ.-99) ALBUMIN_IMP = 30

; UREA_IMP = UREA
; IF (UREA.EQ.-99) UREA_IMP = 5

; CREAT_IMP = CREAT
; IF (CREAT.EQ.-99) CREAT_IMP = 68

; HOSP_EFF = 1
; IF (HOSPITALIZED.EQ.1) HOSP_EFF = THETA(9)

; INDBIO = 1
; IF (FDC.EQ.0) INDBIO = THETA(10)
;-----------------------------------------------------PATIENT EFFECT ON ABSORPTION(KA & MTT)-----------------------------------------------
PATKAMTT = 1 ; FOR OUTPATIENTS
IF (PATIENT.EQ.1)	PATKAMTT=THETA(9) ; SURVIVED INPATIENTS
IF (PATIENT.EQ.2)	PATKAMTT=THETA(10) ; DIED INPATIENTS

;-----------------------------------------------------Log BRC EFFECT ON CLEARANCE----------------------------------------------------------
BRCMED = 6
BRCCL = 1

IF (BRC.NE.-99) THEN
	BRCCL = (BRC/BRCMED)**THETA(11) 
ENDIF

;------------------------------------------------------------------------------------------------------------------------------------------
; ------- BSV
BSVCL   = ETA(1)
BSVV    = ETA(2)
BSVBIO 	= ETA(3)
BSVKA 	= ETA(4)
BSVMTT 	= ETA(5)

; BOV    
BOVBIO 	= ETA(OCC_BIO)
BOVKA   = ETA(OCC_KA)
BOVMTT 	= ETA(OCC_MTT)    
;--------------------------------------------------------------------------------------------------------------------------------------
																																   
;---------Typical values--------------------------------------------------------------------------------------------------------------------------------------------------------
TVCL 	= THETA(1) * ALLMCLFFM * BRCCL
TVV 	= THETA(2) * ALLMVFFM
TVBIO 	= THETA(3) 
TVKA 	= THETA(4) *PATKAMTT
TVMTT 	= THETA(5) /PATKAMTT
TVNN 	= THETA(8)

;-----------Define individual parameters------------------------------------------------------------------------------------------------------------------------------------------
CL  = TVCL*EXP(BSVCL) ; CLEARANCE 
V   = TVV*EXP(BSVV) ; CENTRAL VOL. 
BIO = TVBIO*EXP(BSVBIO + BOVBIO) ; BIOAVAILABILITY
KA  = TVKA*EXP(BSVKA + BOVKA) ; ABS. RATE CONSTANT
MTT = TVMTT*EXP(BSVMTT + BOVMTT) ; LAG TIME
NN 	= TVNN

;--------------------------------------------------------------
; re-parameterization
F1 = BIO
KTR = (NN+1)/MTT
K12	=	KTR  ;Rate between transit CMT
K23	=	KTR   ;Rate between transit CMT
K34	=	KTR   ;Rate between transit CMT
K45	=	KTR   ;Rate between transit CMT
K56	=	KTR   ;Rate between transit CMT
K67	=	KTR   ;Rate between transit CMT
K78	=	KTR   ;Rate between transit CMT
K89	=	KTR   ;Rate between transit CMT
K9T10 =	KTR   ;Rate between transit CMT
K10T11 = KTR   ;Rate between transit CMT
K11T12 = KTR ;Rate between transit CMT
K12T13 = KA
K13T0	=   CL/V
S13  = 	V
;---------------------------------------------------------------------------------------------------------------------------------------------------------
$ERROR                                     
IPRED = A(13)/V

IRES = DV-IPRED


; LLOQ
LLOQ = 0.117

PROP = IPRED*THETA(6)
ADD = 0.2*LLOQ + THETA(7)

; For ADD, in this case we are coding THETA(.) as the additive error on top of 20% of the LLOQ.
; So the lower bound of THETA(.) can be zero. If it goes to zero, we can fix it, and the additive error
; will then be constrained to 20% of LLOQ + the value in THETA(.).
; REMEMBER about this when you report the value of ADD and its uncertainty! NONMEM gives uncertainty on THETA, not ADD
; An alternative approach is to set the lower bound of the THETA for the additive error to 20% of the LLOQ.
; In that case, one does not have to worry about adjusting the precision. On the other hand, this cannot be done if you have different LLOQs within your analysis (e.g. different labs)

; For BLQ==1 (i.e. first BLQ value in a series), we add extra additive error on the concentrations, since the value in DV has been imputed
IF(ICALL.NE.4.AND.RIF_BLQ==1) ADD = ADD + (LLOQ/2)

; For BLQ==2 (i.e. the trailing BLQ values in a series), we don't want these to influence the fit,
; we only want them for simulation-based diagnostics such as the VPC.
; So we define a separate error structure for these points. It has no proportional component
; (PROP = 0, as we would not want these points to affect our estimate of proportional error)
; and a FIXED and HUGE additive component (ADD = 1000000000, large with respect to the readings of concentration),
; so that the values do not affect the fit.
; It's also a good idea to repeat the diagnostic plots without the BLQ=2 points
IF(ICALL.NE.4.AND.RIF_BLQ==2) THEN
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
	RIF_BLQ = 1
ENDIF

; To calculate time after dose. 
IF(AMT.GT.0) THEN
	TIMEDOSE = TIME
	AMOUNTDOSE = AMT
ENDIF

TAD = TIME-TIMEDOSE

VARCL = BSVCL ;+ BOVCL
VARBIO = BSVBIO + BOVBIO
VARAUC = BSVBIO + BOVBIO - BSVCL ;- BOVCL
VARABS = BSVKA + BOVKA - BSVMTT - BOVMTT
VARKA = BSVKA + BOVKA
VARMTT = BSVMTT + BOVMTT

;------------------------------------------RETRIEVE AMOUNT IN EACH COMPARTMENT---------------------------------------------------------------------------------------
AA1 = A(12)
AA2 = A(13)
;--------------------------------------------------------------------------------------------------------------------------------------------------------------------
$THETA  (0,8.8248,20) ; 1 CL [L/h]
$THETA  (0,56.7748,100) ; 2 V [L]
$THETA  1 FIX ; 3 BIO
$THETA  (0,1.3798,5) ; 4 KA [1/h]
$THETA  (0,0.34167,2) ; 5 MTT [h]
$THETA  (0,0.172284,0.5) ; 6 PROP [%]
$THETA  0 FIX ; 7 ADD [mg/L]
$THETA  11 FIX ; 8 NN
$THETA  (-3,0.668258,3) ; 9 SURVIV_ABSORP
$THETA  (-3,0.380226,3) ; 10 DIED_ABSORP
$THETA  (-10,-0.3525,10) ; 11 BRCCL_POWER
;------------------------------------------------------------------------------------------------------------------------------------------------------
$OMEGA  0.179394  ;   1 BSV CL
$OMEGA  0  FIX  ;    2 BSV V
$OMEGA  0  FIX  ;  3 BSV BIO
$OMEGA  0  FIX  ;   4 BSV KA
$OMEGA  0  FIX  ;  5 BSV MTT
;---------------------------------------------------------------------------------------------------------------------------------------------------------------------
$OMEGA  BLOCK(1)
 0.0453815  ; 6,7 BOVBIO
$OMEGA  BLOCK(1) SAME
;---------------------------------------------------------------------------------------------------------------------------------------------------------------------
$OMEGA  BLOCK(1)
 1.41447  ;  8,9 BOVKA
$OMEGA  BLOCK(1) SAME
;---------------------------------------------------------------------------------------------------------------------------------------------------------------------
$OMEGA  BLOCK(1)
 0.879587  ; 10,11 BOVMTT
$OMEGA  BLOCK(1) SAME
;---------------------------------------------------------------------------------------------------------------------------------------------------------------------
$SIGMA  1  FIX
;-------------------------------------------------------------------------------------------------------------------------------------------------------
$ESTIMATION MSFO=RIF.msf MAXEVAL=9999 PRINT=1 METHOD=1
            INTERACTION NOABORT NSIG=3 NONINFETA=1 ETASTYPE=1 SIGL=9
			
$COVARIANCE PRINT=E UNCONDITIONAL MATRIX=R PRECOND=1

