;; 1. Based on: 
; Model desc: Pyrazinamide.Model.KDHTB
; Settings for the memory of NONMEM
$SIZES      PD=-1000 LVR=-150 LTH=-200 MAXFCN=10000000 LNP4=-150000
$PROBLEM    |ADVAN13_TRANS1_PZA_DATA|
$INPUT      ID DAT2=DROP TIME OCC WHAT=DROP 
            EVID_RIF AMT_RIF RIFCONC RIF_BLQ MDV_RIF 
			EVID_INH AMT_INH INHCONC INH_BLQ MDV_INH 
            EVID AMT DV PZA_BLQ MDV 
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
$SUBROUTINE ADVAN5 TRANS1 ;TOL=9 ATOL=9 ; TOL is the precision to solve differential equations

$MODEL      NCOMPS=5 ; NUMBER OF COMPARTMENTS (ABSORPTION COMPATMENT (DEFINED AS FIRST ONE) AND CENTRAL COMPARTMENT DEFIEND AS 2ND COMPARTMENT
            COMP=(TRANSIT1,DEFDOSE) ;1 GUT TRANIST 1 (F1 is associated with first compartment)
            COMP=(TRANSIT2) ;2 GUT TRANIST 2
            COMP=(TRANSIT3) ;3 GUT TRANIST 3
            COMP=(ABS) ;4 GUT ABS
            COMP=("CENTRAL",DEFOBS) ;5 CENTRAL CMT

;----------------------------------------------------------------------------------------------------------------------------------------------
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

;----------------------------Imputing missing values to median--------------------------------------				 					 
; LACTATE_IMP = LACTATE
; IF (LACTATE.EQ.-999) LACTATE_IMP = 1.6

; TIMETOBREAKFAST_IMP = TIMETOBREAKFAST
; IF (TIMETOBREAKFAST.EQ.-999) TIMETOBREAKFAST_IMP = -4.82

; CRP_IMP = CRP
; IF (CRP.EQ.-999) CRP_IMP = 127

; AST_IMP = AST
; IF (AST.EQ.-999) AST_IMP = 37.5

; ALT_IMP = ALT
; IF (ALT.EQ.-999) ALT_IMP = 24

; TOTALPROTEIN_IMP = TOTALPROTEIN
; IF (TOTALPROTEIN.EQ.-999) TOTALPROTEIN_IMP = 80

; ALBUMIN_IMP = ALBUMIN
; IF (ALBUMIN.EQ.-999) ALBUMIN_IMP = 30

; UREA_IMP = UREA
; IF (UREA.EQ.-999) UREA_IMP = 5

; CREATININE_IMP = CREATININE
; IF (CREATININE.EQ.-999) CREATININE_IMP = 68

;--------------------------Hospitalized Effect on Variability in CL---------------------------------
PATETACL = 1
IF(PATIENT.EQ.1) PATETACL = THETA(9)
IF(PATIENT.EQ.2) PATETACL = THETA(10)

;--------------------Typical values (thetas & etas)-------------------------------------------------
;population parameters
TVCL        = THETA(1) * ALLMCLFFM
TVV         = THETA(2) * ALLMVFFM
TVBIO       = THETA(3)
TVKA        = THETA(4) 
TVMTT       = THETA(5) 
TVNN        = THETA(8)

;BSV 
BSVCL       = ETA(1) * PATETACL
BSVV        = ETA(2)
BSVBIO 	= ETA(3)
BSVKA 	= ETA(4)
BSVMTT 	= ETA(5)

;BOV    
BOVBIO 	= ETA(OCC_BIO)
BOVKA       = ETA(OCC_KA)
BOVMTT 	= ETA(OCC_MTT) 
;---------------------------------------------------------------------------------------------------

;-----------Individual parameters-------------------------------------------------------------------
CL    = TVCL*EXP(BSVCL)                   ; Clearance
V     = TVV*EXP(BSVV)                     ; CENTRAL VOL. 
BIO   = TVBIO*EXP(BSVBIO + BOVBIO)        ; BIOAVAILABILITY
KA    = TVKA*EXP(BSVKA + BOVKA) ; ABS. RATE CONSTANT
MTT   = TVMTT*EXP(BSVMTT + BOVMTT) ; MEAN TRANSIT TIME
NN 	= TVNN
;---------------------------------------------------------------------------------------------------
VARCL       = BSVCL ;+ BOVCL
VARBIO      = BSVBIO + BOVBIO
VARAUC      = BSVBIO + BOVBIO - BSVCL ;- BOVCL
VARABS = BSVKA + BOVKA - BSVMTT - BOVMTT
VARKA = BSVKA + BOVKA
VARMTT = BSVMTT + BOVMTT
;===================================================================================================

;------------------------------------Re-parameterization--------------------------------------------
;-----re-parameterization
F1 = BIO
KTR = (NN+1)/MTT
K12	=	KTR  ;Rate between transit CMT
K23	=	KTR   ;Rate between transit CMT
K34	=	KTR   ;Rate between transit CMT

K45 = KA
K50	= CL/V

;----------------------------------------Error model------------------------------------------------
$ERROR                                     
IPRED = A(5)/V
IRES = DV-IPRED

LLOQ = 0.203

PROP = IPRED*THETA(6)
ADD = 0.2*LLOQ + THETA(7)

; For ADD, in this case we are coding THETA(.) as the additive error on top of 20% of the LLOQ.
; So the lower bound of THETA(.) can be zero. If it goes to zero, we can fix it, and the additive error
; will then be constrained to 20% of LLOQ + the value in THETA(.).
; REMEMBER about this when you report the value of ADD and its uncertainty! NONMEM gives uncertainty on THETA, not ADD
; An alternative approach is to set the lower bound of the THETA for the additive error to 20% of the LLOQ.
; In that case, one does not have to worry about adjusting the precision. On the other hand, this cannot be done if you have different LLOQs within your analysis (e.g. different labs)

; For BLQ==1 (i.e. first BLQ value in a series), we add extra additive error on the concentrations, since the value in DV has been imputed
IF(ICALL.NE.4.AND.PZA_BLQ==1) ADD = ADD + (LLOQ/2)

; For BLQ==2 (i.e. the trailing BLQ values in a series), we don't want these to influence the fit,
; we only want them for simulation-based diagnostics such as the VPC.
; So we define a separate error structure for these points. It has no proportional component
; (PROP = 0, as we would not want these points to affect our estimate of proportional error)
; and a FIXED and HUGE additive component (ADD = 1000000000, large with respect to the readings of concentration),
; so that the values do not affect the fit.
; It's also a good idea to repeat the diagnostic plots without the BLQ=2 points
IF(ICALL.NE.4.AND.PZA_BLQ==2) THEN
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
	PZA_BLQ = 1
ENDIF

; To calculate time after dose. 
IF(AMT.GT.0) THEN
	TIMEDOSE = TIME
	AMOUNTDOSE = AMT
ENDIF

TAD = TIME-TIMEDOSE

;------------------------------------------RETRIEVE AMOUNT IN EACH COMPARTMENT---------------------------------------------------------------------------------------
AA1 = A(4)
AA2 = A(5)
;--------------------------------------------------------------------------------------------------------------------------------------------------------------------
$THETA  (0,2.61378,20) ; 1 CL [L/h]
$THETA  (0,35.9666,200) ; 2 V [L]
$THETA  1 FIX ; 3 BIO
$THETA  (0,1.91548,10) ; 4 KA [1/h]
$THETA  (0,0.378744,3) ; 5 MTT [h]
$THETA  (0,0.113882,0.7) ; 6 PROP [%]
$THETA  (0,2.4809,10) ; 7 ADD [mg/L]
$THETA  3 FIX ; 8 NN
$THETA  (0,1.69826) ; 9 SURVIV_PATETACL
$THETA  (0,3.56122) ; 10 DIED_PATETACL
;------------------------------------------------------------------------------------------------------------------------------------------------------
$OMEGA  0.0395359  ;   1 BSV CL
$OMEGA  0  FIX  ;    2 BSV V
$OMEGA  0  FIX  ;  3 BSV BIO
$OMEGA  0  FIX  ;   4 BSV KA
$OMEGA  0  FIX  ;  5 BSV MTT
;---------------------------------------------------------------------------------------------------------------------------------------------------------------------
$OMEGA  BLOCK(1)
 0.011176  ; 6,7 BOVBIO
$OMEGA  BLOCK(1) SAME
;---------------------------------------------------------------------------------------------------------------------------------------------------------------------
$OMEGA  BLOCK(1)
 0.568935  ;  8,9 BOVKA
$OMEGA  BLOCK(1) SAME
;---------------------------------------------------------------------------------------------------------------------------------------------------------------------
$OMEGA  BLOCK(1)
 1.0298  ; 10,11 BOVMTT
$OMEGA  BLOCK(1) SAME
;---------------------------------------------------------------------------------------------------------------------------------------------------------------------
$SIGMA  1  FIX
;-------------------------------------------------------------------------------------------------------------------------------------------------------
$ESTIMATION MSFO=PZA.msf MAXEVAL=9999 PRINT=1 METHOD=1
            INTERACTION NOABORT NSIG=3 NONINFETA=1 ETASTYPE=1 SIGL=9
			
$COVARIANCE PRINT=E UNCONDITIONAL MATRIX=R PRECOND=1
