;--------------------------------
; General
;--------------------------------
OutFile "nemesis-0.9.01.exe"
Name "nemesis. an experimental finite element code"
XPStyle on
InstallDir $PROGRAMFILES\nemesis
LicenseText "nemesis licence"
LicenseData "LICENSE.txt"
BrandingText "nemesis. an experimental finite element code"
;--------------------------------
; Pages
;--------------------------------
Page license
Page directory
Page instfiles
UninstPage uninstConfirm
UninstPage instfiles

;--------------------------------
; The stuff to install
;--------------------------------
Section ""
	; General
	SetOverwrite on

	; Put bin files
	SetOutPath $INSTDIR\bin
	File ..\bin\nemesis.exe
	;File ..\bin\mesh.py
	;File ..\bin\plot.py
	;File ..\bin\formula.py
	
	; Put sci files
	SetOutPath $INSTDIR\sci
	File ..\sci\*.*
 
	; Put rest of the files
	SetOutPath $INSTDIR
	File ..\LICENSE
	File ..\README
	File ..\NEWS
	
	; Put dat files (only well tested!)
	SetOutPath $INSTDIR\dat
	File ..\dat\*.slv
	File ..\dat\*.dat	
	
	; Write uninstaller
	WriteUninstaller "uninstall.exe"

	; Create shortcuts
	CreateDirectory "$SMPROGRAMS\nemesis"
	CreateShortCut  "$SMPROGRAMS\nemesis\nemesis.lnk"   	"$INSTDIR\bin\nemesis.exe"
	CreateShortCut  "$SMPROGRAMS\nemesis\SciTe.lnk"     	"$INSTDIR\sci\SciTe.exe"	
	CreateShortCut  "$SMPROGRAMS\nemesis\README.lnk"     "$INSTDIR\README"	
	CreateShortCut  "$SMPROGRAMS\nemesis\NEWS.lnk"     	"$INSTDIR\NEWS"	
	CreateShortCut  "$SMPROGRAMS\nemesis\LICENSE.lnk"     "$INSTDIR\LICENSE"		
	CreateShortCut  "$SMPROGRAMS\nemesis\Uninstall.lnk" 	"$INSTDIR\uninstall.exe"
SectionEnd

;--------------------------------
; Uninstaller
;--------------------------------
Section "Uninstall"
	; Remove files and uninstaller
	RMDir /r $INSTDIR
	
	; Remove shortcuts
	Delete "$SMPROGRAMS\nemesis\*.*"
	RMDir  "$SMPROGRAMS\nemesis"

SectionEnd
