;--------------------------------
; General
;--------------------------------
OutFile "nemesis-0.9.05.exe"
Name "nemesis"
XPStyle on
InstallDir $PROGRAMFILES\nemesis
LicenseText "nemesis licence"
LicenseData "..\LICENSE.txt"
BrandingText "nemesis"
SetCompressor /SOLID lzma

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
 
	; Put dat files (only well tested!)
	SetOutPath $INSTDIR\dat
	File ..\dat\*.slv
	File ..\dat\*.dat	
	
	; Put msc files
	SetOutPath $INSTDIR\msc
	File ..\msc\nemesis.dxf
	File ..\msc\nemesis.cui
	SetOutPath $INSTDIR\msc\lic
	File ..\msc\*.txt
	

	; Put rest of the files
	SetOutPath $INSTDIR
	File ..\LICENSE.txt
	File ..\CHANGELOG.txt
	;File ..\README.txt

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
