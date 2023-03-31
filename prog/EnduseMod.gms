$TITLE Enduse data modification
$setglobal Maindir
$setglobal CGERepoDir %Maindir%../../AIMCGE
$setglobal IntRepoDir %Maindir%../../../../IntTool
$setglobal outputdir ../

SET
Region,Var,ModName,SCENARIO,Y
Value/Value/
RCGE/
USA	"USA"
XE25	"EU"
XER	"Rest of EU"
TUR	"Turkey"
XOC	"New Zealand and Australia"
CHN	"China"
IND	"India"
JPN	"Japan"
XSE	"Rest of East and South East Asia"
XSA	"Rest of Asia"
CAN	"Canada"
BRA	"Brazil"
XLM	"Rest of Brazil"
CIS	"Former USSR"
XME	"Middle East"
XNF	"North Africa"
XAF	"Sub-Sahara"
World	
"R5OECD90+EU"
R5REF	
R5ASIA	
R5MAF	
R5LAM	
/
Var/
$  include %CGERepoDir%/tools/iiasa_data_submission/define/variables.set
/
;
ALIAS(Var,Var2);

PARAMETER
EnduseCombined(Region,Var,ModName,SCENARIO,Y,*)
EnduseCombined2(RCGE,Var,ModName,SCENARIO,Y,*)
EnduseCombined0(Region,Var,ModName,SCENARIO,Y)
;

$gdxin '%outputdir%/modeloutput/Endusecombine.gdx'
$load Region,ModName,SCENARIO,Y
$load EnduseCombined

SET
MapVar(Var,Var2)/
Fin_Ene_Com	.	Fin_Ene_Oth_Sec
Fin_Ene_Com_Ele	.	Fin_Ene_Oth_Sec_Ele
Fin_Ene_Com_Gas	.	Fin_Ene_Oth_Sec_Gas
Fin_Ene_Com_Heat	.	Fin_Ene_Oth_Sec_Heat
Fin_Ene_Com_Liq	.	Fin_Ene_Oth_Sec_Liq
Fin_Ene_Com_SolidsBio	.	Fin_Ene_Oth_Sec_SolidsBio
Fin_Ene_Com_SolidsCoa	.	Fin_Ene_Oth_Sec_SolidsCoa
Fin_Ene_Com_Liq_and_Gas	.	Fin_Ene_Oth_Sec_Liq
Fin_Ene_Com_Liq_and_Gas	.	Fin_Ene_Oth_Sec_Gas
Fin_Ene_Com_Ele_Heat	.	Fin_Ene_Oth_Sec_Ele
Fin_Ene_Com_Ele_Heat	.	Fin_Ene_Oth_Sec_Heat
Fin_Ene_Com	.	Fin_Ene_Com
Fin_Ene_Com_Ele	.	Fin_Ene_Com_Ele
Fin_Ene_Com_Gas	.	Fin_Ene_Com_Gas
Fin_Ene_Com_Heat	.	Fin_Ene_Com_Heat
Fin_Ene_Com_Liq	.	Fin_Ene_Com_Liq
Fin_Ene_Com_SolidsBio	.	Fin_Ene_Com_SolidsBio
Fin_Ene_Com_SolidsCoa	.	Fin_Ene_Com_SolidsCoa
Fin_Ene_Com_Liq_and_Gas	.	Fin_Ene_Com_Liq_and_Gas
Fin_Ene_Com_Ele_Heat	.	Fin_Ene_Com_Ele_Heat
Sec_Ene_Inp_Coa_Ele_Heat  . Sec_Ene_Inp_Coa_Ele
Sec_Ene_Inp_Coa_Ele_Heat  . Sec_Ene_Inp_Coa_Heat
Sec_Ene_Inp_Gas_Ele_Heat  . Sec_Ene_Inp_Gas_Ele
Sec_Ene_Inp_Gas_Ele_Heat  . Sec_Ene_Inp_Gas_Heat
Sec_Ene_Inp_Oil_Ele_Heat  . Sec_Ene_Inp_Oil_Ele
Sec_Ene_Inp_Oil_Ele_Heat  . Sec_Ene_Inp_Oil_Heat
Sec_Ene_Inp_Bio_Ele_Heat  . Sec_Ene_Inp_Bio_Ele
Sec_Ene_Inp_Bio_Ele_Heat  . Sec_Ene_Inp_Bio_Heat
/
RMAP(Region,RCGE)/
$  include %IntRepoDir%/define/region32.map
World	.	World
"R5OECD90+EU"	.	"R5OECD90+EU"
R5REF	.	R5REF
R5ASIA	.	R5ASIA
R5MAF	.	R5MAF
R5LAM	.	R5LAM
/
VarWtMap_load(var,Var2)/
$  include %CGERepoDir%/tools/iiasa_data_submission/define/Weight.map
/
VarWtMap(var,Var2)
;
VarWtMap(var,Var2)$(VarWtMap_load(var,Var2) AND (NOT SAMEAS(var,var2)))=YES;

EnduseCombined(Region,Var,ModName,SCENARIO,Y,"Value")$(SUM(Var2$(MapVar(Var,Var2)),1))=SUM(Var2$(MapVar(Var,Var2)),EnduseCombined(Region,Var2,ModName,SCENARIO,Y,"Value"));
EnduseCombined2(RCGE,Var,ModName,SCENARIO,Y,"Value")$(NOT SUM(Var2$VarWtMap(var,Var2),1))=SUM(Region$(RMAP(Region,RCGE)),EnduseCombined(Region,Var,ModName,SCENARIO,Y,"Value"));
EnduseCombined2(RCGE,Var,ModName,SCENARIO,Y,"Value")$(SUM(Var2$VarWtMap(var,Var2),1) AND SUM(Var2$VarWtMap(var,Var2),SUM(Region$(RMAP(Region,RCGE)),EnduseCombined(Region,Var2,ModName,SCENARIO,Y,"Value"))))=
  SUM(Var2$VarWtMap(var,Var2),SUM(Region$(RMAP(Region,RCGE)),EnduseCombined(Region,Var,ModName,SCENARIO,Y,"Value")*EnduseCombined(Region,Var2,ModName,SCENARIO,Y,"Value")))/
  SUM(Var2$VarWtMap(var,Var2),SUM(Region$(RMAP(Region,RCGE)),EnduseCombined(Region,Var2,ModName,SCENARIO,Y,"Value")));

execute_unload '%outputdir%modeloutput/EndusecombineMod.gdx'
EnduseCombined2
VarWtMap
Region
;
