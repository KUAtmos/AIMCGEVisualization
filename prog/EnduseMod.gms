$TITLE Enduse data modification
SET
Region,Var,ModName,SCENARIO,Y
Value/Value/
;
ALIAS(Var,Var2);

PARAMETER
EnduseCombined(Region,Var,ModName,SCENARIO,Y,*)
EnduseCombined0(Region,Var,ModName,SCENARIO,Y)
;

$gdxin '../modeloutput/Endusecombine.gdx'
$load Region,Var,ModName,SCENARIO,Y
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
/
;

EnduseCombined(Region,Var,ModName,SCENARIO,Y,"Value")$(SUM(Var2$(MapVar(Var,Var2)),1))=SUM(Var2$(MapVar(Var,Var2)),EnduseCombined(Region,Var2,ModName,SCENARIO,Y,"Value"));
execute_unload '../modeloutput/EndusecombineMod.gdx'
EnduseCombined
;
