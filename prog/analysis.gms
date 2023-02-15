$TITLE Enduse integration analysis
SET
Dum/
"Fin_Ene", "Fin_Ene_Ele_Heat", "Fin_Ene_Liq_and_Gas", "Fin_Ene_SolidsCoa", "Fin_Ene_SolidsBio", "Fin_Ene_Hyd", "Fin_Ene_Ind", "Fin_Ene_Com", "Fin_Ene_Res", "Fin_Ene_Tra"

Sec_Ene_Ele_Heat, Sec_Ene_Ele_Heat_Coa, Sec_Ene_Ele_Oil, Sec_Ene_Ele_Heat_Gas, Sec_Ene_Ele_Heat_Nuc, Sec_Ene_Ele_Hyd, Sec_Ene_Ele_Solar, Sec_Ene_Ele_Win, Sec_Ene_Ele_Heat_Bio, Sec_Ene_Ele_Heat_Geo
Sec_Ene_Hyd, Sec_Ene_Hyd_Gas, Sec_Ene_Hyd_Bio, Sec_Ene_Hyd_Ele

"Fin_Ene_Ind_Ele_Heat", "Fin_Ene_Ind_Liq_and_Gas", "Fin_Ene_Ind_SolidsCoa", "Fin_Ene_Ind_SolidsBio", "Fin_Ene_Ind_Hyd", "Fin_Ene_Com_Ele_Heat", "Fin_Ene_Com_Liq_and_Gas", "Fin_Ene_Com_SolidsCoa", "Fin_Ene_Com_SolidsBio", "Fin_Ene_Com_Hyd", "Fin_Ene_Res_Ele_Heat", "Fin_Ene_Res_Liq_and_Gas", "Fin_Ene_Res_SolidsCoa", "Fin_Ene_Res_SolidsBio", "Fin_Ene_Res_Hyd", "Fin_Ene_Tra_Ele", "Fin_Ene_Tra_Liq_and_Gas", "Fin_Ene_Tra_Liq_Bio", "Fin_Ene_Tra_Liq_Oil", "Fin_Ene_Tra_Hyd"

/
Region,Var,ModName,SCENARIO,Y
Value/Value/
Ite  Iteration/i0*i10/
Indi/
MAE		Mean absolute error
RMSE	Root Mean Square error
RMSEStand	Root Mean Square error standardized
MAEStand	Mean absolute error standardized
MEAN	Mean
Correlation	Correlation coefficient 
/
ModSpe/CGE,Tech/
;
ALIAS(Var,Var2);

PARAMETER
DataLoad(ModName,Region,Var,SCENARIO,Y,Value)
DataLoad1(ModName,Region,Var,SCENARIO,Y)
*DataLoadComp(Ite,ModName,ModName,Region,Var,SCENARIO,Y)
MAE(Ite,SCENARIO,Region,Var)
RMSE(Ite,SCENARIO,Region,Var)
RMSEStand(Ite,SCENARIO,Region,Var)
MAEStand(Ite,SCENARIO,Region,Var)
MEAN(Ite,SCENARIO,Region,Var)
Stat(Indi,Ite,SCENARIO,Region,Var)
ConvStat(Indi,Ite,ModSpe,SCENARIO,Region,Var)
MeanModel(ModName,SCENARIO,Region,Var)
VarModel(ModName,SCENARIO,Region,Var)
;

$gdxin '../output/data/allcombine.gdx'
$load Region,Var,ModName,SCENARIO,Y
$load DataLoad=allmodel

ALIAS(ModName,ModName2,ModName1);

SET
ModCompMap(Ite,ModName,ModName)/
i0	.	CGE0	.	Tech0
i1	.	CGE1	.	Tech0
i2	.	CGE2	.	Tech1
i3	.	CGE3	.	Tech2
i4	.	CGE4	.	Tech3
i5	.	CGE5	.	Tech4
i6	.	CGE6	.	Tech5
/
ModConverMap(ModName1,ModName2,Ite,ModSpe)/
CGE0  . CGE1  . i0  . CGE
CGE1  . CGE2  . i1  . CGE
CGE2  . CGE3  . i2  . CGE
CGE3  . CGE4  . i3  . CGE
CGE4  . CGE5  . i4  . CGE
CGE5  . CGE6  . i5  . CGE
Tech0  . Tech1  . i0 . Tech
Tech1  . Tech2  . i1  . Tech
Tech2  . Tech3  . i2  . Tech
Tech3  . Tech4  . i3  . Tech
Tech4  . Tech5  . i4  . Tech
/

Varmain(Var)/
"Fin_Ene"
"Fin_Ene_Ele_Heat"
"Fin_Ene_Liq_and_Gas"
"Fin_Ene_SolidsCoa"
"Fin_Ene_SolidsBio"
"Fin_Ene_Hyd"
"Fin_Ene_Ind"
"Fin_Ene_Ind_Ele_Heat"
"Fin_Ene_Ind_Liq_and_Gas"
"Fin_Ene_Ind_SolidsCoa"
"Fin_Ene_Ind_SolidsBio"
"Fin_Ene_Ind_Hyd"
"Fin_Ene_Com"
"Fin_Ene_Com_Ele_Heat"
"Fin_Ene_Com_Liq_and_Gas"
"Fin_Ene_Com_SolidsCoa"
"Fin_Ene_Com_SolidsBio"
"Fin_Ene_Com_Hyd"
"Fin_Ene_Res"
"Fin_Ene_Res_Ele_Heat"
"Fin_Ene_Res_Liq_and_Gas"
"Fin_Ene_Res_SolidsCoa"
"Fin_Ene_Res_SolidsBio"
"Fin_Ene_Res_Hyd"
"Fin_Ene_Tra"
"Fin_Ene_Tra_Ele"
"Fin_Ene_Tra_Liq_and_Gas"
"Fin_Ene_Tra_Liq_Bio"
"Fin_Ene_Tra_Liq_Oil"
"Fin_Ene_Tra_Hyd"

Sec_Ene_Ele_Heat
Sec_Ene_Ele_Heat_Coa
Sec_Ene_Ele_Oil
Sec_Ene_Ele_Heat_Gas
Sec_Ene_Ele_Heat_Nuc
Sec_Ene_Ele_Hyd
Sec_Ene_Ele_Solar
Sec_Ene_Ele_Win
Sec_Ene_Ele_Heat_Bio
Sec_Ene_Ele_Heat_Geo

Sec_Ene_Hyd
Sec_Ene_Hyd_Gas
Sec_Ene_Hyd_Bio
Sec_Ene_Hyd_Ele

$ontext
"Fin_Ene_Res_Gas"
"Fin_Ene_Res_Liq"
"Fin_Ene_Ind_Gas"
"Fin_Ene_Ind_Liq"
"Fin_Ene_Gas"
"Fin_Ene_Liq"
"Fin_Ene_Com_Gas"
"Fin_Ene_Com_Liq"
"Fin_Ene_Tra_Gas"
$offtext
/
Y2020Comp(Y)/
2020,2025,2030,2035,2040,2045,2050,2055,2060,2065,2070,2075,2080,2085,2090,2095,2100
/
RegionMain(Region)/
World
R5ASIA
"R5OECD90+EU"
R5REF
R5MAF
R5LAM
/

;
DataLoad1(ModName,Region,Var,SCENARIO,Y)=DataLoad(ModName,Region,Var,SCENARIO,Y,"Value");

Stat("MAE",Ite,SCENARIO,Region,Var)$(Varmain(Var) AND RegionMain(Region) AND SUM((Y,ModName1,ModName2)$(Y2020Comp(Y) AND ModCompMap(Ite,ModName1,ModName2) AND DataLoad1(ModName1,Region,Var,SCENARIO,Y) AND DataLoad1(ModName2,Region,Var,SCENARIO,Y)),1))=
  SUM((Y,ModName1,ModName2)$(Y2020Comp(Y) AND ModCompMap(Ite,ModName1,ModName2) AND DataLoad1(ModName1,Region,Var,SCENARIO,Y) AND DataLoad1(ModName2,Region,Var,SCENARIO,Y)),ABS((DataLoad1(ModName1,Region,Var,SCENARIO,Y)-DataLoad1(ModName2,Region,Var,SCENARIO,Y))))/
  SUM((Y,ModName1,ModName2)$(Y2020Comp(Y) AND ModCompMap(Ite,ModName1,ModName2) AND DataLoad1(ModName1,Region,Var,SCENARIO,Y) AND DataLoad1(ModName2,Region,Var,SCENARIO,Y)),1);

Stat("RMSE",Ite,SCENARIO,Region,Var)$(Varmain(Var) AND RegionMain(Region) AND SUM((Y,ModName1,ModName2)$(Y2020Comp(Y) AND ModCompMap(Ite,ModName1,ModName2) AND DataLoad1(ModName1,Region,Var,SCENARIO,Y) AND DataLoad1(ModName2,Region,Var,SCENARIO,Y)),1))=
  SUM((Y,ModName1,ModName2)$(Y2020Comp(Y) AND ModCompMap(Ite,ModName1,ModName2) AND DataLoad1(ModName1,Region,Var,SCENARIO,Y) AND DataLoad1(ModName1,Region,Var,SCENARIO,Y)),ABS((DataLoad1(ModName1,Region,Var,SCENARIO,Y)-DataLoad1(ModName2,Region,Var,SCENARIO,Y)))**2)/
  SUM((Y,ModName1,ModName2)$(Y2020Comp(Y) AND ModCompMap(Ite,ModName1,ModName2) AND DataLoad1(ModName1,Region,Var,SCENARIO,Y) AND DataLoad1(ModName2,Region,Var,SCENARIO,Y)),1)**2;

Stat("MEAN",Ite,SCENARIO,Region,Var)$(Varmain(Var) AND RegionMain(Region) AND SUM((Y,ModName1,ModName2)$(Y2020Comp(Y) AND ModCompMap(Ite,ModName1,ModName2) AND DataLoad1(ModName1,Region,Var,SCENARIO,Y) AND DataLoad1(ModName2,Region,Var,SCENARIO,Y)),1))=
  SUM((Y,ModName1,ModName2)$(Y2020Comp(Y) AND ModCompMap(Ite,ModName1,ModName2) AND DataLoad1(ModName1,Region,Var,SCENARIO,Y) AND DataLoad1(ModName1,Region,Var,SCENARIO,Y)),ABS((DataLoad1(ModName1,Region,Var,SCENARIO,Y)+DataLoad1(ModName2,Region,Var,SCENARIO,Y))))/
  SUM((Y,ModName1,ModName2)$(Y2020Comp(Y) AND ModCompMap(Ite,ModName1,ModName2) AND DataLoad1(ModName1,Region,Var,SCENARIO,Y) AND DataLoad1(ModName2,Region,Var,SCENARIO,Y)),1);

MeanModel(ModName,SCENARIO,Region,Var)$(Varmain(Var) AND RegionMain(Region) AND SUM(Y$(Y2020Comp(Y) AND DataLoad1(ModName,Region,Var,SCENARIO,Y)),1))=
  SUM(Y$(Y2020Comp(Y) AND DataLoad1(ModName,Region,Var,SCENARIO,Y)),DataLoad1(ModName,Region,Var,SCENARIO,Y))/
  SUM(Y$(Y2020Comp(Y) AND DataLoad1(ModName,Region,Var,SCENARIO,Y)),1);

VarModel(ModName,SCENARIO,Region,Var)$(Varmain(Var) AND RegionMain(Region) AND SUM(Y$(Y2020Comp(Y) AND DataLoad1(ModName,Region,Var,SCENARIO,Y)),1))=
  SUM(Y$(Y2020Comp(Y) AND DataLoad1(ModName,Region,Var,SCENARIO,Y)),ABS(DataLoad1(ModName,Region,Var,SCENARIO,Y)-MeanModel(ModName,SCENARIO,Region,Var))**2)
*  /SUM(Y$(Y2020Comp(Y) AND DataLoad1(ModName,Region,Var,SCENARIO,Y)),1)
  ;

Stat("RMSEStand",Ite,SCENARIO,Region,Var)$(Stat("MEAN",Ite,SCENARIO,Region,Var))=Stat("RMSE",Ite,SCENARIO,Region,Var)/Stat("MEAN",Ite,SCENARIO,Region,Var);
Stat("MAEStand",Ite,SCENARIO,Region,Var)$(Stat("MEAN",Ite,SCENARIO,Region,Var))=Stat("MAE",Ite,SCENARIO,Region,Var)/Stat("MEAN",Ite,SCENARIO,Region,Var);

Stat("Correlation",Ite,SCENARIO,Region,Var)$(SUM((ModName1,ModName2)$ModCompMap(Ite,ModName1,ModName2),VarModel(ModName1,SCENARIO,Region,Var)*VarModel(ModName2,SCENARIO,Region,Var)))=
  SUM(Y$(Y2020Comp(Y)),
    SUM(ModName1$(DataLoad1(ModName1,Region,Var,SCENARIO,Y) AND SUM(ModName2,ModCompMap(Ite,ModName1,ModName2))),DataLoad1(ModName1,Region,Var,SCENARIO,Y)-MeanModel(ModName1,SCENARIO,Region,Var))*
    SUM(ModName2$(DataLoad1(ModName2,Region,Var,SCENARIO,Y) AND SUM(ModName1,ModCompMap(Ite,ModName1,ModName2))),DataLoad1(ModName2,Region,Var,SCENARIO,Y)-MeanModel(ModName2,SCENARIO,Region,Var))
  )/
    SQRT( SUM(ModName1$SUM(ModName2,ModCompMap(Ite,ModName1,ModName2)),VarModel(ModName1,SCENARIO,Region,Var))
	     *SUM(ModName2$SUM(ModName1,ModCompMap(Ite,ModName1,ModName2)),VarModel(ModName2,SCENARIO,Region,Var))
	);

ConvStat("MAE",Ite,ModSpe,SCENARIO,Region,Var)$(Varmain(Var) AND RegionMain(Region) AND SUM((Y,ModName1,ModName2)$(Y2020Comp(Y) AND ModConverMap(ModName1,ModName2,Ite,ModSpe) AND DataLoad1(ModName1,Region,Var,SCENARIO,Y) AND DataLoad1(ModName2,Region,Var,SCENARIO,Y)),1))=
  SUM((Y,ModName1,ModName2)$(Y2020Comp(Y) AND ModConverMap(ModName1,ModName2,Ite,ModSpe) AND DataLoad1(ModName1,Region,Var,SCENARIO,Y) AND DataLoad1(ModName2,Region,Var,SCENARIO,Y)),ABS((DataLoad1(ModName1,Region,Var,SCENARIO,Y)-DataLoad1(ModName2,Region,Var,SCENARIO,Y))))/
  SUM((Y,ModName1,ModName2)$(Y2020Comp(Y) AND ModConverMap(ModName1,ModName2,Ite,ModSpe) AND DataLoad1(ModName1,Region,Var,SCENARIO,Y) AND DataLoad1(ModName2,Region,Var,SCENARIO,Y)),1);
ConvStat("MEAN",Ite,ModSpe,SCENARIO,Region,Var)$(Varmain(Var) AND RegionMain(Region) AND SUM((Y,ModName1,ModName2)$(Y2020Comp(Y) AND ModConverMap(ModName1,ModName2,Ite,ModSpe) AND DataLoad1(ModName1,Region,Var,SCENARIO,Y) AND DataLoad1(ModName2,Region,Var,SCENARIO,Y)),1))=
  SUM((Y,ModName1,ModName2)$(Y2020Comp(Y) AND ModConverMap(ModName1,ModName2,Ite,ModSpe) AND DataLoad1(ModName1,Region,Var,SCENARIO,Y) AND DataLoad1(ModName1,Region,Var,SCENARIO,Y)),ABS((DataLoad1(ModName1,Region,Var,SCENARIO,Y)+DataLoad1(ModName2,Region,Var,SCENARIO,Y))))/
  SUM((Y,ModName1,ModName2)$(Y2020Comp(Y) AND ModConverMap(ModName1,ModName2,Ite,ModSpe) AND DataLoad1(ModName1,Region,Var,SCENARIO,Y) AND DataLoad1(ModName2,Region,Var,SCENARIO,Y)),1);
ConvStat("MAEStand",Ite,ModSpe,SCENARIO,Region,Var)$(ConvStat("MEAN",Ite,ModSpe,SCENARIO,Region,Var))=ConvStat("MAE",Ite,ModSpe,SCENARIO,Region,Var)/ConvStat("MEAN",Ite,ModSpe,SCENARIO,Region,Var);

execute_unload '../output/data/analysis.gdx'
Stat
ConvStat
;
