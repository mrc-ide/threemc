all: Surv_SpaceAgeTime_ByType_withUnknownType_Const_Paed_MMC.so Surv_SpaceAgeTime_ByType_withUnknownType_Const_Paed_MMC_RW.so Surv_SpaceAgeTime.so # Surv_SpaceAgeTime_RW.so  
	# Comment here preserves the prior tab
Surv_SpaceAgeTime_ByType_withUnknownType_Const_Paed_MMC.so: Surv_SpaceAgeTime_ByType_withUnknownType_Const_Paed_MMC.cpp
	"$(R_HOME)/bin/Rscript" --vanilla -e "TMB::compile('Surv_SpaceAgeTime_ByType_withUnknownType_Const_Paed_MMC.cpp','-O0 -g', safeunload = FALSE, safebounds = FALSE)"
	
Surv_SpaceAgeTime_ByType_withUnknownType_Const_Paed_MMC_RW.so: Surv_SpaceAgeTime_ByType_withUnknownType_Const_Paed_MMC_RW.cpp
	"$(R_HOME)/bin/Rscript" --vanilla -e "TMB::compile('Surv_SpaceAgeTime_ByType_withUnknownType_Const_Paed_MMC_RW.cpp','-O0 -g', safeunload = FALSE, safebounds = FALSE)"
	
Surv_SpaceAgeTime_ByType_withUnknownType_No_AgeTime_Interaction.cpp: Surv_SpaceAgeTime_ByType_withUnknownType_No_AgeTime_Interaction.cpp
	"$(R_HOME)/bin/Rscript" --vanilla -e "TMB::compile('Surv_SpaceAgeTime_ByType_withUnknownType_No_AgeTime_Interaction.cpp','-O0 -g', safeunload = FALSE, safebounds = FALSE)"
	
Surv_SpaceAgeTime_ByType_withUnknownType_No_AgeTime_Interaction_RW.cpp: Surv_SpaceAgeTime_ByType_withUnknownType_No_AgeTime_Interaction_RW.cpp
	"$(R_HOME)/bin/Rscript" --vanilla -e "TMB::compile('Surv_SpaceAgeTime_ByType_withUnknownType_No_AgeTime_Interaction_RW.cpp','-O0 -g', safeunload = FALSE, safebounds = FALSE)"

Surv_SpaceAgeTime.so: Surv_SpaceAgeTime.cpp 
	"$(R_HOME)/bin/Rscript" --vanilla -e "TMB::compile('Surv_SpaceAgeTime.cpp','-O0 -g', safeunload = FALSE, safebounds = FALSE)"

Surv_SpaceAgeTime_RW.so: Surv_SpaceAgeTime_RW.cpp 
	"$(R_HOME)/bin/Rscript" --vanilla -e "TMB::compile('Surv_SpaceAgeTime_RW.cpp','-O0 -g', safeunload = FALSE, safebounds = FALSE)"

clean:
	rm -rf *o
