all: Surv_SpaceAgeTime_ByType_withUnknownType.so Surv_SpaceAgeTime_ByType_withUnknownType_No_AgeTime_Interaction.so Surv_SpaceAgeTime.so 
	# Comment here preserves the prior tab
Surv_SpaceAgeTime_ByType_withUnknownType.so: Surv_SpaceAgeTime_ByType_withUnknownType.cpp
	"$(R_HOME)/bin/Rscript" --vanilla -e "TMB::compile('Surv_SpaceAgeTime_ByType_withUnknownType.cpp','-O0 -g', safeunload = FALSE, safebounds = FALSE)"
	
Surv_SpaceAgeTime_ByType_withUnknownType_No_AgeTime_Interaction.so: Surv_SpaceAgeTime_ByType_withUnknownType_No_AgeTime_Interaction.cpp
	"$(R_HOME)/bin/Rscript" --vanilla -e "TMB::compile('Surv_SpaceAgeTime_ByType_withUnknownType_No_AgeTime_Interaction.cpp','-O0 -g', safeunload = FALSE, safebounds = FALSE)"

Surv_SpaceAgeTime.so: Surv_SpaceAgeTime.cpp 
	"$(R_HOME)/bin/Rscript" --vanilla -e "TMB::compile('Surv_SpaceAgeTime.cpp','-O0 -g', safeunload = FALSE, safebounds = FALSE)"

clean:
	rm -rf *o
