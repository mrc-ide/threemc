all: Surv_SpaceAgeTime_ByType_withUnknownType_Const_Paed_MMC.so Surv_SpaceAgeTime_ByType_withUnknownType_Const_Paed_MMC_RW.so Surv_SpaceAgeTime_ByType_withUnknownType_Const_Paed_MMC2.so  Surv_SpaceAgeTime_ByType_withUnknownType.so Surv_SpaceAgeTime_ByType_withUnknownType_RW.so Surv_SpaceAgeTime_ByType_withUnknownType_RW2.so Surv_SpaceAgeTime_ByType_withUnknownType_RW_MMC2.so  Surv_SpaceAgeTime_ByType_withUnknownType2.so Surv_SpaceAgeTime_ByType_withUnknownType_Const_Paed_MMC_RW2.so Surv_SpaceAgeTime_ByType_withUnknownType_Const_Paed_MMC_RW_MMC2.so Surv_SpaceAgeTime.so Surv_SpaceAgeTime_RW.so 

	# Comment here preserves the prior tab
Surv_SpaceAgeTime_ByType_withUnknownType_Const_Paed_MMC.so: Surv_SpaceAgeTime_ByType_withUnknownType_Const_Paed_MMC.cpp
	"$(R_HOME)/bin/Rscript" --vanilla -e "TMB::compile('Surv_SpaceAgeTime_ByType_withUnknownType_Const_Paed_MMC.cpp', framework = 'CppAD', flags = '-O3', safeunload = FALSE, safebounds = FALSE)"

Surv_SpaceAgeTime_ByType_withUnknownType_Const_Paed_MMC_RW.so: Surv_SpaceAgeTime_ByType_withUnknownType_Const_Paed_MMC_RW.cpp
	"$(R_HOME)/bin/Rscript" --vanilla -e "TMB::compile('Surv_SpaceAgeTime_ByType_withUnknownType_Const_Paed_MMC_RW.cpp', framework = 'CppAD', flags = '-O3', safeunload = FALSE, safebounds = FALSE)"
	
Surv_SpaceAgeTime_ByType_withUnknownType_Const_Paed_MMC2.so: Surv_SpaceAgeTime_ByType_withUnknownType_Const_Paed_MMC2.cpp
	"$(R_HOME)/bin/Rscript" --vanilla -e "TMB::compile('Surv_SpaceAgeTime_ByType_withUnknownType_Const_Paed_MMC2.cpp', framework = 'CppAD', flags = '-O3', safeunload = FALSE, safebounds = FALSE)"

Surv_SpaceAgeTime_ByType_withUnknownType.so: Surv_SpaceAgeTime_ByType_withUnknownType.cpp
	"$(R_HOME)/bin/Rscript" --vanilla -e "TMB::compile('Surv_SpaceAgeTime_ByType_withUnknownType.cpp', framework = 'CppAD', flags = '-O3', safeunload = FALSE, safebounds = FALSE)"

Surv_SpaceAgeTime_ByType_withUnknownType_RW.so: Surv_SpaceAgeTime_ByType_withUnknownType_RW.cpp
	"$(R_HOME)/bin/Rscript" --vanilla -e "TMB::compile('Surv_SpaceAgeTime_ByType_withUnknownType_RW.cpp', framework = 'CppAD', flags = '-O3', safeunload = FALSE, safebounds = FALSE)"

Surv_SpaceAgeTime_ByType_withUnknownType_RW2.so: Surv_SpaceAgeTime_ByType_withUnknownType_RW2.cpp
	"$(R_HOME)/bin/Rscript" --vanilla -e "TMB::compile('Surv_SpaceAgeTime_ByType_withUnknownType_RW2.cpp', framework = 'CppAD', flags = '-O3', safeunload = FALSE, safebounds = FALSE)"

Surv_SpaceAgeTime_ByType_withUnknownType_RW_MMC2.so: Surv_SpaceAgeTime_ByType_withUnknownType_RW_MMC2.cpp
	"$(R_HOME)/bin/Rscript" --vanilla -e "TMB::compile('Surv_SpaceAgeTime_ByType_withUnknownType_RW_MMC2.cpp', framework = 'CppAD', flags = '-O3', safeunload = FALSE, safebounds = FALSE)"

Surv_SpaceAgeTime_ByType_withUnknownType2.so: Surv_SpaceAgeTime_ByType_withUnknownType2.cpp
	"$(R_HOME)/bin/Rscript" --vanilla -e "TMB::compile('Surv_SpaceAgeTime_ByType_withUnknownType2.cpp', framework = 'CppAD', flags = '-O3', safeunload = FALSE, safebounds = FALSE)"

Surv_SpaceAgeTime_ByType_withUnknownType_Const_Paed_MMC_RW2.so: Surv_SpaceAgeTime_ByType_withUnknownType_Const_Paed_MMC_RW2.cpp
	"$(R_HOME)/bin/Rscript" --vanilla -e "TMB::compile('Surv_SpaceAgeTime_ByType_withUnknownType_Const_Paed_MMC_RW2.cpp', framework = 'CppAD', flags = '-O3', safeunload = FALSE, safebounds = FALSE)"

Surv_SpaceAgeTime_ByType_withUnknownType_Const_Paed_MMC_RW_MMC2.so: Surv_SpaceAgeTime_ByType_withUnknownType_Const_Paed_MMC_RW_MMC2.cpp
	"$(R_HOME)/bin/Rscript" --vanilla -e "TMB::compile('Surv_SpaceAgeTime_ByType_withUnknownType_Const_Paed_MMC_RW_MMC2.cpp', framework = 'CppAD', flags = '-O3', safeunload = FALSE, safebounds = FALSE)"
	
Surv_SpaceAgeTime.so: Surv_SpaceAgeTime.cpp 
	"$(R_HOME)/bin/Rscript" --vanilla -e "TMB::compile('Surv_SpaceAgeTime.cpp', framework = 'CppAD', flags = '-O3', safeunload = FALSE, safebounds = FALSE)"

Surv_SpaceAgeTime_RW.so: Surv_SpaceAgeTime_RW.cpp 
	"$(R_HOME)/bin/Rscript" --vanilla -e "TMB::compile('Surv_SpaceAgeTime_RW.cpp', framework = 'CppAD', flags = '-O3', safeunload = FALSE, safebounds = FALSE)"
 
clean:
	rm -rf *o
