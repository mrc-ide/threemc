all: threemc.so

	# Comment here preserves the prior tab
threemc.so: threemc.cpp
	"$(R_HOME)/bin/Rscript" --vanilla -e "TMB::compile('threemc.cpp','-I"src/" -O3', framework = 'TMBad', safeunload = FALSE, safebounds = FALSE)"

clean:
	rm -rf *o
