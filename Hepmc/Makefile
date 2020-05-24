include config.in
LOCAL_BIN=bin
LOCAL_EXAMPLE=examples
LOCAL_INCLUDE=include
LOCAL_LIB=lib

#configure: configure.cc
#	g++ configure.cc -o configure
#read: read.cc
#	g++ -I${HEPMC}/${LOCAL_INCLUDE} -Iread.cc -o read
    
    
# FASTJET3 and HEPMC2.
hepqcd : hepqcd.cc
	g++ $< -o $@ -I$(FASTJET3)/$(LOCAL_INCLUDE) -I$(HEPMC2)/$(LOCAL_INCLUDE) \
	 -L$(FASTJET3)/$(LOCAL_LIB) -Wl,-rpath,$(FASTJET3)/$(LOCAL_LIB) -lfastjet\
	 -lQCDAwarePlugin -lm -lfastjettools\
	 -L$(HEPMC2)/$(LOCAL_LIB) -Wl,-rpath,$(HEPMC2)/$(LOCAL_LIB) -lHepMC\


##	  $(GZIP_INC) $(GZIP_FLAGS)
hepreader : hepreader.cc
	g++ $< -o $@ -I$(FASTJET3)/$(LOCAL_INCLUDE) -I$(HEPMC2)/$(LOCAL_INCLUDE) \
	 -L$(FASTJET3)/$(LOCAL_LIB) -Wl,-rpath,$(FASTJET3)/$(LOCAL_LIB) -lfastjet\
	 -L$(HEPMC2)/$(LOCAL_LIB) -Wl,-rpath,$(HEPMC2)/$(LOCAL_LIB) -lHepMC\

##	  $(GZIP_INC) $(GZIP_FLAGS)
hepZ : hepZ.cc 
	g++ $< -o hepZ -I$(FASTJET3)/$(LOCAL_INCLUDE) -I$(HEPMC2)/$(LOCAL_INCLUDE) \
	 -L$(HEPMC2)/$(LOCAL_LIB) -Wl,-rpath,$(HEPMC2)/$(LOCAL_LIB) -lHepMC\
	 -L$(FASTJET3)/$(LOCAL_LIB) -Wl,-rpath,$(FASTJET3)/$(LOCAL_LIB) -lfastjet\
##	  $(GZIP_INC) $(GZIP_FLAGS)



clean:
	rm -f congigure
