# NVIDIA GPU Computing SDK Version 2.30
ifeq ($(emu), 1)
  PROJECTS := $(shell find src -name Makefile | xargs grep -L 'USEDRVAPI')
else
  PROJECTS := $(shell find src -name Makefile)
endif

#%.ph_build : lib/libcutil.so lib/libparamgl.so lib/librendercheckgl.so lib/libdijkstra.so
%.ph_build : lib/libdijkstra.so
	make -C $(dir $*) $(MAKECMDGOALS)

%.ph_clean : 
	make -C $(dir $*) clean $(USE_DEVICE)

%.ph_clobber :
	make -C $(dir $*) clobber $(USE_DEVICE)

all:  $(addsuffix .ph_build,$(PROJECTS))
	@echo "Finished building all"

#lib/libcutil.so:
#	@make -C common
#
#lib/libparamgl.so:
#	@make -C common -f Makefile_paramgl
#
#lib/librendercheckgl.so:
#	@make -C common -f Makefile_rendercheckgl

lib/libdijkstra.so:
	@make -C src/dijkstraLib
tidy:
	@find * | egrep "#" | xargs rm -f
	@find * | egrep "\~" | xargs rm -f

clean: tidy $(addsuffix .ph_clean,$(PROJECTS))
	@make -C common clean

clobber: clean $(addsuffix .ph_clobber,$(PROJECTS))
	@make -C common clobber
