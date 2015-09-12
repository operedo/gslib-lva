default:
	cd Boost_dijkstra; make clean; make; cd ..
	cd CUDA_dijkstra; make clean; make; cd ..
	cd gamv-lva; mkdir bin; cd src; make clean; make; cp GAMV_LVA ../bin; cd ../..
	cd kt3d-lva; mkdir bin; cd src; make clean; make; cp KT3D_LVA ../bin; cd ../..
	cd sgs-lva; mkdir bin; cd src; make clean; make; cp SGS_LVA ../bin; cd ../..

clean:
	cd Boost_dijkstra; make clean; cd ..
	cd CUDA_dijkstra; make clean; cd ..
	cd gamv-lva; cd src; make clean; cd ../..
	cd kt3d-lva; cd src; make clean; cd ../..
	cd sgs-lva; cd src; make clean; cd ../..

install:
	cp Boost_dijkstra/Boost_dijkstra Boost_dijkstra/Boost_dijkstra_openmp $(HOME)/bin	
	cp CUDA_dijkstra/bin/linux/release/CUDA_dijkstra $(HOME)/bin	
	cp gamv-lva/bin/GAMV_LVA $(HOME)/bin	
	cp kt3d-lva/bin/KT3D_LVA $(HOME)/bin	
	cp sgs-lva/bin/SGS_LVA $(HOME)/bin	

test:
	cd gamv-lva; cd test; sh run.sh; cd ../..
	cd kt3d-lva; cd test; sh run.sh; cd ../..
	cd sgs-lva; cd test; sh run.sh; cd ../..

