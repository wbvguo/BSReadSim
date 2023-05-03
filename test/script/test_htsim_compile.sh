# rrcut
c++ -g -O3 -fpermissive   src/rrcut.cpp src/haplo.cpp src/struct.cpp -o rrcut -Bstatic -lz  -I/home/wbguo/iproject/BSReadSim/HTSLIB/htslib/ -L/home/wbguo/iproject/BSReadSim/HTSLIB/ -lhts -Wl,-rpath /home/wbguo/iproject/BSReadSim/HTSLIB -lgsl -lgslcblas

# hstim
c++ -g -O3 -fpermissive   src/htsim_test.cpp src/haplo.cpp src/methdb.cpp src/mode.cpp src/struct.cpp -o htsim -Bstatic -lz  -I/home/wbguo/iproject/BSReadSim/HTSLIB/htslib/ -L/home/wbguo/iproject/BSReadSim/HTSLIB/ -lhts -Wl,-rpath /home/wbguo/iproject/BSReadSim/HTSLIB -lgsl -lgslcblas
