# methdb
c++ -g -O3 -fpermissive   src/methdb.cpp -o methdb -Bstatic -lz  -I/home/wbguo/iproject/BSReadSim/HTSLIB/htslib/ -L/home/wbguo/iproject/BSReadSim/HTSLIB/ -lhts -Wl,-rpath /home/wbguo/iproject/BSReadSim/HTSLIB -lgsl -lgslcblas -Wall


# hstim

c++ -g -O3 -fpermissive   src/methdb.cpp src/htsim.cpp -o htsim -Bstatic -lz  -I/home/wbguo/iproject/BSReadSim/HTSLIB/htslib/ -L/home/wbguo/iproject/BSReadSim/HTSLIB/ -lhts -Wl,-rpath /home/wbguo/iproject/BSReadSim/HTSLIB -lgsl -lgslcblas -Wall