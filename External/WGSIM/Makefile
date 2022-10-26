##/*************************************************************************************
##                           The MIT License
##
##   Permission is hereby granted, free of charge, to any person obtaining
##   a copy of this software and associated documentation files (the
##   "Software"), to deal in the Software without restriction, including
##   without limitation the rights to use, copy, modify, merge, publish,
##   distribute, sublicense, and/or sell copies of the Software, and to
##   permit persons to whom the Software is furnished to do so, subject to
##   the following conditions:
##
##   The above copyright notice and this permission notice shall be
##   included in all copies or substantial portions of the Software.
##
##   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
##   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
##   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
##   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
##   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
##   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
##   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
##   SOFTWARE.
##
##Contacts: Wenbin Guo <wbguo@ucla.edu>
##*****************************************************************************************/
# $(shell export LD_LIBRARY_PATH=$(LD_LIBRARY_PATH):../HTSLIB/)
# $(shell echo $(LD_LIBRARY_PATH))
EXE=		wgsim
LIBS=		-Bstatic -lz 
HTSLIB=		-I../HTSLIB/htslib/ -L../HTSLIB/ -lhts -Wl,-rpath ../HTSLIB
CXX=c++

CXXFLAGS=	-g -O3 -fpermissive $(ARCH_FLAGS) #-Wall ##-xSSE2

wgsim:
	$(CXX) $(CXXFLAGS) src/wgsim.cpp -o $(EXE) $(LIBS) $(HTSLIB) 


clean:
	rm wgsim

depend:
	(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CXXFLAGS) $(CPPFLAGS) -I. -- src/*.cpp)
