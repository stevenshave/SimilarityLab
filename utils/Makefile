CXX = g++
CXXFLAGS = -std=c++0x -Wall -O3

OBJECTS = usrcat_binary_reader_similarity_lab.o

usrcat_binary_reader_similarity_lab: $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^

$(OBJECTS): usrcat_binary_reader_similarity_lab.cpp

clean :
	rm -f usrcat_binary_reader_similarity_lab *.o