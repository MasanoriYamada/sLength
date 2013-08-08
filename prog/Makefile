INCFLAGS = 
LDFLAGS =
LIBLIST = -lgfortran
INCLST = ../include/analys.h ../include/fit.h    ../include/io.h     ../include/jack.h
OPTION = -g

CC = gcc
CXX = g++
GF = gfortran

TARGET_LIST=fit_a.out
############SRC LIST#################
SRC_FIT=fit_a.cpp
#SRC_FIT_FUNCTION=fit_function.cpp
SRC_1=DCOVSRT.f 
SRC_2=DGAULEG.f 
SRC_3=DGAUSSJ.f 
SRC_4=DMRQCOF.f 
SRC_5=DMRQMIN.f
SRC_LIST=$(SRC_1) $(SRC_2) $(SRC_3) $(SRC_4) $(SRC_5) 

#############OBJECT LIST################# 
OBJ_FIT=fit_a.o
#OBJ_FIT_FUNCTION=fit_function.o
OBJ_1=DCOVSRT.o 
OBJ_2=DGAULEG.o
OBJ_3=DGAUSSJ.o
OBJ_4=DMRQCOF.o 
OBJ_5=DMRQMIN.o
OBJ_LIST=$(OBJ_1) $(OBJ_2) $(OBJ_3) $(OBJ_4) $(OBJ_5) 

all:$(TARGET_LIST)

$(TARGET_LIST):$(OBJ_LIST) $(OBJ_FIT) $(OBJ_FIT_FUNCTION)
	$(CXX) $(OPTION) -o $@ $(OBJ_FIT)  $(OBJ_FIT_FUNCTION) $(OBJ_LIST) $(LIBLIST)
$(OBJ_FIT):$(SRC_FIT) $(INCLST)
	$(CXX) $(OPTION) -c -o $@ $(SRC_FIT)
$(OBJ_FIT_FUNCTION):$(SRC_FIT_FUNCTION)
	$(CXX) $(OPTION) -c -o $@ $(SRC_FIT_FUNCTION)
$(OBJ_LIST):$(SRC_LIST)
	$(GF) $(OPTION) -c $(SRC_LIST)


clean: 
	@rm -rf *.out *.o *~ *.bak *.deps
dataclean:
	@rm -rf ../out/*
