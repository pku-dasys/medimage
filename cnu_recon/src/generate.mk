# Makefile template by Guojie Luo

#ROOT_CCFLAGS := $(ROOT_CCFLAGS)

MS3D_OUT = generate
MS3D_TEST = circle sphere ellipse taiji
MS3D_SRC = generate.cpp tracing.cpp utility.cpp ct3d.cpp $(foreach TEST, $(MS3D_TEST), testcase/$(TEST).cpp)

MS3D_LIB =
all :: $(ROOT_BIN)/$(MS3D_OUT)

MS3D_LDFLAGS = $(MS3D_LIB) $(ROOT_LDFLAGS)
MS3D_OBJ = $(MS3D_SRC:%.cpp=$(ROOT_BUILD)/%.o)
MS3D_DEP = $(patsubst %,%.d,$(basename $(MS3D_OBJ)))

$(MS3D_OUT) :: $(MS3D_OBJ)
	$(ROOT_LNK) -o $(ROOT_BIN)/$(MS3D_OUT) $^ $(MS3D_LDFLAGS)

$(ROOT_BIN)/$(MS3D_OUT) :: $(MS3D_OBJ)
	$(ROOT_LNK) -o $@ $^ $(MS3D_LDFLAGS)

clean ::
	rm -f $(MS3D_DEP)
	rm -f $(MS3D_OBJ)
	rm -f $(ROOT_BIN)/$(MS3D_OUT)

-include $(MS3D_DEP)
