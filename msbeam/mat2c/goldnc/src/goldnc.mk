# Makefile template by Guojie Luo

#ROOT_CCFLAGS := $(ROOT_CCFLAGS)

AFFINA_OUT = $(ROOT_NAME)
AFFINA_LIB =
AFFINA_SRC = main.c \
						proj.c

all :: $(ROOT_BIN)/$(AFFINA_OUT)

AFFINA_LDFLAGS = $(AFFINA_LIB) $(ROOT_LDFLAGS)
AFFINA_OBJ = $(AFFINA_SRC:%.c=$(ROOT_BUILD)/%.o)
AFFINA_DEP = $(patsubst %,%.d,$(basename $(AFFINA_OBJ)))

$(AFFINA_OUT) :: $(AFFINA_OBJ)
	$(ROOT_LNK) -o $(ROOT_BIN)/$(AFFINA_OUT) $^ $(AFFINA_LDFLAGS)

$(ROOT_BIN)/$(AFFINA_OUT) :: $(AFFINA_OBJ)
	$(ROOT_LNK) -o $@ $^ $(AFFINA_LDFLAGS)

clean ::
	rm -f $(AFFINA_DEP)
	rm -f $(AFFINA_OBJ)
	rm -f $(ROOT_BIN)/$(AFFINA_OUT)

-include $(AFFINA_DEP)
