
ifdef SingularPath
  singular$O : ExtraCXXFLAGS = -I${SingularPath}/include -DNDEBUG
  LDFLAGS += -L${SingularPath}/lib -Wl,-rpath,${SingularPath}/lib
  LIBS += -lsingular
else
  # disables singular interface unless configured 
  IgnoreModules += singular
endif

