## Declare which search path should be used:
#vpath %.h include 
## After 'include' type ' : otherfoldername : anotherone' etc. if 
## more folders may containn .h files.
#vpath %.cpp src
#vpath lib%.so lib
#vpath lib%.a lib

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#                                                                     #
#                       V A R I A B L E S                             #
#                                                                     #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
SHELL = /bin/sh
CC = gcc
CFLAGS = -Wall -O3
SOURCES = $(wildcard *.c)
OBJECTS = $(SOURCES: .c=.o)
DFILES = $(SOURCES: .c=.d)
LDFLAGS = 
LIBS =  -L/usr/lib64/nvidia -lGL -lGLU -lglut -lm
EXECUTABLE = parts2D

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#                                                                     #
#                        C O M M A N D S                              #
#                                                                     #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(CPPFLAGS)  $^ -o $@ $(LIBS)
	
%.o: %.c
	@echo " CPP   $@"
	$(CC) -c $(CLAGS) $< -o $@

.PHONY: clean	
	
%.o: %.h
	@echo "   CPP   $@"
	$(QUIET)$(CC) $(CFLAGS) -c $< -o $@

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#                                                                     #
#                  C L E A N   C O M M A N D S                        #
#                                                                     #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


clean:
	rm -rf *.o *~ $(EXECUTABLE)
