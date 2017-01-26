# source files 
srcs =
srcs +=constants.cpp
srcs +=SOM.cpp
srcs +=main.cpp

################################################################################

# object files

objs = $(srcs:.cpp=.o)
deps = $(srcs:.cpp=.d)

################################################################################

# list of directories to search
# for c++ source files
srcDirs = 
srcDirs += ./src 
srcDirs += ./include

################################################################################

# header files
headers = 

################################################################################

vpath %.cpp $(srcDirs)

################################################################################

# Program executable file name 
exe = main

################################################################################

#Will not echo the action as its processed. Same as command-line option ``-s''. 
#By default make will echo the action to stdout  prior to invocation.
#.SILENT:

#Ignores all return codes from actions. Same as command-line option ``-l''. 
#By default make will stop processing whenever a non-zero return status is 
#received from an action.
#.IGNORE:

.SUFFIXES:
.SUFFIXES: .cpp .o 

################################################################################

# Program for compiling C++ programs
CXX = g++ 

# Extra flags to give to the C++ compiler. 
# CXXFLAGS = -g3 -Werror -Wall -pedantic-errors -W \
# 	-Wshadow -Wpointer-arith -Wcast-qual -Wcast-align \
# 	-Wwrite-strings -fno-common

#CXXFLAGS = -g3 -Wall -Werror -W -DDEBUG
# CXXFLAGS = -g3 -DDEBUG	
 CXXFLAGS = -g3 

################################################################################

# Extra flags to give to the C preprocessor and programs 
# that use it (the C and Fortran compilers). 
#CPPFLAGS = $(addprefix -I ,$(srcDirs))
CPPFLAGS = $(addprefix -I ,/usr/include/opencv)

################################################################################
# http://www.mail-archive.com/help-make@gnu.org/msg01513.html
# Extra flags to give to compilers when they are supposed to invoke the linker
#LDFLAGS = $(addprefix -L ,~/phdWork/nnlib-1.1/lib/)

LOADLIBES = -lcv -lcxcore -lhighgui -lcvaux -lgsl -lgslcblas

LDLIBS = 

################################################################################

# top-level rule to make everything
all: $(exe)

################################################################################

#
# LINK.cc = $(CXX) $(CXXFLAGS) $(CPPFLAGS) $(LDFLAGS) $(TARGET_ARCH)
# LINK.cpp = $(LINK.cc)
#

# pattern rule to link the program
$(exe): $(objs)
	$(LINK.cpp) $^ $(LOADLIBES) $(LDLIBS) -o $@

################################################################################

#
# COMPILE.cc = $(CXX) $(CXXFLAGS) $(CPPFLAGS) $(TARGET_ARCH) -c
# COMPILE.cpp = $(COMPILE.cc)
#

#
# $^ - stands for all the dependencies of the current rule 
# i.e. all the stuff on right hand side of the colon (:)
#
# $< - the first thing in $^
# $@ - the name of the file to be ``made''
# $? - the set of dependent names that are younger than the target
# $* - the shared prefix of the target and dependent - only for suffix rules
# $$ - escapes macro substitution, returns a single ``$''.
#
# http://www.gnu.org/software/autoconf/manual/make/Automatic-Variables.html
#

# pattern rule to create object files
%.o: %.cpp
	$(COMPILE.cpp) -MMD $< -o $@

-include $(deps)
################################################################################

.PHONY: all clean tags

################################################################################

clean:
	$(RM) *.dat *.out *.jpeg $(objs) $(deps) $(exe) 

################################################################################

tags:
	ctags -e -R --langmap=c++:+.tpp --c++-kinds=+p

################################################################################
