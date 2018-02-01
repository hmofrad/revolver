#!/usr/bin/make
# Makefile: Java makefile
# (c) Mohammad HMofrad, 2017
# (e) mohammad.hmofrad@pitt.edu

WEBGRAPH := /afs/cs.pitt.edu/usr0/hasanzadeh/private/webgraph
JAVAC    := javac
JAVA     := java
CP       := -cp
SOURCE   := LearningAutomata
JFLAGS   := -Xlint:deprecation 

all: install
	$(JAVA) $(CP) .:$(WEBGRAPH)/\* $(SOURCE)
install:
	$(JAVAC) $(CP) $(WEBGRAPH)/\* $(SOURCE).java
clean:
	rm -rf *.class
