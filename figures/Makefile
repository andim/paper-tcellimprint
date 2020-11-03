#### Makefile for Latex file ####
.PHONY: clean

MAKEDIR := $(shell pwd)
IMGDIR = $(MAKEDIR)/pngs/
SKETCHDIR = $(MAKEDIR)


SVGSRC	:= $(notdir $(wildcard $(SKETCHDIR)/*.svg))
SVG	:= $(addprefix $(IMGDIR)/, $(SVGSRC:%.svg=%.png))

RM = /bin/rm -f
OUT = /dev/null 2>&1

all: $(wildcard $(PNGDIR)/*png) $(SVG)
	echo $(SVG) $(PNG)

$(IMGDIR)/%.png : $(SKETCHDIR)/%.svg
	inkscape $< --export-png=$@ --export-area-drawing -d 300

