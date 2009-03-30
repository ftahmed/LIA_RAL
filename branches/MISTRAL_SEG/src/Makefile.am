## Process this file with automake to produce Makefile.in

## Created by Anjuta
##	hmm.cpp hmm.h
##	accumulateStat.cpp accumulateStat.h

LDADD    = -L$(top_builddir)/src -llium_seg @liatools_lib@ @alize_lib@  
INCLUDES = -I$(top_builddir)/src @liatools_include@ @alize_include@

lib_LIBRARIES = liblium_seg.a
liblium_seg_a_SOURCES = cluster.cpp cluster.h \
	clusterServer.cpp clusterServer.h \
	segment.cpp segment.h \
	hclust.cpp hclust.h \
	bichclust.cpp bichclust.h \
	biclclust.cpp biclclust.h \
	clrhclust.cpp clrhclust.h \
	distance.cpp distance.h \
	liaToolsWrapper.cpp liaToolsWrapper.h

bin_PROGRAMS = \
	segInit\
	seg   \
	clust   \
	trainEM   \
	trainInit \
	decode

decode_SOURCES = decode.cpp
decode_DEPENDENCIES = liblium_seg.a

trainInit_SOURCES = trainInit.cpp
trainInit_DEPENDENCIES = liblium_seg.a

trainEM_SOURCES = trainEM.cpp
trainEM_DEPENDENCIES = liblium_seg.a

clust_SOURCES = clust.cpp
clust_DEPENDENCIES = liblium_seg.a

seg_SOURCES = seg.cpp
seg_DEPENDENCIES = liblium_seg.a

segInit_SOURCES = segInit.cpp
segInit_DEPENDENCIES = liblium_seg.a

