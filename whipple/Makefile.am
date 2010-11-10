SUBDIRS = common .

bin_PROGRAMS = whippleeig whipplesim whipplesteady

whippleeig_SOURCES = whippleeig.cpp
whippleeig_CPPFLAGS = -I$(top_srcdir)/whipple/common
whippleeig_LDADD = ./common/libwhipplecommon.a

whipplesim_SOURCES = whipplesim.cpp
whipplesim_CPPFLAGS = -I$(top_srcdir)/whipple/common
whipplesim_LDADD = ./common/libwhipplecommon.a

whipplesteady_SOURCES= whipplesteady.cpp
whipplesteady_CPPFLAGS = -I$(top_srcdir)/whipple/common
whipplesteady_LDADD = ./common/libwhipplecommon.a ./common/libwhipplesteady.a

python_PYTHON = plot.py plotfunctions.py
