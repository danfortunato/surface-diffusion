-include make.inc

default: mex

# Rebuild BandedBatch.m and gateway.cpp via MWrap
mwrap: BandedBatch.m private/gateway.cpp
BandedBatch.m: BandedBatch.mw
	$(MWRAP) -mex gateway -m BandedBatch.m BandedBatch.mw
private/gateway.cpp: BandedBatch.mw
	$(MWRAP) -mex gateway -c private/gateway.cpp BandedBatch.mw

# Rebuild gateway.mex* via MEX
mex: private/gateway.mex*
private/gateway.mex*: private/gateway.cpp
	cd private && \
	$(MEX) gateway.cpp $(INCDIR) $(MEXFLAGS) $(LIBDIR) $(LIB) -output gateway

# Remove the MEX interface
clean:
	cd private && rm -f gateway.mex*

# Remove the MEX interface, MATLAB caller, and MEX file
# Note: You will need MWrap to rebuild the deleted files!
mwrapclean: clean
	rm -f BandedBatch.m private/gateway.cpp

.PHONY: mwrap mex clean mwrapclean
