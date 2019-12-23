MNS_DG=../bin/MnsDg

default: houston.1.bb

houston.1.bb: houston.mesh ../bin/MnsDg
	$(MNS_DG) $< 2 -o $@

clean:
	\rm houston.1.bb houston.1.mesh houston.1.error.bb houston.1.error.mesh

