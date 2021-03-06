# cosp_mmf
Program to run the CFMIP Observation Simulator Package (COSP) on output from SP-CAM

Notes
-----
2015/09/12: The following note is no longer valid. Using the COSP source in
    ${HOME}/cosp/cosp-R83, and this directory only contains the code specific
    to the cosp-mmf driver.

2015/09/12: v6 removed support for passing flags to mimic GCM assumptions.
    These are now handled elsewhere. 
    This made more sense, so that I can actually see what is
    happening to those inputs, and do not have to mess with the code that runs
    COSP. Thus, the inputs are more easily analyzed, rather than trying to
    output them from this driver.

Edits made to COSP code base
----------------------------

cosp-1.4/cosp_types.F90: Roj added code to allow not updating look-up table on each iteration. This was needed for code additions to use adaptive grid version of the MMF (multiple lines).

cosp-1.4/cosp.F90: BRH added code to zero out sgx%frac_out when there is no cloud in CRM mode. sgx%frac_out is not used in CRM mode, so this edit is not really necessary (line 580).

cosp-1.4/quickbeam/radar_simulator.f90: BRH fixed array out of bounds access (line 505).

cmor/COSP_table_2D_roj: BRH commented out valid_min and valid_max from latitude and longitude variable definitions because this caused problems with CMOR (lines 274-275, 297-298).


Other notes
-----------

Added the following F90FLAGs in Makefile:
   -fbounds-check: this revealed the array out of bounds error I was getting in the radar_simulator code
   -ffree-line-length-none: prevents truncation of a few lines in cosp_io.F90
   -ffixed-line-length-none: added just in case.
