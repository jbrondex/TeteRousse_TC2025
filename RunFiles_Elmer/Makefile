SELL = /bin/sh
RM=/bin/rm

TARGETDIR = bin
SOURCEDIR = SRC
ELMERF90 = elmerf90
libs = \
$(TARGETDIR)/USF_TR\
$(TARGETDIR)/My_Grid2DVerticalInterpolator\
$(TARGETDIR)/MyFreeSurfaceSolver\
$(TARGETDIR)/USF_LimitersZeroOne\
$(TARGETDIR)/MyComputeDevStress\
$(TARGETDIR)/My_USF_Damage\
$(TARGETDIR)/ExportSurfaceVelocity\
$(TARGETDIR)/My_ParticleSolver\
$(TARGETDIR)/SaveLine_Fixed


all: $(libs)
clean: 
	$(RM) -f $(libs)
$(TARGETDIR)/USF_TR: $(SOURCEDIR)/USF_TR.F90
	$(ELMERF90)     $(SOURCEDIR)/USF_TR.F90 -o $(TARGETDIR)/USF_TR

$(TARGETDIR)/My_Grid2DVerticalInterpolator: $(SOURCEDIR)/My_Grid2DVerticalInterpolator.F90
	$(ELMERF90)     $(SOURCEDIR)/My_Grid2DVerticalInterpolator.F90 -o $(TARGETDIR)/My_Grid2DVerticalInterpolator

$(TARGETDIR)/MyFreeSurfaceSolver: $(SOURCEDIR)/FreeSurfaceSolver.F90
	$(ELMERF90)     $(SOURCEDIR)/FreeSurfaceSolver.F90 -o $(TARGETDIR)/MyFreeSurfaceSolver

$(TARGETDIR)/USF_LimitersZeroOne: $(SOURCEDIR)/USF_LimitersZeroOne.F90
	$(ELMERF90)     $(SOURCEDIR)/USF_LimitersZeroOne.F90 -o $(TARGETDIR)/USF_LimitersZeroOne

$(TARGETDIR)/MyComputeDevStress: $(SOURCEDIR)/PorousMaterialModels.F90 $(SOURCEDIR)/MyComputeDevStress.F90 
	$(ELMERF90)     $(SOURCEDIR)/PorousMaterialModels.F90 $(SOURCEDIR)/MyComputeDevStress.F90 -o $(TARGETDIR)/MyComputeDevStress

$(TARGETDIR)/My_USF_Damage: $(SOURCEDIR)/My_USF_Damage.F90 
	$(ELMERF90)     $(SOURCEDIR)/My_USF_Damage.F90 -o $(TARGETDIR)/My_USF_Damage

$(TARGETDIR)/ExportSurfaceVelocity: $(SOURCEDIR)/ExportSurfaceVelocity.F90 
	$(ELMERF90)     $(SOURCEDIR)/ExportSurfaceVelocity.F90 -o $(TARGETDIR)/ExportSurfaceVelocity

$(TARGETDIR)/My_ParticleSolver: $(SOURCEDIR)/ParticleUtils_MultiSource.F90 $(SOURCEDIR)/ParticleAdvector_MultiSource.F90
	$(ELMERF90)     $(SOURCEDIR)/ParticleUtils_MultiSource.F90 $(SOURCEDIR)/ParticleAdvector_MultiSource.F90 -o $(TARGETDIR)/My_ParticleSolver

$(TARGETDIR)/SaveLine_Fixed: $(SOURCEDIR)/SaveUtils_Fixed.F90 $(SOURCEDIR)/SaveLine_Fixed.F90
	$(ELMERF90)     $(SOURCEDIR)/SaveUtils_Fixed.F90 $(SOURCEDIR)/SaveLine_Fixed.F90 -o $(TARGETDIR)/SaveLine_Fixed


