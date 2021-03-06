SHELL=cmd.exe 

OBJS=simulator.obj particle.obj voxel.obj field3d.obj scalarfield3d.obj \
     vectorfield3d.obj mtrand.obj util.obj geometry.obj \
     spheremovingobj.obj blendparticles.obj cubesource.obj \
     Vectors.obj CIsoSurface.obj transform.obj sph.obj signeddistancecontainer.obj \
     cylindricalsource.obj cylindricalcontainer.obj cubecontainer.obj cylindricalmovingobj.obj \
     vortexparticles.obj statictarget.obj wallclocktime.obj flip.obj \
     mesh.obj physicsworld.obj configfile.obj aliaswavefrontobj.obj \
     bounding_box_tree.obj mesh_query.obj predicates.obj renderlevelset.obj
     
OBJS0=particle.obj voxel.obj field3d.obj scalarfield3d.obj \
     vectorfield3d.obj mtrand.obj util.obj geometry.obj \
     spheremovingobj.obj blendparticles.obj cubesource.obj cubecontainer.obj \
     Vectors.obj CIsoSurface.obj transform.obj sph.obj vortexparticles.obj \
     statictarget.obj wallclocktime.obj flip.obj \
     mesh.obj physicsworld.obj configfile.obj aliaswavefrontobj.obj \
     bounding_box_tree.obj mesh_query.obj predicates.obj renderlevelset.obj
     
OBJS1=check.obj 

OBJS2=metaball.obj 

OBJS3=prmanimp.obj voxel.obj mtrand.obj
   
OBJS4=reconstruct.obj  

#LIBS=mkl_c.lib libm.lib impi.lib impicxx.lib tbb.lib

#Windows XP x86
#LIBS=mkl_c.lib levelset.lib cuda.lib cudart.lib cutil32.lib

#Windows XP x64 
LIBS=mkl_intel_lp64.lib mkl_intel_thread.lib mkl_core.lib libiomp5md.lib \
     levelset.lib cuda.lib cudart.lib cutil64.lib

EXES1=simulator.exe
EXES4=check.exe
EXES3=metaball.exe
EXES2=reconstruct.exe
EXES5=prmanimp.dll 

all: $(OBJS) $(EXES1) $(EXES2) $(EXES3) $(EXES4) $(EXES5)

$(EXES1):%.exe:%.obj	
#   $(LINK32) /NODEFAULTLIB:libcmt.lib /LARGEADDRESSAWARE /OUT:$@ $(OBJS) $(LIBS)
	$(LINK32) /NODEFAULTLIB:libcmt.lib /OUT:$@ $(OBJS) $(LIBS)
	
$(EXES2):%.exe:%.obj
	$(LINK32) /NODEFAULTLIB:libcmt.lib /OUT:$@ $^ $(OBJS0) $(LIBS)

$(EXES3):%.exe:%.obj
	$(LINK32) /NODEFAULTLIB:libcmt.lib /OUT:$@ $^ $(OBJS0) $(LIBS) 

$(EXES4):%.exe:%.obj
	$(LINK32) /NODEFAULTLIB:libcmt.lib /OUT:$@ $^ $(OBJS0) $(LIBS)
	
$(EXES5):%.dll:%.obj
	$(CPP) /LD $(OBJS3)    	
	
$(OBJS):%.obj:%.cpp
#	$(CPP) /c /Qc99 /DWIN32 /DPOSITIVE_PARTICLES /DAIR_DIV_FREE /w $^
	$(CPP) /c /Qc99 /MD /O3 /DWIN32 /DCUDA /DSPMD /DCORE_SOURCE /DPOSITIVE_PARTICLES /w $^

$(OBJS1):%.obj:%.cpp
#	$(CPP) /c /Qc99 /DPOSITIVE_PARTICLES /DAIR_DIV_FREE /w $^
	$(CPP) /c /Qc99 /DWIN32 /MD /w $^
	
$(OBJS2):%.obj:%.cpp
#	$(CPP) /c /Qc99 /DPOSITIVE_PARTICLES /DAIR_DIV_FREE /w $^
	$(CPP) /c /Qc99 /DWIN32 /MD /w $^

$(OBJS3):%.obj:%.cpp
#	$(CPP) /c /Qc99 /DPOSITIVE_PARTICLES /DAIR_DIV_FREE /w $^
	$(CPP) /c /Qc99 /DWIN32 /MD /O3 /w $^

$(OBJS4):%.obj:%.cpp
#	$(CPP) /c /Qc99 /DPOSITIVE_PARTICLES /DAIR_DIV_FREE /w $^
	$(CPP) /c /Qc99 /DWIN32 /MD /I"C:\Program Files\Pixar\devkit\include" /w $^


	
clean:
	del *.obj
	del *.exe