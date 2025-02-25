cd /home/niraj/1TB-Disc/Workspace-IITG/Packages/Install/autodipole/AutoDipole_V1.2.3/process/Proc_gg_uuxg
g77 -c check.f matrix.f dipole.f reducedm.f cmatrix.f couplings.f collinear.f rambo.f
ar rcs libdipole.a check.o dipole.o cmatrix.o couplings.o reducedm.o rambo.o collinear.o matrix.o
cd /home/niraj/OneDrive/WorkSpaceIITG/Generalfiles/integrators/DY1/ggaa_working_lib/dipole_standalone 
make clean && make && ./runDipole
