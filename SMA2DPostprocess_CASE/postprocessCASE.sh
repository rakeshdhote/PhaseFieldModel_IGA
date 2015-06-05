clear
date

echo "--Post-processing MPI code--"
rm -rf postp
mkdir postp

d1=$PWD
d2=$PWD/postp/

ln -s $d1/restart.* $d2
ln -s $d1/SMA2DPostprocess_CASE/* $d2
ln -s $d1/step.dat $d2
ln -s $d1/scalegeom.dat $d2
ln -s $d1/meshu*.dat $d2

cd $d2

vi step.dat
vi scalegeom.dat

COMPILER=ifort
#COMPILER=gfortran -fbounds-check

$COMPILER -o out aAdjKeep.f basisfuns.f input_3D.f database_3D.f eval_SHAPE_3D.f dersbasisfuns.f results_3D.f VizGMV_3D.f genIEN_INN_3D.f

echo "Running program"
./out

vi SMA2DF.case

mkdir postptransfer
cd postptransfer

mv ../*.res .
mv ../*.case .
mv ../*.geo .

echo "Completed postprocessing.Current directory: "
echo $PWD
echo "== Open paraview in Visualization System == "
paraview




