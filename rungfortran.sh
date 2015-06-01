clear
rm -f a.out output.txt *.csv *.mod
gfortran -g -ffpe-trap=invalid,zero,overflow common.f aAdjKeep.f maindef.f gen_FACE_IEN_3D.f e3LHS_3D.f genGP_GW.f basisfuns.f e3Rhs_3D.f genIEN_INN_3D.f SparseMatLoc_3D.f BCLhs_3D.f genSparStruc_3D.f SparsProd_3D.f input_3D.f eval_SHAPE_3D.f eval_FACE_3D.f IntElmAss_3D.f dersbasisfuns.f LocaltoGlobal_3D.f FillSparseMat_3D.f gen_EL_CON_3D_p.f gen_FACE_CON_3D_p.f driver.f solflow.f e3int.f e3STAB_3D.f SparseGMRES_3D.f constraint.f gral.f processors.f
#./a.out | tee SMA2DSerial.txt
gdb a.out
