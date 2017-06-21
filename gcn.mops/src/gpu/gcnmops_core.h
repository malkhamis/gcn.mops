extern "C" void gcnmops_core(unsigned int grid_size, unsigned int block_size, struct X_norm d_Xchunk, struct gcnmops_args d_args, struct
	gcnmops_tmpVars d_workspace, struct gcnmops_buf d_buf);
extern "C" void cnmops_leqMRC(struct X_norm h_Xchunk, struct gcnmops_args h_args, struct gcnmops_buf h_buf, bool returnPosterior);

extern "C" void cpyArgsToDevSym(struct gcnmops_args h_args);
extern "C" void cpyDataToDev(struct X_norm d_Xchunk, struct X_norm h_Xchunk, unsigned int nRanges_tot);
extern "C" void cpyResToHost(struct gcnmops_buf h_buf, struct gcnmops_buf d_buf);
extern "C" void devSync();
