
#ifndef _REINIT_CUDA_H_
#define _REINIT_CUDA_H_

extern "C"
void ReinitializeCUDA(float *phi, const char *obj_cpu, const char *movobj_cpu, const char *source,
					  bool init,  int iterations, int I, int J, int K,
					   float delta, float dtau, float eps, float limit,
					   int dim[]);
extern "C"
void ReinitializeCUDANeeded(float *phi, const char *obj_cpu, const char *movobj_cpu, const char *need,
							int iterations, int I, int J, int K,
 							float delta, float dtau, int dim[]);

extern "C"
void ReinitializeCUDANeededRLE(float *phi,  const char *need,
							int iterations, int I, int J, int K,
 							float delta, float dtau, int dim[]);

extern "C"
void ExtrapolatePhiCUDA(float *phi, const float *phi_obj, const char *obj_cpu, 
					    int I, int J, int K,
 				      float delta, float dtau, float inf, float limit,
					  int dim[]);

extern "C"
void ExtrapolateVelCUDA(float *vel, const float *phi, const char *needed_cpu, const char *valid_cpu, 
					  int I, int J, int K,
 				      float delta, float dtau, float inf, float limit, 
					  int dim[]);

#endif // #ifndef _REINIT_CUDA_H_