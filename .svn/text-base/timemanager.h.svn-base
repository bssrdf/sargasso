#ifndef TIMEMANAGER_H_
#define TIMEMANAGER_H_

#define FRAMERATE  24
#define EPS  1.e-5

class TimeManager{
public:
	TimeManager(float st, float en, float slr, float rs, int si, int fix=0, float fixt=0.f)
	: startTime(st), endTime(en), slToEulerRatio(slr), reSeedRatio(rs),
	  subCyclingIters(si), subCyclingFixedIters(si), fixedSteps(fix) {
		frameDt = 1.f / FRAMERATE;
		dt = frameDt;
		frameStart  = st;
		frameTime   = st;
		simulTime   = st;
		reSeedTime  = st;
		n = 0;
		if(fix > 0){
			dt = fixt;
		}
	}

	float GetDt() const{
		return dt;
	}

	float GetSubcyclingDt() const{
		return dt / subCyclingIters;
	}
	
	float GetSubcyclingIterations() const{
		return subCyclingIters;
	}

	bool TimeToOutputFrame() const{
		printf("dt = %f, simultime = %f, frametime = %f \n", dt, simulTime, frameTime);
		if( simulTime-frameTime+EPS > 0.f )
			return true;
		else
			return false;
	}

	bool TimeToReseed() const{
		if( simulTime > EPS && simulTime-reSeedTime+EPS > 0.f  )
			return true;
		else
			return false;
	}

	void AdvanceFrameTime(){
		frameTime += frameDt;
	}

	void AdvanceReSeedTime(){
		reSeedTime += reSeedRatio*frameDt;
		printf(" Next particle reseeding time = %f \n", reSeedTime);
	}

	void AdvanceSimulationTime(float dtt){
		if(fixedSteps > 0){
			++n;
		}
		if(n >= fixedSteps){
			dt = min(frameDt, slToEulerRatio*dtt);
//			printf("A. dtt = %f, dt = %f simultime = %f frametime = %f \n", dtt, dt, simulTime, frameTime);
			if( simulTime+dt-frameTime > EPS ){
				dt = frameTime - simulTime;
				dt = min(dt, min(frameDt, slToEulerRatio*dtt));				
			}
			if(dt <= dtt)
				subCyclingIters = 1;
			else
				subCyclingIters = subCyclingFixedIters;
		}
		simulTime += dt;
//		printf("B. dtt = %f, dt = %f simultime = %f frametime = %f \n", dtt, dt, simulTime, frameTime);

	}

	bool Stop() const{
		if(simulTime < endTime)
			return false;
		else
			return true;
	}


	void Print() const{
		printf(" ********************************************************* \n ");
		printf(" Simulation starts at %f, and will stop at %f \n", startTime, endTime);
		printf(" Current simulation time step = %f, simulation time = %f \n ", dt, simulTime);
		printf(" Current frame time step = %f, next output time = %f \n ", frameDt, frameTime);
		printf(" Current particle reseeding time step = %f, next particle reseeding time = %f \n ", reSeedRatio*frameDt, reSeedTime);
		printf(" ********************************************************* \n ");
	}
	float simulTime;
private:
	float startTime, endTime;
	float dt;
	float slToEulerRatio; // typically use 4.9
	float frameDt;
	float frameStart,frameTime;
	float reSeedRatio;
	float reSeedTime;
	int subCyclingIters;
	int subCyclingFixedIters;	
	int fixedSteps;	
	int n;
	

};


#endif /*TIMEMANAGER_H_*/
