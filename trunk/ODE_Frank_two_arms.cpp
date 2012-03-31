/*
ODEFrank Copyright (c) 2012, Alessio Rocchi
All rights reserved. Email: rocchi.alessio@gmail.com 
					 Web:https://sites.google.com/site/rocchialessio

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of Alessio Rocchi nor the
      names of his contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL Alessio Rocchi BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/


#define S_FUNCTION_LEVEL 2
#define S_FUNCTION_NAME ODE_Frank_two_arms

#define NUM_INPUTS          1
/* Input Port  0 */
#define IN_PORT_0_NAME      tau
#define INPUT_0_WIDTH       8
#define INPUT_DIMS_0_COL    1
#define INPUT_0_DTYPE       real_T
#define INPUT_0_COMPLEX     COMPLEX_NO
#define IN_0_FRAME_BASED    FRAME_NO
#define IN_0_BUS_BASED      0
#define IN_0_BUS_NAME
#define IN_0_DIMS           1-D
#define INPUT_0_FEEDTHROUGH 0
#define IN_0_ISSIGNED        0
#define IN_0_WORDLENGTH      8
#define IN_0_FIXPOINTSCALING 1
#define IN_0_FRACTIONLENGTH  9
#define IN_0_BIAS            0
#define IN_0_SLOPE           0.125

#define NUM_OUTPUTS          3
/* Output Port  0, joint positions q */
#define OUT_PORT_0_NAME      q
#define OUTPUT_0_WIDTH       8
#define OUTPUT_DIMS_0_COL    1
#define OUTPUT_0_DTYPE       real_T
#define OUTPUT_0_COMPLEX     COMPLEX_NO
#define OUT_0_FRAME_BASED    FRAME_NO
#define OUT_0_BUS_BASED      0
#define OUT_0_BUS_NAME
#define OUT_0_DIMS           1-D
#define OUT_0_ISSIGNED        1
#define OUT_0_WORDLENGTH      8
#define OUT_0_FIXPOINTSCALING 1
#define OUT_0_FRACTIONLENGTH  3
#define OUT_0_BIAS            0
#define OUT_0_SLOPE           0.125


/* Output Port  1, joint velocities dq*/
#define OUT_PORT_1_NAME      dq
#define OUTPUT_1_WIDTH       8
#define OUTPUT_DIMS_1_COL    1
#define OUTPUT_1_DTYPE       real_T
#define OUTPUT_1_COMPLEX     COMPLEX_NO
#define OUT_1_FRAME_BASED    FRAME_NO
#define OUT_1_BUS_BASED      0
#define OUT_1_BUS_NAME
#define OUT_1_DIMS           1-D
#define OUT_1_ISSIGNED        1
#define OUT_1_WORDLENGTH      8
#define OUT_1_FIXPOINTSCALING 1
#define OUT_1_FRACTIONLENGTH  3
#define OUT_1_BIAS            0
#define OUT_1_SLOPE           0.125

/* Output Port  2, joints accelerations ddq */
#define OUT_PORT_2_NAME      ddq
#define OUTPUT_2_WIDTH       8
#define OUTPUT_DIMS_2_COL    1
#define OUTPUT_2_DTYPE       real_T
#define OUTPUT_2_COMPLEX     COMPLEX_NO
#define OUT_2_FRAME_BASED    FRAME_NO
#define OUT_2_BUS_BASED      0
#define OUT_2_BUS_NAME
#define OUT_2_DIMS           1-D
#define OUT_2_ISSIGNED        1
#define OUT_2_WORDLENGTH      8
#define OUT_2_FIXPOINTSCALING 1
#define OUT_2_FRACTIONLENGTH  3
#define OUT_2_BIAS            0
#define OUT_2_SLOPE           0.125

#define NPARAMS              0

#define SAMPLE_TIME_0        INHERITED_SAMPLE_TIME
#define NUM_DISC_STATES      1
#define DISC_STATES_IC       [0]
#define NUM_CONT_STATES      0
#define CONT_STATES_IC       [0]

#define SFUNWIZ_GENERATE_TLC 1
#define SOURCEFILES "__SFB__"
#define PANELINDEX           6
#define USE_SIMSTRUCT        0
#define SHOW_COMPILE_STEPS   0
#define CREATE_DEBUG_MEXFILE 0
#define SAVE_CODE_ONLY       0
#define SFUNWIZ_REVISION     3.0

#define NRWORK               8

#include "simstruc.h"
#include <math.h>

#include <assert.h>
#include <ode/ode.h>
#include <drawstuff/drawstuff.h>
#include <trimesh/TriMesh.h>
#include <time.h>
#include <sys/time.h>
#include <pthread.h>
#include <errno.h>

/* implements an always stable friction mechanism requiring joint speeds to be 0, and setting a maximum torque to reach that speed */
#define FRICTION_HACK

/* if not defined, the peg and hole are not attached to the end effectors */
#define ATTACH_PEG_AND_HOLE

/** 
 * Trough a simple switch system you can use a different set of collision parameters
 * USE_ORIGINAL_COLLISION imposes the use of ad-hoc shosen collision parameters
 * USE_BASKET_COLLISION uses the parameters from the basket demo
 * USE_BUNNY_COLLISION uses the parameters form the bunny demo
 **/
/*#define USE_ORIGINAL_COLLISION*/
#define USE_BASKET_COLLISION
/*#define USE_BUNNY_COLLISION*/


/* 
 * USE_CONTACT_FILTERING enables a simple filtering algorithm that discards contct points
 * which are too near between them
 */
#define USE_CONTACT_FILTERING
/*#define USE_CONTACT_HISTORY*/

/*
 * kinematic and dynamic parameters
 */
#define NUM   8         /* Number of links; */
#define shoulder_height		1
#define shoulder_separation 0.297
#define motor_size			0.055
#define motor_mass			0.260
#define hole_mass			0.001
#define peg_mass			0.001
#define length_biceps		0.16
#define length_forearm		0.115

#define peg_radius			0.014
#define peg_height         0.06
#define peg_offset_z		0.055           /* from ee to clamps */
#define peg_offset_x       -0.0125			/* motor_size/2 - hole_radius */

#define	 hole_offset_y		0.0325-0.02		/* from ee to hole surface */
#define	 hole_offset_x		-0.0075         /* from ee to hole y-center */
#define hole_2_ee   		0.02            /* from ee to hole center */
#define hole_radius			0.04            /* this is the external radius; hole internal radius is .015*/
#define hole_height			0.07

#define time_step			0.0315/3

dsFunctions fn;
dWorldID    world;
dSpaceID	space;

dGeomID		ground;
dGeomID 	hole_geom;
dGeomID		peg_geom;
dJointGroupID contactgroup;

dJointID	peg_to_hand;
dJointID	hole_to_hand;
dJointID	body_to_world;
dJointID	left_arm_to_body;
dJointID	right_arm_to_body;

dBodyID		hole;
dBodyID		body;
dBodyID		peg;
dBodyID     links[NUM];
dJointID    joint[NUM];

static dReal* vertices;
static int numv;
static dTriIndex* indices;
static int numf;
static dTriMeshDataID hole_data;
static TriMesh *hole_trimesh;

static double dim_x[NUM]  = 	{motor_size, 								motor_size,										motor_size,									motor_size,														motor_size, 								motor_size,										motor_size,									motor_size};													/*x size of links; right arm, left arm*/
static double dim_y[NUM]  = 	{motor_size,								motor_size,										motor_size,									motor_size,														motor_size,									motor_size,										motor_size,									motor_size};													/*y size of links; right arm, left arm*/
static double dim_z[NUM]  = 	{motor_size,								length_biceps,									motor_size,									length_forearm,													motor_size,									length_biceps,									motor_size,									length_forearm};												/*z size of links; right arm, left arm*/
static double link_x[NUM] = 	{0.00,										0.00,											0.00,										0.00,															0.00,										0.00,											0.00,										0.00};															/* center of gravity positions */
static double link_y[NUM] = 	{shoulder_separation/2,                     shoulder_separation/2,                          shoulder_separation/2,                      shoulder_separation/2,                                          -(shoulder_separation/2),                   -(shoulder_separation/2),                       -(shoulder_separation/2),                   -(shoulder_separation/2)};
static double link_z[NUM] = 	{shoulder_height,							shoulder_height-length_biceps/2+motor_size/2,	shoulder_height-length_biceps,				shoulder_height-length_biceps-length_forearm/2+motor_size/2,	shoulder_height,							shoulder_height-length_biceps/2+motor_size/2,	shoulder_height-length_biceps,				shoulder_height-length_biceps-length_forearm/2+motor_size/2};

static double x[NUM] = 		{0.00,										0.00,											0.00,										0.00,															0.00,										0.00,											0.00,										0.00};															/* center of gravity positions */
static double y[NUM] = 		{shoulder_separation/2,                     shoulder_separation/2,                          shoulder_separation/2,                      shoulder_separation/2,                                          -(shoulder_separation/2),                   -(shoulder_separation/2),                       -(shoulder_separation/2),                   -(shoulder_separation/2)};
static double z[NUM] = 		{shoulder_height,							shoulder_height-length_biceps+motor_size,		shoulder_height-length_biceps,				shoulder_height-length_biceps-length_forearm+motor_size,		shoulder_height,							shoulder_height-length_biceps+motor_size,		shoulder_height-length_biceps,				shoulder_height-length_biceps-length_forearm+motor_size};

static double m[NUM] = 		{motor_mass,								motor_mass,										motor_mass,									motor_mass,														motor_mass,									motor_mass,										motor_mass,									motor_mass};           										/* mass for every motor */

static double anchor_x[NUM] =	{0.00,										0.00,											0.00,										0.00,															0.00,										0.00,											0.00,										0.00}; /* anchor points for hinges: left arm, right arm*/
static double anchor_y[NUM] =	{shoulder_separation/2,                     shoulder_separation/2,                          shoulder_separation/2,                      shoulder_separation/2,                                          -(shoulder_separation/2),                   -(shoulder_separation/2),                       -(shoulder_separation/2),                   -(shoulder_separation/2)};
static double anchor_z[NUM] =	{shoulder_height,							shoulder_height,								shoulder_height-length_biceps+motor_size/2,	shoulder_height-length_biceps,									shoulder_height,							shoulder_height,								shoulder_height-length_biceps+motor_size/2,	shoulder_height-length_biceps};

static double axis_x[NUM]  = 	{0.00, 1.00, 0.00, 0.00,		0.00, 1.00, 0.00, 0.00};  	/* rotation axis, left arm, right arm */
static double axis_y[NUM]  = 	{1.00, 0.00, 0.00, -1.00,		1.00, 0.00, 0.00, -1.00};
static double axis_z[NUM]  = 	{0.00, 0.00, -1.00, 0.00,		0.00, 0.00, 1.00, 0.00};

pthread_t 		sim_loop_thread;
pthread_mutex_t sync_mutex			= PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t  data_ready_cond 	= PTHREAD_COND_INITIALIZER;
pthread_cond_t  data_consumed_cond	= PTHREAD_COND_INITIALIZER;

static int				consumed			= 0;
static int				ready				= 0;
static int				drawstuff_run		= 1;
static int				simulink_run		= 1;

static real_T			dq_old[NUM]			= {0, 0, 0, 0, 0, 0, 0, 0};

/** caltulates magnitude of a dVector3 vector 
 */
double vectorMagnitude(dVector3 vector) {
	return dSqrt(vector[0]*vector[0] + vector[1]*vector[1] + vector[2]*vector[2]);
}

/** gets called everrytime a collision is detected
 */
void nearCallback(void *data, dGeomID o1, dGeomID o2)
{
  int i,j;
  const int N = 10; /*for the 32 faces hole*/
  dContact contact[N];

  int nAdded = 0;
  int toAdd = 1;
  dContact* contactAdded[N];


  int n =  dCollide(o1,o2,N,&contact[0].geom,sizeof(dContact));

  if (n > 0)  {
	  nAdded = 0;
	  for (i = 0; i < n; i++) {
		  toAdd = 1;
#ifdef USE_CONTACT_FILTERING
		  if(nAdded > 0) {
			  for(j = 0; j < nAdded; j++) {
				  if(dFabs(vectorMagnitude(contact[i].geom.pos) - vectorMagnitude((*contactAdded[j]).geom.pos)) < 0.0005) {
					  toAdd = 0;
					  break;
				  }
			  }
		  }
#ifdef USE_CONTACT_FILTERING
		  if(toAdd == 0) {
			  continue;
		  }
#endif
		  contactAdded[nAdded] = &contact[i];
		  nAdded++;
#endif

#ifdef USE_ORIGINAL_COLLISION
		  contact[i].surface.mode = dContactBounce;
          contact[i].surface.mu   = 0.15; /*plastic on plastic*/
		  contact[i].surface.bounce     = 0.2;
		  contact[i].surface.bounce_vel = 0.0;
#endif
#ifdef USE_BASKET_COLLISION
		  contact[i].surface.slip1 = 0.15;
		  contact[i].surface.slip2 = 0.15;
		  contact[i].surface.mode = dContactSoftERP | dContactSoftCFM | dContactApprox1 | dContactSlip1 | dContactSlip2;
		  contact[i].surface.mu = 50;
		  contact[i].surface.soft_erp = 0.96;
		  contact[i].surface.soft_cfm = 0.04;
#endif
#ifdef USE_BUNNY_COLLISION
		  contact[i].surface.mode = dContactBounce | dContactSoftCFM;
		  contact[i].surface.mu = dInfinity;
		  contact[i].surface.mu2 = 0;
		  contact[i].surface.bounce = 0.1;
		  contact[i].surface.bounce_vel = 0.1;
		  contact[i].surface.soft_cfm = 0.01;
#endif

	      dJointID c = dJointCreateContact (world,contactgroup,&contact[i]);
	      dJointAttach (c,
			    dGeomGetBody(contact[i].geom.g1),
			    dGeomGetBody(contact[i].geom.g2));
    }
  }
}

/** initializes drawstuff visualizer
 */
void start()
{
    float xyz[3] = {  -.75, 0, 1 };
    float hpr[3] = { 0.0, -4.50, 0.00};
    dsSetViewpoint(xyz,hpr);
}

/** stub to insert keybindings during the visualization using drawstuff
 */
void command(int cmd)   /*  key control function */
{
    switch (cmd)
    {
		default:
			break;
    }
}

/** executed during each step of the simulation by the dsSimulationLoop function
 * In this function, one major integration step inside ODE is computed,
 * collisions are computed and drawing instructions are performed. */
void simLoop(int pause)
{
	int i = 0;
	int err;
	double* v0;
	double* v1;
	double* v2;
	int i0;
	int i1;
	int i2;

	timespec abstime;
	timeval now;

	/* simulation thread waits for input data */
    pthread_mutex_lock( &sync_mutex );
    if(!ready) {
    	gettimeofday(&now, NULL);
    	abstime.tv_sec = now.tv_sec+10;
    	abstime.tv_nsec = now.tv_usec*1000;
    	err = pthread_cond_timedwait(&data_ready_cond, &sync_mutex, &abstime);
    	if(err == ETIMEDOUT || !simulink_run) {
    		dsStop();
    		pthread_mutex_unlock( &sync_mutex );
    		return;
    	}
    }
    ready = 0;

    /* collision detection */
    dSpaceCollide(space,0,&nearCallback);
    dWorldStep(world, time_step);
    dJointGroupEmpty(contactgroup);

    /* drawing code */
    for (i = 0; i <NUM; i++ )
    {
    	double link_size[3] = {dim_x[i],dim_y[i],dim_z[i]};
    	double mass_size[3] = {motor_size,motor_size,motor_size};
    	dVector3 box_position;
    	dBodyGetRelPointPos(links[i],link_x[i]-x[i],link_y[i]-y[i],link_z[i]-z[i],box_position); /* dBodyGetPosition refers to Center of Mass; for links that is translated to joint end */

    	dsDrawBoxD(box_position,dBodyGetRotation(links[i]),link_size);
    	dsSetColor(1.0,1.0,0.0);
    	dsSetDrawMode(1);
    	dsDrawBoxD(dBodyGetPosition(links[i]),dBodyGetRotation(links[i]),mass_size);
    	dsSetColor(1.0,1.0,1.0);
    	dsSetDrawMode(0);
    }
    	dsDrawCylinderD(dBodyGetPosition(peg), dBodyGetRotation(peg),peg_height+peg_offset_z,peg_radius);
    if(true)
    {
    	dVector3 small_cylinder_position;
    	dBodyGetRelPointPos(peg,0,0,+peg_height/2,small_cylinder_position);
    	dsSetColor(1.0,1.0,0.0);
    	dsSetDrawMode(1);
    	dsDrawCylinderD(small_cylinder_position, dBodyGetRotation(peg),peg_offset_z,peg_radius);
    	dsSetColor(1.0,1.0,1.0);
    	dsSetDrawMode(0);
    }

    	for (int i=0; i<numf;i++)
    	{
    	  i0 = indices[i*3+0];
    	  i1 = indices[i*3+1];
    	  i2 = indices[i*3+2];
    	  v0 = vertices+i0*3;
    	  v1 = vertices+i1*3;
    	  v2 = vertices+i2*3;
    	  dsDrawTriangleD(dBodyGetPosition(hole), dBodyGetRotation(hole), v0,v1,v2, true);
    	}

    	/* alerting the main thread the simulation and drawing code has been executed */
    	consumed = 1;
    	pthread_cond_signal(&data_consumed_cond);
    	pthread_mutex_unlock( &sync_mutex );
}

/** executes the visualization code in another thread
 */
void *simLoopThread( void *data )
{
	dAllocateODEDataForThread(dAllocateMaskAll);
	dsFunctions* fn = (dsFunctions*)data;
	dsSimulationLoop(0, NULL, 640, 570, fn);
	pthread_exit(NULL);
    return NULL;
}


/** Output functions
 */
void ODEFrank_RR_Outputs_wrapper(real_T *q,
								 real_T *dq,
								 real_T *ddq,
                          const real_T *xC,
                          const real_T *rWork)
{
		int i = 0;
		int err;

		timespec abstime;
		timeval now;
		
		/* reading inputs */
        for(i = 0; i < NUM; i++)
        {
        	dJointAddHingeTorque(joint[i],rWork[i]);
        }

        /* signaling simulation thread that inputs are ready */
		pthread_mutex_lock( &sync_mutex );
		ready = 1;
		pthread_cond_signal(&data_ready_cond);

		/* waiting for simulation thread to execute */
		if(!consumed) {
			gettimeofday(&now, NULL);
			abstime.tv_sec = now.tv_sec+1;
			abstime.tv_nsec = now.tv_usec*1000;
			err = pthread_cond_timedwait(&data_consumed_cond, &sync_mutex, &abstime);
			if(err == ETIMEDOUT) {
				drawstuff_run = 0;
				pthread_mutex_unlock( &sync_mutex );
				return;
			}
		}

		/* writing outputs */
		for (i = 0; i <NUM; i++ )
		{
			q[i] = dJointGetHingeAngle(joint[i]);
			dq[i] = dJointGetHingeAngleRate(joint[i]);
			ddq[i] = (dq[i] - dq_old[i])/time_step;
			dq_old[i] = dq[i];
		}

		consumed = 0;
		pthread_mutex_unlock( &sync_mutex );

}

/**  update function which gets called for every simulation step.
 * 	 Note how, in order to obtain consistent results, simulink simulation time
 * 	 and integration methods must the the same as the one used inside ODE. 
 */
void ODEFrank_RR_mdlUpdate_wrapper(const real_T *tau,
                          const real_T *q,
                          const real_T *dq,
                          const real_T *ddq,
                          real_T *xC,
                          real_T *rWork)
{
	int i = 0;
	xC[0] = 0;
	for(i = 0; i < NUM; ++i)
		rWork[i] = tau[i];
}

/*====================*
 * S-function methods *
 *====================*/

/* Function: mdlInitializeSizes ===============================================
 * Abstract:
 *   Setup sizes of the various vectors.
 */
static void mdlInitializeSizes(SimStruct *S)
{

    DECL_AND_INIT_DIMSINFO(inputDimsInfo);
    DECL_AND_INIT_DIMSINFO(outputDimsInfo);
    ssSetNumSFcnParams(S, NPARAMS);
     if (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S)) {
	 return; /* Parameter mismatch will be reported by Simulink */
     }

    ssSetNumContStates(S, NUM_CONT_STATES);
    ssSetNumDiscStates(S, NUM_DISC_STATES);

    if (!ssSetNumInputPorts(S, NUM_INPUTS)) return;
    inputDimsInfo.width = INPUT_0_WIDTH;
    ssSetInputPortDimensionInfo(S, 0, &inputDimsInfo);
    ssSetInputPortFrameData(S, 0, IN_0_FRAME_BASED);
    ssSetInputPortDataType(S, 0, SS_DOUBLE);
    ssSetInputPortComplexSignal(S, 0, INPUT_0_COMPLEX);
    ssSetInputPortDirectFeedThrough(S, 0, INPUT_0_FEEDTHROUGH);
    ssSetInputPortRequiredContiguous(S, 0, 1); /*direct input signal access*/

    if (!ssSetNumOutputPorts(S, NUM_OUTPUTS)) return;
    outputDimsInfo.width = OUTPUT_0_WIDTH;
    ssSetOutputPortDimensionInfo(S, 0, &outputDimsInfo);
    ssSetOutputPortFrameData(S, 0, OUT_0_FRAME_BASED);
    ssSetOutputPortDataType(S, 0, SS_DOUBLE);
    ssSetOutputPortComplexSignal(S, 0, OUTPUT_0_COMPLEX);

    outputDimsInfo.width = OUTPUT_1_WIDTH;
    ssSetOutputPortDimensionInfo(S, 1, &outputDimsInfo);
    ssSetOutputPortFrameData(S, 1, OUT_1_FRAME_BASED);
    ssSetOutputPortDataType(S, 1, SS_DOUBLE);
    ssSetOutputPortComplexSignal(S, 1, OUTPUT_1_COMPLEX);

    outputDimsInfo.width = OUTPUT_2_WIDTH;
    ssSetOutputPortDimensionInfo(S, 2, &outputDimsInfo);
    ssSetOutputPortFrameData(S, 2, OUT_2_FRAME_BASED);
    ssSetOutputPortDataType(S, 2, SS_DOUBLE);
    ssSetOutputPortComplexSignal(S, 2, OUTPUT_2_COMPLEX);

    ssSetNumSampleTimes(S, 1);
    ssSetNumRWork(S, NUM);
    ssSetNumIWork(S, 0);
    ssSetNumPWork(S, 0);
    ssSetNumModes(S, 0);
    ssSetNumNonsampledZCs(S, 0);

    /* Take care when specifying exception free code - see sfuntmpl_doc.c */
    ssSetOptions(S, (SS_OPTION_EXCEPTION_FREE_CODE |
                     SS_OPTION_USE_TLC_WITH_ACCELERATOR |
		     SS_OPTION_WORKS_WITH_CODE_REUSE));
}

#if defined(MATLAB_MEX_FILE)
#define MDL_SET_INPUT_PORT_DIMENSION_INFO
static void mdlSetInputPortDimensionInfo(SimStruct        *S,
                                         int_T            port,
                                         const DimsInfo_T *dimsInfo)
{
    if(!ssSetInputPortDimensionInfo(S, port, dimsInfo)) return;
}
#endif

#define MDL_SET_OUTPUT_PORT_DIMENSION_INFO
#if defined(MDL_SET_OUTPUT_PORT_DIMENSION_INFO)
static void mdlSetOutputPortDimensionInfo(SimStruct        *S,
                                          int_T            port,
                                          const DimsInfo_T *dimsInfo)
{
 if (!ssSetOutputPortDimensionInfo(S, port, dimsInfo)) return;
}
#endif
# define MDL_SET_INPUT_PORT_FRAME_DATA
static void mdlSetInputPortFrameData(SimStruct  *S,
                                     int_T      port,
                                     Frame_T    frameData)
{
    ssSetInputPortFrameData(S, port, frameData);
}

/* Function: mdlInitializeSampleTimes =========================================
 * Abstract:
 *    Specifiy  the sample time.
 */
static void mdlInitializeSampleTimes(SimStruct *S)
{
    ssSetSampleTime(S, 0, time_step);
    ssSetOffsetTime(S, 0, 0.0);
}

#define MDL_INITIALIZE_CONDITIONS
 /* Function: mdlInitializeConditions ========================================
  * Abstract:
  *    Initialize the states
  */
 static void mdlInitializeConditions(SimStruct *S)
 {
	 int i;
	 
	 real_T *xC   = ssGetDiscStates(S);
	 xC[0] =  0;

	 for(i=0; i < NUM; ++i)
		 dq_old[i]		= 0;
 }

#define MDL_START  /* Change to #undef to remove function */
#if defined(MDL_START)
  /** mdlStart is called once at start of model execution
   * Here ODE and drawstuff are initialized, the simulation thread is started,
   * and physics objects (bodies) and geometric objects (geoms for collision detection)
   * are created. 
   */
  static void mdlStart(SimStruct *S)
  {
	  	int i,j;
	    dMass mass;
	    dMatrix3 rotation;

	    fn.version = DS_VERSION;
	    fn.start   = &start;
	    fn.step    = &simLoop;
	    fn.command = &command;
#ifdef WIN32
	    fn.path_to_textures = "C:\\odewin\\share\\drawstuff\\textures";
#else
        fn.path_to_textures = "/usr/local/share/drawstuff/textures";
#endif

	    dInitODE();

	    consumed = 0;
	    ready = 0;
	    simulink_run = 1;
	    drawstuff_run = 1;
	    for(i=0;i<NUM;i++) {
	    	dq_old[i]=0;
	    }

	    world = dWorldCreate();
	    dWorldSetGravity(world, 0, 0, -9.81);
	    /*
	    dWorldSetERP(world, 0.2);
	    dWorldSetCFM(world, 1e-4);
	    */

	    space = dHashSpaceCreate(0);
	    ground = dCreatePlane(space,0,0,1,0);
	    contactgroup = dJointGroupCreate(0);

	    for (i = 0; i <NUM; i++)
	    {
	        links[i] = dBodyCreate(world);
	        dBodySetPosition(links[i], x[i], y[i], z[i]);
	        dMassSetZero(&mass);
	        dMassSetBoxTotal(&mass,	m[i],	motor_size,	motor_size,	motor_size);
	        dBodySetMass(links[i], 	&mass);
	    }
#ifdef WIN32
        hole_trimesh = TriMesh::read("C:\\odewin\\share\\ode\\hand_hole.obj");
#else
	    hole_trimesh = TriMesh::read("/usr/local/share/ode/hand_hole.obj");
#endif
	    hole_trimesh->need_normals();
	    hole_trimesh->need_faces();
	    numv = hole_trimesh->vertices.size();
	    numf = hole_trimesh->faces.size();
	    vertices = (dReal*)malloc(numv*3*sizeof(dReal));
	    indices = (dTriIndex*)malloc(numf*3*sizeof(dTriIndex));
	    for(i = 0; i<numv; i++) {
	    	vertices[i*3] =	hole_trimesh->vertices[i][0];
	      	vertices[i*3+1] = hole_trimesh->vertices[i][1];
	      	vertices[i*3+2] = hole_trimesh->vertices[i][2];
	    }
	    for(i = 0; i<numf; i++) {
	    	indices[i*3] = 	hole_trimesh->faces[i][0];
	    	indices[i*3+1] =	hole_trimesh->faces[i][1];
	    	indices[i*3+2] =	hole_trimesh->faces[i][2];
	    }

	    hole_data = dGeomTriMeshDataCreate();
	    dGeomTriMeshDataBuildDouble
	    (
	    		hole_data,
	    		vertices,
	    		3 * sizeof(dReal),
	    		numv,
	    		indices,
	    		numf*3,
	    		3 * sizeof(dTriIndex)
	    );


	    hole_geom = dCreateTriMesh(space, hole_data, 0, 0, 0);
	    dGeomSetData(hole_geom,hole_data);
	    hole = dBodyCreate(world);
	    /* dMassSetCylinderTotal directions : 1 == x, 2 == y, 3 == z */

	    dMassSetZero(&mass);
	    dMassSetCylinderTotal(&mass,hole_mass,2,hole_radius,hole_height+hole_offset_y);
	    dBodySetMass(hole, 	&mass);


#ifndef USE_CONTACT_HISTORY
	    dGeomTriMeshEnableTC(hole_geom, dCylinderClass, false);
#endif

	    dGeomSetBody(hole_geom,hole);
	    dBodySetPosition(hole,x[NUM-1]+hole_offset_x,y[NUM-1]+hole_offset_y,z[NUM-1]-hole_2_ee-motor_size/2);
	    /*dRFromEulerAngles(rotation,0,M_PI_2,M_PI_2);*/
	    /*dBodySetRotation(hole,rotation);*/
#ifdef  ATTACH_PEG_AND_HOLE
	    hole_to_hand = dJointCreateFixed(world, 0);
	    dJointAttach(hole_to_hand, hole, links[NUM-1]);
	    dJointSetFixed(hole_to_hand);
#endif

	    body = dBodyCreate(world);
	    body_to_world = dJointCreateFixed(world, 0);
	    dJointAttach(body_to_world, body, 0);
	    dJointSetFixed(body_to_world);

	    for (j = 0; j <NUM; j++)
	    {
	        joint[j] = dJointCreateHinge(world, 0);
	        if(j != 0 && j != 4)
	        	dJointAttach(joint[j], links[j-1], links[j]);
	        else
	        	dJointAttach(joint[j], body, links[j]);
	        dJointSetHingeAnchor(joint[j], anchor_x[j], anchor_y[j],anchor_z[j]);
	        dJointSetHingeAxis(joint[j], axis_x[j], axis_y[j], axis_z[j]);
	        dJointSetHingeParam (joint[j], dParamLoStop, -3.14*10/9);
	        dJointSetHingeParam (joint[j], dParamHiStop, 3.14*10/9);
#ifdef FRICTION_HACK
	        dJointSetHingeParam (joint[j], dParamVel, 0);
	        dJointSetHingeParam (joint[j], dParamFMax, 0.2);
#endif
	    }

	    peg_geom = dCreateCylinder(space, peg_radius, peg_height+peg_offset_z);
	    peg = dBodyCreate(world);
	    /* dMassSetCylinderTotal directions : 1 == x, 2 == y, 3 == z */

	    dMassSetZero(&mass);
	    dMassSetCylinderTotal(&mass,peg_mass,3,peg_radius,peg_height+peg_offset_z);
	    dBodySetMass(peg, 	&mass);

	    dGeomSetBody(peg_geom,peg);
	    dBodySetPosition(peg, x[NUM/2-1]+peg_offset_x, y[NUM/2-1], z[NUM/2-1]-(peg_height+peg_offset_z)/2-motor_size/2);
#ifdef	ATTACH_PEG_AND_HOLE
	    peg_to_hand = dJointCreateFixed(world, 0);
	    dJointAttach(peg_to_hand, peg, links[NUM/2-1]);
	    dJointSetFixed(peg_to_hand);
#endif

	    pthread_create( &sim_loop_thread, NULL, &simLoopThread, (void*) &fn);
  }
#endif /*  MDL_START */

#define MDL_SET_INPUT_PORT_DATA_TYPE
static void mdlSetInputPortDataType(SimStruct *S, int port, DTypeId dType)
{
    ssSetInputPortDataType( S, 0, dType);
}
#define MDL_SET_OUTPUT_PORT_DATA_TYPE
static void mdlSetOutputPortDataType(SimStruct *S, int port, DTypeId dType)
{
    ssSetOutputPortDataType(S, 0, dType);
    ssSetOutputPortDataType(S, 1, dType);
    ssSetOutputPortDataType(S, 2, dType);
}

#define MDL_SET_DEFAULT_PORT_DATA_TYPES
static void mdlSetDefaultPortDataTypes(SimStruct *S)
{
  ssSetInputPortDataType( S, 0, SS_DOUBLE);
  ssSetOutputPortDataType(S, 0, SS_DOUBLE);
  ssSetOutputPortDataType(S, 1, SS_DOUBLE);
  ssSetOutputPortDataType(S, 2, SS_DOUBLE);
}




/* Function: mdlOutputs =======================================================
 *
*/
static void mdlOutputs(SimStruct *S, int_T tid)
{
    real_T        *q  = (real_T *)ssGetOutputPortRealSignal(S,0);
    real_T        *dq  = (real_T *)ssGetOutputPortRealSignal(S,1);
    real_T		  *ddq = (real_T *)ssGetOutputPortRealSignal(S,2);
    const real_T   *xC = ssGetDiscStates(S);
    const real_T   *rWork = ssGetRWork(S);

    ODEFrank_RR_Outputs_wrapper(q, dq, ddq, xC, rWork);
    if(!drawstuff_run) ssSetErrorStatus(S,"error encountered while processing Output. Did you close the simulation window?");
}

#define MDL_UPDATE	// Change to #undef to remove function
#if defined(MDL_UPDATE)
static void mdlUpdate( SimStruct *S, int_T tid )
{
    const real_T   *tau  = (const real_T*) ssGetInputPortSignal(S,0);
    real_T         *xC  = ssGetDiscStates(S);
    real_T        *q  = (real_T *) ssGetOutputPortRealSignal(S,0);
    real_T        *dq  = (real_T *) ssGetOutputPortRealSignal(S,1);
    real_T        *ddq  = (real_T *) ssGetOutputPortRealSignal(S,2);
    real_T        *rWork = (real_T *) ssGetRWork(S);

    ODEFrank_RR_mdlUpdate_wrapper(tau,q,dq,ddq,xC,rWork);
}
#endif /* MDL_UPDATE */

/* Function: mdlTerminate =====================================================
 * Abstract:
 *    In this function, you should perform any actions that are necessary
 *    at the termination of a simulation.  For example, if memory was
 *    allocated in mdlStart, this is the place to free it.
 */
static void mdlTerminate(SimStruct *S)
{
	simulink_run = 0;
	pthread_cond_signal(&data_ready_cond);

	pthread_join(sim_loop_thread,NULL);
	dJointGroupDestroy(contactgroup);
	dSpaceDestroy(space);
	dWorldDestroy(world);
	dCloseODE();
    free(vertices);
    free(indices);
}
#ifdef  MATLAB_MEX_FILE    /* Is this file being compiled as a MEX-file? */
#include "simulink.c"      /* MEX-file interface mechanism */
#else
#include "cg_sfun.h"       /* Code generation registration function */
#endif


