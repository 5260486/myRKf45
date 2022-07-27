#include "mws_ivp_solver.h"
#include "mws_ls_solver.h"
#include "mws_nls_solver.h"

#include "my_euler.c"    /*自定义算法头文件*/

void MwsRegisterUserAlgorithm1(void* mdl_data)
{
	/* Register user defined LS and NLS algorithm. 线性或非线性算法*/
}

void MwsRegisterUserAlgorithm2(void* sim_data)
{
	/* Register user defined IVP algorithm.积分算法 */
    MwsIVPSolverProp ivp_prop;
    MwsIVPSolverFcns ivp_fcns;

    ivp_prop.m_name = "myeuler";                        /*名称，不区分大小写，保证唯一*/
    ivp_prop.m_desc = "MYEULER";                        /*算法描述*/
    ivp_prop.m_fixedStep = moTrue;                      /*是否为定步长，moTrue或moFalse*/
    ivp_prop.m_ivpType = MWS_IVP_ODE;                   /*积分算法类型，ODE为MWS_IVP_ODE，DAE为MWS_IVP_DAE*/
    ivp_fcns.m_createPtr = &myEulerCreate;              /*创建算法对象函数的指针*/
    ivp_fcns.m_createPBPtr = &myEulerProblemCreate;     /*创建问题对象函数的指针*/
    ivp_fcns.m_destroyPBPtr = &myEulerProblemDestroy;   /*销毁问题对象函数的指针*/
    ivp_fcns.m_destroyPtr = &myEulerDestroy;            /*销毁算法对象函数的指针*/
    ivp_fcns.m_initPtr = &myEulerInit;                  /*初始化函数指针*/
    ivp_fcns.m_interpolatePtr = &myEulerInterpolate;    /*插值函数函数*/
    ivp_fcns.m_solvePtr = &myEulerSolve;                /*求解函数指针*/
    isimRegisterIVPSolver(sim_data, &ivp_prop, &ivp_fcns);
}

void MwsUnregisterUserAlgorithm1(void* mdl_data)
{
	/* Unregister user defined LS and NLS algorithm. */
}

void MwsUnregisterUserAlgorithm2(void* sim_data)
{
	/* Unregister user defined IVP algorithm. */
    isimUnregisterIVPSolver(sim_data, "myeuler");
}