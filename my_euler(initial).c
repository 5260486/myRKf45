/***************************************************************************
///
/// Copyright (c) 2020, 苏州同元软控信息技术有限公司
/// All rights reserved.
///
/// @file           my_euler.c
/// @brief          显式欧拉积分算法示例
///
/// @version        v1.0
/// @author         田显钊
/// @date           2020/07/07
///
***************************************************************************/

#include "mo_types.h"
#include "mws_ivp_solver.h"

#include <memory.h>

#ifdef __cplusplus
extern "C"{
#endif

/* 算法对象 */
typedef struct  
{
    MwsIVPUtilFcns	m_utils;
    void*           m_userData;

} MyEuler;

/* 求解问题的数据，用数据内存分配和释放函数 */
typedef struct  
{
    MoReal *m_preY;
    MoReal *m_curY;  

    MoReal m_curTime;
    MoReal m_initialStep;
} MyEulerProblemData;

/* 求解问题对象 */
typedef struct  
{
    MoSize          m_nStates;

    MwsIVPOptions	m_opt;
    MwsIVPCallback	m_callback;
    void*			m_userData;

    MyEulerProblemData* m_data;
    MyEuler* m_solverWork;

} MyEulerProblem;

void myEulerProblemDestroy(MwsIVPSolverObj solver, MwsIVPObj ivp);
void myEulerDestroy(MwsIVPSolverObj solver);

MwsIVPSolverObj myEulerCreate(MwsIVPUtilFcns* util_fcns, void* user_data)
{
    MyEuler* sw = (MyEuler*)util_fcns->m_allocMemory(user_data, 1, sizeof(MyEuler));

    if (sw)
    {
        memset(sw, 0, sizeof(*sw));
        sw->m_utils = *util_fcns;
        sw->m_userData = user_data;
    }

    return sw;
}


MwsIVPObj myEulerProblemCreate(MwsIVPSolverObj solver, MwsSize n, MwsIVPCallback* call_back, MwsIVPOptions* opt, void* ivp_user_data)
{
    MyEuler* sw = (MyEuler*)solver;
    MyEulerProblem* spw = (MyEulerProblem*)sw->m_utils.m_allocMemory(sw->m_userData, 1, sizeof(MyEulerProblem));

    if (spw)
    {
        MyEulerProblemData* ds = (MyEulerProblemData*)sw->m_utils.m_allocDataMemory(sw->m_userData, 1, sizeof(MyEulerProblemData));

        if (ds == mwsNullPtr)
        {
            sw->m_utils.m_freeMemory(sw->m_userData, spw);
            return mwsNullPtr;
        }

        memset(spw, 0, sizeof(*spw));
        memset(ds, 0, sizeof(*ds));

        spw->m_callback = *call_back;
        spw->m_userData = ivp_user_data;
        spw->m_nStates = n;
        spw->m_data = ds;
        spw->m_solverWork = sw;

        spw->m_data->m_curTime = 0;
        spw->m_data->m_initialStep = 0.002;

        if (spw->m_nStates > 0)
        {
            spw->m_data->m_preY = (MoReal *)sw->m_utils.m_allocDataMemory(sw->m_userData, 1, n*sizeof(MoReal));
            spw->m_data->m_curY = (MoReal *)sw->m_utils.m_allocDataMemory(sw->m_userData, 1, n*sizeof(MoReal));

            if (!spw->m_data->m_preY || !spw->m_data->m_curY)
            {
                myEulerProblemDestroy(sw, spw);
                spw = MWnullptr;
            }
            else
            {
                MoSize index;
                for (index=0; index<n; ++index)
                {
                    spw->m_data->m_preY[index] = 0;
                    spw->m_data->m_curY[index] = 0;
                }
            }
        }
    }

    return spw;
}

MwsInteger myEulerInit(MwsIVPSolverObj solver, MwsIVPObj ivp, MwsReal t0, const MwsReal* y0, 
    const MwsReal* yp0, MwsBoolean is_reinit, void* reserve)
{
    /* nothing to do */
    return MWS_IVP_SUCCESS;
}

MwsInteger myEulerSolve(MwsIVPSolverObj solver, MwsIVPObj ivp, MwsReal step_size, MwsReal t, 
    MwsReal tout, MwsReal* tret, MwsReal* yret, MwsReal* ypret, void* reserve)
{
    MyEuler* sw = (MyEuler*)solver;
    MyEulerProblem* spw = (MyEulerProblem*)ivp;

    MoSize index = 0;
    MoSize nState = spw->m_nStates;

    MoReal* preY = spw->m_data->m_preY;
    MoReal* curY = spw->m_data->m_curY;

    spw->m_data->m_curTime = t;
    spw->m_data->m_initialStep = step_size;

    for (index = 0; index<nState; ++index)
    {
        preY[index] = yret[index];
        yret[index] += spw->m_data->m_initialStep * ypret[index];
        curY[index] = yret[index];
    }

    /* 更新当前积分时间 */
    *tret += spw->m_data->m_initialStep;

    return MWS_IVP_SUCCESS;
}


MwsInteger myEulerInterpolate(MwsIVPSolverObj solver, MwsIVPObj ivp, MwsReal tout, MwsReal* yret, void* reserve)
{
    /* 线性插值 */
    MyEuler* sw = (MyEuler*)solver;
    MyEulerProblem* spw = (MyEulerProblem*)ivp;

    MoSize nStates = spw->m_nStates;
    MoSize index;

    for (index=0; index<nStates; ++index)
    {
        /* 线性插值: 
                    f(x2)-f(x1)
            L(x) = ------------ * (x - x1) + f(x1)
                       x2-x1 
        */
        yret[index] = (spw->m_data->m_curY[index] - spw->m_data->m_preY[index]) * (
            tout - spw->m_data->m_curTime + spw->m_data->m_initialStep) / spw->m_data->m_initialStep + spw->m_data->m_preY[index]; 
    }

    return MWS_IVP_SUCCESS;
}

void myEulerProblemDestroy(MwsIVPSolverObj solver, MwsIVPObj ivp)
{
    MyEuler* sw = (MyEuler*)solver;
    MyEulerProblem* spw = (MyEulerProblem*)ivp;

    if (spw)
    {
        if (spw->m_data->m_preY)
        {
            (*sw->m_utils.m_freeDataMemory)(sw->m_userData, spw->m_data->m_preY);
        }

        if (spw->m_data->m_curY)
        {
            (*sw->m_utils.m_freeDataMemory)(sw->m_userData, spw->m_data->m_curY);
        }

        if (spw->m_data)
        {
            (*sw->m_utils.m_freeDataMemory)(sw->m_userData, spw->m_data);
        }

        (*sw->m_utils.m_freeMemory)(sw->m_userData, spw);
    }
}

void myEulerDestroy(MwsIVPSolverObj solver)
{
    MyEuler* sw = (MyEuler*)solver;

    if (sw)
    {
        (*sw->m_utils.m_freeMemory)(sw->m_userData, sw);
    }
}

#ifdef __cplusplus
}
#endif


/***************************************************************************
//   end of file
***************************************************************************/

