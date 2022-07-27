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
    MoReal *m_preYp;

    MoReal m_curTime;
    MoReal m_initialStep;
} MyEulerProblemData;

/* 求解问题对象 */
typedef struct  
{
    MoSize          m_nStates;

    MwsIVPOptions	m_opt;
    MwsIVPCallback  m_callback;
    void*			m_userData;

    MyEulerProblemData* m_data;
    MyEuler* m_solverWork;

} MyEulerProblem;

void myEulerProblemDestroy(MwsIVPSolverObj solver, MwsIVPObj ivp);
void myEulerDestroy(MwsIVPSolverObj solver);

/// <summary>
/// 创建算法
/// </summary>
/// <param name="util_fcns">工具函数（由求解器提供，用户调用）</param>
/// <param name="user_data">用户数据（求解器内部数据，传递给工具函数util_fcns，自定义算法无需关心其内容）</param>
/// <returns></returns>
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

/// <summary>
/// 创建问题
/// </summary>
/// <param name="solver">积分算法对象</param>
/// <param name="n">问题规模，即状态变量数量==微分方程阶数</param>
/// <param name="call_back">回调函数</param>
/// <param name="opt">积分器选项</param>
/// <param name="ivp_user_data">用户数据(求解器内部数据，传递给回调函数call_back，算法不需关心)</param>
/// <returns></returns>
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
            spw->m_data->m_preY = (MoReal *)sw->m_utils.m_allocDataMemory(sw->m_userData, 1, n*sizeof(MoReal));  //(MoReal *)强制转换
            spw->m_data->m_curY = (MoReal *)sw->m_utils.m_allocDataMemory(sw->m_userData, 1, n*sizeof(MoReal));
            spw->m_data->m_preYp = (MoReal*)sw->m_utils.m_allocDataMemory(sw->m_userData, 1, n * sizeof(MoReal));

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
                    spw->m_data->m_preYp[index] = 0;
                }
            }
        }
    }

    return spw;   //返回：问题对象(自定义算法内部数据，求解器不关心其内容)，作为其它接口函数的第二个参数，即下边的ivp
}

/// <summary>
/// 初始化
/// </summary>
/// <param name="solver">积分算法对象</param>
/// <param name="ivp">问题对象</param>
/// <param name="t0">起始时间</param>
/// <param name="y0">y的初始值（必定不为空）</param>
/// <param name="yp0">y’的初始值（DAE）</param>
/// <param name="is_reinit">是否重新初始化（重启）</param>
/// <param name="reserve">保留参数，暂不使用</param>
/// <returns></returns>
MwsInteger myEulerInit(MwsIVPSolverObj solver, MwsIVPObj ivp, MwsReal t0, const MwsReal* y0, 
    const MwsReal* yp0, MwsBoolean is_reinit, void* reserve)
{
    /* nothing to do */
    return MWS_IVP_SUCCESS;  //返回状态，取MwsIVPStatus的值
}

/// <summary>
/// 求解
/// </summary>
/// 输入：
/// <param name="solver">积分算法对象</param>
/// <param name="ivp">问题对象</param>
/// <param name="step_size">步长（定步长，或初始化积分步长）</param>
/// <param name="t">当前时间</param>
/// <param name="tout">期望输出时间</param>
/// 输出：
/// <param name="tret">求解实际达到的时间</param>
/// <param name="yret">y的结果值</param>
/// <param name="ypret">y’的结果值（DAE）</param>
/// <param name="reserve">保留参数，暂不使用</param>
/// <returns></returns>
MwsInteger myEulerSolve(MwsIVPSolverObj solver, MwsIVPObj ivp, MwsReal step_size, MwsReal t, 
    MwsReal tout, MwsReal* tret, MwsReal* yret, MwsReal* ypret, void* reserve)
{
    MyEuler* sw = (MyEuler*)solver;
    MyEulerProblem* spw = (MyEulerProblem*)ivp;

    MoSize index = 0;
    MoSize nState = spw->m_nStates;

    spw->m_data->m_curTime = t;
    spw->m_data->m_initialStep = step_size;

    MoReal* preY = spw->m_data->m_preY;                 //上个y
    MoReal* curY = spw->m_data->m_curY;                 //当前y
    MoReal* preYp = spw->m_data->m_preYp;               //上个y'
   
    for (index = 0; index<nState; ++index)              
    {
        preY[index] = yret[index];
        preYp[index] = ypret[index];

        MoReal* k1 = preYp[index];
        MoReal* k2 = spw->m_callback.m_rshFunction(spw->m_userData, spw->m_data->m_curTime + spw->m_data->m_initialStep, 
                                                    yret[index] + spw->m_data->m_initialStep * k1, ypret[index]);
                                                                //k2=f(xn+h,yn+h*k1)
        yret[index] += spw->m_data->m_initialStep / 2 * (k1 + k2);      //梯形公式
        
        curY[index] = yret[index];
    }
    

    /* 更新当前积分时间 */
    *tret += spw->m_data->m_initialStep;                //积分时间增加一步，说明此处是每一步内的处理，外部有对算法的循环引用

    return MWS_IVP_SUCCESS;  //返回状态，取MwsIVPStatus的值
}

/// <summary>
/// 生成插值函数（过yn点的曲线）
/// </summary>
/// 输入：
/// <param name="solver">积分算法对象</param>
/// <param name="ivp">问题对象</param>
/// <param name="tout">期望的输出时间</param>
/// 输出：
/// <param name="yret">y的结果值</param>
/// <param name="reserve"></param>
/// <returns></returns>
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

/// <summary>
/// 销毁问题
/// </summary>
/// <param name="solver"></param>
/// <param name="ivp"></param>
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

/// <summary>
/// 销毁积分算法
/// </summary>
/// <param name="solver"></param>
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

