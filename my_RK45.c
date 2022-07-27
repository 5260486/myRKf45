/***************************************************************************
///
/// Copyright (c) 2020, 苏州同元软控信息技术有限公司
/// All rights reserved.
///
/// @file           my_RK45.c
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
#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

    /* 算法对象 */
    typedef struct
    {
        MwsIVPUtilFcns	m_utils;
        void* m_userData;

    } MyRK45;

    /* 求解问题的数据，用数据内存分配和释放函数 */
    typedef struct
    {
        MoReal* m_preY;
        MoReal* m_curY;

        MoReal* m_preYp;
        MoReal* m_curYp;

        MoReal* k1;
        MoReal* k2;
        MoReal* k3;
        MoReal* k4;
        MoReal* k5;
        MoReal* k6;

        MoReal* k2y;
        MoReal* k3y;
        MoReal* k4y;
        MoReal* k5y;
        MoReal* k6y;

        MoReal m_curTime;
        MoReal m_initialStep;
        MoReal m_h;
        
        MoReal* m_D;
        MoReal m_Q;
    } MyRK45ProblemData;

    /* 求解问题对象 */
    typedef struct
    {
        MoSize          m_nStates;

        MwsIVPOptions	m_opt;
        MwsIVPCallback  m_callback;
        void* m_userData;

        MyRK45ProblemData* m_data;
        MyRK45* m_solverWork;

    } MyRK45Problem;

    void myRK45ProblemDestroy(MwsIVPSolverObj solver, MwsIVPObj ivp);
    void myRK45Destroy(MwsIVPSolverObj solver);

    /// <summary>
    /// 创建算法
    /// </summary>
    /// <param name="util_fcns">工具函数（由求解器提供，用户调用）</param>
    /// <param name="user_data">用户数据（求解器内部数据，传递给工具函数util_fcns，自定义算法无需关心其内容）</param>
    /// <returns></returns>
    MwsIVPSolverObj myRK45Create(MwsIVPUtilFcns* util_fcns, void* user_data)
    {
        MyRK45* sw = (MyRK45*)util_fcns->m_allocMemory(user_data, 1, sizeof(MyRK45));

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
    MwsIVPObj myRK45ProblemCreate(MwsIVPSolverObj solver, MwsSize n, MwsIVPCallback* call_back, MwsIVPOptions* opt, void* ivp_user_data)
    {
        MyRK45* sw = (MyRK45*)solver;
        MyRK45Problem* spw = (MyRK45Problem*)sw->m_utils.m_allocMemory(sw->m_userData, 1, sizeof(MyRK45Problem));

        if (spw)
        {
            MyRK45ProblemData* ds = (MyRK45ProblemData*)sw->m_utils.m_allocDataMemory(sw->m_userData, 1, sizeof(MyRK45ProblemData));

            if (ds == mwsNullPtr)
            {
                sw->m_utils.m_freeMemory(sw->m_userData, spw);
                return mwsNullPtr;
            }

            memset(spw, 0, sizeof(*spw));
            memset(ds, 0, sizeof(*ds));

            spw->m_callback = *call_back;
            spw->m_userData = ivp_user_data;
            spw->m_opt = *opt;
            spw->m_nStates = n;
            spw->m_data = ds;
            spw->m_solverWork = sw;

            spw->m_data->m_curTime = 0;
            spw->m_data->m_initialStep = spw->m_opt.m_stopTime/1000;          //默认初始积分步长default
            spw->m_data->m_Q = 0;
            spw->m_data->m_h = 0;

            spw->m_opt.m_maxStepSize = spw->m_opt.m_stopTime;
            
            if (spw->m_nStates > 0)
            {
                spw->m_data->m_preY = (MoReal*)sw->m_utils.m_allocDataMemory(sw->m_userData, 1, n * sizeof(MoReal));  //(MoReal *)强制转换
                spw->m_data->m_curY = (MoReal*)sw->m_utils.m_allocDataMemory(sw->m_userData, 1, n * sizeof(MoReal));

                spw->m_data->m_preYp = (MoReal*)sw->m_utils.m_allocDataMemory(sw->m_userData, 1, n * sizeof(MoReal));
                spw->m_data->m_curYp = (MoReal*)sw->m_utils.m_allocDataMemory(sw->m_userData, 1, n * sizeof(MoReal));

                spw->m_data->k1 = (MoReal*)sw->m_utils.m_allocDataMemory(sw->m_userData, 1, n * sizeof(MoReal));
                spw->m_data->k2 = (MoReal*)sw->m_utils.m_allocDataMemory(sw->m_userData, 1, n * sizeof(MoReal));
                spw->m_data->k3 = (MoReal*)sw->m_utils.m_allocDataMemory(sw->m_userData, 1, n * sizeof(MoReal));
                spw->m_data->k4 = (MoReal*)sw->m_utils.m_allocDataMemory(sw->m_userData, 1, n * sizeof(MoReal));
                spw->m_data->k5 = (MoReal*)sw->m_utils.m_allocDataMemory(sw->m_userData, 1, n * sizeof(MoReal));
                spw->m_data->k6 = (MoReal*)sw->m_utils.m_allocDataMemory(sw->m_userData, 1, n * sizeof(MoReal));

                spw->m_data->k2y = (MoReal*)sw->m_utils.m_allocDataMemory(sw->m_userData, 1, n * sizeof(MoReal));
                spw->m_data->k3y = (MoReal*)sw->m_utils.m_allocDataMemory(sw->m_userData, 1, n * sizeof(MoReal));
                spw->m_data->k4y = (MoReal*)sw->m_utils.m_allocDataMemory(sw->m_userData, 1, n * sizeof(MoReal));
                spw->m_data->k5y = (MoReal*)sw->m_utils.m_allocDataMemory(sw->m_userData, 1, n * sizeof(MoReal));
                spw->m_data->k6y = (MoReal*)sw->m_utils.m_allocDataMemory(sw->m_userData, 1, n * sizeof(MoReal));
                spw->m_data->m_D = (MoReal*)sw->m_utils.m_allocDataMemory(sw->m_userData, 1, n * sizeof(MoReal));
                if (!spw->m_data->m_preY || !spw->m_data->m_curY)
                {
                    myRK45ProblemDestroy(sw, spw);
                    spw = MWnullptr;
                }
                else
                {
                    MoSize index;
                    for (index = 0; index < n; ++index)
                    {
                        spw->m_data->m_preY[index] = 0;
                        spw->m_data->m_curY[index] = 0;

                        spw->m_data->m_preYp[index] = 0;
                        spw->m_data->m_curYp[index] = 0;

                        spw->m_data->m_D[index] = 0;

                        spw->m_data->k1[index] = 0;
                        spw->m_data->k2[index] = 0;
                        spw->m_data->k3[index] = 0;
                        spw->m_data->k4[index] = 0;
                        spw->m_data->k5[index] = 0;
                        spw->m_data->k6[index] = 0;

                        spw->m_data->k2y[index] = 0;
                        spw->m_data->k3y[index] = 0;
                        spw->m_data->k4y[index] = 0;
                        spw->m_data->k5y[index] = 0;
                        spw->m_data->k5y[index] = 0;

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
    MwsInteger myRK45Init(MwsIVPSolverObj solver, MwsIVPObj ivp, MwsReal t0, const MwsReal* y0,
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
    MwsInteger myRK45Solve(MwsIVPSolverObj solver, MwsIVPObj ivp, MwsReal step_size, MwsReal t,
        MwsReal tout, MwsReal* tret, MwsReal* yret, MwsReal* ypret, void* reserve)
    {
        MyRK45* sw = (MyRK45*)solver;
        MyRK45Problem* spw = (MyRK45Problem*)ivp;

        MoSize index = 0;
        MoSize nState = spw->m_nStates;

        spw->m_data->m_curTime = t;
        spw->m_data->m_initialStep = step_size;             //初始化积分步长

        MoReal* D = spw->m_data->m_D;                       //D=w(i+1)-y(i+1)
        MoReal Q = spw->m_data->m_Q;
        MoReal h = spw->m_data->m_initialStep;
        
        MoReal* preY = spw->m_data->m_preY;                 //上个y
        MoReal* curY = spw->m_data->m_curY;                //当前y

        MoReal* preYp = spw->m_data->m_preYp;               //上个y'
        MoReal* curYp = spw->m_data->m_curYp;               //当前y’

        MoReal* K1 = spw->m_data->k1;
        MoReal* K2 = spw->m_data->k2;
        MoReal* K3 = spw->m_data->k3;
        MoReal* K4 = spw->m_data->k4;
        MoReal* K5 = spw->m_data->k5;
        MoReal* K6 = spw->m_data->k6;

        MoReal* K2y = spw->m_data->k2y;
        MoReal* K3y = spw->m_data->k3y;
        MoReal* K4y = spw->m_data->k4y;
        MoReal* K5y = spw->m_data->k5y;
        MoReal* K6y = spw->m_data->k6y;
        
        MoBoolean Flag = moTrue;
        while (Flag)
        {
            for (index = 0; index < nState; ++index)
            {
                preY[index] = yret[index];
                preYp[index] = ypret[index];

                //公式主体

                K1[index] = preYp[index];
                K2y[index] = preY[index] + 1 / 4 * K1[index];
                spw->m_callback.m_rshFunction(spw->m_userData, spw->m_data->m_curTime + h / 4,
                    K2y, K2);
                K3y[index] = preY[index] + 3 / 32 * K1[index] + 9 / 32 * K2[index];
                spw->m_callback.m_rshFunction(spw->m_userData, spw->m_data->m_curTime + h * 3 / 8,
                    K3y, K3);
                K4y[index] = preY[index] + 1932 / 2197 * K1[index] - 7200 / 2197 * K2[index] + 7296 / 2197 * K3[index];
                spw->m_callback.m_rshFunction(spw->m_userData, spw->m_data->m_curTime + h * 12 / 13,
                    K4y, K4);
                K5y[index] = preY[index] + 439 / 216 * K1[index] - 8 * K2[index] + 3680 / 513 * K3[index]
                    - 845 / 4104 * K4[index];
                spw->m_callback.m_rshFunction(spw->m_userData, spw->m_data->m_curTime + h,
                    K5y, K5);
                K6y[index] = preY[index] - 8 / 27 * K1[index] + 2 * K2[index] - 3544 / 2565 * K3[index]
                    + 1859 / 4104 * K4[index] - 11 / 40 * K5[index];
                spw->m_callback.m_rshFunction(spw->m_userData, spw->m_data->m_curTime + h / 2,
                    K6y, K6);

                D[index] = K1[index] / 360 - 128 / 4275 * K3[index] - 2197 / 75240 * K4[index]
                    + K5[index] / 50 + 2 / 55 * K6[index];

                ypret[index] = D[index];

            }
            //MoReal sum = 0;
            //MoReal normD = 0;

            //for (index = 0; index < nState; ++index)
            //{
            //    sum = sum + fabs(D[index]) *fabs( D[index]);
            //}   
            //normD = sqrt(sum);
            //MoReal R = normD / h; 
            MoReal smallestD = D[0];
            for (index = 0; index < nState; ++index)
            {
                if (fabs(D[index])<fabs(smallestD))
                {
                    smallestD =fabs(D[index]);
                }
            }

            MoReal R = smallestD / h;

            if (R < spw->m_opt.m_absoluteTolerance[0])            //误差容限满足精度
            {
                *tret += h;
                for (index = 0; index < nState; ++index)
                {
                    //improved RK45
                    yret[index] = preY[index] + 25 / 216 * K1[index] + 1408 / 2565 * K3[index]
                        + 2197 / 4104 * K4[index] - 1 / 5 * K5[index] + D[index];
                
                    curY[index] = yret[index];

                    spw->m_callback.m_rshFunction(spw->m_userData, spw->m_data->m_curTime + spw->m_data->m_initialStep, curY, preYp);

                    ypret[index] = preYp[index];
                    curYp[index] = ypret[index];

                }
                
                Flag = moFalse;
            }
            
                Q = 0.84*pow(spw->m_opt.m_absoluteTolerance[0] / R, 1 / 4);
                if (Q <= 0.1)
                {
                    h = 0.2 * h;
                }
                else if (Q >= 4)
                {
                    h = 4 * h;
                }
                else
                {
                    h = h * Q;
                }

                if (h > spw->m_opt.m_maxStepSize)      //步长大于设定的最大步长，则h=hmax
                {
                    h = spw->m_opt.m_maxStepSize;
                }

                if (spw->m_data->m_curTime >= spw->m_opt.m_stopTime)       //当前仿真时间大于设定的停止仿真时间
                {
                    *tret = tout;
                    Flag = moFalse;

                   // return MWS_IVP_TSTOP_RETURN;        //到达终止时间而返回
                }
                else if (spw->m_data->m_curTime + h > spw->m_opt.m_stopTime)     //ti+qh>L,即超出边界，则h=L-ti
                {
                    h = spw->m_opt.m_stopTime - spw->m_data->m_curTime;
                }
                if (h < 0)
                {
                    Flag = moFalse;

                    return MWS_IVP_FAIL;
                }
            

        }

            
                /* 更新当前积分时间 */
      // *tret += h;                //积分时间增加一步。外部进行循环                 
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
    MwsInteger myRK45Interpolate(MwsIVPSolverObj solver, MwsIVPObj ivp, MwsReal tout, MwsReal* yret, void* reserve)
    {
        /* 线性插值 */
        MyRK45* sw = (MyRK45*)solver;
        MyRK45Problem* spw = (MyRK45Problem*)ivp;

        MoSize nStates = spw->m_nStates;
        MoSize index;

        for (index = 0; index < nStates; ++index)
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
    void myRK45ProblemDestroy(MwsIVPSolverObj solver, MwsIVPObj ivp)
    {
        MyRK45* sw = (MyRK45*)solver;
        MyRK45Problem* spw = (MyRK45Problem*)ivp;

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
    void myRK45Destroy(MwsIVPSolverObj solver)
    {
        MyRK45* sw = (MyRK45*)solver;

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

