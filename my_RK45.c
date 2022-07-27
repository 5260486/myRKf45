/***************************************************************************
///
/// Copyright (c) 2020, ����ͬԪ�����Ϣ�������޹�˾
/// All rights reserved.
///
/// @file           my_RK45.c
/// @brief          ��ʽŷ�������㷨ʾ��
///
/// @version        v1.0
/// @author         ������
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

    /* �㷨���� */
    typedef struct
    {
        MwsIVPUtilFcns	m_utils;
        void* m_userData;

    } MyRK45;

    /* �����������ݣ��������ڴ������ͷź��� */
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

    /* ���������� */
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
    /// �����㷨
    /// </summary>
    /// <param name="util_fcns">���ߺ�������������ṩ���û����ã�</param>
    /// <param name="user_data">�û����ݣ�������ڲ����ݣ����ݸ����ߺ���util_fcns���Զ����㷨������������ݣ�</param>
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
    /// ��������
    /// </summary>
    /// <param name="solver">�����㷨����</param>
    /// <param name="n">�����ģ����״̬��������==΢�ַ��̽���</param>
    /// <param name="call_back">�ص�����</param>
    /// <param name="opt">������ѡ��</param>
    /// <param name="ivp_user_data">�û�����(������ڲ����ݣ����ݸ��ص�����call_back���㷨�������)</param>
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
            spw->m_data->m_initialStep = spw->m_opt.m_stopTime/1000;          //Ĭ�ϳ�ʼ���ֲ���default
            spw->m_data->m_Q = 0;
            spw->m_data->m_h = 0;

            spw->m_opt.m_maxStepSize = spw->m_opt.m_stopTime;
            
            if (spw->m_nStates > 0)
            {
                spw->m_data->m_preY = (MoReal*)sw->m_utils.m_allocDataMemory(sw->m_userData, 1, n * sizeof(MoReal));  //(MoReal *)ǿ��ת��
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

        return spw;   //���أ��������(�Զ����㷨�ڲ����ݣ������������������)����Ϊ�����ӿں����ĵڶ������������±ߵ�ivp
    }

    /// <summary>
    /// ��ʼ��
    /// </summary>
    /// <param name="solver">�����㷨����</param>
    /// <param name="ivp">�������</param>
    /// <param name="t0">��ʼʱ��</param>
    /// <param name="y0">y�ĳ�ʼֵ���ض���Ϊ�գ�</param>
    /// <param name="yp0">y���ĳ�ʼֵ��DAE��</param>
    /// <param name="is_reinit">�Ƿ����³�ʼ����������</param>
    /// <param name="reserve">�����������ݲ�ʹ��</param>
    /// <returns></returns>
    MwsInteger myRK45Init(MwsIVPSolverObj solver, MwsIVPObj ivp, MwsReal t0, const MwsReal* y0,
        const MwsReal* yp0, MwsBoolean is_reinit, void* reserve)
    {
        /* nothing to do */
        return MWS_IVP_SUCCESS;  //����״̬��ȡMwsIVPStatus��ֵ
    }

    /// <summary>
    /// ���
    /// </summary>
    /// ���룺
    /// <param name="solver">�����㷨����</param>
    /// <param name="ivp">�������</param>
    /// <param name="step_size">�����������������ʼ�����ֲ�����</param>
    /// <param name="t">��ǰʱ��</param>
    /// <param name="tout">�������ʱ��</param>
    /// �����
    /// <param name="tret">���ʵ�ʴﵽ��ʱ��</param>
    /// <param name="yret">y�Ľ��ֵ</param>
    /// <param name="ypret">y���Ľ��ֵ��DAE��</param>
    /// <param name="reserve">�����������ݲ�ʹ��</param>
    /// <returns></returns>
    MwsInteger myRK45Solve(MwsIVPSolverObj solver, MwsIVPObj ivp, MwsReal step_size, MwsReal t,
        MwsReal tout, MwsReal* tret, MwsReal* yret, MwsReal* ypret, void* reserve)
    {
        MyRK45* sw = (MyRK45*)solver;
        MyRK45Problem* spw = (MyRK45Problem*)ivp;

        MoSize index = 0;
        MoSize nState = spw->m_nStates;

        spw->m_data->m_curTime = t;
        spw->m_data->m_initialStep = step_size;             //��ʼ�����ֲ���

        MoReal* D = spw->m_data->m_D;                       //D=w(i+1)-y(i+1)
        MoReal Q = spw->m_data->m_Q;
        MoReal h = spw->m_data->m_initialStep;
        
        MoReal* preY = spw->m_data->m_preY;                 //�ϸ�y
        MoReal* curY = spw->m_data->m_curY;                //��ǰy

        MoReal* preYp = spw->m_data->m_preYp;               //�ϸ�y'
        MoReal* curYp = spw->m_data->m_curYp;               //��ǰy��

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

                //��ʽ����

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

            if (R < spw->m_opt.m_absoluteTolerance[0])            //����������㾫��
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

                if (h > spw->m_opt.m_maxStepSize)      //���������趨����󲽳�����h=hmax
                {
                    h = spw->m_opt.m_maxStepSize;
                }

                if (spw->m_data->m_curTime >= spw->m_opt.m_stopTime)       //��ǰ����ʱ������趨��ֹͣ����ʱ��
                {
                    *tret = tout;
                    Flag = moFalse;

                   // return MWS_IVP_TSTOP_RETURN;        //������ֹʱ�������
                }
                else if (spw->m_data->m_curTime + h > spw->m_opt.m_stopTime)     //ti+qh>L,�������߽磬��h=L-ti
                {
                    h = spw->m_opt.m_stopTime - spw->m_data->m_curTime;
                }
                if (h < 0)
                {
                    Flag = moFalse;

                    return MWS_IVP_FAIL;
                }
            

        }

            
                /* ���µ�ǰ����ʱ�� */
      // *tret += h;                //����ʱ������һ�����ⲿ����ѭ��                 
       return MWS_IVP_SUCCESS;  //����״̬��ȡMwsIVPStatus��ֵ
    }

    /// <summary>
    /// ���ɲ�ֵ��������yn������ߣ�
    /// </summary>
    /// ���룺
    /// <param name="solver">�����㷨����</param>
    /// <param name="ivp">�������</param>
    /// <param name="tout">���������ʱ��</param>
    /// �����
    /// <param name="yret">y�Ľ��ֵ</param>
    /// <param name="reserve"></param>
    /// <returns></returns>
    MwsInteger myRK45Interpolate(MwsIVPSolverObj solver, MwsIVPObj ivp, MwsReal tout, MwsReal* yret, void* reserve)
    {
        /* ���Բ�ֵ */
        MyRK45* sw = (MyRK45*)solver;
        MyRK45Problem* spw = (MyRK45Problem*)ivp;

        MoSize nStates = spw->m_nStates;
        MoSize index;

        for (index = 0; index < nStates; ++index)
        {
            /* ���Բ�ֵ:
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
    /// ��������
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
    /// ���ٻ����㷨
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

