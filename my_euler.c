/***************************************************************************
///
/// Copyright (c) 2020, ����ͬԪ�����Ϣ�������޹�˾
/// All rights reserved.
///
/// @file           my_euler.c
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

#ifdef __cplusplus
extern "C"{
#endif

/* �㷨���� */
typedef struct  
{
    MwsIVPUtilFcns	m_utils;
    void*           m_userData;

} MyEuler;

/* �����������ݣ��������ڴ������ͷź��� */
typedef struct  
{
    MoReal *m_preY;
    MoReal *m_curY;  
    MoReal *m_preYp;

    MoReal m_curTime;
    MoReal m_initialStep;
} MyEulerProblemData;

/* ���������� */
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
/// �����㷨
/// </summary>
/// <param name="util_fcns">���ߺ�������������ṩ���û����ã�</param>
/// <param name="user_data">�û����ݣ�������ڲ����ݣ����ݸ����ߺ���util_fcns���Զ����㷨������������ݣ�</param>
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
/// ��������
/// </summary>
/// <param name="solver">�����㷨����</param>
/// <param name="n">�����ģ����״̬��������==΢�ַ��̽���</param>
/// <param name="call_back">�ص�����</param>
/// <param name="opt">������ѡ��</param>
/// <param name="ivp_user_data">�û�����(������ڲ����ݣ����ݸ��ص�����call_back���㷨�������)</param>
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
            spw->m_data->m_preY = (MoReal *)sw->m_utils.m_allocDataMemory(sw->m_userData, 1, n*sizeof(MoReal));  //(MoReal *)ǿ��ת��
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
MwsInteger myEulerInit(MwsIVPSolverObj solver, MwsIVPObj ivp, MwsReal t0, const MwsReal* y0, 
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
MwsInteger myEulerSolve(MwsIVPSolverObj solver, MwsIVPObj ivp, MwsReal step_size, MwsReal t, 
    MwsReal tout, MwsReal* tret, MwsReal* yret, MwsReal* ypret, void* reserve)
{
    MyEuler* sw = (MyEuler*)solver;
    MyEulerProblem* spw = (MyEulerProblem*)ivp;

    MoSize index = 0;
    MoSize nState = spw->m_nStates;

    spw->m_data->m_curTime = t;
    spw->m_data->m_initialStep = step_size;

    MoReal* preY = spw->m_data->m_preY;                 //�ϸ�y
    MoReal* curY = spw->m_data->m_curY;                 //��ǰy
    MoReal* preYp = spw->m_data->m_preYp;               //�ϸ�y'
   
    for (index = 0; index<nState; ++index)              
    {
        preY[index] = yret[index];
        preYp[index] = ypret[index];

        MoReal* k1 = preYp[index];
        MoReal* k2 = spw->m_callback.m_rshFunction(spw->m_userData, spw->m_data->m_curTime + spw->m_data->m_initialStep, 
                                                    yret[index] + spw->m_data->m_initialStep * k1, ypret[index]);
                                                                //k2=f(xn+h,yn+h*k1)
        yret[index] += spw->m_data->m_initialStep / 2 * (k1 + k2);      //���ι�ʽ
        
        curY[index] = yret[index];
    }
    

    /* ���µ�ǰ����ʱ�� */
    *tret += spw->m_data->m_initialStep;                //����ʱ������һ����˵���˴���ÿһ���ڵĴ����ⲿ�ж��㷨��ѭ������

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
MwsInteger myEulerInterpolate(MwsIVPSolverObj solver, MwsIVPObj ivp, MwsReal tout, MwsReal* yret, void* reserve)
{
    /* ���Բ�ֵ */
    MyEuler* sw = (MyEuler*)solver;
    MyEulerProblem* spw = (MyEulerProblem*)ivp;

    MoSize nStates = spw->m_nStates;
    MoSize index;

    for (index=0; index<nStates; ++index)
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
/// ���ٻ����㷨
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

