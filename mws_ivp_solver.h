/***************************************************************************
///
/// Copyright (c) 2020, ����ͬԪ�����Ϣ�������޹�˾
/// All rights reserved.
///
/// @file           mws_ivp_solver.h
/// @brief          �����㷨
///
/// @version        v1.0
/// @author         ������
/// @date           2020/05/25
///
***************************************************************************/

#ifndef MWS_IVP_SOLVER_H
#define MWS_IVP_SOLVER_H

#include "mws_types.h"

MW_EXTERN_C_BEGIN

/* DAE����ʽΪ:
 *  F(t,y,y')=0
 *  y(t0)=y0
 *  y'(t0)=y'0
 * ODE����ʽΪ:
 *  y'=f(t,y)
 *  y(t0)=y0
 */

typedef void* MwsIVPSolverObj;  /* �����㷨���� */
typedef void* MwsIVPObj;        /* ���������� */

/**
 * @brief �����㷨��������״̬
 */
typedef enum
{
    MWS_IVP_SUCCESS = 0,    /* �ɹ� */
    MWS_IVP_WARNING = 1,    /* �о��� */
    MWS_IVP_TSTOP_RETURN,   /* ������ֹʱ������� */
    MWS_IVP_MEM_FAIL,       /* �ڴ����ʧ�� */
    MWS_IVP_FAIL,           /* ʧ�ܣ�ԭ��δ֪�� */
    MWS_IVP_INVALID_INPUT,  /* ������Ч */
    MWS_IVP_RHSFN_FAIL,     /* ODEͨ����ʽ���Ҷ˺������� */
    MWS_IVP_RESFN_FAIL,     /* DAE���ຯ������ */
}MwsIVPStatus;

/**
 * @brief �㷨�����������
 */
typedef enum
{
    MWS_IVP_ODE,    /* ODE���� */
    MWS_IVP_DAE,    /* DAE���� */
}MwsIVPType;

/**
 * @brief ���ߺ���(ƽ̨�ṩ�������㷨ʹ��)
 */
typedef struct
{
    /*
     * @breief ��־����
     * @param[in] user_data     �û����ݣ���������ӿڵ����һ����������ͬ
     * @param[in] error_code    �����
     * @param[in] where         �������Ժδ�
     * @param[in] msg           ������Ϣ
     */
    void (*m_logger)(void* user_data, MwsInteger error_code, MwsString where, MwsString msg);

    /*
     * @breief �ڴ���亯��
     * @param[in] user_data     �û�����
     * @param[in] nobj          ��������
     * @param[in] size          �����С���ֽ�����
     * @return  �ɹ��򷵻ط�����׵�ַ�����򷵻ؿ�
     */
    void* (*m_allocMemory)(void* user_data, MwsSize nobj, MwsSize size);

    /*
     * @breief �ڴ��ͷź���
     * @param[in] user_data     �û�����
     * @param[in] p             �ڴ��׵�ַ
     */
    void (*m_freeMemory)(void* user_data, void* p);

    /*
     * @breief �����ڴ���亯��
     * @param[in] user_data     �û�����
     * @param[in] nobj          ��������
     * @param[in] size          �����С���ֽ�����
     * @return  �ɹ��򷵻ط�����׵�ַ�����򷵻ؿ�
     */
    void* (*m_allocDataMemory)(void* user_data, MwsSize nobj, MwsSize size);

    /*
     * @breief �����ڴ��ͷź���
     * @param[in] user_data     �û�����
     * @param[in] p             �ڴ��׵�ַ
     */
    void (*m_freeDataMemory)(void* user_data, void* p);

}MwsIVPUtilFcns;

/**
 * @brief �����㷨����
 */
typedef struct MwsIVPSolverProp
{
    MwsString m_name;           /* ���ƣ������ִ�Сд������֤Ψһ */
    MwsString m_desc;           /* �㷨���� */
    MwsBoolean m_fixedStep;     /* �Ƿ񶨲��� */
    MwsIVPType m_ivpType;       /* �㷨����(ODE��DAE) */
}MwsIVPSolverProp;

/**
 * @brief �����㷨ѡ��
 */
typedef struct MwsIVPOptions
{
    MwsBoolean m_stopTimeDefined;       /* ��ֹʱ���Ƿ��� */
    MwsReal m_stopTime;                 /* ��ֹʱ�� */
    MwsBoolean m_toleranceDefined;      /* ��������Ƿ��� */
    MwsReal* m_relativeTolerance;       /* ������������� */
    MwsReal* m_absoluteTolerance;       /* �������������� */
    MwsBoolean m_maxStepSizeDefined;    /* ��󲽳��Ƿ��� */
    MwsReal m_maxStepSize;              /* ��󲽳� */

}MwsIVPOptions;

/*
 * �ص�������ƽ̨�ṩ�������㷨ʹ�ã�
 */

/*
 * @breief ODEͨ����ʽ�Ҷ˺��� y'=f(t,y)
 * @param[in]user_data �û����ݣ������������ӿڵ����һ����������ͬ
 * @param[in]t          ʱ��
 * @param[in]y          y��ֵ
 * @param[out]yp        y'��ֵ
 * @return ״̬��MwsIVPStatus���͵�ֵ
 */
typedef MwsInteger (*MwsIVPRshFcnPtr)(void* user_data, MwsReal t, const MwsReal* y, MwsReal* yp);

/*
 * @breief DAE����������ص�����
 * @param[in]user_data �û�����
 * @param[in]t          ʱ��
 * @param[in]y          y��ֵ
 * @param[in]yp         y'��ֵ
 * @param[out]delta     ������
 * @return ״̬��MwsIVPStatus���͵�ֵ
 */
typedef MwsInteger (*MwsIVPResFcnPtr)(void* user_data, MwsReal t, const MwsReal* y, 
    const MwsReal* yp, MwsReal* delta);

/**
 * @brief Jacobian����ص�����
 * @param[in]user_data  �û�����
 * @param[in]t          ʱ��      
 * @param[in]y          y��ֵ
 * @param[in]yp         y'��ֵ��DAE��
 * @param[in]cj         ������Jacobian�ı�����DAE��
 * @param[out]pd        ƫ����pd=dF/dy+cj*dF/dy'��pd=df/dy��    
 * @return ״̬��MwsIVPStatus���͵�ֵ 
 * @note 
 */
typedef MwsInteger(*MwsIVPJacFcnPtr)(void* user_data, MwsReal t, const MwsReal* y, const MwsReal* yp, 
    MwsReal cj, MwsReal* pd);

/*
 * @breief ���ֲ���ɻص�����
 * @param[in]user_data �û�����
 * @param[in]t          ʱ��
 * @param[in]y          y��ֵ
 * @return ״̬��MwsIVPStatus���͵�ֵ
 */
typedef MwsInteger (*MwsIVPStepFinishedPtr)(void* user_data, MwsReal t, const MwsReal* y);

/**
 * @brief �����㷨�ص�����
 */
typedef struct  
{
    MwsIVPRshFcnPtr m_rshFunction;          /* ODE�Ҷ˺��� */
    MwsIVPResFcnPtr m_resFunction;          /* DAE���ຯ������m_rshFunction���⣩*/
    MwsIVPJacFcnPtr m_jacFunction;          /* Jacobian���㺯�� */
    MwsIVPStepFinishedPtr m_stepFinished;   /* ���ֲ���ɻص����� */
}MwsIVPCallback;
 
/*
 * �����㷨�ӿں����������㷨�ṩ��ƽ̨ʹ�ã�
 */

/**
  * @breief �����㷨���󴴽�����
  * @param[in]util_fcns     ���ߺ���      
  * @param[in]user_data     �û�����(���ݸ����ߺ���)
  * @return �����㷨����
  * @note  ���ߺ����ض����ṩ
  */
typedef MwsIVPSolverObj (*MwsIVPSolverCreatePtr)(MwsIVPUtilFcns* util_fcns, void* user_data);

/**
 * @brief ���������󴴽�����
 * @param[in]solver         �����㷨����
 * @param[in]n              ��ģ��״̬��������      
 * @param[in]call_back      �ص�����
 * @param[in]opt            �����㷨ѡ��
 * @param[in]ivp_user_data  �û�����(���ݸ��ص�����)
 * @return ����������
 * @note  �ص������ض����ṩ 
 */
typedef MwsIVPObj (*MwsIVPCreatePtr)( MwsIVPSolverObj solver, MwsSize n, MwsIVPCallback* call_back,
    MwsIVPOptions* opt, void* ivp_user_data);

/*
 * @breief ��ʼ������
 * @param[in]solver     �����㷨����
 * @param[in]ivp        ����������
 * @param[in]t0         ��ʼʱ��
 * @param[in]y0         y�ĳ�ʼֵ���ض���Ϊ�գ�
 * @param[in]yp0        y'�ĳ�ʼֵ��DAE��
 * @param[in]is_reinit  �Ƿ����³�ʼ����������
 * @param[in]reserve    �����������ݲ�ʹ��
 * @return ״̬��ȡMwsIVPStatus��ֵ
 */
typedef MwsInteger (*MwsIVPInitPtr)(MwsIVPSolverObj solver, MwsIVPObj ivp, MwsReal t0, 
    const MwsReal* y0, const MwsReal* yp0, MwsBoolean is_reinit, void* reserve);

/*
 * @breief ��⺯��
 * @param[in]solver     �����㷨����
 * @param[in]ivp        ����������
 * @param[in]step_size  ������������/��ʼ�����ֲ�����
 * @param[in]t          ��ǰʱ��
 * @param[in]tout       ���������ʱ��   
 * @param[out]tret      ���ʵ�ʵ����ʱ��   
 * @param[out]yret      y�Ľ��ֵ
 * @param[out]ypret     y'�Ľ��ֵ��DAE��
 * @param[in]reserve    �����������ݲ�ʹ��
 * @return ״̬��ȡMwsIVPStatus��ֵ
 */
typedef MwsInteger (*MwsIVPSolvePtr)(MwsIVPSolverObj solver, MwsIVPObj ivp, MwsReal step_size, MwsReal t, 
    MwsReal tout, MwsReal* tret, MwsReal* yret, MwsReal* ypret, void* reserve);

/*
 * @breief ��ֵ����
 * @param[in]solver     �����㷨����
 * @param[in]ivp        ����������
 * @param[in]tout       ���������ʱ��   
 * @param[out]yret      y�Ľ��ֵ
 * @param[in]reserve    �����������ݲ�ʹ��
 * @return ״̬��ȡMwsIVPStatus��ֵ
 */
typedef MwsInteger (*MwsIVPInterpolatePtr)(MwsIVPSolverObj solver, MwsIVPObj ivp, MwsReal tout, 
    MwsReal* yret, void* reserve);
/*
 * @breief �������������ٺ���
 * @param[in]solver     �����㷨����
 * @param[in]ivp        ����������
 */
typedef void (*MwsIVPDestroyPtr)(MwsIVPSolverObj solver, MwsIVPObj ivp);

/*
 * @breief �����㷨�������ٺ���
 * @param[in]solver     �����㷨����
 */
typedef void (*MwsIVPSolverDestroyPtr)(MwsIVPSolverObj solver);

/**
 * @brief �����㷨�ӿں���
 */
typedef struct MwsIVPSolverFcns
{
    MwsIVPSolverCreatePtr m_createPtr;      /* ���󴴽��ĺ���ָ�� */
    MwsIVPCreatePtr       m_createPBPtr;    /* ������󴴽��ĺ���ָ�� */
    MwsIVPInitPtr m_initPtr;                /* ��ʼ���ĺ���ָ�� */
    MwsIVPSolvePtr m_solvePtr;              /* ���ĺ���ָ�� */
    MwsIVPInterpolatePtr m_interpolatePtr;  /* ��ֵ�ĺ���ָ�� */
    MwsIVPDestroyPtr m_destroyPBPtr;        /* ����������ٵĺ���ָ�� */
    MwsIVPSolverDestroyPtr m_destroyPtr;    /* �������ٵĺ���ָ�� */
}MwsIVPSolverFcns;


/*
 * @breief ע������㷨
 * @param[in] sim_data ����������
 * @param[in] prop     ���������
 * @param[in] fcns     ������ӿں���
 */
MoBoolean isimRegisterIVPSolver(void* sim_data, const MwsIVPSolverProp* prop, const MwsIVPSolverFcns* fcns);

/*
 * @breief �Ƴ�ע��Ļ����㷨
 * @param[in] sim_data  ����������
 * @param[in] name      �����㷨����
 */
MoBoolean isimUnregisterIVPSolver(void* sim_data, MwsString name);

struct MwsLSSolverProp;
struct MwsLSSolverFcns;
/*
 * @breief ע�����Է���������㷨
 * @param[in] sim_data ����������
 * @param[in] prop     ���������
 * @param[in] fcns     ������ӿں���
 */
MoBoolean isimRegisterLSSolver(void* sim_data, const struct MwsLSSolverProp* prop, const struct MwsLSSolverFcns* fcns);

/*
 * @breief �Ƴ�ע������Է���������㷨
 * @param[in] sim_data  ����������
 * @param[in] name      �����㷨����
 */
MoBoolean isimUnregisterLSSolver(void* sim_data, MwsString name);

struct MwsNLSSolverProp;
struct MwsNLSSolverFcns;
/*
 * @breief ע������Է���������㷨
 * @param[in] sim_data ����������
 * @param[in] prop     ���������
 * @param[in] fcns     ������ӿں���
 */
MoBoolean isimRegisterNLSSolver(void* sim_data, const struct MwsNLSSolverProp* prop, const struct MwsNLSSolverFcns* fcns);

/*
 * @breief �Ƴ�ע��ķ����Է���������㷨
 * @param[in] sim_data  ����������
 * @param[in] name      �����㷨����
 */
MoBoolean isimUnregisterNLSSolver(void* sim_data, MwsString name);

MW_EXTERN_C_END

#endif /* !MWS_IVP_SOLVER_H */

/***************************************************************************
//   end of file
***************************************************************************/

